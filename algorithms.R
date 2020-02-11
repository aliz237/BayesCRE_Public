## algorithms.R
##
## This file contains the functions necessary for pre-processing the KB, runing the inference, 
## and post-processing the output results for BayesCRE
## 
##
## Computational Sceinces CoE
## Pfizer Worldwide Research and Development
##
## Date: April 16, 2013
##
## Last Updated: Oct 10, 2015
##
## Author: Kourosh Zarringhalam
###############################################################################

## processKB:
##
## This function cleans up the KB. Reads the KB in and only consideres the network of Proteins, 
## Compounds and mRNAs.
##
## Arguments:
##
## ents.all.file     -  entry files from the knowledge base
##
## rels.all.file     -  relation files from the knowledge base
##
## L                 -  the returned list of pruned entries and relations
processKB = function(ents.all.file, rels.all.file ){
  ## Note: There are ' and # characters in KB that mess up the tables. The following will clean them up.
  ents.all = read.table(ents.all.file, header = T, stringsAsFactors = F, strip.white=T, sep = '\t', 
                        quote = NULL, comment.char = '')
  colnames(ents.all) = c('uid', 'name', 'id', 'type')
  rels.all = read.table(rels.all.file, header = T, stringsAsFactors = F, strip.white=T, sep = '\t', 
                        quote = NULL, comment.char = '')
  colnames(rels.all) = c('uid', 'srcuid', 'trguid', 'type', 'pmids', 'nls')
  
  ## Remove the ", ', # characters
  ents.all = data.frame(apply(ents.all, 2, function(x) gsub('\"', '', x)), stringsAsFactors = F)
  ents.all = data.frame(apply(ents.all, 2, function(x) gsub('\'', 'p', x)), stringsAsFactors = F)
  ents.all = data.frame(apply(ents.all, 2, function(x) gsub('#', '_', x)), stringsAsFactors = F)
  
  rels.all = data.frame(apply(rels.all, 2, function(x) gsub('\"', '', x)), stringsAsFactors = F)
  rels.all = data.frame(apply(rels.all, 2, function(x) gsub('\'', 'p', x)), stringsAsFactors = F)
  rels.all = data.frame(apply(rels.all, 2, function(x) gsub('#', '_', x)), stringsAsFactors = F)
  
  ## Get rid of ambiguities (edges with conflicts in signs)
  conf.ind = which(rels.all$type == "conflict")
  print(paste('total number of edges:', length(rels.all$type)))
  print(paste('number of ambiguous edges:', length(conf.ind), ' (', 
              round(100 * length(conf.ind)/length(rels.all$type), 3), 'percent.)'))
  print('removing ambiguous edges.')
  rels.all = rels.all[which(rels.all$type != "conflict"), ]
  rownames(rels.all) = 1:nrow(rels.all)
  
  ##Map the ids to integeres
  uid.orig = ents.all$uid
  uid.new = seq(1, length(ents.all$uid))
  id.map = data.frame(uid.orig = uid.orig, uid.new = uid.new, stringsAsFactors = F)
  ents.all$uid = uid.new
  rels.all$srcuid = uid.new[match(rels.all$srcuid, uid.orig)]
  rels.all$trguid = uid.new[match(rels.all$trguid, uid.orig)]
  
  L = list(ents.all = ents.all, rels.all = rels.all, id.map = id.map)
  
  L
}

## genOneLevel:
##
## Given a knowledge-base, this function generates a one-level causal network from the protein/compound
## /mRNA entries. The hypothesis layer consists of protiens and chemicals while the evidence layer 
## consistes of mRNAs only. The i-->i and cyles are removed from the causal network as well.
##
## Arguments:
##
## ents.all.file     -  entry files from the knowledge base
##
## rels.all.file     -  relation files from the knowledge base
##
## L                 -  the returned list of pruned entries and relations
genOneLevel = function(ents.all.file, rels.all.file ){
  ## pre-process KB
  print('processing the network...') 
  L = processKB(ents.all.file, rels.all.file)
  ents.all = L$ents.all
  rels.all = L$rels.all
  
  ## Protein, Compound or mRNA entries
  ents = unique(ents.all[which(ents.all$type %in% c('mRNA', 'Protein', 'Compound')),])
  
  ## (unique) Protein/Compound entries
  ents.pc = unique(ents[which(ents$type %in% c('Protein', 'Compound')),]) 
  ## (unique) mRNA entries
  ents.mRNA = unique(ents[which(ents$type == 'mRNA'),]) 
  
  ## src has to be protiens or a compunds
  rels = rels.all[which(rels.all$srcuid %in% ents.pc$uid),]
  rownames(rels) = 1:nrow(rels)
  ## target has to be mRNA
  rels = rels[which(rels$trguid %in% ents.mRNA$uid),]
  rownames(rels) = 1:nrow(rels)
  
  ## Only consider entries that are participating in valid relations
  ents = ents[which(ents$uid %in% c(rels$srcuid, rels$trguid)), ]
  rownames(ents) = 1:nrow(ents)
  
  rels = rels[!duplicated(rels[,c(2,3)]),]
  rownames(rels) = 1:nrow(rels)
   
  
  print('processed network dimenssions:')
  print(paste('ents:', dim(ents)[1]))
  print(paste('rels:', dim(rels)[1]))
    
  L = list(ents = ents, rels = rels)
 
  L
}


## genContextGraph:
## Given a network of values and relationships, this function generates a Bayesian network by first, 
## removing the nodes (relationships) that are more than one level away from the targets, removing 
## i->i relations, performing enrichment analysis on MeSH terms and adding the context nodes as 
## well as the applicability nodes per edge. 
##
## Arguments
##
## ents.all.file     -       entries
##
## rels.all.file     -       relationships
##
## evidence.file     -       evidence file (gene expression data)
##
## MeSH.file         -       MeSH file
##
## cutoff.q          -       the cutoff q-value for MeSH enrichment analysis
##
## Max.MeSH          -       Max number of edges that a MeSH term can be connect to.
##                           (This is to reduce the inference time)
##
## L                 -       output list of context specific entiries and relations
genContextGraph = function(ents.all.file, rels.all.file, evidence.file, MeSH.file, 
                           cutoff.q, Max.MeSH, minscore ){
  
  ## Generating the one level-network form the knowledge base
  print('pre-processing...')
  L = genOneLevel(ents.all.file, rels.all.file)
  ents = L$ents
  rels = L$rels
  ents.mRNA = ents[which(ents$type == 'mRNA'),]
  
  ##repeat rows per pmid
  rels.pmid = lapply(rels$pmids, function(x) strsplit(x, '\\|')[[1]])
  pmid.nums = unlist(lapply(rels.pmid, length))
  rel.names = colnames(rels)
  rels = data.frame(rels[rep(1:nrow(rels), pmid.nums), c(1:4)], unlist(rels.pmid), stringsAsFactors = F)
  colnames(rels) = rel.names[1:5]
  rownames(rels) = 1:nrow(rels)
  
  ## Read in the evidence. Evidence has to be a valid file. Either 3 columns or 2 columns
  ## Check for validity.
  evidence = read.table(evidence.file, header = T, sep = '\t', stringsAsFactors = F)
  if(ncol(evidence) == 3){
    pval.ind = grep('pval|p.val|p-val|P-val', colnames(evidence), ignore.case = T) ## Looking for p-value
    fc.ind = grep('fc|FC|fold', colnames(evidence), ignore.case = T) ## Looking for fold change
    id.ind = grep('id|entr', colnames(evidence), ignore.case = T) ## Looking for Entrez ids
    if(length(id.ind) == 0 | length(fc.ind) == 0 | length(pval.ind) == 0){
      print('Please make sure the expression files column names are labled as entrez, fc, pvalue')
      quit(save = "no", status = 1, runLast = FALSE)
    }
    
    p.thresh = 0.01 ## Cutoff threshold for p-value and fold change
    isLog = grep('log', colnames(evidence)[fc.ind], ignore.case = T)
    if(length(isLog) > 0){
      fc.thresh = log2(1.3)
    }else{
      fc.thresh = 1.3
    }
    
    ## Generate evidence file (Up or Down regulated)
    evidence = evidence[which(abs(evidence[,fc.ind]) >= fc.thresh & 
                                evidence[,pval.ind] <= p.thresh ),c(id.ind, fc.ind)]
    evidence[,2] = ifelse(evidence[,2] > 0, 1, -1)
    
    evidence = data.frame(id = suppressWarnings(as.numeric(evidence[,1])), val = as.numeric(evidence[,2]),
                          stringsAsFactors = F)
    evidence = evidence[!is.na(evidence[,1]),]
    ## The file maybe already in the right format, Entrez ID, Up/Down regulation
  }else if(ncol(evidence) == 2){
    id.ind = grep('id|entr|Entr', colnames(evidence))
    val.ind = grep('val', colnames(evidence))
    if(length(id.ind) == 0 | length(val.ind) == 0){
      print('Please provide a valid expression data file:')
      print('1: a file containing columns: entrez, fc, pvalue')
      print('2: a file containing columns: entrez, val, whre val is 1, -1, or 0')
      quit(save = "no", status = 1, runLast = FALSE)
    }
    colnames(evidence) = c('id', 'val')
  }else{
    print('Please provide a valid expression data file:')
    print('1: a file containing columns: entrez, fc, pvalue')
    print('2: a file containing columns: entrez, val, whre val is 1, -1, or 0')
    quit(save = "no", status = 1, runLast = FALSE)
  }
  ## Remove duplacated genes
  evidence = evidence[!duplicated(evidence[,1]), ]

  n.e1 = nrow(evidence)
  ## Make sure evidence id is in the ents. Same gene with multiple ids may exist
  ## in the original ents. 
  evidence = evidence[which(evidence[,1] %in% ents.mRNA$id),]
  n.e2 = nrow(evidence)
  print(paste((n.e1-n.e2), "evidence removed!"))
  
  ##Change id back to uid
  evidence.tmp = merge(evidence, ents.mRNA, by.x = 1, by.y = 3)
  if(nrow(evidence.tmp) > nrow(evidence)){
    print('Warning')
    print('Entrez ids in evidence file mapped to multiple uids in KB')
    print(paste('total number of duplicates:' , (nrow(evidence.tmp) - nrow(evidence))))
    print('selecting one uid at random')
  }
  evidence = data.frame(uid = evidence.tmp$uid, val = evidence.tmp$val, stringsAsFactors = F)
  evidence = evidence[!duplicated(evidence),]

  ## Read in the MeSH terms
  MeSH = read.table(MeSH.file, header = F, stringsAsFactors = F, sep = '\t')
  colnames(MeSH) = c('pmids', 'MeSH')
  
  ## Generating the in context network.
  print('Performing MeSH term enrichment analysis.....')
  
  ## zero and nonzero networks
  targets         = unique(rels[,3])
  zero.targets    = targets[which(!(targets %in% evidence[,1]))]
  nonzero.targets = targets[which((targets %in% evidence[,1]))]
  
  ## Get rid of hyps that do not make at least min.score correct predictions
  all.hyps = unique(rels$srcuid)
  get.num.correct = function(x){
    xx = merge(unique(rels[which(rels$srcuid == x), c('trguid', 'type')]), evidence, by.x = 1, by.y = 1)
    up = length(which(xx[,2] == 'increase' & xx[,3] == 1)) + length(which(xx[,2] == 'decrease' & xx[,3] == -1))
    down = length(which(xx[,2] == 'increase' & xx[,3] == -1)) + length(which(xx[,2] == 'decrease' & xx[,3] == 1))
    num.correct = max(up, down)
  }
  hyp.stat = unlist(lapply(all.hyps, function(x) get.num.correct(x)))
  nonzero.hyps = all.hyps[which(hyp.stat >= minscore)]
  ## Check to see if any are left. If not lower the min.score
  if(length(nonzero.hyps) < 10){
    nonzero.hyps = all.hyps[which(hyp.stat >= floor(minscore/2))]
  }
  if(length(nonzero.hyps) < 10){
    nonzero.hyps = all.hyps[which(hyp.stat > 0)]
  }
  
  ## Update the rels  
  rels = rels[which(rels$srcuid %in% nonzero.hyps), ]

  
  zero.rels = merge(rels, zero.targets, by.x = 3, by.y = 1)
  zero.rels = zero.rels[,c('uid', 'srcuid', 'trguid', 'type', 'pmids')]
  
  nonzero.rels = merge(rels, nonzero.targets, by.x = 3, by.y = 1)
  nonzero.rels = nonzero.rels[,c('uid', 'srcuid', 'trguid', 'type', 'pmids')]
  
  ## update the rels: This is because some pmids in KB do not have any match
  ## in the MeSH file. This is done for zero.rels only as nonzero.rels will
  ## be annotated to an artificail MeSH term
  zero.rels = zero.rels[which(zero.rels$pmids %in% MeSH$pmids), ] ## KZ: March 10, 2015
  

  rels = rbind(nonzero.rels, zero.rels)
  

  ## MeSH terms assocaited with nonzero network
  nonzero.MeSH = unique(MeSH$MeSH[which(MeSH$pmids %in% nonzero.rels$pmids)])
  ##nonzero.MeSH = unique(merge(MeSH, nonzero.rels, by.x = 1, by.y = 5)[,2])
  
  ## All pmids annotated to non.zero MeSH tersm
  pmid.nonzero.MeSH = merge(MeSH, nonzero.MeSH, by.x = 2, by.y = 1)
  colnames(pmid.nonzero.MeSH) = c('MeSH', 'pmids')
  
  ## uids annotated with nonzero MeSH terms in non-zero network
  nonzero.uid.nonzero.MeSH = merge(pmid.nonzero.MeSH, nonzero.rels, by.x = 2, by.y = 5)
  nonzero.uid.nonzero.MeSH = nonzero.uid.nonzero.MeSH[,c('uid', 'srcuid', 'trguid', 'type', 'pmids', 'MeSH')]
  
  ## uids annotated with nonzero MeSH terms in zero network
  zero.uid.nonzero.MeSH = merge(pmid.nonzero.MeSH, zero.rels, by.x = 2, by.y = 5)
  zero.uid.nonzero.MeSH = zero.uid.nonzero.MeSH[,c('uid', 'srcuid', 'trguid', 'type', 'pmids', 'MeSH')]
  
  
  ## uids relations annotated to nonzero Mesh terms in all network
  all.uid.nonzero.MeSH = merge(pmid.nonzero.MeSH, rels, by.x = 2, by.y = 5)
  all.uid.nonzero.MeSH = all.uid.nonzero.MeSH[,c('uid', 'srcuid', 'trguid', 'type', 'pmids', 'MeSH')]
  
  ## Counting the frquencies of the nonzero MeSH terms in non-zero network
  ## (prevent double counting because of multiple pmids)
  tmp.nonzero.uid.nonzero.MeSH = nonzero.uid.nonzero.MeSH[!duplicated(nonzero.uid.nonzero.MeSH[,c(1,2,3,4,6)]),]
  ## Counting the frquencies of the nonzero MeSH terms in non-zero network
  rle.nonzero.uid.nonzero.MeSH = rle(sort(tmp.nonzero.uid.nonzero.MeSH[,'MeSH']))
  terms.nonzero.uid.nonzero.MeSH = rle.nonzero.uid.nonzero.MeSH$values
  freqs.nonzero.uid.nonzero.MeSH = rle.nonzero.uid.nonzero.MeSH$lengths
  
  ## compute the frequecy of the nonzero MeSH terms in all network
  tmp.all.uid.nonzero.MeSH = all.uid.nonzero.MeSH[!duplicated(all.uid.nonzero.MeSH[,c(1,2,3,4,6)]),]
  rle.all.uid.nonzero.MeSH = rle(sort(tmp.all.uid.nonzero.MeSH[,'MeSH']))
  terms.all.uid.nonzero.MeSH = rle.all.uid.nonzero.MeSH$values ## Should be the same as nonzero
  freqs.all.uid.nonzero.MeSH = rle.all.uid.nonzero.MeSH$lengths
  terms.freqs.uid.nonzero.MeSH = cbind(terms.all.uid.nonzero.MeSH, 
                                       freqs.nonzero.uid.nonzero.MeSH, 
                                       freqs.all.uid.nonzero.MeSH)
  
  
  
  ## Testing enrichment of MeSH terms with hypergeometric distribution
  ## q - total number of causal relations in the non-zero network annotated with the MeSH term
  ## m - total number of causal relations in the network annotated to the MeSH term
  ## n - total number of causal relations in the network Not annotated with the MeSH term
  ## k - total number of non-zero causal relations.
  ## N - total number of causal relations
  
  N = dim(rels)[1]
  k = dim(nonzero.rels)[1]
  
  pvals = {}
  for(t in 1:dim(terms.freqs.uid.nonzero.MeSH)[1]){
    MeSH.term = terms.freqs.uid.nonzero.MeSH[t,1]
    q = as.numeric(terms.freqs.uid.nonzero.MeSH[t,2])
    m = as.numeric(terms.freqs.uid.nonzero.MeSH[t,3])
    n = N - m
    pval = phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
    pvals = c(pvals, pval)
  }
  
  ## Adjusted p-values with fdr
  fdr = fdrtool(pvals, statistic = 'pvalue', plot=F, verbose=F)
  terms.freqs.uid.nonzero.MeSH = cbind(terms.freqs.uid.nonzero.MeSH, pvals, fdr$qval)
  colnames(terms.freqs.uid.nonzero.MeSH) = c('MeSH', 'in-net', 'all', 'pval', 'qval')

  ##Note: qval of 0.05 means that 5% of the "significant" discoveries are false
  ind = which((as.numeric(terms.freqs.uid.nonzero.MeSH[,5])) <= cutoff.q)
  if(length(ind) == 0){
    print('lowering significnace threshold for MeSH terms')
    ind = which((as.numeric(terms.freqs.uid.nonzero.MeSH[,5])) <= 0.05)
  }
  
  ## If the MeSH term is connected to more than Max.MeSH edges remove it.
  ## This is to speed up the inference.
  ind.NN = which((as.numeric(terms.freqs.uid.nonzero.MeSH[ind,3])) <= Max.MeSH)
  M = data.frame(terms.freqs.uid.nonzero.MeSH[ind[ind.NN],1], 
                 as.numeric(terms.freqs.uid.nonzero.MeSH[ind[ind.NN],2]) , 
                 as.numeric(terms.freqs.uid.nonzero.MeSH[ind[ind.NN],3]) , 
                 terms.freqs.uid.nonzero.MeSH[ind[ind.NN],4],
                 terms.freqs.uid.nonzero.MeSH[ind[ind.NN],5], stringsAsFactors = F)
  
  colnames(M) = c('MeSH', 'in-net-edges', 'overal-edges', 'p-val', 'q-val')
  
  ind = sort(as.numeric(M[,2]), index.return = T)$ix
  M.sorted = M[ind,]
  
  ## uids in the entire network annotated with at least 1 top MeSH term
  all.uid.top.MeSH = merge(all.uid.nonzero.MeSH, M.sorted, by.x = 6, by.y = 1)
  all.uid.top.MeSH[,7:10] = apply(all.uid.top.MeSH[,7:10], 2, as.numeric)
  
  ## KZ: March 10, 2015 concatinate the pmids back
  rep.uids = unique(sort(all.uid.top.MeSH$uid))
  rep.uids.pmids = lapply(rep.uids, function(x) unique(sort(all.uid.top.MeSH$pmids[which(all.uid.top.MeSH$uid == x)])))
  rep.uids.pmids = unlist(lapply(rep.uids.pmids, function(x) paste(x, collapse = '|')))
  all.uid.top.MeSH = all.uid.top.MeSH[!duplicated(all.uid.top.MeSH[,c(1,2,3,4,5)]), ]
  all.uid.top.MeSH = merge(all.uid.top.MeSH, data.frame(rep = rep.uids, rep.pmids = rep.uids.pmids), by.x = 2, by.y = 1)
  all.uid.top.MeSH = all.uid.top.MeSH[, c('MeSH', 'uid', 'srcuid', 'trguid', 'type','rep.pmids',
                                          'in-net-edges', 'overal-edges','p-val', 'q-val')]
  colnames(all.uid.top.MeSH)[6] = 'pmids' 
  ## Add an additional MeSH term to represent the non-zero edges. This MeSH term is
  ## Added to keep all the nonzero edges in the network
  tmp.nonzero.rels = nonzero.rels[!duplicated(nonzero.rels[,1:4]), ]
  nonzero.rels.extra.MeSH = data.frame(rep('nonzeroMeSH', dim(tmp.nonzero.rels)[1]), tmp.nonzero.rels, 
                                       rep(dim(tmp.nonzero.rels)[1], dim(tmp.nonzero.rels)[1]), 
                                       rep(dim(tmp.nonzero.rels)[1], dim(tmp.nonzero.rels)[1]), 
                                       rep(0, dim(tmp.nonzero.rels)[1]), rep(0, dim(tmp.nonzero.rels)[1]), 
                                       stringsAsFactors = F)
  
  colnames(nonzero.rels.extra.MeSH) = colnames(all.uid.top.MeSH)
  
  ##update all.uids
  all.uid.top.MeSH = rbind(all.uid.top.MeSH, nonzero.rels.extra.MeSH)
  
  ## sort by uid
  ind = sort(as.character(all.uid.top.MeSH[,2]), index.return = T)$ix
  rels.top.MeSH = all.uid.top.MeSH[ind,c(2,3,4,5,6,1)]
  rownames(rels.top.MeSH) = 1:dim(rels.top.MeSH)[1]
  ents.top.MeSH = ents[which(ents[,1] %in% c(rels.top.MeSH[,2], rels.top.MeSH[,3])),]
  rownames(ents.top.MeSH) = 1:dim(ents.top.MeSH)[1]
  
  print('MeSH enriched network dimenssions:')
  print(dim(ents.top.MeSH))
  print(dim(rels.top.MeSH))
  
  ## Adding the MeSH term nodes to the network
  names = unique(as.character(rels.top.MeSH[,6]))
  ##uid = (max(as.numeric(ents.top.MeSH$uid)) + 1):(max(as.numeric(ents.top.MeSH$uid)) + length(names))
  uid = paste('MeSHuid',1:length(names), sep = '')
  type = rep('MeSH', length(names))
  MeSH.data.frame = data.frame(cbind(uid = uid, name = names, id = names, type = type), stringsAsFactors = F)
  ents = rbind(ents.top.MeSH, MeSH.data.frame)
  
  print('Adding Applicability and Context nodes to the network')
  ## Adding the Applicability nodes to the network
  urels = unique(rels.top.MeSH[,2:3])
  names = paste(urels[,1], urels[,2], sep = '_')
  uid = paste('APPuid',1:length(names), sep = '')
  ##uid = (max(as.numeric(ents$uid)) + 1):(max(as.numeric(ents$uid)) + length(names))
  type = rep('Applicability', length(names))
  App.data.frame = data.frame(cbind(uid = uid, name = names, id = names, type = type), stringsAsFactors = F)
  ents = rbind(ents, App.data.frame)
  
  ## Adding the Mesh term ids to rels
  rels.top.MeSH.MeSH.ids = merge(rels.top.MeSH, MeSH.data.frame, by.x = 6, by.y = 2)
  
  ## Adding the Applicability ids to rels
  App.names = paste(rels.top.MeSH.MeSH.ids[,3], rels.top.MeSH.MeSH.ids[,4], sep = '_')
  rels.top.MeSH.MeSH.ids.App.names = data.frame(rels.top.MeSH.MeSH.ids, App.names, stringsAsFactors = F)
  
  rels.top.MeSH.App.ids = merge(rels.top.MeSH.MeSH.ids.App.names, App.data.frame, by.x = 10, by.y = 2)
  
  rels = rels.top.MeSH.App.ids[,c(3,4,5,6,7,8,11)]
  colnames(rels) = c('uid', 'srcuid', 'trguid', 'type', 'pmids', 'meshid', 'appid')
  
  ## Adding the index of the src, trg, and mesh in the ents to rels data frame
  src.ind  = lapply(rels[,2], function(x) which(ents[,1] == x))
  trg.ind  = lapply(rels[,3], function(x) which(ents[,1] == x))
  mesh.ind = lapply(rels[,6], function(x) which(ents[,1] == x))
  app.ind  = lapply(rels[,7], function(x) which(ents[,1] == x))
  
  rels = data.frame(rels, unlist(src.ind), unlist(trg.ind), unlist(mesh.ind), unlist(app.ind), stringsAsFactors = F)
  colnames(rels) = c('uid', 'srcuid', 'trguid', 'type', 'pmids', 'meshid', 'appid', 'srcind', 'trgind', 'meshind', 'appind')
  
  
  ## Update the evidence: Some may have been removed from ents because of restriction
  ## to hyps with at least score n0
  ents.mRNA = ents[which(ents$type == 'mRNA'), ]
  evidence = evidence[which(evidence$uid %in% ents.mRNA$uid),]
  
  L = list(ents = ents, rels = rels, evidence = evidence)
  
  L
}  


## The following few functions are auxilary function for speeding up the inference.
########### Start of auxiliary functions
addToChain = function(...) {
  if (LOG_CHAIN == TRUE) {
    numberchain <<- c(numberchain, ...);
  }
}

which_central = function(x) { which(x); }

which_1 = function(x) { which_central(x); }
which_2 = function(x) { which_central(x); }
which_3 = function(x) { which_central(x); }
which_4 = function(x) { which_central(x); }
which_5 = function(x) { which_central(x); }
which_6 = function(x) { which_central(x); }
which_7 = function(x) { which_central(x); }
which_8 = function(x) { which_central(x); }
which_9 = function(x) { which_central(x); }
which_10 = function(x) { which_central(x); }
which_11 = function(x) { which_central(x); }
which_12 = function(x) { which_central(x); }
which_13 = function(x) { which_central(x); }
which_14 = function(x) { x; } 
which_15 = function(x) { which_central(x); }
which_16 = function(x) { which_central(x); }
which_17 = function(x) { which_central(x); }
which_18 = function(x) { which_central(x); }
which_19 = function(x) { which_central(x); }

loc_1 = function(x) { x };
loc_2 = function(x) { x };
loc_3 = function(x) { x };
loc_4 = function(x) { x };
loc_5 = function(x) { x };
loc_6 = function(x) { x };
loc_7 = function(x) { x };
loc_8 = function(x) { x };
loc_9 = function(x) { x };
loc_10 = function(x) { x };
loc_11 = function(x) { x };

unique_1 = function(x) { unique(x) };
unique_2 = function(x) { unique(x) };
unique_3 = function(x) { unique(x) };
unique_4 = function(x) { unique(x) };
unique_5 = function(x) { unique(x) };
######## End of auxiliary functions

## This function computes the the combinitorial part of the sum in calculating P(H | X)
comb = function(n.p, n.m, k, l, p.c){
  n = n.p + n.m
  x = choose(n.p, k) * choose(n.m, l) * (((1-p.c)^(k+l))*(p.c)^(n - k - l))
  y = abs(k - l) / ((k+l)^2)
  B = x * k * y
  C = x * l * y
  BC = list(B=B,C=C)
  BC
}

## This function computes P(H | X)
compProbH = function(n.m, n.p, p.c, p.m, p.z, p.p){
  n = n.m + n.p
  
  if(n > 0){
    phm = (p.c^n) * p.m
    ph0 = (p.c^n) * p.z
    php = (p.c^n) * p.p
    pha = 0
    if(n.m >= 1){
      phm = phm + p.c^n.p - p.c^n
    }
    if(n.p >= 1){
      php = php + p.c^n.m - p.c^n
    }
    B = 0
    C = 0
    if(n.p >= 1 && n.m >= 1){
      for(k in 1:n.p){
        for(l in 1:n.m){
          x = comb(n.p, n.m, k, l, p.c)
          B = B + x$B
          C = C + x$C
        }
      }
      A = B + C
      phm = phm + C
      php = php + B
      pha = (1 - (p.c^n.m + p.c^n.p - p.c^n)) - A
    }
    ph = c(phm,ph0,php,pha)
  }else{## All parents are zero or ambigues
    ph = c(p.m, p.z, p.p, 0)
  }
  
  ph
}

## This function computes P(Z | H) 
pZ = function(Z, alpha, beta, p.m, p.z, p.p){
  if(Z == 0){
    #pz = c(beta,1-2*alpha,beta,p.z)
    pz = c(beta,1-2*alpha,beta,(1/3))
  }else if(Z == -1){
    #pz = c(1-2*beta,alpha,beta,p.m)
    pz = c(1-2*beta,alpha,beta,(1/3))
  }else if(Z == 1){
    #pz = c(beta,alpha,1-2*beta,p.p)
    pz = c(beta,alpha,1-2*beta,(1/3))
  }
  pz
}

## gibbsSampler
## This is the main function that performs the gibbs sampling inference:
## Given a Bayesian Network and evidence values, this function generates samples from the 
## joint probability distribution of the BN using Gibbs sampling. The hidden state of the 
## transcripts are implicitly added to the network and the FP and FN rates are factored in
## by merging the hidden and the observed states. The context nodes are added according to 
## the enriched MeSH terms and applicability nodes are added per edge. Hypothesis are selected 
## in an iterative fashion. In the second iteration, the selected hypothesis of the first iteration 
## are set to zero as evidence and so on. The function has been optimized to do more pre-processing 
## of lookup tables and requires the plyr package.
##
## Arguments:
##
## chain     -   Number of chains to be run (in parallel)
##
## ents      -   Matrix of the entries
##
## rels      -   Matrix of relationships
##
## evidence  -   Value of evidence nodes (mRNA, must be in the same order as in ents) 
##
## sim.num   -   Number of simulations
##
## burn.in   -   Burn in samples.
##
## iter      -   number of iteration
##
## p.c       -   probability of the edge being wrong for nonzero edges
##
## p.a       -   probability of the edge being wrong for zero edges
##
## p.m       -   prior probability of the ture state of the genes being -1
##
## p.z       -   prior probability of the ture state of the genes being 0
##
## p.p       -   prior probability of the ture state of the genes being 1
##
## h.z       -   prior probability of the hypothesis being off
##
## alpha     -   false positive rate
##
## beta      -   false negative rate
##
## w         -   weight of the sigmoid
##
## hyps.ind  -   index of the selected hypothes in the previous iterations
##
## results   -   list of marginal probability and values of non-evidence nodes.

gibbsSampler = function(chain, sim.num, burn.in, iter.num, p.c, p.a, p.m, 
                                     p.z, p.p, h.z, c.p, alpha, beta, w, hyps.ind){
  ## Initialize the non-evidence variables
  fix.nonzero.MeSH.ind = which(ents[,2] == 'nonzeroMeSH')
  pcma.ind = which(ents[,'type'] %in% c('Protein', 'Compound', 'MeSH', 'Applicability'))
  pcma.ind = pcma.ind[-which(pcma.ind == fix.nonzero.MeSH.ind)]
  pc.ind   = which(ents[,'type'] %in% c('Protein', 'Compound'))
  ma.ind   = which(ents[,'type'] %in% c('MeSH', 'Applicability'))
  ma.ind   = ma.ind[-which(ma.ind == fix.nonzero.MeSH.ind)]
  m.ind    = which(ents[,'type'] == 'MeSH')
  m.ind    = m.ind[-which(m.ind == fix.nonzero.MeSH.ind)]
  a.ind    = which(ents[,'type'] == 'Applicability')
  
  if(iter.num > 1){
    if(length(hyps.ind) == 0){
      print(paste('no hypothesis selected at iteration ', (iter-1)))
      results = {}
      results
    }
  }
  

  num.pcma = length(pcma.ind)
  num.pc   = length(pc.ind)
  num.ma   = length(ma.ind)
  num.a    = length(a.ind)
  num.m    = length(m.ind)
  
  ## X holds the current node values in the same order as in ents
  X = rep(0, length(ents[,1])) ## Current sample
  mrna.ind = which(ents[,'type'] == 'mRNA')
  num.mrna = length(mrna.ind)
  
  initialize.pc = rep(0, num.pc) 
  ##initialize.pc = sample(c(-1,0,1), num.pc, replace = T) 
  initialize.ma  = rep(0, num.ma)
  ##initialize.ma  = rep(1, num.ma)
  ##initialize.ma = sample(c(0,1), num.ma, replace = T) 
  
  ## Set the value of evidence
  X[match(evidence[,'uid'], ents[,'uid'])] = as.numeric(as.character(evidence[,'val']))
  
  ## nonzero context nodes set to 1
  X[fix.nonzero.MeSH.ind] = 1 
  
  X[pc.ind] = initialize.pc
  X[ma.ind] = initialize.ma
  X[hyps.ind] = 0 ## From previous iteration
  
  ## Marginals
  marg.probs.pc = matrix(0, nrow = 3, ncol = num.pc)
  marg.probs.a  = matrix(0, nrow = 2, ncol = num.a)
  marg.probs.m  = matrix(0, nrow = 2, ncol = num.m)
   
  samp.num.pc = 0
  samp.num.m  = 0
  samp.num.a  = 0
  
  ents[,2:4] = 0;
  rels[,"type"] = ifelse(rels[,"type"] == 'increase', 1, -1)
  rels[,"type_srcind"] = rels[,"type"] * as.double(rels[,"srcind"])

  ### Prepare relationships lookup tables
  ### Each table is *stratified* by a variable for quick lookup.
  ### Result is a list that contains data.frames under the name of the lookup variable.
  ### needs plyr package to work
  
  relsBySrcUid = dlply(rels, .(srcuid))
  relsByTrgUid = dlply(rels, .(trguid))
  relsByMeSHId = dlply(rels, .(meshid))
  relsByAppId = dlply(rels, .(appid))
  
  ### To make subselection fast in some cases, two lookup "data.frames" are converted to matrices!
  ##relsByTrgUid = lapply(relsByTrgUid, data.matrix)
  ##relsByAppId = lapply(relsByAppId, data.matrix)
  
  print("Sampling starting...")
  for(t in 1:(sim.num + burn.in - 1)){
    
    if ((t %% 10000) == 0) {
      print(paste("Iteration ", t, " / ", (sim.num + burn.in - 1)));
    }
    
    i.pc = pc.ind[(t %% length(pc.ind) + 1)]
    if(i.pc %in% hyps.ind){
      X[i.pc] = 0
    }else{
      addToChain( i.pc);
      
      ## Determining the network of the childrens of the hypothesis (Markov Blanket)
      i.pc.ch = unique(relsBySrcUid[[ents[i.pc,1]]][,3])
      
     
      ## Prior probability of hypothesis
      i.pc.prior.probs = c(((1-h.z) / 2.0),h.z,((1-h.z) / 2.0)) 
      ## Log(P(ch(x) | Parents(ch(x))))
      ## Due to underflow, computations need to be done in log scale
      i.pc.log.ch.probs = c(0,0,0) ## holds the product of childs, components are different current Parent values
      
      for(ch in i.pc.ch){
        addToChain( ch)
        ch.net = relsByTrgUid[[as.character(ch)]];
        
        
        ## Value of the current child
        ch.val = X[ch.net[1,'trgind']] ## index of the current child
        ## Find the applicable edges
        
        ### DZ I added a number of ", drop=FALSE" to the code as these are now matrices and they need to
        ### stay matrices and not become vectors if only one row gets returned.
        
        ch.net.app = ch.net[which_3(X[ch.net[,"appind"]] == 1),, drop=FALSE]
        
        ## No applicable edges?
        if (dim(ch.net.app)[1] == 0){
          H.p = compProbH(0, 0,  p.a * (1 - abs(ch.val)) + p.c * abs(ch.val) , p.m, p.z, p.p)
          Z.p = pZ(ch.val, alpha, beta, p.m, p.z, p.p)
          ch.p = rep(sum(Z.p * H.p), 3) ## Pr(Ch = val | Par = .) does not have to add up to 1
          i.pc.log.ch.probs = i.pc.log.ch.probs + log(ch.p)
          next
        }
        
        ## unique applicable network
        pa.ind.lab = unique_1(ch.net.app[,"type_srcind", drop=FALSE])           
        ##pa.ind = as.double(abs(pa.ind.lab)) ## changed by KZ on March 5 2014
        pa.ind = as.double(abs(pa.ind.lab[[1]])) ## changed by KZ on March 5 2014
        ##pa.lab = as.double(sign(pa.ind.lab))  ## changed by KZ on March 5 2014
        pa.lab = as.double(sign(pa.ind.lab[[1]]))  ## changed by KZ on March 5 2014
        
        Val = X[pa.ind]
        ch.p = {} ## prob of child for diffeent values of Xi
        ## the child edge is not applicable, hence different values of X make no difference
        if(!(ents[i.pc,1] %in% unique_2(ch.net.app[,'srcuid']))){
          pred.vals = Val * pa.lab
          n.m = sum(pred.vals == -1)
          n.p = sum(pred.vals == 1)
          n.all = length(pred.vals)
          H.p = compProbH(n.m, n.p,  p.a * (1 - abs(ch.val)) + p.c * abs(ch.val) , p.m, p.z, p.p)
          Z.p = pZ(ch.val, alpha, beta, p.m, p.z, p.p)
          ch.p = rep(sum(Z.p * H.p), 3)
        }else{
          ## find the parent currently being updated
          cur.par.ind = which_5(pa.ind == i.pc)
          for(k in c(-1,0,1)){
            ## Compute Pr(H = h | Pa(H) = x) for h = -1, 0, and 1
            Val[cur.par.ind] = k ## Value of the node currently being updated
            pred.vals = Val * pa.lab
            n.m = sum(pred.vals == -1)
            n.p = sum(pred.vals == 1)
            n.all = length(pred.vals)
            H.p = compProbH(n.m, n.p,  p.a * (1 - abs(ch.val)) + p.c * abs(ch.val) , p.m, p.z, p.p)
            Z.p = pZ(ch.val, alpha, beta, p.m, p.z, p.p)
            ch.p = c(ch.p, sum(Z.p * H.p))
          }
        }
        i.pc.log.ch.probs = i.pc.log.ch.probs + log(ch.p)
      }

     
      i.pc.logX = log(i.pc.prior.probs) + i.pc.log.ch.probs
      i.pc.probs = c(1 / (1 + exp(i.pc.logX[2] - i.pc.logX[1]) + exp(i.pc.logX[3] - i.pc.logX[1])),
                     1 / (1 + exp(i.pc.logX[1] - i.pc.logX[2]) + exp(i.pc.logX[3] - i.pc.logX[2])),
                     1 / (1 + exp(i.pc.logX[1] - i.pc.logX[3]) + exp(i.pc.logX[2] - i.pc.logX[3])))
      
      ## computint p(X=.) = pa.probs * ch.probs / sum(pa.probs * ch.probs)
      i.pc.samp = sample(c(-1, 0, 1), size = 1, prob = i.pc.probs)
      X[i.pc] = i.pc.samp ##update i-th node
      
    }
   

    if(t > burn.in){
      x.pc.ind = X[pc.ind];
      marg.idx = which_14(x.pc.ind == -1);
      marg.probs.pc[1, marg.idx] =  marg.probs.pc[1, marg.idx] + 1
      marg.idx = which_14(x.pc.ind == 0);
      marg.probs.pc[2, marg.idx]  =  marg.probs.pc[2, marg.idx] + 1
      marg.idx = which_14(x.pc.ind == 1);
      marg.probs.pc[3, marg.idx]  =  marg.probs.pc[3, marg.idx] + 1
      
      samp.num.pc = samp.num.pc + 1
    }
    
    ##Updating the applicability nodes.
    i.a = a.ind[(t %% length(a.ind) + 1)]
    if(i.a %in% hyps.ind){
      X[i.a] = 0
    }else{
      ## Determining the network of current applicability
      i.a.pa.ind = relsByAppId[[as.character(ents[i.a, 1])]][, 'meshind', drop=FALSE]
      
      ##Val = X[i.a.pa.ind] ## changed by KZ MArch 5th 2014
      Val = X[i.a.pa.ind[[1]]] ## changed by KZ MArch 5th 2014
      ## Pr(A=a | Markov(A)) = Pr(A = a | Pa(A)) * Pr(Z | A = a) / sum{Pr(A = . | Pa(A)) * Pr(H | A = .)}
      z = 1 - (1 - w)^sum(Val != 0)
      i.a.pa.probs = c(1 - z, z)
      
      ## Child network (there is only one child)
      i.a.ch = unique_3(relsByAppId[[as.character(ents[i.a,1])]][, "trguid", drop=FALSE])
      
      addToChain( i.a.ch)
      
      i.a.ch.net = relsByTrgUid[[as.character(i.a.ch)]];
      
      ## Value of the child
      ##i.a.ch.val = X[i.a.ch.net[1,'trgind']] ## changed by KZ MArch 5th 2014
      i.a.ch.val = X[i.a.ch.net[1,'trgind'][[1]]] ## changed by KZ MArch 5th 2014
      ## Temporarly set the curret applicability node to 1 so it would be included in the network
      X[i.a] = 1
      ## Find the applicable edges
      ##ch.net.app.with.A = i.a.ch.net[which_11(X[i.a.ch.net[, "appind", drop = FALSE]] == 1),, drop = FALSE] ## changed by KZ MArch 5th 2014
      ch.net.app.with.A = i.a.ch.net[which_11(X[i.a.ch.net[, "appind", drop = FALSE][[1]]] == 1),, drop = FALSE] ## changed by KZ MArch 5th 2014
      ## Temporarly set the curret applicability node to 0 so it would be removed from the network
      X[i.a] = 0
      ## Find the applicable edges
      ##ch.net.app.without.A = i.a.ch.net[which_12(X[i.a.ch.net[, "appind", drop = FALSE]] == 1),, drop = FALSE] ## changed by KZ MArch 5th 2014
      ch.net.app.without.A = i.a.ch.net[which_12(X[i.a.ch.net[, "appind", drop = FALSE][[1]]] == 1),, drop = FALSE] ## changed by KZ MArch 5th 2014
      
      ## In each case compute Pr(Z | pa(Z))
      i.a.ch.p = {} ## prob of child for diffeent values of A
      
      ## When removing A are there any applicable edges? 
      if (dim(ch.net.app.without.A)[1] == 0){
        H.p = compProbH(0, 0,  p.c * (1 - abs(i.a.ch.val)) + p.c * abs(i.a.ch.val) , p.m, p.z, p.p)
        Z.p = pZ(i.a.ch.val, alpha, beta, p.m, p.z, p.p)
      }else{      
        pa.ind.lab = unique_4(ch.net.app.without.A[,"type_srcind", drop=FALSE])           
        ##pa.ind = as.double(abs(pa.ind.lab)) ## changed by KZ MArch 5th 2014
        pa.ind = as.double(abs(pa.ind.lab[[1]])) ## changed by KZ MArch 5th 2014
        ##pa.lab = as.double(sign(pa.ind.lab)) ## changed by KZ MArch 5th 2014
        pa.lab = as.double(sign(pa.ind.lab[[1]])) ## changed by KZ MArch 5th 2014
        
        addToChain( pa.ind, pa.lab)
        
        Val = X[pa.ind]
        
        pred.vals = Val * pa.lab
        n.m = sum(pred.vals == -1)
        n.p = sum(pred.vals == 1)
        H.p = compProbH(n.m, n.p,  p.c * (1 - abs(i.a.ch.val)) + p.c * abs(i.a.ch.val) , p.m, p.z, p.p)
        Z.p = pZ(i.a.ch.val, alpha, beta, p.m, p.z, p.p)
      }
      
      i.a.ch.p = c(i.a.ch.p, sum(Z.p * H.p))
      
      ## unique applicable network with A
      
      pa.ind.lab = unique_5(ch.net.app.with.A[,"type_srcind", drop=FALSE])           
      ##pa.ind = as.double(abs(pa.ind.lab[[1]])) ## changed by KZ MArch 5th 2014
      pa.ind = as.double(abs(pa.ind.lab[[1]])) ## changed by KZ MArch 5th 2014
      ##pa.lab = as.double(sign(pa.ind.lab)) ## changed by KZ MArch 5th 2014
      pa.lab = as.double(sign(pa.ind.lab[[1]])) ## changed by KZ MArch 5th 2014
      
      Val = X[pa.ind]
      pred.vals = Val * pa.lab
      n.m = sum(pred.vals == -1)
      n.p = sum(pred.vals == 1)
      
      H.p = compProbH(n.m, n.p,  p.c * (1 - abs(i.a.ch.val)) + p.c * abs(i.a.ch.val) , p.m, p.z, p.p)
      Z.p = pZ(i.a.ch.val, alpha, beta, p.m, p.z, p.p)
      
      i.a.ch.p = c(i.a.ch.p, sum(Z.p * H.p))
      
      ## computint p(A=.) = pa.probs * ch.probs / sum(pa.probs * ch.probs)
      i.a.probs = i.a.pa.probs * i.a.ch.p / sum(i.a.pa.probs * i.a.ch.p)
      i.a.samp = sample(c(0, 1), size = 1, prob = i.a.probs)
      X[i.a] = i.a.samp ##update A
      
    }
   

    if(t > burn.in){
      x.a.ind = X[a.ind]
      marg.idx = which_14(x.a.ind == 0);
      marg.probs.a[1, marg.idx]  =  marg.probs.a[1, marg.idx] + 1
      marg.idx = which_14(x.a.ind == 1);
      marg.probs.a[2, marg.idx]  =  marg.probs.a[2, marg.idx] + 1
      
      samp.num.a = samp.num.a + 1
    }
    
    
    ##Updating the context nodes.
    i.m = m.ind[(t %% length(m.ind) + 1)]
    if(i.m %in% hyps.ind){
      X[i.m] = 0
    }else{
      ## Determining the network of current applicability
      i.m.prior.probs = c(c.p, (1 - c.p))
      
      ## Determining the network of the childrens of the context (Markov Blanket)
      i.m.ch = unique(relsByMeSHId[[as.character(ents[i.m,1])]][,"appid"])
      
      ## Log(P(ch(x) | Parents(ch(x))))
      ## Due to underflow, computations need to be done in log scale
      i.m.log.ch.probs = c(0,0) ## holds the product of childs, components are different current Parent values
      
      for(ch in i.m.ch){
        ch.net = relsByAppId[[as.character(ch)]];
        
        ## Value of the current child
        ch.val = X[ch.net[1,'appind']]
        pa.ind = ch.net[,'meshind']
        Val = X[pa.ind]
        cur.par.ind = which_16(pa.ind == i.m)
        ch.p = {}
        for(k in c(0,1)){
          Val[cur.par.ind] = k ## Value of the node currently being updated
          z = 1 - (1 - w)^sum(Val != 0)
          if(ch.val == 1){
            ch.p = c(ch.p, z)
          }else{
            ch.p = c(ch.p, (1 - z))
          }
        }
        i.m.log.ch.probs = i.m.log.ch.probs + log(ch.p)
      }
      
      i.m.logX = log(i.m.prior.probs) + i.m.log.ch.probs
      i.m.probs = c(1 / (1 + exp(i.m.logX[2] - i.m.logX[1])),
                    1 / (1 + exp(i.m.logX[1] - i.m.logX[2])))
      
      ## computint p(X=.) = pa.probs * ch.probs / sum(pa.probs * ch.probs)
      i.m.samp = sample(c(0, 1), size = 1, prob = i.m.probs)
      X[i.m] = i.m.samp ##update i-th node
      
    }
   
    
    if(t > burn.in){
      x.m.ind = X[m.ind]
      marg.idx = which_14(x.m.ind == 0);
      marg.probs.m[1, marg.idx]  =  marg.probs.m[1, marg.idx] + 1
      marg.idx = which_14(x.m.ind == 1);
      marg.probs.m[2, marg.idx]  =  marg.probs.m[2, marg.idx] + 1
      
      samp.num.m = samp.num.m + 1
    }
    
  }
  
  marg.probs.pc = marg.probs.pc / samp.num.pc
  marg.probs.m  = marg.probs.m / samp.num.m
  marg.probs.a  = marg.probs.a / samp.num.a 
  
  results.pc = cbind(ents[pc.ind,1], t(marg.probs.pc))  
  results.m  = rbind(c(ents[fix.nonzero.MeSH.ind,1], 0, 1), cbind(ents[m.ind,1], t(marg.probs.m)))
  results.a  = cbind(ents[a.ind,1], t(marg.probs.a))
  
  results = list(results.pc =  results.pc, results.m =  results.m, results.a = results.a)
  
  results
}


### The following functions perform some statistics and post-processing:

########## Stat, pos-processing Start
## This function returns the network directly connected to a given node.
##
## Arguments:
##
## node     -    id of the node
##
## ents     -    entries of the network
##
## rels     -    relationships
##
## levels   -    boolean (include other level for proteins)
nodeNet = function(node, ents, rels, levels = F){
  node.ents = {}
  node.rels = {}
  
  ## Find the targets of the node
  targ = rels[which(rels[,2] == node),]
  node.ents = rbind(node.ents, ents[which(ents[,1] %in% targ[,3]), ])
  node.rels = rbind(node.rels, targ)
  
  ## Find the sources of the node
  src = rels[which(rels[,3] == node),]
  node.ents = rbind(node.ents, ents[which(ents[,1] %in% src[,2]), ])
  node.rels = rbind(node.rels, src)
  
  ## add the node
  node.ents = rbind(node.ents, ents[which(ents[,1] == node), ])
  
  ## consider other level
  if(levels){
    if(ents[which(ents[,1] == node), 4] == 'protein'){
      if(substr(node,nchar(node), nchar(node)) == '1'){
        node2 = paste(substr(node,1, (nchar(node) - 1)), '0', sep = '')
      }else{
        node2 = paste(substr(node,1, (nchar(node) - 1)), '1', sep = '')
      }
    }
    ## Find the targets of the node
    targ = rels[which(rels[,2] == node2),]
    node.ents = rbind(node.ents, ents[which(ents[,1] %in% targ[,3]), ])
    node.rels = rbind(node.rels, targ)
    
    ## Find the sources of the node
    src = rels[which(rels[,3] == node2),]
    node.ents = rbind(node.ents, ents[which(ents[,1] %in% src[,2]), ])
    node.rels = rbind(node.rels, src)
    
    ## add the node
    node.ents = rbind(node.ents, ents[which(ents[,1] == node2), ])
    
  }
  
  L = list(ents = unique(node.ents), rels = unique(node.rels))
  
}

## This function returns the network directly connected to a list of given node.
##
## Arguments:
##
## nodes    -    list of id of the nodes
##
## ents     -    entries of the network
##
## rels     -    relationships
##
## levels   -    boolean (include other level for proteins)
nodeNetList = function(nodes, ents, rels, levels = F){
  node.ents = {}
  node.rels = {}
  
  for(node in nodes){
    ## Find the targets of the node
    targ = rels[which(rels[,2] == node),]
    node.ents = rbind(node.ents, ents[which(ents[,1] %in% targ[,3]), ])
    node.rels = rbind(node.rels, targ)
    
    ## Find the sources of the node
    src = rels[which(rels[,3] == node),]
    node.ents = rbind(node.ents, ents[which(ents[,1] %in% src[,2]), ])
    node.rels = rbind(node.rels, src)
    
    ## add the node
    node.ents = rbind(node.ents, ents[which(ents[,1] == node), ])
    
    ## consider other level
    if(levels){
      if(ents[which(ents[,1] == node), 4] == 'protein'){
        if(substr(node,nchar(node), nchar(node)) == '1'){
          node2 = paste(substr(node,1, (nchar(node) - 1)), '0', sep = '')
        }else{
          node2 = paste(substr(node,1, (nchar(node) - 1)), '1', sep = '')
        }
      }
      ## Find the targets of the node
      targ = rels[which(rels[,2] == node2),]
      node.ents = rbind(node.ents, ents[which(ents[,1] %in% targ[,3]), ])
      node.rels = rbind(node.rels, targ)
      
      ## Find the sources of the node
      src = rels[which(rels[,3] == node2),]
      node.ents = rbind(node.ents, ents[which(ents[,1] %in% src[,2]), ])
      node.rels = rbind(node.rels, src)
      
      ## add the node
      node.ents = rbind(node.ents, ents[which(ents[,1] == node2), ])
      
    }
    
  }
  
  L = list(ents = unique(node.ents), rels = unique(node.rels))
  
}

## This function prints the stat of a node
##
## Arguments:
##
## reg      -    direction of regulation
##
## node     -    id of the node
##
## ents     -    entries of the network
##
## rels     -    relationships
##
## evidence -    gene values
node.stat = function(reg, node, ents, rels, evidence){
  
  targ = rels[which(rels[,2] == node),3]
  dirs = ifelse(rels[which(rels[,2] == node), 4] == 'increase', 1, -1)
  
  targ.vals = rep(0, length(targ))
  for(t in targ[which(targ %in% evidence[,1])]){
    targ.vals[which(targ == t)] = as.numeric(as.character(evidence[which(evidence[,1] == t),2]))
  }
  
  prediction = reg * as.numeric(dirs)
  r = prediction - targ.vals
  c = length(which(r == 0))
  z = length(which(targ.vals == 0))
  i = length(targ.vals) - (z + c)
  
  
  L = list(c = c, i = i, z = z)
  L
}

## This function prints the stat of a list of nodes
##
## Arguments:
##
## nodes      -    list of nodes
##
## node.vals  -    node values (direction of regulation)
##
## ents       -    entries of the network
##
## rels       -    relationships
##
## evidence   -    gene values

node.stat.list = function(nodes, node.vals, ents, rels, evidence){
  
  L = nodeNetList(nodes, ents, rels, levels = F)
  on.trans = merge(L$ents[which(L$ents[,4] == 'mRNA'),1, drop = F], evidence, by.x = 1, by.y = 1)
  no.pred.num = dim(evidence)[1] - dim(on.trans)[1]
  
  correct.num = 0 #explained by at least one node
  incorrect.num = 0 #wrong prediction by all of the nodes
  for(t in 1:dim(on.trans)[1]){
    tr = on.trans[t,1]
    tr.val = on.trans[t,2]
    tr.par = L$rels[which(L$rels[,3] == tr),2]
    tr.dir = ifelse(L$rels[which(L$rels[,3] == tr),4]  == 'increase', 1, -1)
    par.vals = node.vals[which(nodes %in% tr.par)]
    pred.vals = par.vals * tr.dir
    if(any(pred.vals == tr.val)){
      correct.num = correct.num + 1
    }else{
      incorrect.num = incorrect.num + 1
    }
  }
  
  L = list(correct = correct.num, incorrect = incorrect.num, zero = no.pred.num)
}
########## Stat, pos-processing End

########## Functions for generating evidence data and simulating data:

###############################################################
## Given Micro Array data, this function extracts the significant genes and maps the entrez gene ids to
## uid of the network. 
## Given a (one-level) network, this function simulates the values using the
## given parameters.
##
## Arguments:
##
## ents    -    one-level network ents (generated from a given KB using genOneLevel)
##
## rels    -    one-level network rels (generated from a given KB using genOneLevel)
##
## nodes   -    list of uid of selected hypothesis
##
## values  -    values of selected hypothesis
##
## p.c     -   probability of the edge being wrong for nonzero edges
##
## p.m     -   prior probability of the ture state of the genes being -1
##
## p.z     -   prior probability of the ture state of the genes being 0
##
## p.p     -   prior probability of the ture state of the genes being 1
##
## h.z     -   prior probability of the hypothesis being off
##
## alpha   -   false positive rate
##
## beta    -   false negative rate
simValues = function(ents, rels, nodes, values, p.c, p.m, p.z, p.p, alpha, beta){
  ## Simulate values.
  mrna.ind = which(ents[,'type'] == 'mRNA')
  num.mrna = length(mrna.ind)
  sim.data = cbind(ents[mrna.ind, 'uid'], rep(0,num.mrna))
  colnames(sim.data) = c('id', 'val')
  
  ## simulate the value of each mRNA
  for(i in 1:dim(sim.data)[1]){
    i.ind = which(rels[,'trguid'] == sim.data[i,1])
    direction = ifelse(rels[i.ind,'type'] == 'increase', 1, -1)
    
    src.ind = which(nodes %in% rels[i.ind,'srcuid'])
    pred.val = sum(values[src.ind] * direction)
    
    if(pred.val >= 1){
      ph = c(p.c * p.m, p.c * p.z, 1 - p.c + p.c * p.p, 0)
    }else if(pred.val <= -1){
      ph = c(1 - p.c + p.c * p.m, p.c * p.z, p.c * p.p, 0)
    }else{ ## Conflict case with equal 1s and -1s: Flip a coin.
      z = runif(1)
      if(z > 0.5){
        ph = c(p.c * p.m, p.c * p.z, 1 - p.c + p.c * p.p, 0)
      }else{
        ph = c(1 - p.c + p.c * p.m, p.c * p.z, p.c * p.p, 0)
      }
    }
    
    pzm1 = c(1-2*beta,alpha,beta,(1/3))
    pz0  = c(beta,1-2*alpha,beta,(1/3))
    pzp1 = c(beta,alpha,1-2*beta,(1/3))
    
    sim.prob = c(sum(pzm1*ph), sum(pz0*ph), sum(pzp1*ph))
    
    
    sim.data[i,2] = sample(c(-1,0,1), 1, prob=sim.prob)
  }
  
  ## change uid to entrez id
  sim.data[,1] = ents$id[match(sim.data[,1], ents$uid)]
  
  sim.data
}
##############################

######
