#! /usr/bin/Rscript
##
## This file is a run script for BayesCRE Gibbs Sampler.
##
## Computational Sceinces CoE
## Pfizer Worldwide Research and Development
##
## Date: April 16, 2013
##
## Author: Kourosh Zarringhalam
###############################################################################

## Required libraires
library('plyr')
library('parallel')
library('optparse')
library('fdrtool')

suppressPackageStartupMessages(library("optparse"))

## Parsing the commandline input
option_list = list(
  make_option(c("-e", "--entries"),
              help="path to the KB entries File"),
  make_option(c("-r", "--relations"),
              help="path to KB relations File"),
  make_option(c("-d", "--dat"),
              help="path to evidence File (differentially expressed genes with ENTREZ IDs)"),
  make_option(c("-m", "--mesh"),
              help="path to the MeSH file (MeSH.txt)"),
  make_option(c("-a", "--alpha"), default = 0.05,
              help="FP rate: Default Value: 0.05"),
  make_option(c("-b", "--beta"), default = 0.1,
              help="FN rate: Default Value: 0.1"),
  make_option(c("-j", "--mincorrect"), default = 1,
              help="Exclude hypothesis that do not have this much correct predictions. Default value 1"),
  make_option(c("-c", "--pc"), default = 0.1,
              help="probability of edge being inapplicable for non-zero network: Default Value: 0.1"),
  make_option(c("-g", "--pa"), default = 0.9,
              help="probability of edge being inapplicable for zero network: Default Value: 0.9"),
  make_option(c("-z", "--hz"), default = 0.9,
              help="prior probability of zero for hypothesis: Default Value: 0.9"),
  make_option(c("-p", "--cp"), default = 0.7,
              help="prior probability of zero for context (MeSH): Default Value: 0.7"),
  make_option(c("-w", "--weight"), default = 0.3,
              help="weight for the applicability sigmoid: Default Value: 0.3"),
  make_option(c("-n", "--samples"), default = 2000000,
              help="total number of samples: Default Value: 1000000"),
  make_option(c("-u", "--burnin"), default = 1000000,
              help="Number of burnin samples for MCMC: Default Value 200000"),
  make_option(c("-l", "--iter"), default = 1,
              help="number of iterations: Default Value: 1"),
  make_option(c("-s", "--chains"), default = 1,
              help="number of chains: Default Value: 1"),
  make_option(c("-q", "--cutoffq"), default = 1e-4,
              help="cutoff for fdr q value: Default Value: 1e-4"),
  make_option(c("-f", "--cutoffh"), default = 0.4,
              help="cutoff for top hypothesis selection: Default Value: 0.4"),
  make_option(c("-x", "--mn"), default = 600,
              help="max number of edges that MeSH can be annotated to: Default Value: 600"),
  make_option(c("-k", "--meshprob"),
              help="path to the output joint prob File"),
  make_option(c("-y", "--marg"),
              help="path to the output marginal prob File"),
  make_option(c("-t", "--tophyp"),
              help="path to the top hypothesis output"),
  make_option(c("-o", "--outmeshapp"),
              help="path to the output Applicability MeSH File"))

## Parsing the input argument
opt = parse_args(OptionParser(option_list=option_list))

entFile         = opt$entries
relFile         = opt$relations
dataFile        = opt$dat
meshFile        = opt$mesh
alpha           = as.numeric(opt$alpha)
beta            = as.numeric(opt$beta)
minscore        = as.numeric(opt$minscore)
p.c             = as.numeric(opt$pc)
p.a             = as.numeric(opt$pa)
h.z             = as.numeric(opt$hz)
c.p             = as.numeric(opt$cp)
w               = as.numeric(opt$weight)
N               = as.numeric(opt$samples)
burn.in         = as.numeric(opt$burnin)
iter.num        = as.numeric(opt$iter)
chain.num       = as.numeric(opt$chains)
cutoff.q        = as.numeric(opt$cutoffq)
cutoff.h        = as.numeric(opt$cutoffh)
Mn              = as.numeric(opt$mn)
meshProbFile    = opt$meshprob
margFile        = opt$marg
hypFile         = opt$tophyp
outMeSHAppFile  = opt$outmeshapp


setwd('./') ## Change this as appropriate
source('./algorithms.R')

## Generaing the contex network
print('Generaing the contex network')
L = genContextGraph(entFile, relFile, dataFile, meshFile, cutoff.q , Mn, minscore)

ents     = L$ents
rels     = L$rels
evidence = L$evidence


print('dimensions of the context network:')
print(dim(ents))
print(dim(rels))


## Read in and process all ents and rels for annotations
LL = processKB(entFile, relFile)
ents.all = LL$ents.all
rels.all = LL$rels.all
id.map  = LL$id.map

mRNA.ents = ents[which(ents[,'type'] == 'mRNA'), 1]

total.mRNA = length(mRNA.ents)
up.mRNA    = length(which(evidence[,2] == 1))
down.mRNA  = length(which(evidence[,2] == -1))
zero.mRNA  = total.mRNA - (up.mRNA + down.mRNA)

## Estimate background probabilities
p.z = zero.mRNA / total.mRNA
p.m = down.mRNA / total.mRNA
p.p = up.mRNA / total.mRNA

print('Estimate of the background probabilities:')
print(c(p.m, p.z, p.p))

## Unique relations
rels.unique = rels[!duplicated(rels[,c(2,3)]),]
#rels.unique = rels[match(unique(rels[,1]), rels[,1]), ]

LOG_CHAIN = FALSE
hyps.ind = {}
up.down.pred = {}
for(iter in 1:iter.num){
  print(paste('Running the Gibbs Sampler. Iteration', iter))
  
  tt <- system.time(R <- mclapply(1:chain.num, gibbsSampler, sim.num = N, 
                                  burn.in = burn.in, iter.num = iter, p.c = p.c, 
                                  p.a = p.a, p.m = p.m, p.z = p.z, p.p = p.p, h.z = h.z, c.p = c.p, 
                                  alpha = alpha, beta = beta, w = w, hyps.ind = hyps.ind))
  
  print(paste('Gibbs sampler run time:', tt[1]))
  
  print('Generating marginal probabilities')
  marg.probs.pc = matrix(0, nrow = dim(R[1][[1]][[1]])[1], ncol = 3)
  for(ch in 1:chain.num){
    marg.probs.pc = marg.probs.pc + apply(R[ch][[1]][[1]][,2:4],2,as.numeric)
  }
  
  marg.probs.pc = marg.probs.pc / chain.num
  marg.probs.pc = data.frame(cbind(R[1][[1]][[1]][,1], marg.probs.pc), stringsAsFactors = F)
  colnames(marg.probs.pc) = c('uid', 'prob_down', 'prob_zero', 'prob_up')
  
  ## Applicable edges
  print(paste('Determining applicable edges for iteration', iter))
  marg.probs.a = matrix(0, nrow = dim(R[1][[1]][[3]])[1], ncol = 2)
  for(ch in 1:chain.num){
    marg.probs.a = marg.probs.a + as.numeric(R[ch][[1]][[3]][,2:3])
  }
  
  marg.probs.a = marg.probs.a / chain.num
  marg.probs.a = data.frame(cbind(R[1][[1]][[3]][,1], marg.probs.a), stringsAsFactors = F)
  colnames(marg.probs.a) = c('uid', 'prob_zero', 'prob_one')
  App.edge = merge(rels, marg.probs.a, by.x = 7, by.y = 1)[,c(2,3,4,5,6,1,7,12,13)]
 
  ## write the top MeSH file
  print(paste('Determining MeSH terms for iteration', iter))
  marg.probs.m = matrix(0, nrow = dim(R[1][[1]][[2]])[1], ncol = 2)
  for(ch in 1:chain.num){
    marg.probs.m = marg.probs.m + as.numeric(R[ch][[1]][[2]][,2:3])
  }
  
  marg.probs.m = marg.probs.m / chain.num
  marg.probs.m = data.frame(cbind(R[1][[1]][[2]][,1], marg.probs.m), stringsAsFactors = F)
  colnames(marg.probs.m) = c('uid', 'prob_zero', 'prob_one')
  sig.MeSH = merge(ents, marg.probs.m, by.x = 1, by.y = 1)
  
  App.MeSH = merge(App.edge, sig.MeSH, by.x = 7, by.y = 1)[,c(2,3,4,5,6,7,9,10, 14)]

  ## Add gene values
  colnames(App.MeSH) = c('uid','srcuid','trguid', 'type','pmids','appid', 'edge_prob','MeSH','MeSH_prob')
  App.MeSH = merge(App.MeSH, evidence, by.x = 3, by.y = 1, all.x=T)[,c(2,3,1,10,4,6,7,8,9,5)]
  App.MeSH$val[is.na(App.MeSH$val)] = 0

  ## Add srcname and trgname
  tmp.frame = merge(App.MeSH, ents.all, by.x = 2, by.y = 1)
  App.MeSH = merge(tmp.frame, ents.all, by.x = 3, by.y = 1)[,c(3,2,1,11,13,14,16,4,5,6,7,8,9,10)]
  colnames(App.MeSH) = c('uid','srcuid','trguid','srcname','srctype','trgname','trgtype','trgval', 
                         'type','appid', 'edge_prob','MeSH','MeSH_prob', 'pmids')
  ## Add nls
  App.MeSH = merge(App.MeSH, rels.all, by.x = 1, by.y = 1)[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,19)]
  colnames(App.MeSH) = c('uid','srcuid','trguid','srcname','srctype','trgname','trgtype','trgval', 
                         'type','appid', 'edge_prob','MeSH','MeSH_prob', 'pmids', 'nls')
  
  ## sort by src id
  ind = sort(as.character(App.MeSH[,2]), index.return = T)$ix
  App.MeSH.sort = App.MeSH[ind,]
  rownames(App.MeSH.sort) = 1:dim(App.MeSH.sort)[1]
  
  iter.MeSHAppFile = paste(dirname(outMeSHAppFile), '/iteration_', iter, '_', basename(outMeSHAppFile), sep = '')
  write.table(App.MeSH.sort, iter.MeSHAppFile, row.names = F, sep = '\t', quote = F)
  
  select.down.ind1 = which(as.numeric(marg.probs.pc[,2]) >= cutoff.h )
  select.down.ind2 = which(as.numeric(marg.probs.pc[select.down.ind1,2]) - 
                             as.numeric(marg.probs.pc[select.down.ind1,3]) > 0 )
  select.up.ind1 = which(as.numeric(marg.probs.pc[,4]) >= cutoff.h )
  select.up.ind2 = which(as.numeric(marg.probs.pc[select.up.ind1,4]) - 
                             as.numeric(marg.probs.pc[select.up.ind1,3]) > 0 )
  select.sig.ind = unique(c(select.up.ind1[select.up.ind2], select.down.ind1[select.down.ind2]))
  
  sig.hyp = marg.probs.pc[select.sig.ind, , drop = F]
  ITE = rep(paste('iter', iter, sep = ''), dim(sig.hyp)[1])
  
  colnames(sig.hyp)[1] = 'srcuid'
  tmp.Tab = merge(sig.hyp, App.MeSH.sort, by.x = 1, by.y = 2)
  colnames(sig.hyp)[1] = 'uid'
  hyp.MeSH.list = data.frame(uid= character(0), MeSH=character(0))
  for(hyp in unique(tmp.Tab$srcuid)){
    hyp.MeSH = unique(tmp.Tab[tmp.Tab$srcuid == hyp,][,c('MeSH', 'MeSH_prob')])
    sig.hype.MeSH = c(hyp.MeSH$MeSH[hyp.MeSH$MeSH_prob >= 0.5])
    hyp.MeSH.list = rbind(hyp.MeSH.list, cbind(hyp, paste(sig.hype.MeSH, ' ', sep = ',',collapse='')))
  }
  
  colnames(hyp.MeSH.list) = c('uid', 'MeSH')  
  up.down.pred = rbind(up.down.pred, merge(cbind(sig.hyp, ITE), hyp.MeSH.list, by.x = 1, by.y = 1))
  
  m = merge(ents, up.down.pred, by.x = 1, by.y = 1)
  hyps.id = m[,1]
  hyps.ind = which(ents[,1] %in% hyps.id)
  
  ##### KZ: June 9, 2015: These will not be upadted in second iteration
  rm.apps.id = unique(rels$appid[which(rels$srcuid %in% hyps.id)])
  rm.apps.ind = which(ents[,1] %in% rm.apps.id)
  ## May have to remove these if no other app is connected to them
  rm.hyps.MeSH = unique(rels$meshid[which(rels$srcuid %in% hyps.id)]) 
  rm.MeSH = rm.hyps.MeSH[which(!(rm.hyps.MeSH %in% rels$meshid[-which(rels$srcuid %in% hyps.id)]))]
  if(length(rm.MeSH) > 0){
    rm.MeSH.ind = which(ents$uid %in% rm.MeSH)
    hyps.ind = c(hyps.ind, rm.apps.ind, rm.MeSH.ind)
  }else{
    hyps.ind = c(hyps.ind, rm.apps.ind)
  }
  ##### KZ: June 9, 2015: These will not be upadted in second iteration
  
  print(paste('Generating output of marginals for hypothesis for iteration', iter))
  m = merge(ents, marg.probs.pc, by.x = 1, by.y = 1)
  Tab = matrix(0, nrow = (dim(m)[1]*3), ncol = 10)
  row.count = 1
 
  for(node in 1:dim(m)[1]){
    L = nodeNet(m[node,1], ents, rels.unique, levels = F)
    s.n = node.stat(-1, m[node,1], L$ents, L$rels, evidence)
    E = c(uid = m[node,1], regulation = -1, name = m[node,2], id = m[node,3], type = m[node,4], 
          prob = m[node,5], correct = s.n$c, 
          incorrect = s.n$i, zero = s.n$z, score = (s.n$c - s.n$i))
    Tab[((node - 1) * 3 + 1), ] = t(as.matrix(E))
    
    s.n = list(c = 0,i = 0,z = 0)
    E = c(uid = m[node,1], regulation = 0, name = m[node,2], id = m[node,3], type = m[node,4], 
          prob = m[node,6],  correct = s.n$c, 
          incorrect = s.n$i, zero = s.n$z, score = (s.n$c - s.n$i))
    Tab[((node - 1) * 3 + 2), ] = t(as.matrix(E))
    
    
    s.n = node.stat(1, m[node,1], L$ents, L$rels, evidence)
    E = c(uid = m[node,1], regulation = 1, name = m[node,2], id = m[node,3], type = m[node,4], 
          prob = m[node,7],  correct = s.n$c, 
          incorrect = s.n$i, zero = s.n$z, score = (s.n$c - s.n$i))
    Tab[((node - 1) * 3 + 3), ] = t(as.matrix(E))
  }
  
  colnames(Tab) = c('uid', 'regulation', 'name', 'id', 'type', 'prob',  'correct', 'incorrect', 'zero', 'score')
  rownames(Tab) = 1:dim(Tab)[1]
  
  iter.margFile = paste(dirname(margFile), '/iteration_', iter, '_', basename(margFile), sep = '')
  write.table(Tab, iter.margFile, row.names = F, sep = '\t', quote = F)
  
  ## MeSH probabilities
  MeSH.probs = unique(App.MeSH.sort[,c('MeSH', 'MeSH_prob')])
  rownames(MeSH.probs) = 1:dim(MeSH.probs)[1]
  iter.meshProbFile = paste(dirname(meshProbFile), '/iteration_', iter, '_', basename(meshProbFile), sep = '')
  write.table(MeSH.probs, iter.meshProbFile, row.names = F, sep = '\t', quote = F)
  
  ## Put the original ids back and generate a uniqe relation network
  plotNW = App.MeSH.sort[,c("uid", "srcuid", 'trguid', "srcname", "srctype", "trgname", "trgtype", "trgval", "type", 
                   "pmids", "nls")]
  plotNW = plotNW[!duplicated(plotNW), ]
  plotNW$srcuid = id.map$uid.orig[match(plotNW$srcuid, id.map$uid.new)]
  plotNW$trguid = id.map$uid.orig[match(plotNW$trguid, id.map$uid.new)]
  
  iter.plotNW = paste(dirname(hypFile), '/iteration_', iter, '_', 
                      paste(gsub('.txt', '', gsub('hyp', '', basename(hypFile))), 'network.txt', sep = '_')
                      , sep = '')
  write.table(plotNW, iter.plotNW, row.names = F, sep = '\t', quote = F)
    
}


print('Writing predicted up/down regulated hypothesis')

m = merge(ents, up.down.pred, by.x = 1, by.y = 1)

up.down.stat = {}
up.down.dir = {}
if(dim(m)[1] >= 1){
  for(node in 1:dim(m)[1]){
    ## entire network
    L = nodeNet(m[node,1], ents, rels.unique, levels = F)
    if(as.numeric(as.character(m[node,5])) >= as.numeric(as.character(m[node,7]))){
      up.down.dir = c(up.down.dir, -1)
      p.d = -1
    }else{
      p.d = 1
      up.down.dir = c(up.down.dir, 1)
    }
    s.n = node.stat(p.d, m[node,1], L$ents, L$rels, evidence)
    hyp.stat = as.data.frame(s.n)
    colnames(hyp.stat) = c('c', 'i', 'z')
    up.down.stat = rbind(up.down.stat, hyp.stat)
  }
  
  score = up.down.stat$c - up.down.stat$i
  E = cbind(up.down.dir,up.down.stat[,c(1,2,3)], score)
  colnames(E) = c('regulation', 'correct', 'incorrect', 'zero', 'score')
  
  Tab = cbind(m[,1:(ncol(m)-1)], E,m[,ncol(m)])
  colnames(Tab)[ncol(Tab)] = 'MeSH'
  rownames(Tab) = 1:dim(Tab)[1]
  
  nodes = as.character(Tab$uid)
  node.vals = Tab$regulation
  S = node.stat.list(nodes, node.vals, ents, rels.unique, evidence)
  
  print('number of genes correctly predicted by at least one hyp in the list')
  print(S$correct)
  
  print('number of genes incorrectly predicted by all hyps in the list')
  print(S$incorrect)
  
  print('number of genes not in the network of selected hyps')
  print(S$zero)
  
  print('percentage of the observation explained by the selected hypes')
  pr = S$correct / (S$correct + S$incorrect + S$zero)
  print(pr)
  write.table(Tab, hypFile, row.names = F, sep = '\t', quote = F)
  
}else{
  print('No Significant Hypothesis selected!')
}

