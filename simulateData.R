##### Simulating Data
source('./algorithms.R')


## Read in the BEL KB and generate a one-level network:
entFile = './data/KB/BEL_LargeCorpus.ents'
relFile = './data/KB/BEL_LargeCorpus.rels'
L = genOneLevel(entFile, relFile)

ents = L$ents
rels = L$rels

## Need to pick a few hypothesis. Selecting hypothesis with reasonable number
## of downstream genes (N)

N = 50
pc.ents = ents[which(ents[,'type'] %in% c('Protein', 'Compound')),]

pc.list = {}
for(item in pc.ents[,1]){
  L = nodeNet(node=item,ents,rels)
  if(dim(L$rels)[1] > N){
    pc.list = c(pc.list, item)
  }
}

## Candidate hypothesis.
print(ents[which(ents[,1] %in% pc.list), ])

## Select a few nodes, for example:
## 1) HRAS
## 2) FOXO3
## 3) KLF1
## nodes contains the uid of the selected nodes
nodes = c('1117', '5459', '8376')
values = c(1, -1, 1)

## Generate the network of the selected hypothesis and downstream genes.
L = nodeNetList(nodes,ents,rels)

## Parameter Values
alpha     = 0.05 ## False positive rate
beta      = 0.1  ## False negative rate
p.c       = 0.1  ## Probability of edge being inapplicable
p.m       = 1/3  ## backgorund probability of -1
p.z       = 1/3  ## backgorund probability of  0
p.p       = 1/3  ## backgorund probability of +1


## Simulate the data:
sim.data = simValues(L$ents, L$rels, nodes, values, p.c, p.m, p.z, p.p, alpha, beta)

## Save the evidence.
write.table(sim.data, './data/Simulations/simulated_evidence.txt', row.names = F, quote = F, sep = '\t')

