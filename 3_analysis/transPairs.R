# Original code by Xavier Didelot: https://github.com/xavierdidelot/TransPairs

library(graph)
library(RBGL)
library(igraph)
library(ape)
library(BactDating)

source('likelihoodSEIR.R')

load('~/data/timeTrees_outbreaks.RData')


for (id in names(res_full)){
  
  tree = res_full[[id]]
  
  n = length(tree$tip.label)
  res = matrix(n,n, data = 0)

  if (rep==1){allres=res}
  for (i in 1:n) {
    for (j in 1:n) {
      if (i==j) next
      # Consider transmission from i to j
      tree2=drop.tip(tree, tree$tip.label[!tree$tip.label %in% tree$tip.label[c(i,j)]])
      
      if (length(tree2$tip.label) != 2) {stop('Error: duplicate tree labels')}
      
      d=adephylo::distRoot(tree2)
      
      ti = d[names(d) == tree$tip.label[i]]
      tj = d[names(d) == tree$tip.label[j]]
      tij=0
      if (TRUE) {
        res[i,j] = likelihoodSEIR(tij,ti,tj) # Calculate pairwise likelihood using the likelihoodSEIR function
      } else {
        # This is a previous attempt at calculating the pairwise likelihood by averaging the integral over several points
        npoints=10;
        for (si in 1:npoints){
          s=tij+(min(ti,tj)-tij)*si/(npoints+1);
          res(i,j)=res(i,j)+exppdf(abs(tij-s),1/neg)*exppdf(abs(tj-s),1/gamma);
        }
        res[i,j]=res[i,j]*(1-expcdf(abs(tij-ti),1/gamma));
      }
    }
  }
  allres=allres+res/nrep
  colnames(allres) = tree$tip.label
  rownames(allres) = tree$tip.label
  
}


allres = log10(allres)



test1 = allres
test1[test1==-Inf]<-0;test1[test1==0]<-NaN;test2=exp(test1);test2[test2=='NaN']<-0

colnames(test2) = sapply(colnames(test2), function(x) strsplit(x, '[.]')[[1]][1])
rownames(test2) = sapply(rownames(test2), function(x) strsplit(x, '[.]')[[1]][1])


source('transPairs_network.R')
source('transPairs_heatmap.R')

