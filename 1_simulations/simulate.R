library(ape)
library(phangorn)

rm(list=ls())

for (seed in 1:100){
# seed=1
set.seed(seed)


dir_prefix = '~/all_sims/'
dir.create(dir_prefix, showWarnings = FALSE)
setwd(dir_prefix)


l=30000#number of sites
t=rtree(100)
t$edge.length=t$edge.length/sum(t$edge.length)*0.05 # normalise tree branch lengths so that about x% of sites will mutate
#plot(t);axisPhylo(1)
write.tree(t,sprintf('tree%d.nwk',seed))

levels=c( 1 , 2 , 3 , 4 , 5  , 6  , 7  , 8  , 9   ,10  ,11,  12 , 13 , 14 , 15 , 16 )#State indexes
levfig=c('A','C','G','T','Ac','Ag','At','Cg','Ct','Gt','Ca','Ga','Ta','Gc','Tc','Tg')#State names used in figure
levels=c('A','C','G','T','M', 'R', 'W', 'S', 'Y', 'K', 'N', 'D', 'Q', 'E', 'H', 'I' )#State names taken from AA alphabet

p1=1    #Rate at which minor variant evolves [keep equal to 1 as reference rate]
p2=100   #Rate at which minor variants are lost  [relative to p1]
p3=200  #Rate at which minor/major variants are switched [relative to p1]

Q=matrix(0,length(levels),length(levels)) 
colnames(Q)<-levfig
rownames(Q)<-levfig
Q[1,c(5,6,7)]=p1            #A evolves minor variant
Q[2,c(8,9,11)]=p1           #C evolves minor variant
Q[3,c(10,12,14)]=p1         #G evolves minor variant
Q[4,c(13,15,16)]=p1         #T evolves minor variant
Q[c(5,6,7),1]=p2            #loss of minor variant
Q[c(8,9,11),2]=p2           #loss of minor variant
Q[c(10,12,14),3]=p2         #loss of minor variant
Q[c(13,15,16),4]=p2         #loss of minor variant
for (i in 5:10) {           #minor/major switch
  Q[i,i+6]=p3
  Q[i+6,i]=p3
}
diag(Q)=0;diag(Q)=-rowSums(Q) # Replace diagonal with minus sum of rows

#Calculate steady state
A=Q
A[,1]=1
bf=solve(t(A),c(1,rep(0,nrow(A)-1)))
print('Steady state');print(round(bf,3))
print(paste('Expect non-hets to be more likely than hets by factor p2/p1=',p2/p1)) # Note this fully determines p2 given p1=1 and empirical bf

#Check time-reversibility
for (i in 1:nrow(Q)) for (j in 1:nrow(Q)) A[i,j]=bf[i]*Q[i,j]-bf[j]*Q[j,i]
print(paste('Checking time-reversibility:',max(abs(A))))

pamlorder=c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')
Q2=matrix(0,length(pamlorder),length(pamlorder))
bf2=rep(0,length(pamlorder))
for (i in 1:length(levels)) for (j in 1:length(levels)) {
  wi=which(pamlorder==levels[i])
  wj=which(pamlorder==levels[j])
  Q2[wi,wj]=Q[i,j]
  bf2[wi]=bf[i]
}
s=simSeq(t,l,type='AA',Q=Q2[lower.tri(Q2)],bf=bf2)

# Write file containing full information about major/minor variants
write.phyDat(s,format='fasta',file = sprintf('simuFull%d.fas',seed)) 

# Write file without information on major/minor variants
s2=as.character(s)
cons=c('A','C','G','T','M', 'R', 'W', 'S', 'Y', 'K','M', 'R', 'W', 'S', 'Y', 'K')
for (i in 1:16) s2[which(s2==levels[i])]=cons[i]
s2=phyDat(s2)
write.phyDat(s2,format='fasta',file = sprintf('simuNoMaj%d.fas',seed)) 

# Write file with consensus sequences only
s3=as.character(s)
cons=c('A','C','G','T','A','A','A','C','C','G','C','G','T','G','T','T')
for (i in 1:16) s3[which(s3==levels[i])]=cons[i]
s3=phyDat(s3)
write.phyDat(s3,format='fasta',file=sprintf('simuCons%d.fas',seed)) 

# Write rates.txt file
if (file.exists(sprintf('rates%d.txt', seed))) o=file.remove(sprintf('rates%d.txt', seed))
for (i in 2:nrow(Q2)) {
  todo=Q2[i,1:(i-1),drop=F]
  write.table(todo,file = sprintf('rates%d.txt', seed),row.names = F,col.names=F,append = T)
}
cat('\n',file=sprintf('rates%d.txt', seed),append = T)
write.table(t(bf2),file = sprintf('rates%d.txt', seed),row.names = F,col.names=F,append = T)

# Compute some stats to see if the alignment is similar to real data

t=as.character(s)
nsnps=0
for (i in 1:ncol(t)) if (length(unique(t[,i]))>1) nsnps=nsnps+1
print(paste('Number of SNPs in full data:',nsnps))

t=as.character(s3)
nsnps=0
for (i in 1:ncol(t)) if (length(unique(t[,i]))>1) nsnps=nsnps+1
print(paste('Number of SNPs in consensus:',nsnps))

t=as.character(s2)
nsnps=0
nhet=0
for (i in 1:ncol(t)) if (length(unique(t[,i]))>1) {
  nsnps=nsnps+1
  if (length(setdiff(t[,i],c('a','c','g','t')))>0) nhet=nhet+1
  }
print(paste('Number of SNPs when using hets:',nsnps))
print(paste('Proportion of these SNPs that have hets:',nhet/nsnps))
}
