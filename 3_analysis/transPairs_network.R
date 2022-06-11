
###################
## Create graphs ##
###################

sample_meta = read.csv('~/data/metadata/all_metadata.csv', stringsAsFactors = F,na.strings = c(NA,''))
rownames(sample_meta) = sample_meta$id


rownames(test2) = sapply(rownames(test2), function(x) strsplit(x, '-')[[1]][2])
colnames(test2) = sapply(colnames(test2), function(x) strsplit(x, '-')[[1]][2])


## Threshold to show 
test2[test2<lt] = 0

# Create a nel.graph object from a matrix
nel.graph <-as(test2, Class="graphNEL")

# Find an optimum using Edmonds algorithm
search <- edmondsOptimumBranching(nel.graph)

#Output transmission tree
mat=cbind(search$edgeList[1,],search$edgeList[2,],t(search$weights))
mat = as.data.frame(mat)
mat[,3] = as.numeric(mat[,3])

links = as.data.frame(mat)
colnames(links) = c('from', 'to', 'weight')


# Graph with all possible links
net = graph_from_graphnel(nel.graph, name = TRUE, weight = TRUE, unlist.attrs = TRUE)



###################
##  PLOT graphs  ##
###################

# Only one out of two curved

ENDS = ends(net, E(net))
Curv = rep(F, nrow(ENDS))
skip = c()
for(i in 1:nrow(ENDS)) {
  if(i %in% skip){next()}
  Curv[i] = are.connected(net, ENDS[i,2], ENDS[i,1])
  skip = c(skip, intersect(which(ENDS[,1] == ENDS[i,2]), which(ENDS[,2] == ENDS[i,1])))
}

E(net)$curved = Curv
E(net)$curved[E(net)$curved] = 1

edge.col = rep('gray85', length(E(net)))
edge.col[apply(links, 1, function(x) which(attr(E(net), 'vnames') == paste(x[1:2], collapse = '|')))] = '#E69F00'


pdf(paste('~/data/figures/transmission/',id,'_network.pdf', sep = ''), width = 5.5, height = 5.5, useDingbats = F)
par(mar = c(7,7,7,7))
plot(net, edge.width=E(net)$weight, edge.curved=Curv, edge.color=edge.col, vertex.shape="square",
     vertex.label.font=2, vertex.label.color="gray40", vertex.size = 30, vertex.color = 'white',
     vertex.label.cex=.70)
dev.off()


