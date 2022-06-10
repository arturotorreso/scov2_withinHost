rm(list=ls())


##### MEASURE TEMPORAL SIGNAL OF CLUSTERS #####

library(ape)
library(RColorBrewer)
library(BactDating)
library(lubridate)
library(dplyr)


dir_prefix = "~/data/"

#### 1. PREPARE DATA ####

# Metadata
sample_meta = read.csv(paste(dir_prefix, 'metadata/all_metadata.csv', sep = ''), stringsAsFactors = F,na.strings = c(NA,''))
epi_meta = read.csv(paste(dir_prefix, 'metadata/tracing_meta.csv', sep = ''),
                    stringsAsFactors = F, na.strings = c('',NA))

seq_info = read.csv(paste(dir_prefix, 'metadata/sequence_info.csv', sep = ''), stringsAsFactors = F)
seq_info$id2 = paste(seq_info$id, seq_info$seq_id, sep='.')

seq_info = merge(x = seq_info[,c("id", "seq_id", "id2")],
                 y = sample_meta[,c("id", "GOSH._Staff", "Contact", "Patient", "External", "CT", "CT2", "Date_Sampled", 
                                    "Date_of_onset", "GOSH_Directorate", "Other_Directorate", "Code", "Outbreak_code", 
                                    "sequence", "pass")],
                 by = 'id', all.x = T, all.y= F)


seq_info$External = trimws(seq_info$External)
seq_info$GOSH_Directorate = trimws(seq_info$GOSH_Directorate)

seq_info = seq_info[!is.na(seq_info$id2),]

seq_info$External = trimws(seq_info$External)
seq_info$Code = trimws(seq_info$Code)


colnames(seq_info) = c("id", "seq_id", "id2", "GOSH_Staff", "Contact", "Patient", 
                       "External", "CT", "CT2", "Date_Sampled", "Date_of_onset", "GOSH_Directorate", 
                       "Other_Directorate", "Code", "Outbreak_code", "sequence", "pass")



#### 2. Get Trees ####

cons_tree = read.tree(paste(dir_prefix, 'phylogenetics/ml_tree/scov2.cons.tree', sep = ''))
full_tree = read.tree(paste(dir_prefix, 'phylogenetics/ml_tree/scov2.full.tree', sep = ''))

cons_tree = root(cons_tree, "NC_045512.2")
cons_tree = drop.tip(cons_tree, "NC_045512.2")

full_tree = root(full_tree, "NC_045512.2")
full_tree = drop.tip(full_tree, "NC_045512.2")


cons_tree$tip.label = unname(sapply(cons_tree$tip.label, function(x) strsplit(x, '[.]')[[1]][1]))
full_tree$tip.label = unname(sapply(full_tree$tip.label, function(x) strsplit(x, '[.]')[[1]][1]))

cons_tree = di2multi(cons_tree,tol = 1e-08)
full_tree = di2multi(full_tree,tol = 1e-08)

cons_tree$edge.length = cons_tree$edge.length*29903
full_tree$edge.length = full_tree$edge.length*29903


dates = setNames(decimal_date(mdy(sample_meta$Date_Sampled)), sample_meta$id)
dates = dates[full_tree$tip.label]

dna = read.dna(paste(dir_prefix, "/alignment/scov2.consensus.pass.fa", sep = ''), as.matrix = T, as.character = T, format='fasta')
colnames(dna) = 1:ncol(dna)
rownames(dna) = sapply(rownames(dna), function(x) strsplit(x, '[.]')[[1]][1])


## GET DISTANCES

cons_dists = cophenetic(cons_tree)
cons_hc = hclust(as.dist(cons_dists))
cons_clusts <- cutree(cons_hc, h=0.1)

dist_zero = names(table(cons_clusts)[table(cons_clusts) > 1])
df.clusts = data.frame('tips'=names(cons_clusts[cons_clusts %in% dist_zero]), 'cluster'=cons_clusts[cons_clusts %in% dist_zero])

df.clusts$dates = dates[as.character(df.clusts$tips)]
df.clusts$root2tip = setNames(leafDates(full_tree), full_tree$tip.label)[as.character(df.clusts$tips)]

df.clusts = df.clusts %>% group_by(cluster) %>%
  mutate(slope = coef(lm(root2tip~dates))[2],
         r2 = summary(lm(root2tip~dates))$r.squared,
         pval = ifelse(length(summary(lm(root2tip~dates))$coefficients[-1,4])>0,
                       summary(lm(root2tip~dates))$coefficients[-1,4],NA)) %>%
  as.data.frame()


df.clusts$slope[df.clusts$slope == 0] = NA
df.clusts = df.clusts[!df.clusts$cluster %in% df.clusts$cluster[is.na(df.clusts$dates)],]

cls = names(table(df.clusts$cluster)[table(df.clusts$cluster)>2])
df.clusts = df.clusts[df.clusts$cluster %in% cls,]



#### PLOT ####

pdf(paste(dir_prefix, 'figures/fig4_v2.pdf',sep = ''), height = 6.3, width = 8, useDingbats = F)
par(mfrow = c(3,4), oma = c(5,5,3,3))
np = 1
for (cls in unique(df.clusts$cluster)){
  df = df.clusts[df.clusts$cluster %in% cls,]
  
  if (is.na(unique(df$slope))){
    col = 'darkgrey'
    xlim = c(min(df$dates)-0.01, max(df$dates) + 0.01)
    ylim = c(min(df$root2tip - 5), max(df$root2tip + 5))
  } else if (unique(df$slope) > 0){
    col = 'dodgerblue3'
    xlim = c(min(df$dates)-0.01, max(df$dates) + 0.01)
    ylim = c(min(df$root2tip - 5), max(df$root2tip + 5))
  } else {
    col = 'red'
    xlim = c(min(df$dates)-0.01, max(df$dates) + 0.01)
    ylim = c(min(df$root2tip - 5), max(df$root2tip + 5))
  }
  
  par(mar = c(2, 2, 1.5, 1.5))
  
  plot(lubridate::date_decimal(df$dates), df$root2tip, type = 'n',
       xlim  = lubridate::date_decimal(xlim), ylim = ylim, yaxt = 'n', bty="n", ylab = '', xlab = '', cex.axis = 0.9)
  
  grid(lwd = 1.3)
  axis(2, las = 2)
  box(lwd = 1.3)

  
  points(lubridate::date_decimal(df$dates), df$root2tip, 
         pch = 19, col = adjustcolor(col, 0.6), cex = 1.2)
  
  if (!is.na(unique(df$slope))){
    m1 = lm(root2tip~lubridate::date_decimal(dates), data = df)
    abline(m1, lwd = 1.5, col = col)
  }
  
  if (np %in% c(1, 5, 9)){
    mtext(text = 'Root to tip\ndistance', side = 2, line = 3, outer = NA, cex = 0.8, font = 2)
  }
  if (np %in% c(9:12)){
    mtext(text = 'Collection date', side = 1, line = 2.8, outer = NA, cex = 0.8, font = 2)
  }
  
  
  np = np+1
  
}


dev.off()

