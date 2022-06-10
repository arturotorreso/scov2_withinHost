rm(list = ls())

library(ape)
library(vcfR)
library(dplyr)


dir_prefix = "~/data/"
meta_file = paste(dir_prefix, 'metadata/sequence_info.csv', sep = '')
mask_file = paste(dir_prefix, 'reference/scov2_masking.txt', sep = '')

# LOAD METADATA

sample_meta = read.csv(meta_file, stringsAsFactors = F)
mask = read.table(mask_file)$V1

# LOAD DNA
dna = read.dna(paste(dir_prefix, "illumina/alignment/scov2.consensus.pass.fa", sep = ''),
               as.matrix = T, as.character = T, format='fasta')

colnames(dna) = 1:ncol(dna)
dna = dna[-1,] # Remove reference genome


# FREQUENCIES OF DUPLICATES

# List of technical duplicates
dup_seqs =  read.table(paste(dir_prefix, '/illumina/dups_ids.txt', header = F))$V1

sample_meta = sample_meta[!is.na(sample_meta$id),]
sample_meta = sample_meta[paste(sample_meta$id, sample_meta$seq_id, sep='.') %in% dup_seqs,]
dups_ids = unique(sample_meta$id[duplicated(sample_meta$id)])
dups_ids = dups_ids[!is.na(dups_ids)]
dups_meta = sample_meta[sample_meta$id %in% dups_ids,]

out = data.frame('pos' = numeric(),
                 'id' = character(),
                 's1' = character(),
                 's2' = character(),
                 'af1' = numeric(),
                 'af2' = numeric(),
                 'alt1' = character(),
                 'alt2' = character(),
                 'gt1' = character(),
                 'gt2' = character(),
                 'pass1' = logical(),
                 'pass2' = logical(),
                 stringsAsFactors = F)

pb <- utils::txtProgressBar(min = 1, max = length(dups_ids), style = 3)
n = 0
for (id in dups_ids){
  n = n + 1
  utils::setTxtProgressBar(pb, n)
  sample_ids = dups_meta[dups_meta$id == id,]
  if (nrow(sample_ids) < 2){next()}
  
  for (i in 1:(nrow(sample_ids)-1)){
    for (j in (i+1):length(rownames(sample_ids))) {
      
      s1 = sample_ids$seq_id[i]
      s2 = sample_ids$seq_id[j]
      
      s1.code = paste(id,'.', s1, sep = '')
      s2.code = paste(id,'.', s2, sep = '')
      
      vcf1 = read.vcfR(paste(dir_prefix, 'illumina/vcf_freqs/',s1.code,'.hets.vcf.gz', sep = ''),verbose = F)
      vcf2 = read.vcfR(paste(dir_prefix, 'illumina/vcf_freqs/',s2.code,'.hets.vcf.gz', sep = ''),verbose = F)
      
      if (any(c(nrow(vcf1@fix), nrow(vcf2@fix))==0)){next()}
      # Remove indels
      vcf1 = vcf1[!grepl('INDEL', vcf1@fix[,8])]
      vcf2 = vcf2[!grepl('INDEL', vcf2@fix[,8])]
      
      # Filter variants
      vcf1 = vcf1[apply(vcf1@gt, 1, function(x) {
          length(strsplit(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'AD')], ',')[[1]]) > 1 &
          as.numeric(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'DP')]) >= 100 &
          as.numeric(strsplit(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'AD')], ',')[[1]][2]) >= 20 &
          as.numeric(strsplit(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'ADR')], ',')[[1]][2]) >= 5 &
          as.numeric(strsplit(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'ADF')], ',')[[1]][2]) >= 5 &
          as.numeric(strsplit(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'AD')], ',')[[1]][2]) / sum(as.numeric(strsplit(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'AD')], ',')[[1]])) < 0.95
      })]
      vcf2 = vcf2[apply(vcf2@gt, 1, function(x) {
          length(strsplit(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'AD')], ',')[[1]]) > 1 &
          as.numeric(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'DP')]) >= 50 &
          as.numeric(strsplit(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'AD')], ',')[[1]][2]) >= 20 &
          as.numeric(strsplit(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'ADR')], ',')[[1]][2]) >= 5 &
          as.numeric(strsplit(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'ADF')], ',')[[1]][2]) >= 5 &
          as.numeric(strsplit(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'AD')], ',')[[1]][2]) / sum(as.numeric(strsplit(strsplit(x[2], ':')[[1]][which(strsplit(x[1],':')[[1]] == 'AD')], ',')[[1]])) < 0.95
      })]
      
      
      # Remove missing genotypes (keep reference to compare allele frequencies)
      vcf1.pos = vcf1[!sapply(vcf1@gt[,2], function(x) strsplit(x, ':')[[1]][1] %in% c('./.'))]
      vcf2.pos = vcf2[!sapply(vcf2@gt[,2], function(x) strsplit(x, ':')[[1]][1] %in% c('./.'))]
      
      # Get all positions to be considered
      vcf.pos = sort(as.integer(unique(c(vcf1.pos@fix[,'POS'], vcf2.pos@fix[,'POS']))))
      
      vcf1 = vcf1[vcf1@fix[,'POS'] %in% vcf.pos]
      vcf2 = vcf2[vcf2@fix[,'POS'] %in% vcf.pos]
      
      pos = sort(as.integer(unique(c(vcf1@fix[,'POS'],vcf2@fix[,'POS']))))
      
      out_id = data.frame('pos' = numeric(),
                          'id' = character(),
                          's1' = character(),
                          's2' = character(),
                          'af1' = numeric(),
                          'af2' = numeric(),
                          'alt1' = character(),
                          'alt2' = character(),
                          'gt1' = character(),
                          'gt2' = character(),
                          'pass1' = logical(),
                          'pass2' = logical(),
                          stringsAsFactors = F)
      
      for (p in pos){
        
        vcf1p = vcf1[vcf1@fix[,'POS'] == p]
        vcf2p = vcf2[vcf2@fix[,'POS'] == p]
        
        # If the position is not in the vcf, then is very low depth
        if (!p %in% vcf1@fix[,'POS']){
          gt1 = '0'
          af1p = 0
          alt1p = unname(vcf2p@fix[,'ALT'])
          pass1 = TRUE
        } else {
          ad.pos = which(strsplit(vcf1p@gt[,1],':')[[1]] == 'AD') 
          dp.pos = which(strsplit(vcf1p@gt[,1],':')[[1]] == 'DP')
          dp = as.numeric(strsplit(vcf1p@gt, ':')[[2]][dp.pos])
          ad = as.numeric(strsplit(strsplit(vcf1p@gt, ':')[[2]][ad.pos], ',')[[1]][2])
          af1p =  ad / dp
          alt1p = unname(vcf1p@fix[,'ALT'])
          if (is.na(af1p)) af1p = 0
          if (is.na(alt1p)) alt1p = unname(vcf2p@fix[,'ALT'])
          pass1 = ifelse(strsplit(vcf1p@gt, ':')[[2]][length(strsplit(vcf1p@gt, ':')[[2]])] == 1, TRUE, FALSE)
          if (p %in% mask){pass1 = FALSE}
          gt1 = ifelse(af1p <= 0.95, '0/1', '1')
          if (gt1 == '0/1' &  dp < 100){pass1 = FALSE}
        }
        
        if (!p %in% vcf2@fix[,'POS']){
          gt2 = '0'
          af2p = 0
          alt2p = unname(vcf1p@fix[,'ALT'])
          pass2 = TRUE
        } else {
          ad.pos = which(strsplit(vcf2p@gt[,1],':')[[1]] == 'AD') 
          dp.pos = which(strsplit(vcf2p@gt[,1],':')[[1]] == 'DP')
          dp = as.numeric(strsplit(vcf2p@gt, ':')[[2]][dp.pos])
          ad = as.numeric(strsplit(strsplit(vcf2p@gt, ':')[[2]][ad.pos], ',')[[1]][2])
          af2p = ad / dp 
          alt2p = unname(vcf2p@fix[,'ALT'])
          if (is.na(af2p)) af2p = 0
          if (is.na(alt2p)) alt2p = unname(vcf1p@fix[,'ALT'])
          pass2 = ifelse(strsplit(vcf2p@gt, ':')[[2]][length(strsplit(vcf2p@gt, ':')[[2]])] == 1, TRUE, FALSE)
          if (p %in% mask){pass2 = FALSE}
          gt2 = ifelse(af2p <= 0.95, '0/1', '1')
          if (gt2 == '0/1' & dp < 100){pass2 = FALSE}
        }
        
        
        out_id = rbind(out_id, 
                       data.frame('pos' = p,
                                  'id' = id,
                                  's1' = s1.code,
                                  's2' = s2.code,
                                  'af1' = af1p,
                                  'af2' = af2p,
                                  'alt1' = alt1p,
                                  'alt2' = alt2p,
                                  'gt1' = gt1,
                                  'gt2' = gt2,
                                  'pass1' = pass1,
                                  'pass2' = pass2,
                                  stringsAsFactors = F))
      }
      
      out = rbind(out, out_id)

    }
  }
}


out2 = out

out2 = out2[!(out2$gt1 == '1' & out2$gt2 == '1'),]
out2 = out2[!(out2$gt1 == '0' & out2$gt2 == '0'),]

out2$af1[out2$gt1 == "0/1" & out2$af1 > 0.5] = 1 - out2$af1[out2$gt1 == "0/1" & out2$af1 > 0.5]
out2$af2[out2$gt2 == "0/1" & out2$af2 > 0.5] = 1 - out2$af2[out2$gt2 == "0/1" & out2$af2 > 0.5]


out2 = out2[apply(out2, 1, function(x) all(x['pass1'] == TRUE, x['pass2'] == TRUE)),]


sample_meta = read.csv(paste(dir_prefix, 'metadata/all_metadata.csv', sep = ''), stringsAsFactors = F,na.strings = c(NA,''))
dp = read.delim(paste(dir_prefix, 'illumina/all_samples.DP.txt', sep = ''), stringsAsFactors = F,na.strings = c(NA,''), header = T, check.names = F)
bq = read.delim(paste(dir_prefix, 'illumina/all_samples.meanBQ.txt', sep = ''), stringsAsFactors = F,na.strings = c(NA,''), header = T, check.names = F)


out2$bq1 = apply(out2, 1, function(x) bq[as.numeric(x['pos']), as.character(x['s1'])])
out2$bq2 = apply(out2, 1, function(x) bq[as.numeric(x['pos']), as.character(x['s2'])])
out2$dp1 = apply(out2, 1, function(x) dp[as.numeric(x['pos']), as.character(x['s1'])])
out2$dp2 = apply(out2, 1, function(x) dp[as.numeric(x['pos']), as.character(x['s2'])])

out2$bq = apply(out2[,c('bq1','bq2')], 1, min)
out2$dp = apply(out2[,c('dp1','dp2')], 1, min)

out2$shared = as.numeric(out2$gt1 == "0/1" & out2$gt1 == "0/1")
out2 = merge(x = out2, y = sample_meta[,c('id','CT')], all.x = T, all.y = F, by = 'id')



###################
## Ct value plot ##
###################

het_n = out2 %>% group_by(id) %>% mutate(nhets = n()) %>% distinct(id,nhets) %>% as.data.frame()
het_n = merge(x = het_n, y = sample_meta[,c('id','CT')], all.x = T, all.y = F, by = 'id')



het_n$phets = 0
for (i in 1:nrow(het_n)){
  id = het_n$id[i]
  id_hets = out2[out2$id == id,]
  het_n$phets[i] = nrow(id_hets[id_hets$gt1 == '0/1' & id_hets$gt2 == '0/1',])/nrow(id_hets)
}


pdf(paste(dir_prefix, 'figures/fig3b.pdf', sep =''), height = 5, width = 5)

par(mar = c(7,7,6,6))
plot(het_n$CT, het_n$phets, type = 'n', xaxt='n', yaxt='n', xlab = '',ylab = '', ann=T, axes=F)
grid(lwd = 1.7)
box(lwd = 2)

axis(1)
mtext('Ct value', 1, 3, font = 2, cex = 1.2)
axis(2, las = 2)
mtext('Proportion of shared\nmixed variants', 2, 3.5, font = 2, cex = 1.2)



m1 = glm(as.numeric(phets) ~ as.numeric(CT), data = het_n, family = quasibinomial(link = 'logit'))
# m1 = glm(as.numeric(phets) ~ as.numeric(CT), data = het_n)
preds = predict(m1, newdata = data.frame('CT'=seq(20,40,0.1)),se.fit = T, type = 'response')

critval <- 1.96 ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

lwr[lwr < 0] = 0
upr[upr > 1] = 1

polygon(x = c(seq(20,40,0.1), rev(seq(20,40,0.1))), y = c(lwr, rev(upr)), col = adjustcolor('steelblue', 0.3), border = F)
# polygon(x = c(seq(20,40,0.5), rev(seq(20,40,0.5))), y = c(preds$fit - preds$se.fit, rev(preds$fit + preds$se.fit)), col = adjustcolor('steelblue', 0.5), border = F)

lines(seq(20,40,0.1), fit, lwd = 2)
points(het_n$CT, het_n$phets, cex = 1.1, pch = 21, bg = 'grey', lwd = 1.6)


dev.off()



## PER SAMPLE 2
library(colorspace)
library(randomcoloR)

het_n = out2 %>%
  group_by(s1) %>% mutate(nhets1 = sum(gt1 == '0/1')) %>% 
  group_by(s2) %>% mutate(nhets2 = sum(gt2 == '0/1')) %>%
  distinct(id,nhets1,nhets2) %>% as.data.frame()

het_n = merge(x = het_n, y = sample_meta[,c('id','CT')], all.x = T, all.y = F, by = 'id')

cols1 = distinctColorPalette(nrow(het_n))
cols2 = lighten(cols1, amount = 0.5)

pdf(paste(dir_prefix, 'figures/suppfig2.pdf', sep =''), height = 5, width = 5)

par(mar = c(7,7,6,6))

plot(c(het_n$CT,het_n$CT), c(het_n$nhets1, het_n$nhets2) + 10, type = 'n', xaxt='n', yaxt='n', xlab = '',ylab = '', ann=T, axes=F)
grid(lwd = 1.7)
box(lwd = 2)

axis(1)
mtext('Ct value', 1, 3, font = 2, cex = 1.2)
axis(2, las = 2)
mtext('Number of mixed variants', 2, 3.5, font = 2, cex = 1.2)

apply(het_n, 1, function(x) lines(x = c(x['CT'],x['CT']), y = c(x['nhets1'], x['nhets2']), lwd = 1.3))
points(het_n$CT, het_n$nhets1, pch = 21, bg = adjustcolor(cols1, 1), lwd = 1.3, cex = .9)
points(het_n$CT, het_n$nhets2, pch = 21, bg = adjustcolor(cols2, 1),lwd = 1.3, cex = .9)

dev.off()



####################
## Freq1 vs Freq2 ##
####################


####  1. Per Sample ####

cols = apply(out2[, c('gt1', 'gt2')], 1, function (x)
               ifelse(any(x == '0' | x == '1'), '#939597',
                      '#F5DF4D'))

plot_id = cbind(out2, cols)


####  2. All together + Color by CT ####

plot_id2 = merge(x = plot_id,y = sample_meta[,c('id','CT')],
                 by = 'id', all.x = T, all.y = F)
colnames(plot_id2)[colnames(plot_id2) == 'CT'] = 'ct'
plot_id2$ct = as.numeric(plot_id2$ct)

range.i <- range(plot_id2$ct[!is.na(plot_id2$ct)])
range.i[1] = range.i[1] - 0.1
range.i[2] = range.i[2] + 0.1
extr <- min(abs(range.i))
# extr  <- 3



scale01 <- function(x, low = min(x), high = max(x)) {
  x <- (x - low)/(high - low)
  x
}


pdf(paste(dir_prefix, 'figures/fig3a.pdf', sep =''), height = 5, width = 5)

ncols = 50
par(mar=c(0,0,0,0))
layout(matrix(c(1,2,3,2), nrow=2, byrow=F), heights=c(0.15, 0.8, 0.15), widths=c(0.4,0.8,0.6))

#1 -scale
par(mar=c(2,2,2,0))
z <- seq(range.i[1],range.i[2], length = ncols+1)

image(z = matrix(z), col = colorRampPalette(rev(brewer.pal(n = 7, 
                                                           name = "RdBu")))(ncols), breaks = z, xaxt = "n", yaxt = "n", axes=F)
box(lwd=1)
lv <- pretty(z)
xv <- scale01(as.numeric(lv),range.i[1], range.i[2])
axis(1, at = xv, labels = lv, las=1, font=1, lwd = 1, cex=0.8, cex.axis=0.7, padj=-1.5)
mtext(1, text = 'Ct value',1.5, font = 2)

#2 - plot
# par(mar = c(2,7,4,4))

plot_df = plot_id2
plot_df$ct_bin = cut(plot_id2$ct, breaks = z)

cols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(ncols)
cols = setNames(cols, levels(plot_df$ct_bin))
plot_df$cols = cols[plot_df$ct_bin]

par(mar = c(6,7,4,4))

plot(plot_id2$af1[order(plot_id2$ct,decreasing = T)], plot_id2$af2[order(plot_id2$ct,decreasing = T)],
     bg = plot_df$cols[order(plot_id2$ct, decreasing = T)],
     cex = 1,pch=21,
     ylab = "", xlab = "", axes = FALSE, xlim = c(0,0.5), ylim = c(0,0.5), xpd = NA)
box(lwd = 2)
axis(1); mtext('Allele frequency of duplicate 1', 1, 3.5,font = 2)
axis(2, las = 1); mtext('Allele frequency of duplicate 2', 2, 4,font=2)


par(mar = c(0,0,0,0))
plot.new()

dev.off()




####  2. All together + Color by CT + zoom-in  ####

pdf(paste(dir_prefix, 'figures/fig3a.pdf', sep =''), height = 5, width = 8)

ncols = 50
par(mar=c(0,0,0,0))
layout(matrix(c(1,3,4,2), nrow=2, byrow=F), heights=c(0.2, 0.8, 0.2, 0.8), widths=c(0.5,0.5,0.5,0.5))

#1 -scale
par(mar=c(1.5,7,2,6))
z <- seq(range.i[1],range.i[2], length = ncols+1)

image(z = matrix(z), col = colorRampPalette(rev(brewer.pal(n = 7, 
                                                           name = "RdBu")))(ncols), breaks = z, xaxt = "n", yaxt = "n", axes=F)
box(lwd=1)
lv <- pretty(z)
xv <- scale01(as.numeric(lv),range.i[1], range.i[2])
axis(1, at = xv, labels = lv, las=1, font=1, lwd = 1, cex=0.8, cex.axis=0.7, padj=-1.5)
mtext(1, text = 'Ct value',1.5, font = 2)

#2 - plot
# par(mar = c(2,7,4,4))

plot_df = plot_id2
plot_df$ct_bin = cut(plot_id2$ct, breaks = z)

cols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(ncols)
cols = setNames(cols, levels(plot_df$ct_bin))
plot_df$cols = cols[plot_df$ct_bin]


par(mar = c(6,4,4,6))

plot(plot_id2$af1[order(plot_id2$ct,decreasing = T)], plot_id2$af2[order(plot_id2$ct,decreasing = T)],
     bg = plot_df$cols[order(plot_id2$ct, decreasing = T)],
     cex = 1,pch=21, type = 'n',
     ylab = "", xlab = "", axes = FALSE, xlim = c(0,0.5), ylim = c(0,0.5), xpd = NA)


polygon(x = c(0.1, 0.1, par('usr')[1], -0.29, -0.83, -0.83, par('usr')[1]), 
        y = c(par('usr')[1], 0.1, 0.1, 0.52, 0.52, par('usr')[1], par('usr')[1]),
        
        col = adjustcolor('grey', 0.5), xpd = NA, border = F)


points(plot_id2$af1[order(plot_id2$ct,decreasing = T)], plot_id2$af2[order(plot_id2$ct,decreasing = T)],
       bg = plot_df$cols[order(plot_id2$ct, decreasing = T)],cex = 1,pch=21)
box(lwd = 2)
axis(1); mtext('Allele frequency of duplicate 1', 1, 3.5,font = 2)
axis(2, las = 1)

par(mar = c(6,7,4,3))

plot(plot_id2$af1[order(plot_id2$ct,decreasing = T)], plot_id2$af2[order(plot_id2$ct,decreasing = T)],
     bg = plot_df$cols[order(plot_id2$ct, decreasing = T)],
     cex = 1,pch=21,
     ylab = "", xlab = "", axes = FALSE, xlim = c(0,0.1), ylim = c(0,0.1))
box(lwd = 2)
axis(1); mtext('Allele frequency of duplicate 1', 1, 3.5,font = 2)
axis(2, las = 1); mtext('Allele frequency of duplicate 2', 2, 4,font=2)


par(mar = c(0,0,0,0))
plot.new()

dev.off()





################
## DP/BQ plot ##
################


out2$bq1 = apply(out2, 1, function(x) bq[as.numeric(x['pos']), as.character(x['s1'])])
out2$bq2 = apply(out2, 1, function(x) bq[as.numeric(x['pos']), as.character(x['s2'])])

out2$dp1 = apply(out2, 1, function(x) dp[as.numeric(x['pos']), as.character(x['s1'])])
out2$dp2 = apply(out2, 1, function(x) dp[as.numeric(x['pos']), as.character(x['s2'])])


pdf(paste(dir_prefix, 'figures/techRep_bq_dp.pdf', sep =''), height = 6, width = 18)

par(mfrow = c(1,2), cex.lab = 1.3, font.lab = 2, oma = c(1,1,1,1))
plot(out2$dp1, out2$dp2, xlab = 'Depth in replicate 1', ylab = 'Depth in replicate 2')
plot(out2$bq1[out2$bq2 > 20], out2$bq2[out2$bq2 > 20],
     xlab = 'Base Quality in replicate 1', ylab = 'Base Quality in replicate 2')

text(x = 13, y = 52,'Depth and base quality in technical replicates',
     cex = 1.5, font = 2,xpd = NA, adj = 0)
dev.off()


pdf(paste(dir_prefix, 'figures/techRep_bq_dp2.pdf', sep =''), height = 5, width = 10)
par(mfrow = c(1,2), cex.lab = 1.3, font.lab = 2, oma = c(1,1,1,1), mar = c(4,6,6,4))

with(out2, boxplot(minDP ~ shared, names = c('Not shared', 'Shared'),
     xlab = '', ylab = 'Minimum depth\nof the technical pair'))
with(out2, boxplot(minBQ ~ shared, names = c('Not shared', 'Shared'),
                   xlab = '', ylab = 'Minimum base quality\nof the technical pair'))
text(x = -4, y = 55,
'Differences in depth and base
quality between shared and non-shared
subconsensus calls in technical replicates',
     cex = 1.5, font = 2,xpd = NA, adj = 0)
dev.off()



##############################
## Share hets by AF dotplot ##
##############################

het_n = merge(x = out2, y = sample_meta[,c('id','CT')], all.x = T, all.y = F, by = 'id')
het_n = het_n[!duplicated(het_n$id),c('id','CT')]

for (maf in c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2)){
  het_n$phets = 0
  het_n$nhets = 0
  for (i in 1:nrow(het_n)){
    id = het_n$id[i]
    id_hets = out2[out2$id == id,]
    id_hets$maf = apply(id_hets, 1, function(x) ifelse(all(c(as.numeric(x['af1']),as.numeric(x['af2'])) != 0),
                                                       min(c(as.numeric(x['af1']),as.numeric(x['af2']))),
                                                       c(as.numeric(x['af1']),as.numeric(x['af2']))[which(c(as.numeric(x['af1']),as.numeric(x['af2'])) != 0)]))
    id_hets = id_hets[id_hets$maf >= maf,]
    het_n$phets[i] = nrow(id_hets[id_hets$gt1 == '0/1' & id_hets$gt2 == '0/1',])/nrow(id_hets)
    het_n$nhets[i] = nrow(id_hets)
  }
  colnames(het_n)[colnames(het_n) == 'phets'] = paste('phets_',maf,sep='')
  colnames(het_n)[colnames(het_n) == 'nhets'] = paste('nhets_',maf,sep='')
}

het_n = het_n[order(het_n$CT),]
freqs = as.factor(c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2))

ct=as.numeric(het_n$CT)
ii <- cut(ct, breaks = seq(floor(min(ct)), ceiling(max(ct)), by = 2), 
          include.lowest = TRUE)

cols <- colorRampPalette(c('blue', 'white','red'))(length(levels(ii)))[ii]


pdf(paste(dir_prefix, 'figures/suppfig3.pdf', sep =''), height = 6, width = 8)

par(mfrow=c(5,4), oma = c(6,7,2,2))

# Color legend
par(mar = c(6,0,2,1))
image(z = matrix(c(19,21,23,25,27,29,31)),
      col = colorRampPalette(c('blue', 'white','red'))(length(levels(ii))), xaxt = "n", yaxt = "n", axes=F)

axis(1, at = seq(0, 1, length = 7), labels = c(19,21,23,25,27,29,31), las=1,
     font=1, lwd = 1, cex=0.8, cex.axis=1, tck= -0.2, padj=0)
box()
mtext(text = 'Ct value', side = 1, line = 3, cex = 0.9, font=2)

# Size legend
pt_sizes = c(1, 5, 10, 20, 100, 200)
par(mar = c(6,2,2,4))
plot(1,1,type='n', xaxt='n', yaxt='n', ann=FALSE, axes=F, xlim=c(0,1), ylim=c(0,1))
legend(-0.3,1, cex = 1.2, legend = pt_sizes, pch=21, pt.cex = log10(pt_sizes) + 0.5, ncol = 3,bty = "n",xpd=NA, adj = 0)
mtext(text = 'Number of\nmixed variants', side = 1, line = 5.5, cex = 0.9, font=2, xpd=NA,
      padj = -1)


for(i in 1:nrow(het_n)){
  par(mar=c(1,0.2,1.2,0.2))
  freqs = as.factor(c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2))
  plot(freqs,
       rep(-10, length(freqs)), ylim=c(0,1.1),
       type = 'n', xaxt='n', yaxt='n', xlab = '',ylab = '', ann=T, axes=F)
  mtext(text = het_n$id[i], side = 3, line = 0.2, font = 2, cex = 0.8)
  box(lwd=1.5)
  grid(lwd=1.5)
  if (i %in% c(14,15,16,17)){
    axis(side = 1, at = freqs, labels = as.numeric(as.vector(freqs))*100,
         cex.axis = 1.1,  gap.axis=0)
  } else {
    axis(side = 1, at = freqs, labels = F,lwd.tick=0)
  }
  
  if (i %in% c(3,7,11,15)){
    axis(side = 2, las = 2, cex.axis = 1.2)
  } else {
    axis(side = 2, labels = F, lwd.tick=0)
  }
  
  het_sample = het_n[i,]
  sample_pts = het_sample[,c("phets_0", "phets_0.01", 
                             "phets_0.02", "phets_0.03", "phets_0.04", 
                             "phets_0.05", "phets_0.1", "phets_0.2")]
  
  sample_ns = het_sample[,c("nhets_0", "nhets_0.01", 
                            "nhets_0.02", "nhets_0.03", "nhets_0.04", 
                            "nhets_0.05", "nhets_0.1", "nhets_0.2")]
  
  sample_col = cols[i]
  sample_sizes = log10(sample_ns) + 0.5
  
  lines(freqs, sample_pts, lwd = 1.3, type = 'c')
  points(freqs, sample_pts, cex = as.numeric(sample_sizes), lwd = 1.3, pch = 21, bg = sample_col, xpd = NA)
  
}

mtext("Proportion of shared mixed variants", side = 2, line = 4,outer = TRUE, cex = 1.1, font=2)
mtext("Minimum allele frequency (%)", side = 1, line = 3,outer = TRUE, cex = 1.1, font=2)


dev.off()

