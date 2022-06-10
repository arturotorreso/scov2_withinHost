rm(list = ls())

library(ape)
library(vcfR)
library(dplyr)

dir_prefix = "~/data/"

getRelation = function(s1,s2,meta){
  
  if (!any(c(s1,s2) %in% meta$id)){
    return(NA)
  }
  # Default option: there is no relationship
  rel = 'none'
  
  # Option 1: Same external hospital
  if (!any(is.na(meta$External[meta$id %in% c(s1,s2)]))){
    if (length(unique(meta$External[meta$id %in% c(s1,s2)])) == 1){
      rel = 'hosp'
    }
  }
  
  # Option 1.2: GOSH
  if (!any(is.na(meta$GOSH._Staff[meta$id %in% c(s1,s2)]))){
    rel = 'gosh'
  }
  
  # Option 2: Same department within GOSH
  if (all(c(s1,s2) %in% meta$id)){
    if (!any(is.na(meta$GOSH_Directorate[meta$id %in% c(s1,s2)]))){
      if (length(unique(meta$GOSH_Directorate[meta$id %in% c(s1,s2)])) == 1){
        rel = 'dept'
      }
    }
  }
  
  # Option 3: Samples belong to same outbreak
  if (!any(is.na(meta$Outbreak_code[meta$id %in% c(s1,s2)]))){
    if (length(unique(meta$Outbreak_code[meta$id %in% c(s1,s2)])) == 1){
      rel = 'epi'
    }
  }
  
  # Option 4: Sample it's a duplicate
  if (!any(is.na(meta$Code[meta$id %in% c(s1,s2)]))){
    if (length(unique(meta$Code[meta$id %in% c(s1,s2)])) == 1){
      rel = 'dup'
    }
  }
  
  # Option 5: It's a technical replicate
  if (s1 == s2){
    rel = 'tech_dup'
  }
  
  return(rel)
}


# LOAD METADATA
mask_file = paste(dir_prefix, 'reference/scov2_masking.txt', sep = '')
mask = read.table(mask_file)$V1

seq_info = read.csv(paste(dir_prefix, 'metadata/sequence_info.csv',sep = ''), stringsAsFactors = F)
sample_meta = read.csv(paste(dir_prefix, '/metadata/all_metadata.csv', sep = ''), stringsAsFactors = F,na.strings = c(NA,''))


# LOAD DNA
full_dna = as.matrix(seqinr::read.alignment(paste(dir_prefix, "alignment/scov2.full.pass.fa", sep = ''),format='fasta'))
colnames(full_dna) = 1:ncol(full_dna)
full_dna = full_dna[-1,]

cons_dna = as.matrix(seqinr::read.alignment(paste(dir_prefix, "alignment/scov2.consensus.pass.fa", sep = ''),format='fasta'))
colnames(cons_dna) = 1:ncol(cons_dna)
cons_dna = cons_dna[-1,]

table(apply(full_dna, 2, function(x) length(unique(x[! x %in% '-']))>1))
table(apply(cons_dna, 2, function(x) length(unique(x[! x %in% 'n']))>1))

hets = full_dna[,apply(full_dna, 2, function(x) any(!x %in% c('a', 'c', 'g', 't', '-')))]
hets2 = apply(hets, 1, function(x) table(x)[!names(table(x)) %in% c('a', 'c', 'g', 't', '-')])
table(unlist(lapply(hets2, length)))


# LOAD TREE
cons_tree = read.tree(paste(dir_prefix, 'phylogenetics/ml_tree/scov2.cons.tree', sep = ''))
full_tree = read.tree(paste(dir_prefix, 'phylogenetics/ml_tree/scov2.full.tree', sep = ''))

cons_tree = root(cons_tree, "NC_045512.2")
cons_tree = drop.tip(cons_tree, "NC_045512.2")

full_tree = root(full_tree, "NC_045512.2")
full_tree = drop.tip(full_tree, "NC_045512.2")

full_tree$edge.length = full_tree$edge.length * 29903
cons_tree$edge.length = cons_tree$edge.length * 29903


# FRECUENCIES OF DUPLICATES

seq_info = seq_info[!is.na(seq_info$id),]
seq_info = seq_info[paste(seq_info$id, seq_info$seq_id, sep='.') %in% rownames(full_dna),]

sample_meta = sample_meta[sample_meta$id %in% seq_info$id,]



vcf_dir = paste(dir_prefix, 'vcf/', sep ='')

vcfs_list = parallel::mclapply(1:nrow(seq_info), function(i){
  s1 = seq_info$id[i]
  s1.seq = seq_info$seq_id[i]
  s1.code = paste(s1, '.', s1.seq, sep = '')
  vcf1 = read.vcfR(paste(vcf_dir, s1.code, '.hets.vcf.gz', sep = ''),verbose = F)
  vcf1 = vcf1[!grepl('INDEL', vcf1@fix[,8])]
  
  vcf1 = vcf1[!vcf1@fix[,'POS'] %in% mask]
  pos = names(full_dna[s1.code,][!full_dna[s1.code,] %in% c('a','c','g','t','-')])
  vcf1 = vcf1[vcf1@fix[,'POS'] %in% pos]
  return(vcf1)
}, mc.cores = 6)


names(vcfs_list) = unlist(lapply(vcfs_list, function(x) colnames(x@gt)[2]))



########################
### Get frequencies ###
#######################

seq_info2 = merge(x = seq_info, y = sample_meta[,c('id', 'CT')], all.x = T, all.y = F, by = 'id')

seq_info2 = seq_info2[paste(seq_info2$id, '.', seq_info2$seq_id,sep = '') %in% rownames(full_dna),]
out = setNames(data.frame(matrix(ncol = 11, nrow = 0)), c('id1', 'id2',"s1", "s2", "rel", "af1", 'af2', 'ct1', 'ct2', 'is_shared', 'pos'))

pb <- utils::txtProgressBar(min = 1, max = nrow(seq_info2), style = 3)
for (i in 1:(nrow(seq_info2) - 1)) {
  utils::setTxtProgressBar(pb, i)
  for (j in (i + 1):nrow(seq_info2)) {
    s1 = seq_info2$id[i]
    s2 = seq_info2$id[j]
    
    s1.seq = seq_info2$seq_id[i]
    s2.seq = seq_info2$seq_id[j]
    
    s1.code = paste(s1, '.', s1.seq, sep = '')
    s2.code = paste(s2, '.', s2.seq, sep = '')
    
    vcf1 = vcfs_list[[s1.code]]
    vcf2 = vcfs_list[[s2.code]]
    
    rel = getRelation(s1,s2,sample_meta)
    ct1 = as.numeric(sample_meta$CT[sample_meta$id == s1])
    ct2 = as.numeric(sample_meta$CT[sample_meta$id == s2])
    
    if (length(ct1) == 0) {ct1 = NA}
    if (length(ct2) == 0) {ct2 = NA}
    
    het_pos = sort(as.integer(unique(c(vcf1@fix[,'POS'],vcf2@fix[,'POS']))))
    if(length(het_pos) == 0){next()}
    hets.common = as.integer(Reduce(intersect, list(vcf1@fix[,'POS'],vcf2@fix[,'POS'])))
    is_shared = if (length(hets.common) > 0) (het_pos %in% hets.common) else {rep(FALSE, length(het_pos))}
    
    af1 = setNames(apply(vcf1@gt, 1, function (x) strsplit(x[2], ':')[[1]][which(strsplit(x[1], ':')[[1]] == 'AF')]),
                   vcf1@fix[,'POS'])
    af2 = setNames(apply(vcf2@gt, 1, function (x) strsplit(x[2], ':')[[1]][which(strsplit(x[1], ':')[[1]] == 'AF')]),
                   vcf2@fix[,'POS'])
    
    out = rbind(out, data.frame('id1' = s1, 'id2' = s2, 's1'  =  s1.code, 's2'  =  s2.code,
                                'rel' = rel,
                                'af1' = unname(af1[as.character(het_pos)]), 'af2' = unname(af2[as.character(het_pos)]),
                                'ct1' = ct1, 'ct2' = ct2, 'is_shared' = is_shared, pos = het_pos))
  }
}


out[out$rel == 'gosh','rel'] = 'hosp'

out2 = out
out2$af1 = as.numeric(out2$af1)
out2$af2 = as.numeric(out2$af2)


out2$freq = as.numeric(apply(out2, 1, function(x) ifelse(!any(is.na(c(x['af1'], x['af2']))), min(c(x['af1'], x['af2'])),
                                                         c(x['af1'], x['af2'])[which(!is.na(c(x['af1'], x['af2'])))])))

out2$freq = ifelse(out2$freq > 0.5, 1 - out2$freq, out2$freq)


out2$freq2 = cut(out2$freq, breaks = c(0, 0.05, 0.1, 0.5))
out2$pair = apply(out2, 1, function(x) paste(x['s1'], '_', x['s2'], sep = ''))


####################
###   Logistic  ####
####################

shared_vars = out2[,c("id1", "id2", "s1", "s2", "rel", "af1", "af2", "ct1", "ct2", 
                      "is_shared", "freq","pos", "pair")]

shared_vars$is_shared = as.numeric(shared_vars$is_shared)
shared_vars$rel = factor(shared_vars$rel, levels = c("none", "hosp", "dept", "epi", "dup", "tech_dup"))

m1 = glm(is_shared ~ freq + rel, data = shared_vars, family = binomial(link=logit))

intercept = coef(m1)[1]
coefs = confint(m1)
coefs = cbind(coef(m1),coefs)
colnames(coefs) = c('fit',  "2.5 %",  "97.5 %")

exp(intercept + coefs)/(1+(exp(intercept+coefs)))
                                     
                                     
exp(summary(m1)$coefficients[,1])/(1+exp(summary(m1)$coefficients[,1]))

new_data = data.frame('freq' = rep(seq(0,0.5,0.01), 6),
                      'rel' = as.vector(sapply(c('none', 'hosp', 'dept', 'epi', 'dup', 'tech_dup'), function(x) rep(x, length(seq(0,0.5,0.01))))))

## grad the inverse link function
ilink <- family(m1)$linkinv
## add fit and se.fit on the **link** scale
new_data <- bind_cols(new_data, setNames(as_tibble(predict(m1, new_data, se.fit = TRUE)[1:2]),
                                   c('fit_link','se_link')))
## create the interval and backtransform
new_data <- mutate(new_data,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))


new_data$fit_resp = log10(new_data$fit_resp)
new_data$right_upr = log10(new_data$right_upr)
new_data$right_lwr = log10(new_data$right_lwr)


pdf(paste(dir_prefix, 'figures/suppFig_logProb.pdf', sep = ''), height = 5.60, width = 6.89, useDingbats = F)

par(mar = c(6,8,4,8))
ys = c(0.001, 0.01, 0.1, 1)

cols = brewer.pal(9, "Blues")[4:9]
with(new_data[new_data$rel == 'none',],plot(freq, fit_resp,bty="n",
                                            col = cols[1], lwd = 2, type = 'l',ylim = c(min(log10(ys)),log10(1)),
                                            xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', panel.first = grid(lwd = 1.3, col = 'grey80',ny = NA)))
# box(lwd=1.3)
abline(h=log10(ys), col="grey80", lty=3, lwd = 1.3)

with(new_data[new_data$rel == 'none',],polygon(c(freq,rev(freq)), c(right_lwr, rev(right_upr)),
                                               col = adjustcolor(cols[1], 0.4),border = F))

axis(1, cex.axis = 1.1)
axis(2, at = log10(ys), labels = ys, cex.axis = 1.1, las = 2)
for (i in log10(ys)){
  axis(side=2, at=log10(seq(2,9))+ i,
       labels=NA, tck = -0.01)
}

mtext('Allele frequency', 1, 3, font = 2, cex = 1.1)
mtext('Probability a\nmixed variant is shared', 2, 4, font = 2, cex = 1.1)



with(new_data[new_data$rel == 'hosp',],lines(freq, fit_resp, col = cols[2], lwd = 2))
with(new_data[new_data$rel == 'hosp',],polygon(c(freq,rev(freq)), c(right_lwr, rev(right_upr)),
                                               col = adjustcolor(cols[2], 0.4), border = F))
with(new_data[new_data$rel == 'dept',],lines(freq, fit_resp, col = cols[3], lwd = 2))
with(new_data[new_data$rel == 'dept',],polygon(c(freq,rev(freq)), c(right_lwr, rev(right_upr)),
                                               col = adjustcolor(cols[3], 0.4), border = F))
with(new_data[new_data$rel == 'epi',],lines(freq, fit_resp, col = cols[4], lwd = 2))
with(new_data[new_data$rel == 'epi',],polygon(c(freq,rev(freq)), c(right_lwr, rev(right_upr)),
                                               col = adjustcolor(cols[4], 0.4), border = F))
with(new_data[new_data$rel == 'dup',],lines(freq, fit_resp, col = cols[5], lwd = 2))
with(new_data[new_data$rel == 'dup',],polygon(c(freq,rev(freq)), c(right_lwr, rev(right_upr)),
                                              col = adjustcolor(cols[5], 0.4), border = F))
with(new_data[new_data$rel == 'tech_dup',],lines(freq, fit_resp, col = cols[6], lwd = 2))
with(new_data[new_data$rel == 'tech_dup',],polygon(c(freq,rev(freq)), c(right_lwr, rev(right_upr)),
                                              col = adjustcolor(cols[6], 0.4), border = F))

text(x = 0.53, y = new_data$fit_resp[new_data$rel == 'none'][nrow(new_data[new_data$rel == 'none',])]-0.1,
     'None', col = cols[1], lwd = 1.3, cex = 1, xpd = NA, adj = 0)
text(x = 0.53, y = new_data$fit_resp[new_data$rel == 'hosp'][nrow(new_data[new_data$rel == 'hosp',])]+0.1,
     'Hospital', col = cols[2], lwd = 1.3, cex = 1, xpd = NA, adj = 0)
text(x = 0.53, y = new_data$fit_resp[new_data$rel == 'dept'][nrow(new_data[new_data$rel == 'dept',])],
     'Department', col = cols[3], lwd = 1.3, cex = 1, xpd = NA, adj = 0)
text(x = 0.53, y = new_data$fit_resp[new_data$rel == 'epi'][nrow(new_data[new_data$rel == 'epi',])],
     'Epidemiological', col = cols[4], lwd = 1.3, cex = 1, xpd = NA, adj = 0)
text(x = 0.53, y = new_data$fit_resp[new_data$rel == 'dup'][nrow(new_data[new_data$rel == 'dup',])],
     'Replicate', col = cols[5], lwd = 1.3, cex = 1, xpd = NA, adj = 0)
text(x = 0.53, y = new_data$fit_resp[new_data$rel == 'tech_dup'][nrow(new_data[new_data$rel == 'tech_dup',])],
     'Tech. Replicate', col = cols[6], lwd = 1.3, cex = 1, xpd = NA, adj = 0)

dev.off()


#####################################################################################################


###############
## Distance  ##
###############

out_prop_inds = out2 %>% group_by(pair, rel) %>% summarise(prop = mean(is_shared)) %>% as.data.frame()
out_prop_inds$s1 = sapply(out_prop_inds$pair, function(x) paste(strsplit(x, '_')[[1]][1:4], collapse = '_'))
out_prop_inds$s2 = sapply(out_prop_inds$pair, function(x) paste(strsplit(x, '_')[[1]][5:8], collapse = '_'))


dists_matrix = as.matrix(cophenetic.phylo(cons_tree))
out_prop_inds$dist = apply(out_prop_inds, 1, function(x) dists_matrix[as.character(x['s1']),as.character(x['s2'])])

# out_prop_inds = merge(x = out_prop_inds, y = out2[,c('pair', 'rel')], all.x = T, all.y = F, by = 'pair')
out_prop_inds$rel = factor(out_prop_inds$rel, levels = c("none", "hosp", "dept", "epi", "dup", "tech_dup"))


# 1. Distances vs relation
m_dists = glm(dist~rel, data = out_prop_inds)

