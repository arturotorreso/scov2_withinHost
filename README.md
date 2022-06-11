Code for Within-host diversity improves phylogenetic and transmission reconstruction of SARS-CoV-2 outbreaks
https://www.biorxiv.org/content/10.1101/2022.06.07.495142v1 (not peer-reviewed)

1. Model for within-host diversity (1_simulations)
> Sequence simulations using a model with within-host diversity: simulate.R <br/>
> Tree comparison: compareTrees.R<br/>

2. Phylogenetic inference (2_genomics)
> Maximum likelihood tree inference: raxml_ng.sh <br/>

3. Whole-genome sequence analysis (2_genomics)
> Raw reads mapping: align.sh <br/>
> Variant calling: varCalling.sh, softFilt_vcf_by_bed.py, addFT_pileup.py <br/>
> Pseudosequence: pseudoseq_pileup.sh, vcf2pseudoseq.py, add_maj_min.py <br/>

4. Analysis (3_analysis)
> Comparison of technical replicates: 1_techReps.R <br/>
> Pairwise comparison of samples: 2_epiLinks.R <br/>
> Temporal signal in clusters: 3_clust_signal.R <br/>
> Transmission inference: transPairs.R, likelihoodSEIR.R, transPairs_heatmap.R, transPairs_network.R <br/>
