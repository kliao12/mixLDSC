# mixLDSC
mixLDSC uses a bayesian mixture of cross trait LD score regression to follow up regions of significant local genetic correlation across two complex traits. mixLDSC infers latent groups of genetically correlated and uncorrelated SNP sets to identify SNPs and genes driving signals of local genetic correlation

Input: 
1) GWAS summary statistics across two traits (z scores)
2) LD scores for samples matching ancestry of GWAS
   
Output: 
1) Per-snp posterior probability of belonging to genetically correlated group. 
