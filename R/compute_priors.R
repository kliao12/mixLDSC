#' Get prior means and variances using empirical bayes approach
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param inFile Path to the input file
#' @return A dataframe of summary statistics with columns chr, pos, a0, a1, beta
#' @export
compute_priors <- function(df, bigSNP_file){
  #1) run original cross trait ldsc with all data
  base_mod <- summary(lm(data = df, z1_z2 ~ ld, weights=weight))
  b0_se_prec <- 1/base_mod$coefficients[1,2]; b1_se_prec <- 1/base_mod$coefficients[2,2] 
  b0_est <- base_mod$coefficients[1,1]; b1_est <- base_mod$coefficients[2,1]
  
  #2) run cross trait ldsc with clumped data
  if(b1_est > 0){
    clumped_snps <- snp_clumping(bigSNP_file$genotypes, infos.chr = rep(1, length(bigSNP_file$map$physical.pos)), 
                                 S = df$z1_z2,
                                 infos.pos = bigSNP_file$map$physical.pos)
  } else {
    clumped_snps <- snp_clumping(bigSNP_file$genotypes, infos.chr = rep(1, length(bigSNP_file$map$physical.pos)), 
                                 S = abs(df$z1_z2),
                                 infos.pos = bigSNP_file$map$physical.pos)
  }
  
  ### Get priors values via empirical bayes approach
  base_mod_clumped <- summary(lm(data = df[clumped_snps, ], z1_z2 ~ ld, weights=weight))
  b1_est_clumped <- base_mod_clumped$coefficients[2,1]
  
  avg_b1 <- (b1_est + b1_est_clumped)/2
  
  prior_list <- list("b0_est" = b0_est, "b0_se_prec" = b0_se_prec, 
                     "b1_est" = b1_est, "b1_se_prec" = b1_se_prec, 
                     "avg_b1" = avg_b1)
  return(prior_list) 
}
