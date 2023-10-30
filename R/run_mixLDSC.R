#' Get prior means and variances using empirical bayes approach
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param inFile Path to the input file
#' @return A dataframe of summary statistics with columns 1) z1_z2  2)ld 3) weight corresponding to product of z-scores, ld scores, and weights for each SNP
#' @export
run_mixLDSC <- function(df, p_list, N1, N2){
  new_df = df 
  
  #Define jags models
  pos_jags_model3 = paste0(" 
  model {
    #likelihood
    for (i in 1:n) {
      #y[i] ~ dnorm(beta0[zeta[i]] + beta1_sorted[zeta[i]]*x[i], tau[zeta[i]]) 
      y[i] ~ dnorm(beta0[zeta[i]] + beta1_sorted[zeta[i]]*x[i], tau[zeta[i]]/weight[i])
      #y[i] ~ dnorm(beta1_sorted[zeta[i]]*x[i], tau[zeta[i]]*weight[i])
      zeta[i] ~ dcat(pi[])
    }
  
    #priors
    beta0[1] ~ dnorm(", p_list$b0_est, ", ", sqrt(N1*N2)*p_list$b0_se_prec, ")
    beta0[2] ~ dnorm(", p_list$b0_est, ", ", sqrt(N1*N2)*p_list$b0_se_prec, ")
    
    beta1[1] ~ dnorm(0, 100)
    beta1[2] ~ dnorm(", p_list$avg_b1, ", ", sqrt(N1*N2)*p_list$b1_se_prec, ") 
    
    beta1_sorted[1:2] <- sort(beta1) 
    
    for (i in 1:H) {
      tau[i] ~ dgamma(1,1) 
      sigma[i] <- 1/sqrt(tau[i])
    }
  
    pi ~ ddirich(a)
    #sorted_pi[1:2] <- sort(pi)
  }
  ")
  
  neg_jags_model3 = paste0(" 
  model {
    #likelihood
    for (i in 1:n) {
      #y[i] ~ dnorm(beta0[zeta[i]] + beta1_sorted[zeta[i]]*x[i], tau[zeta[i]]) 
      y[i] ~ dnorm(beta0[zeta[i]] + beta1_sorted[zeta[i]]*x[i], tau[zeta[i]]/weight[i])
      #y[i] ~ dnorm(beta1_sorted[zeta[i]]*x[i], tau[zeta[i]]*weight[i])
      zeta[i] ~ dcat(pi[])
    }
  
    #priors
    beta0[1] ~ dnorm(", p_list$b0_est, ", ", sqrt(N1*N2)*p_list$b0_se_prec, ")
    beta0[2] ~ dnorm(", p_list$b0_est, ", ", sqrt(N1*N2)*p_list$b0_se_prec, ")
    
    beta1[1] ~ dnorm(0, 100)
    beta1[2] ~ dnorm(", p_list$b1_est, ", ", sqrt(N1*N2)*p_list$b1_se_prec, ") 
    
    beta1_sorted[1:2] <- beta1 
    
    for (i in 1:H) {
      tau[i] ~ dgamma(1,1) 
      sigma[i] <- 1/sqrt(tau[i])
    }
  
    pi ~ ddirich(a)
    #sorted_pi[1:2] <- sort(pi)
  }
  ")
  
  ##### Run jags
  dat = list(n=nrow(new_df), H=2, y=new_df$z1_z2, x=new_df$ld, a=c(99,1), weight=new_df$weight)
  
  model.inits <- list(beta0=c(p_list$b0_est,p_list$b0_est), beta1=c(0,p_list$avg_b1), pi=c(0.99, 0.01))
  #model.inits <- list(beta1=c(0,orig_beta), pi=c(0.7, 0.3))
  
  if(p_list$avg_b1 >= 0){
    print(paste0("Empirical bayes b1 est: ", p_list$avg_b1))
    print("using positive jags model")
    jm = jags.model(textConnection(pos_jags_model3), data = dat, n.chains = 1, inits = model.inits) 
  } else {
    print(paste0("Empirical bayes b1 est: ", p_list$avg_b1))
    print("using negative jags model")
    jm = jags.model(textConnection(neg_jags_model3), data = dat, n.chains = 1, inits = model.inits) 
  }
  
  update(jm, n.iter = 1000)
  r = coda.samples(jm, c('beta0','beta1', 'sigma','pi', 'zeta'), n.iter=10000, thin = 2)
  r_sum <- summary(r)
  
  beta0_1 <- r_sum$quantiles[1,3]; beta0_2 <- r_sum$quantiles[2,3]
  beta1_1 <- r_sum$quantiles[3,3]; beta1_2 <- r_sum$quantiles[4,3]

  pi_1 <- r_sum$quantiles[5,3]; pi_2 <- r_sum$quantiles[6,3]

  #Get counts of class assignments. Use most common across posterior samples
  col_end <- dim(r[[1]])[2]
  zetas <- r[[1]][, 9:col_end]

  res <- apply(zetas,2,function(x) names(which.max(table(x))))
  M1 = table(res)[1]; M2 = table(res)[2]

  ### Compute genetic corrs by inferred class
  class_vec <- res
  class_vec[class_vec == '1'] <- 'neutral'; class_vec[class_vec == '2'] <- 'pos'
  
  new_df$inf_class <- class_vec
  
  ### Look at posterior probability 
  col_end <- dim(r[[1]])[2]
  zetas <- r[[1]][, 9:col_end]
  zetas_minus_one <- zetas-1
  
  new_df$post_prob_pos <- apply(zetas_minus_one, 2, mean)
  output_list <- list("new_df" = new_df, "beta0_1" = beta0_1, 
                      "beta0_2" = beta0_2, "beta1_1" = beta1_1, 
                      "beta1_2" = beta1_2)
  return(output_list) 
}
