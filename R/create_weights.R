#' Create weights for LDSC
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param inFile Path to the input file
#' @return A dataframe of summary statistics with columns chr, pos, a0, a1, beta
#' @export
create_weights <- function(df){
  new_df = df
  M <- nrow(new_df) #number of sites 
  N1 = N2 = Ns = 25000 #number of samples for traits 1 and 2
  
  rho_g = sum(new_df$z1_z2)/sqrt(N1*N2)
  new_df$cond_var_base <- ((N1*new_df$ld)/M + 1)*((N2*new_df$ld)/M + 1) + 
    ((sqrt(N1*N2)*rho_g*new_df$ld)/M)^2
  
  base_mod <- summary(lm(data = new_df, z1_z2 ~ ld, weights=1/cond_var_base))
  int <- base_mod$coefficients[1,1]; slope <- base_mod$coefficients[2,1]
  rho_g_hat <- slope*M/sqrt(N1*N2); rho_Ns_hat <- int*sqrt(N1*N2)
  
  #Use estimated rho_g_hat and rho_Ns_hat in final weights 
  new_df$cond_var_final <- ((N1*new_df$ld)/M + 1)*((N2*new_df$ld)/M + 1) + 
    ((sqrt(N1*N2)*rho_g_hat*new_df$ld)/M + rho_Ns_hat/sqrt(N1*N2))^2
  
  #clumped_data$weight <- clumped_data$cond_var_final*(1/clumped_data$ld)
  new_df$weight <- 1/new_df$cond_var_final*(1/new_df$ld)
  return(new_df)
}
