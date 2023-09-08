#' Load summary statistics 
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param inFile Path to the input file
#' @return A dataframe of summary statistics with columns chr, pos, a0, a1, beta
#' @export
load_sumStats <- function(inFile){
  sumStats <- bigreadr::fread2(inFile)
  colnames(sumStats) <- c('chr', 'pos',"a0", "a1", 'beta')
  return(sumStats)
}

#' Load plink file using bigSNP package  
#'
#' This function loads a plink file into R. Uses bigSNP package to first check if .rds file exists. If not will create from .bed file and then import 
#'
#' @param inFile Path to plink bed file including .bed extension 
#' @return A matrix of the infile
#' @export
load_plinkFile <- function(inFile){
  rds_path <- paste0(substr(inFile, 1, nchar(inFile) - 4), ".rds")
  if(file.exists(rds_path) == FALSE){
    snp_readBed(inFile)
  }
  obj.bigSNP <- snp_attach(rds_path)
  return(obj.bigSNP)
}