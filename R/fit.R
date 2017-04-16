#' Fit parameters to a dataset
#'
#' \code{fit_parameters} find parameters from data to simulate new dataset
#'
#'
#' @param data matrix; gene by cell dge matrix of counts
#' 
#' @param plot logical; if TRUE the real vs theoretical distribution plots qqplot etc. are printed
#' 
#' @return An object of class sinsynthr_parameters containing fitted parameters for the mean and library size lognormal distributions
#'
#' @examples
#' # new_parameters <- fit_parameters(your_data_matrix)
#'
#' @export
#' @import data.table fitdistrplus

fit_parameters <- function(data, plot=TRUE, method="mge", ...){
  # Dropseq normalise to estimate mean dist params
  lib.sizes <- colSums(data)
  norm.counts <- sweep(data, 2, lib.sizes / median(sqrt(lib.sizes)), "/")
  norm.counts <- norm.counts[rowSums(norm.counts > 0) >= 1, ]
  
  # Fit gene mean parameters
  mean_var_n <- summariseDGE(norm.counts)
  
  lnorm_parameters <- fitdistrplus::fitdist(mean_var_n$Mean, "lnorm", method = "mge", gof="CvM")
  
  plot(lnorm_parameters)
  mtext("Gene Mean fits", outer = TRUE, cex = 1.5)
  
  # Fit cell library parameters
  library_parameters <- fitdistrplus::fitdist(lib.sizes, "lnorm", method = "mge", gof="CvM")
  
  plot(library_parameters)
  mtext("Library Size fits", outer = TRUE, cex = 1.5)
  
  new_parameters <- new("sinsynthr_parameters",
                        n_genes = nrow(data),
                        n_cells = ncol(data),
                        gene_meanlog = lnorm_parameters$estimate["meanlog"],
                        gene_sdlog = lnorm_parameters$estimate["sdlog"],
                        library_meanlog = library_parameters$estimate["meanlog"],
                        library_sdlog = library_parameters$estimate["sdlog"],
                        groups = data.table(
                          scale = numeric(),
                          cells = list(),
                          names = character()
                        )
  )
  
  return(new_parameters)
  
}

