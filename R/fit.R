#' Fit parameters to a dataset
#'
#' \code{fit_parameters} find parameters from data to simulate new dataset
#'
#'
#' @param data matrix; gene by cell dge matrix of counts
#'
#' @param plot logical; if TRUE the real vs theoretical distribution plots qqplot etc. are printed
#'
#' @return An object of class dropsim_parameters containing fitted parameters for the mean and library size lognormal distributions
#'
#' @examples
#' # new_parameters <- fit_parameters(your_data_matrix)
#'
#' @export
#' @import data.table fitdistrplus

fit_parameters <- function(data, plot=TRUE, method="mge", ...){
  # Dropseq normalise to estimate mean dist params
  library_sizes <- colSums(data)
  normalised_counts <- sweep(data, 2, library_sizes , "/") #/ median(sqrt(library_sizes))
  normalised_counts <- normalised_counts[rowSums(normalised_counts > 0) >= 1, ]

  # Fit gene mean parameters
  summarised_counts <- summariseDGE(normalised_counts)

  lnorm_parameters <- fitdistrplus::fitdist(summarised_counts$Mean, "lnorm", method = "mge", gof="CvM")

  gene_means <- data.table(
    Empirical_data = summarised_counts$Mean,
    lognormal = rlnorm(n=nrow(summarised_counts),
                       meanlog=lnorm_parameters$estimate["meanlog"],
                       sdlog=lnorm_parameters$estimate["sdlog"])
  )

  print(compare_distributions(gene_means, title="Gene Means"))

  # Fit cell library parameters
  library_parameters <- fitdistrplus::fitdist(library_sizes, "lnorm", method = "mge", gof="CvM")

  library_dists <- data.table(
    Empirical_data = library_sizes,
    lognormal = rlnorm(n=length(library_sizes),
                       meanlog=library_parameters$estimate["meanlog"],
                       sdlog=library_parameters$estimate["sdlog"])
  )

  print(compare_distributions(library_dists, title="Library Sizes"))

  new_parameters <- new("dropsim_parameters",
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

