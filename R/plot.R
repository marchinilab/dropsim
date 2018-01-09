
#' Calculate gene wise summaries
#'
#' \code{summariseDGE} Calculate gene wise summaries from a Digital Gene Expression matrix
#'
#' @param data matrix; A digital gene expression matrix, e.g. output from simulateDGE()
#'
#' @param name character; A descriptive name for the dataset
#'
#' @param rm.0 logical; if TRUE genes with a mean of 0 will be removed
#'
#' @return A data.table of summary statistics of each gene
#'
#' @export
#' @import data.table

summariseDGE <- function(data, name=NA, rm.0=TRUE){

  summary <- data.table(Mean = rowMeans(data),
                        Variance = apply(data, 1, var),
                        Dropout = apply(data, 1, function(x){sum(x==0)}) / ncol(data),
                        Name = name,
                        Gene = rownames(data)
                        )

  summary$Dispersion <- (summary$Variance - summary$Mean) / summary$Mean^2

  if(rm.0){
    summary <- summary[Mean!=0]
  }

  return(summary)

}

# #number of non-zero genes
# print(data[, .(non_zero_genes = .N), by = type][order(-non_zero_genes)])

#' Plot Mean Histogram
#'
#' \code{histogramDGE} Plot a histogram of the gene means
#'
#'
#' @param data data.table; the output of summariseDGE
#'
#' @return A ggplot2 object
#'
#' @examples
#' # see vignette
#'
#' @export
#' @import ggplot2
histogramDGE <- function(data){
  ggplot(data, aes(Mean)) +
    geom_histogram(bins=300) +
    scale_x_log10() +
    ggtitle("Genewise Mean Histogram")
}

#' Plot Mean Variance Relationship
#'
#' \code{histogramDGE} Scatter plot of Variance by Mean
#'
#' @param data data.table; the output of summariseDGE
#'
#' @return A ggplot2 object
#'
#' @examples
#' # see vignette
#'
#' @export
#' @import ggplot2
mean_varianceDGE <- function(data){
  ggplot(data, aes(Mean, Variance)) +
    geom_point(alpha=0.2, size=0.5) +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(intercept=0,slope=1, col="red") +
    ggtitle("Mean-Variance Relationship")
}

#' Plot Normalised Variance
#'
#' \code{histogramDGE} Scatter plot of Variance/Mean by Mean
#'
#' @param data data.table; the output of summariseDGE
#'
#' @return A ggplot2 object
#'
#' @examples
#' # see vignette
#'
#' @export
#' @import ggplot2
normalised_varianceDGE <- function(data){
  ggplot(data, aes(Mean, Variance/Mean)) +
    geom_point(alpha=0.2, size=0.5) +
    scale_x_log10() +
    scale_y_log10() +
    ggtitle("Mean-Dispersion Relationship")
}


#' Plot Dispersion
#'
#' \code{dispersionDGE} Scatter plot of dispersion by gene mean
#'
#' @param data data.table; the output of summariseDGE
#'
#' @return A ggplot2 object
#'
#' @examples
#' # see vignette
#'
#' @export
#' @import ggplot2
  # plot dispersion
dispersionDGE <- function(data, ...){
  ggplot(data, aes(Mean, Dispersion)) +
    geom_point(alpha=0.2, size=0.5) +
    scale_x_log10() +
    scale_y_log10() +
    ggtitle("Mean-Dispersion Relationship")
}


#' Plot Dropout Rate
#'
#' \code{dropoutsDGE} Plot dropout rate by gene mean
#'
#' @param data data.table; the output of summariseDGE
#'
#' @param K numeric; Km parameter of the michaelis menten equation shown in red in the output
#'
#' @return A ggplot2 object
#'
#' @examples
#' # see vignette
#'
#' @export
#' @import ggplot2
dropoutsDGE <- function(data, K=1.056){
  mean_temp <- seq(from=0, to=13, by=0.01)
  mean_var_mm <- data.table(Mean = mean_temp, Dropout = 1 - (mean_temp / (K + mean_temp)) )
  #mean_var_zifa <- data.table(Mean = mean_temp, Dropout = exp(-5 * mean_temp^2) )

  ggplot(data, aes(Mean, Dropout)) +
    geom_point(alpha = 0.1, size=0.5) +
    scale_x_log10() +
    ggtitle("Mean-Dropout Relationship") +
    geom_line(data = mean_var_mm, col="red")
}

#' Plot Summary Grid
#'
#' \code{plot_summaryDGE} Plot multiple gene level summaries
#'
#' @param data data.table; the output of summariseDGE
#'
#' @return A ggplot2 object
#'
#' @examples
#' # see vignette
#'
#' @export
#' @importFrom cowplot plot_grid
plot_summaryDGE <- function(data){
plot_grid(histogramDGE(data),
          mean_varianceDGE(data),
          normalised_varianceDGE(data),
          dispersionDGE(data),
          dropoutsDGE(data), ncol = 3)
}

#' Compare Distributions
#'
#' \code{compare_distributions} Compare an empirical data distribution to fitted distributions using pdf, cdf, qq, and pp plots.
#'
#' @param data data.frame or matrix; The first column is the empirical data, the others being the data from the theoretical distribution e.g. gene means or library sizes
#'
#' @param quantity string; Name of thing being compared - used in the axis labels
#'
#' @param limits; passed to scale_x_log10 "A numeric vector of length two providing limits of the scale.". Defaults to empirical range of the data.
#'
#' @return A ggplot2 object of pdf, cdf, qq, and pp plots.
#'
#' @examples
#' fake_data <- data.frame(data1 = rlnorm(n = 500, meanlog = 0, sdlog = 1), # would usually empirical data
#'                         data2 = rlnorm(n = 500, meanlog = 0, sdlog = 1.5),
#'                         data3 = rlnorm(n = 500, meanlog = 0.5, sdlog = 1))
#' compare_distributions(fake_data)
#'
#' @export
#' @importFrom mltools empirical_cdf
#' @importFrom cowplot plot_grid
compare_distributions <- function(data, limits=NULL, quantity="Data", title=NULL){

  comp <- data

  comp <- data.table(apply(comp,2,sort,decreasing=F))

  empirical_variable <- names(comp)[1]

  colours <- c("black",brewer.pal(9, "Set1")[1:length(names(comp))-1])

  names(colours) <- names(comp)

  if (is.null(limits)){
    empirical_range <- range(comp[,1])
  }else{
    empirical_range <- limits
  }

  library(mltools)
  comp_tidy <- melt(comp, id.vars = empirical_variable, variable.name = "Model")
  ecdf_mltools <- apply(comp, 2, function(x) mltools::empirical_cdf(x, ubounds=comp[,1][[1]]))

  ecdf_ml <- data.table(upperbound = ecdf_mltools[[1]]$UpperBound)
  for (i in 1:length(ecdf_mltools)){
    ecdf_ml <- cbind(ecdf_ml, ecdf_mltools[[i]]$CDF)
  }
  setnames(ecdf_ml, c("upperbound", names(ecdf_mltools)))

  ecdf_melt <- melt(ecdf_ml, id.vars = names(ecdf_ml)[1:2], variable.name = "Model")

  # Cumulative distribution
  pp <- ggplot(ecdf_melt, aes(value, get(empirical_variable), colour=Model)) +
    geom_line(size=1) +
    geom_abline(slope=1, intercept=0) +
    scale_color_manual(values=colours) +
    theme(legend.position="bottom") +
    labs(title="PP Plot", x="Theoretical Probabilities", y="True Probabilities")

  # ggplot(ecdf_melt, aes(x=upperbound)) +
  #   geom_line(aes(y=value, colour=Model), size=1) +
  #   geom_line(aes(y=empirical, colour="Empirical"), size=1) +
  #   xlim(0,1) +
  #   scale_color_manual(values=c("black", brewer.pal(4, "Set1")[c(3,4,2,1)])) +
  #   theme(legend.position="bottom") +
  #   labs(title="CDF", x="Mean", y="CDF")

  cdf <- ggplot(ecdf_melt, aes(x=upperbound)) +
    geom_line(aes(y=value, colour=Model), size=1) +
    geom_line(aes(y=get(empirical_variable), colour=empirical_variable), size=1) +
    scale_x_log10(limits=empirical_range) +
    scale_color_manual(values=colours) +
    theme(legend.position="bottom") +
    labs(title="CDF (Log)", x=quantity, y="CDF")
  #c(exp(-8),exp(5))

  # Probability Density
  # Histogram
  pdf <- ggplot(comp_tidy, aes(get(empirical_variable))) +
    geom_histogram(aes(y=..density..), bins=100) +
    geom_density(aes(value, colour=Model), size=1) +
    scale_x_log10(limits=empirical_range) +
    scale_color_manual(values=colours) +
    theme(legend.position="bottom") +
    labs(title="PDF (Log)", x=quantity, y="Density")
  #c(exp(-8),exp(6))

  # QQ Log
  qq <- ggplot(comp_tidy, aes(value, get(empirical_variable), colour=Model)) +
    geom_line(size=1) +
    geom_abline(slope=1, intercept=0) +
    scale_x_log10(limits=empirical_range) + scale_y_log10() +
    scale_color_manual(values=colours) +
    theme(legend.position="bottom") +
    labs(title="QQ (Log)", x=paste0("Fitted ",quantity), y=paste0("True ",quantity))
  #c(exp(-10),exp(6))

  # QQ
  # ggplot(comp_tidy, aes(value, Empirical_data, colour=Model)) +
  #   geom_point() +
  #   geom_abline(slope=1, intercept=0) +
  #   ylim(0,200) + xlim(0,200) +
  #   scale_color_brewer(palette="Set1") +
  #   theme(legend.position="bottom") +
  #   labs(title="QQ Plot", x="Fitted Means", y="True Means")

  p <- plot_grid(pdf, qq, cdf, pp,
            ncol = 2)

  if(is.null(title)){
    return(p)
  }else{
    title <- ggdraw() + draw_label(title, fontface='bold')
    plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
  }
}
