
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
