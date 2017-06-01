#' An S4 class to hold simulation parameters.
#'
#' @slot n_genes integer; number of genes to simulate
#' 
#' @slot n_cells integer; number of cells to simulate
#' 
#' @slot gene_meanlog numeric; meanlog of lognormal distribution to simulate means from, see ?rlnorm for more details
#' 
#' @slot gene_sdlog numeric; sdlog of lognormal distribution to simulate means from, see ?rlnorm for more details
#' 
#' @slot library_meanlog numeric; meanlog of lognormal distribution to simulate library sizes from, see ?rlnorm for more details
#' 
#' @slot library_sdlog numeric; sdlog of lognormal distribution to simulate library sizes from, see ?rlnorm for more details
#' 
#' @slot groups data.table with variables 'scale', 'cells' and 'names'
#' scale is a numeric vector of the scale parameter for the logistic distribution to simulate the fold changes from
#' cells is a list of integer vector giving the indexes of cells each differential expression should be applied to, cells can be a member of more than one group
#' names is a character vector of cell group names

sinsynthr_parameters <- setClass("sinsynthr_parameters",
                                 
                                 slots = list(n_genes = "integer",
                                              n_cells = "integer",
                                              gene_meanlog = "numeric",
                                              gene_sdlog = "numeric",
                                              library_meanlog = "numeric",
                                              library_sdlog = "numeric",
                                              groups = "data.table"),
                                 
                                 prototype = list(n_genes = 10000L,
                                                  n_cells = 500L,
                                                  gene_meanlog = -2.5,
                                                  gene_sdlog = 2.25,
                                                  library_meanlog = 9.5,
                                                  library_sdlog = 0.35,
                                                  groups = data.table(
                                                    scale = c(0.06),
                                                    cells = list(1:250),
                                                    names = "B"
                                                  )
                                 )
)


#' Simulate Single Cell RNAseq count data
#'
#' \code{simulateDGE} simulate digital gene expression matrix containing count data from given parameters
#'
#'
#' @param parameters object of class sinsynthr containing parameters
#'
#' @param sparse logical; if true the counts are returned as a sparse matrix
#' 
#' @param dge logical; if false counts a returned as a matrix cells by genes rather than genes by cells
#' (If sparse=FALSE, setting dge=FALSE can save time on larger datasets.)
#' 
#' @param cell_prefix character; the default group name for cells
#' 
#' @return A matrix of read counts
#'
#' @examples
#' # new_parameters <- sinsynthr_parameters()
#' # dge <- simulate(new_parameters)
#' @export
#' @import data.table Matrix

simulateDGE <- function(parameters, sparse=TRUE, cell_prefix = "cell", dge=TRUE){
  
  if (class(parameters)!="sinsynthr_parameters") stop("parameters should be a 'sinsynthr_parameters' object")
  
  true_means <- matrix(rep(rlnorm(n = parameters@n_genes,
                                  meanlog = parameters@gene_meanlog,
                                  sdlog = parameters@gene_sdlog),
                           times = parameters@n_cells),
                       ncol = parameters@n_cells)
  
  library_sizes <- rlnorm(n = parameters@n_cells,
                          meanlog = parameters@library_meanlog,
                          sdlog = parameters@library_sdlog)
  
  colnames(true_means) <- rep(cell_prefix, parameters@n_cells)
  
  for (group in seq_along(parameters@groups$scale)){
    
    cell_indexes <- parameters@groups$cells[[group]]
    
    # simulate fold changes for the groupings
    fold_change <- exp(rlogis(n = parameters@n_genes,
                              location = 0.0,
                              scale = parameters@groups$scale[group]))
    
    # apply differential expression
    true_means[, cell_indexes] <- true_means[, cell_indexes] * (matrix(rep(fold_change, times = length(cell_indexes)),
                                                                     ncol = length(cell_indexes)))
    
    colnames(true_means)[cell_indexes] <- paste0(colnames(true_means)[cell_indexes],
                                                 parameters@groups$names[group])
    
  }
  
  # Convert means to counts
  
  true_means_t <- t(true_means)
  
  corrected_means <- true_means_t * (library_sizes / rowSums(true_means_t))
  
  counts <- Matrix(rpois(parameters@n_cells * parameters@n_genes, lambda = corrected_means),
                   nrow = parameters@n_cells, ncol = parameters@n_genes, sparse = sparse)
  
  # add labels to counts
  colnames(counts) <- 1:parameters@n_genes
  rownames(counts) <- rownames(true_means_t)
  
  if (dge == TRUE){  
    counts <- t(counts)
  }
  
  return(counts)
}


#' Normalise Digital Gene Expression Matrix
#'
#' \code{normaliseDGE} simulate digital gene expression matrix containing count data from given parameters
#'
#'
#' @param dge matrix; a digital gene expression matrix containing count data
#'
#' @param center logical; if true genes are centered to mean of 0
#' 
#' @param scale logical; if TRUE genes are scaled to unit variance
#' (If sparse=FALSE, setting dge=FALSE can save time on larger datasets.)
#' 
#' @param threshold integer; any values larger than this will be rounded down
#' 
#' @param verbose logical; if TRUE the function will print diagnostic graphs along the way
#' 
#' @return A matrix of read counts
#'
#' @examples
#' # normalised_counts <- normalise_dge(dge)
#' @export
#' @import data.table Matrix
#' 
normaliseDGE <- function(dge, verbose=FALSE, center=TRUE, scale=TRUE, threshold = 5, min_library_size=1){
  
  library_size <- colSums(dge)
  dge <- sweep(dge, 2, library_size / median(sqrt(library_size)), "/")
  
  # subset cells and genes
  cell_subset <- library_size > min_library_size
  gene_mean <- rowMeans(dge[, cell_subset])
  gene_subset <- data.table(gene_mean, names(gene_mean))[order(-gene_mean)][1:round(length(gene_mean) / 2)]$V2
  
  dge <- dge[gene_subset, cell_subset]
  
  # multiply & sqrt
  dge <- sqrt(10000 * dge)
  
  # transpose so cells are rows
  dge <- t(dge)
  
  # check everything looks ok
  if(verbose){
    str(dge)
    
    hist(dge, breaks = 2000, ylim = c(0, 3e4), xlim=c(0,26), xlab="Expression", main="Histogram of data post cell-wise (& sqrt) normalisation")
  }
  
  # scale to unit variance all columns (genes) - don't center so 0 values remain exactly 0 and all data >= 0
  dge <- scale(dge, center = center, scale = scale)
  dge <- dge[, !is.na(colSums(dge))] # remove NA genes (sd calculation fails if all 0s)
  
  if(verbose){hist(dge, breaks = 2000, ylim = c(0, 3e4), main="Histogram of data post gene-wise normalisation")}
  
  # make sure no NA values left
  stopifnot(sum(is.na(dge)) == 0)
  
  # assess degree of outliers
  if(verbose){
    max_by_cell <- apply(dge, 1, max)
    max_by_gene <- apply(dge, 2, max)
    plot(max_by_gene, main="Maximum expression value per gene", ylab="Maximum Expression", xlab="Gene Index (by mean expression)")
    abline(h=10, col="red")
    plot(max_by_cell[], main="Maximum expression value per cell", ylab="Maximum Expression", xlab="Cell Index (by library size)")
    abline(h=10, col="red")
  }
  
  # round down large values which could be due to small library size and normalisation OR low gene wise expression and scaling
  sum(dge > threshold)
  dge[dge > threshold] <- threshold
  dge[dge < (-threshold)] <- (-threshold)
  
  if(verbose){
    hist(dge, breaks = 2000, ylim = c(0, 3e4), main="Histogram of final data", xlab="Expression")
    
    # check everything looks ok
    str(dge)
  }
  
  return(dge)
}
