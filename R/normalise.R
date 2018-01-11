#' Normalise Digital Gene Expression Matrix
#'
#' \code{normaliseDGE} normalise a gene by cell matrix
#'
#'
#' @param dge matrix; a digital gene expression matrix containing count data (columns = cells, rows = genes), can be a sparce Matrix object.
#'
#' @param center logical; if true genes are centered to mean of 0
#'
#' @param scale logical; if TRUE genes are scaled to unit variance
#'
#' @param threshold integer; any normalised values larger than this will be rounded down
#'
#' @param gene_subset real; what fraction of the genes do you want to keep (by mean) e.g. 0.5 is top half most expressed
#'
#' @param min_library_size real; only cells with library size equal or greater to this will be kept
#'
#' @param verbose logical; if TRUE the function will print diagnostic graphs along the way
#'
#' @param outliers logical; calculate max value per gene and cell to assess threshold position
#'
#' @param transformation character; either sqrt (default) or asinh
#'
#' @return A matrix of read counts
#'
#' @examples
#' # normalised_counts <- normalise_dge(dge)
#' @export
#' @import data.table Matrix
#'
normaliseDGE <- function(dge, verbose=FALSE, center=TRUE, scale=TRUE, outliers = FALSE, threshold = 5, min_library_size=200, gene_subset=0.5, transformation="sqrt"){

  library_size <- colSums(dge)
  lib_size_correction_factor <- library_size / median(sqrt(library_size))
  dge <- t(t(dge) / lib_size_correction_factor) # sweep(dge, 2, , "/")

  # subset cells and genes
  cell_subset <- library_size >= min_library_size
  gene_mean <- rowMeans(dge[, cell_subset])
  gene_subset <- data.table(gene_mean, names(gene_mean))[order(-gene_mean)][1:round(length(gene_mean) * gene_subset)]$V2

  dge <- dge[gene_subset, cell_subset]

  # multiply & sqrt
  if (transformation=="asinh"){
  	  dge <- asinh(10000 * dge)
  } else{
  	  dge <- sqrt(10000 * dge)
  }


  # transpose so cells are rows
  dge <- t(dge)

  # check everything looks ok
  if(verbose){
    str(dge)

    if(prod(dim(dge))>1e6){
      tmp <- sample(dge, 1e6)
    } else if(prod(dim(dge))<1e6 & class(dge)=="dgCMatrix"){
      tmp <- as.matrix(dge)
    } else{
      tmp <- dge
    }

    hist(tmp, breaks = 2000, ylim = c(0, 1e3), xlim=c(0,50), xlab="Expression", main="Histogram of data post cell-wise (& sqrt) normalisation")
  }

  # scale to unit variance all columns (genes) - don't center so 0 values remain exactly 0 and all data >= 0
  if (center==TRUE){

    dge <- scale(dge, center = TRUE, scale = TRUE)

  } else{
  colSdColMeans <- function(x, na.rm=TRUE) {
    if (na.rm) {
      n <- colSums(!is.na(x)) # thanks @flodel
    } else {
      n <- nrow(x)
    }
    colVar <- colMeans(x*x, na.rm=na.rm) - (colMeans(x, na.rm=na.rm))^2
    return(sqrt(colVar * n/(n-1)))
  }

  dge <- t(t(dge) / colSdColMeans(dge)) #, na.rm=TRUE

  }

  dge <- dge[, !is.na(colSums(dge))] # remove NA genes (sd calculation fails if all 0s)

  if(verbose){
    if(prod(dim(dge))>1e6){
      tmp <- sample(dge, 1e6)
    }else if(prod(dim(dge))<1e6 & class(dge)=="dgCMatrix"){
      tmp <- as.matrix(dge)
    }else{
      tmp <- dge
    }

    hist(tmp, breaks = 2000, ylim = c(0, 1e4), main="Histogram of data post gene-wise normalisation")
    }

  # make sure no NA values left
  stopifnot(sum(is.na(dge)) == 0)

  # assess degree of outliers
  if(outliers){
    max_by_cell <- apply(dge, 1, max)
    max_by_gene <- apply(dge, 2, max)
    plot(max_by_gene, main="Maximum expression value per gene", ylab="Maximum Expression", xlab="Gene Index (by mean expression)")
    abline(h=threshold, col="red")
    plot(max_by_cell[], main="Maximum expression value per cell", ylab="Maximum Expression", xlab="Cell Index (by library size)")
    abline(h=threshold, col="red")
  }

  # round down large values which could be due to small library size and normalisation OR low gene wise expression and scaling
  sum(dge > threshold)
  dge[dge > threshold] <- threshold
  dge[dge < (-threshold)] <- (-threshold)

  if(verbose){
    if(prod(dim(dge)) > 1e6){
      tmp <- sample(dge, 1e6)
    }else if(prod(dim(dge)) < 1e6 & class(dge)=="dgCMatrix"){
      tmp <- as.matrix(dge)
    }else{
      tmp <- dge
    }

    hist(tmp, breaks = 500, ylim = c(0, 1e3), main="Histogram of final data", xlab="Expression")

    # check everything looks ok
    str(dge)
  }

  return(dge)
}
