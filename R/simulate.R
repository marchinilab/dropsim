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
#'
#' @export

dropsim_parameters <- setClass("dropsim_parameters",

                                 slots = list(n_genes = "integer",
                                              n_cells = "integer",
                                              gene_meanlog = "numeric",
                                              gene_sdlog = "numeric",
                                              library_meanlog = "numeric",
                                              library_sdlog = "numeric",
                                              groups = "data.table"),

                                 prototype = list(n_genes = 10000L,
                                                  n_cells = 500L,
                                                  gene_meanlog = -12.5,
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
#' @param parameters object of class dropsim containing parameters
#'
#' @param sparse logical; if true the counts are returned as a sparse matrix
#'
#' @param dge logical; if false counts a returned as a matrix cells by genes rather than genes by cells
#' (If sparse=FALSE, setting dge=FALSE can save time on larger datasets.)
#'
#' @param cell_prefix character; the default group name for cells
#'
#' @param seed integer; seed for random number generation, set this for repoduciable simulations.
#'
#' @return A list containing a matrix of read counts, true underlying expression values,
#' fold change in expression for the groups simulated, and the seed used for random number generation.
#'
#' @examples
#' # new_parameters <- dropsim_parameters()
#' # dge <- simulate(new_parameters)
#' @export
#' @import data.table Matrix

simulateDGE <- function(parameters, sparse=TRUE, cell_prefix = "cell", dge=TRUE, seed=NULL){

  if (class(parameters)!="dropsim_parameters") stop("parameters should be a 'dropsim_parameters' object")

  if (is.null(seed)){
    seed <- sample.int(2^20, 1)
  }
  set.seed(seed)

  true_expression_vector <- signif(rlnorm(n = parameters@n_genes,
                                   meanlog = parameters@gene_meanlog,
                                   sdlog = parameters@gene_sdlog)) # rounding for cross-OS reproducibility

  true_expression <- matrix(rep(true_expression_vector,
                                times = parameters@n_cells),
                            ncol = parameters@n_cells)

  library_sizes <- as.integer(rlnorm(n = parameters@n_cells,
                          meanlog = parameters@library_meanlog,
                          sdlog = parameters@library_sdlog))

  colnames(true_expression) <- rep(cell_prefix, parameters@n_cells)

  fold_change_matrix <- matrix(nrow=parameters@n_genes, ncol=length(parameters@groups$scale))
  colnames(fold_change_matrix) <- parameters@groups$names
  names(dimnames(fold_change_matrix)) <- c("gene", "group")

  for (group in seq_along(parameters@groups$scale)){

    cell_indexes <- parameters@groups$cells[[group]]

    # simulate fold changes for the groupings
    fold_change_matrix[,group] <- signif(exp(rlogis(n = parameters@n_genes,
                                             location = 0.0,
                                             scale = parameters@groups$scale[group]))) # rounding for cross-OS reproducibility

    # apply differential expression
    true_expression[, cell_indexes] <- true_expression[, cell_indexes] * (matrix(rep(fold_change_matrix[,group], times = length(cell_indexes)),
                                                                     ncol = length(cell_indexes)))

    colnames(true_expression)[cell_indexes] <- paste0(colnames(true_expression)[cell_indexes],
                                                 parameters@groups$names[group])

  }

  # Convert means to counts

  true_expression_t <- t(true_expression)

  corrected_means <- true_expression_t * (library_sizes / rowSums(true_expression_t))

  counts <- Matrix(rpois(parameters@n_cells * parameters@n_genes, lambda = corrected_means),
                   nrow = parameters@n_cells, ncol = parameters@n_genes, sparse = sparse)

  # add labels to counts
  colnames(counts) <- 1:parameters@n_genes
  rownames(counts) <- rownames(true_expression_t)
  names(dimnames(counts)) <- c("cell","gene")

  names(true_expression_vector) <- 1:parameters@n_genes
  rownames(fold_change_matrix) <- 1:parameters@n_genes

  if (dge == TRUE){
    counts <- t(counts)
  }

  return(list(counts = counts,
              true_expression = true_expression_vector,
              fold_change = fold_change_matrix,
              library_sizes = library_sizes,
              seed = seed))
}
