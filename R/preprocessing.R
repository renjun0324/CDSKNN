#' Marchenko-Pastur Significant PCs (URD)
#'
#' This function comes from URD, please check [github] (https://github.com/farrellja/URD) for more function details
#'
#' @param M (Numeric) Number of rows in input data, pc.gene
#' @param N (Numeric) Number of columns in input data, cell number
#' @param pca.sdev (Numeric vector) Standard deviations for each principal component
#' @param factor (Numeric) Factor to multiply eigenvalue null upper bound before determining significance.
#' @param do.print (Logical) Whether to report the results
#' @return Logical vector of whether each PC is significant.
#'
#' @export
#'
pcaMarchenkoPastur <- function(M, N, pca.sdev, factor = 2, do.print=T) {
  pca.eigenvalue <- (pca.sdev)^2
  marchenko.pastur.max <- (1+sqrt(M/N))^2
  pca.sig <- pca.eigenvalue > (marchenko.pastur.max * factor)
  if (do.print) {
    print(paste("Marchenko-Pastur eigenvalue null upper bound:", marchenko.pastur.max))
    if (factor != 1) {
      print(paste(length(which(pca.sig)), "PCs have eigenvalues larger than", factor, "times null upper bound."))
    } else {
      print(paste(length(which(pca.eigenvalue > marchenko.pastur.max)), "PCs have larger eigenvalues."))
    }
  }
  return(pca.sig)
}
