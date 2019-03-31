#' @title Index of inteacting loci
#' @description Index of inteacting loci
#' 
#' @param snp An object of S4 class \linkS4class{snpMatrix}. 
#' @param minDist minimal physical distence
#' 
#' @importFrom qtcat snpInfo
#' @export
iqtcatInx <- function(snp, minDist = 10000) {
  stopifnot(is(snp, "snpMatrix"))
  out <- interactionIndices(snp@snpData, 
                            as.character(snpInfo(snp)[, 1L]), 
                            snpInfo(snp)[, 2L], 
                            minDist)
  class(out) <- c("iqtcatInx", class(out))
  out
}


#' @title Test index of inteacting markers
#' @description Test if inteacting markers have all four homozygous combination, 
#' if not, they are removed.
#' 
#' @param snp An object of S4 class \linkS4class{snpMatrix}.
#' @param iInx An object of class \code{\link{iqtcatInx}}.
#' @param minDist minimal physical distence
#' 
#' @importFrom qtcat snpInfo
#' @export
test.iInx <- function(snp, iInx, minDist = 10000) {
  stopifnot(is(snp, "snpMatrix"))
  if (missing(iInx))
    stop("'iInx' must be specifid")
  if (any(colnames(iInx) != c("inx1", "inx2")))
    stop("Column of 'iInx' have to be named 'inx1' and 'inx2'")
  out <- testIndices(snp@snpData, iInx[, 1L], iInx[, 2L], 
                     as.character(snpInfo(snp)[, 1L]), 
                     snpInfo(snp)[, 2L], 
                     minDist)
  class(out) <- c("iqtcatInx", class(out))
  out
}
