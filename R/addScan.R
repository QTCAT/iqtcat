#' @title Additive ANOVA genome scan
#' @description Additive ANOVA genome scan, scanning each marker.
#' 
#' @param pheno An object of class \code{\link{qtcatPheno}}.
#' @param snp An object of S4 class \linkS4class{snpMatrix}.
#' @param mc.cores Number of cores for parallelising. Theoretical maximum is
#' \code{'B'}. For details see \code{\link[parallel]{mclapply}}.
#' 
#' @importFrom parallel mclapply
#' @export
addScan <- function(pheno, snp, mc.cores = 1L) {
  stopifnot(is(pheno, "qtcatPheno"))
  stopifnot(is(snp, "snpMatrix"))
  id <- intersect(pheno$names, rownames(snp))
  if (!length(id))
    stop("The ID intersect of 'pheno' and 'snp' is emty")
  if (length(id.uniqueGeno <- setdiff(rownames(snp), id)))
    cat("The following individuals are part of 'snp' but not of 'pheno':\n",
        paste(id.uniqueGeno, collapse = " "), "\n")
  if (length(id.uniquePheno <- setdiff(pheno$names, id)))
    cat("The following individuals are part of 'pheno' but not of 'snp':\n",
        paste(id.uniquePheno, collapse = " "), "\n")
  phenoInx <- which(pheno$names %in% id)
  genoInx <- match(pheno$names[phenoInx], rownames(snp))
  if (ncol(pheno$covariates) == 0L) {
    snp <- snp[genoInx, ]
    x.des <- NULL
  } else {
    snp <- snp[genoInx, ]
    x.des <- pheno$covariates[phenoInx, ]
  }
  y <- pheno$pheno[phenoInx]
  family <- pheno$family
  pval <- unlist(mclapply(1:ncol(snp), singleTests, y, snp, x.des, family,
                          mc.cores = mc.cores, mc.preschedule = TRUE,
                          mc.cleanup = TRUE))
  out <- cbind(snpInfo(snp)[, 1L:2L], pValues = pval)
  class(out) <- c("addScan", class(out))
  out
} 


#' @title Single marker ANOVA
#' @description Internal singele marker ANOVA
#' 
#' @param i ...
#' @param y ...
#' @param snp An object of S4 class \linkS4class{snpMatrix}.
#' @param x.des ...
#' @param family ...
#' @param test ...
#'
#' @importFrom hit fast.anova
#' @keywords internal
singleTests <- function(i, y, snp, x.des, family, test = c("LRT", "F")) {
  if (is.null(x.des)) {
    x <- cbind(rep(1, nrow(snp)), as.matrix(snp[, i]))
    assign <- 0L:1L
    p <- fast.anova(x, y, assign, family, test)
  } else {
    x <- cbind(rep(1, nrow(snp)), x.des, as.matrix(snp[, i]))
    assign <- c(0L, rep(1L, ncol(x.des)), 2L)
    p <- fast.anova(x, y, assign, family, test)[2L]
  }
  p
}


#' @title Plot of significance of marker-phenotype association.
#'
#' @description Plot of significance of marker-phenotype association at their
#' genome positions.
#'
#' @param x A data.frame object of class \code{\link{addScan}}.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param col.axis Colors for axis line, tick marks, and title respectively.
#' @param ... Other graphical parameters may also be passed as arguments to this function.
#'
#'
#' @importFrom graphics plot.default axis mtext
#' @export
plot.addScan <- function(x, xlab = "Chromosomes", 
                              ylab = expression(-log[10](italic(p))), 
                              col.axis = NULL, ...) {
  # make positions linear with gaps between chr's
  stopifnot(colnames(x) == c("chr", "pos", "pValues"))
  x[, 3L] <- -log10(x[, 3L])
  chrminmax <- vapply(split(x[, 2L], x[, 1L]), function(x) c(min(x), max(x)), c(1, 2))
  chrgap <- round(sum(chrminmax[2L, ] - chrminmax[1L, ]) * .01, 0)
  chrsize <- cumsum(c(0, chrminmax[2L, -ncol(chrminmax)])) -
    cumsum(chrminmax[1L, ]) +
    cumsum(c(0, rep(chrgap, ncol(chrminmax) - 1L)))
  chrstartend <- chrminmax + rbind(chrsize, chrsize) + chrgap / 2
  chr <- sort(unique(x[, 1L]))
  for (i in seq_along(chr)) {
    inx <- which(x[, 1L] == chr[i])
    if (length(inx))
      x[inx, 2L] <- x[inx, 2L] + chrsize[i] + chrgap / 2
  }
  xlim <- c(chrstartend[1L, 1L] - chrgap, chrstartend[2L, ncol(chrstartend)] + chrgap)
  # plot
  plot.default(x[, c(2L, 3L)], xlim = xlim, axes = FALSE, xlab = "", ylab = "", ...)
  # x
  for (i in seq_along(chr))
    axis(1, labels = FALSE, at = c(chrstartend[1, i], chrstartend[2, i]), col = col.axis)
  axis(1, at = colMeans(chrstartend), labels = chr, col = NA, col.axis = col.axis)
  mtext(xlab, 1, 2.5, col = col.axis)
  # y
  axis(2, col = col.axis, col.axis = col.axis)
  mtext(ylab, 2, 2.5, col = col.axis)
}
