#' @title Epistatic ANOVA genome scan
#' @description  Epistatic ANOVA genome scan, scanning all two marker interactions.
#' 
#' @param pheno An object of class \code{\link{qtcatPheno}}.
#' @param snp An object of S4 class \linkS4class{snpMatrix}.
#' @param iInx An object of class \code{\link{iqtcatInx}}.
#' @param mc.cores Number of cores for parallelising. Theoretical maximum is
#' \code{'B'}. For details see \code{\link[parallel]{mclapply}}.
#' 
#' @importFrom parallel mclapply
#' @importFrom hit fast.anova
#' @export
epiScan <- function(pheno, snp, iInx, mc.cores = 1L) {
  stopifnot(is(pheno, "qtcatPheno"))
  stopifnot(is(snp, "snpMatrix"))
  if (missing(iInx))
    stop("'iInx' must be specifid")
  if (any(colnames(iInx) != c("inx1", "inx2")))
    stop("Column of 'iInx' have to be named 'inx1' and 'inx2'")
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
  # Running multiple ANOVAs
  pVal <- mclapply(1:nrow(iInx), 
                   iSinglTests, 
                   y, snp, iInx, x.des, family, test = c("LRT", "F"), 
                   mc.cores = mc.cores)
  out <- data.frame(chr1 = snpInfo(snp)[iInx$inx1, 1L],
                    pos1 = snpInfo(snp)[iInx$inx1, 2L],
                    chr2 = snpInfo(snp)[iInx$inx2, 1L],
                    pos2 = snpInfo(snp)[iInx$inx2, 2L],
                    pValues = unlist(pVal),
                    row.names = paste(colnames(snp)[iInx$inx1], 
                                      colnames(snp)[iInx$inx2], 
                                      sep = ":"),
                    stringsAsFactors = FALSE)
  class(out) <- c("epiScan", class(out))
  out
} 


#' @title A two marker interaction  ANOVA
#' @description Internal two marker interaction  ANOVA
#' 
#' @param i ...
#' @param y ...
#' @param snp An object of S4 class \linkS4class{snpMatrix}.
#' @param iInx An object of class \code{\link{iqtcatInx}}.
#' @param x.des ...
#' 
#' @importFrom hit fast.anova
#' @importFrom qtcat as.matrix
#' @keywords internal
iSinglTests <- function(i, y, snp, iInx, x.des, family, test = c("LRT", "F")) {
  l <- iInx$inx1[i]
  k <- iInx$inx2[i]
  x.snp <- as.matrix(snp[, c(l, k)])
  if (is.null(x.des)) {
    x <- cbind(rep(1, nrow(x.snp)),
               x.snp,
               x.snp[, 1L] * x.snp[, 2L])
    assign <- 0L:3L
    p <- fast.anova(x, y, assign)[3L]
  } else {
    x <- cbind(rep(1, nrow(x.snp)),
               x.des,
               x.snp,
               x.snp[, 1L] * x.snp[, 2L])
    assign <- c(0L, rep(1L, ncol(x.des)), 2L:4L)
    p <- fast.anova(x, y, assign)[4L]
  }
  p
}


#' @title Plot of significance of (marker by marker)-phenotype association.
#'
#' @description Plot of significance of (marker by marker)-phenotype association at their
#' genome positions.
#'
#' @param x a data.frame object of class \code{\link{epiScan}}.
#' @param threshold a significance threshold, used for coloring of plot.
#' @param col.axis the color to be used for axis annotation. Defaults to "black".
#' @param col.lab the color to be used for x and y labels. Defaults to "black".
#' @param ... other graphical parameters may also be passed as arguments to this function.
#'
#'
#' @importFrom graphics plot lines polygon rect text
#' @importFrom grDevices rainbow
#' @importFrom stats p.adjust p.adjust.methods
#' @export
plot.epiScan <- function(x, threshold = 0.05, col.axis = "black", col.lab = "black", ...) {
  # make positions linear with gaps between chr's
  stopifnot(colnames(x) == c("chr1", "pos1", "chr2", "pos2", "pValues"))
  chr <- unique(c(x[, 1L], x[, 3L]))
  chrminmax <- vapply(split(c(x[, 2L], x[, 4L]), c(x[, 1L], x[, 3L])), 
                      function(x) c(min(x), max(x)), c(1, 2))
  chrgap <- round(sum(chrminmax[2L, ] - chrminmax[1L, ]) * .01, 0)
  chrsize <- cumsum(c(0, chrminmax[2L, -ncol(chrminmax)])) -
    cumsum(chrminmax[1L, ]) +
    cumsum(c(0, rep(chrgap, ncol(chrminmax) - 1L)))
  chrstartend <- chrminmax + rbind(chrsize, chrsize) + chrgap / 2
  xthr <- x[x[, 5L] <= threshold, , drop = FALSE]
  for (i in seq_along(chr)) {
    inx1 <- which(xthr[, 1L] == chr[i])
    if (length(inx1))
      xthr[inx1, 2L] <- xthr[inx1, 2L] + chrsize[i] + chrgap / 2
    inx3 <- which(xthr[, 3L] == chr[i])
    if (length(inx3))
      xthr[inx3, 4L] <- xthr[inx3, 4L] + chrsize[i] + chrgap / 2
  }
  a <- pi * 2 / (chrstartend[2L, ncol(chrstartend)] + chrgap / 2)
  xthr[, 5L] <- -log10(xthr[, 5L])
  thr <- -log10(threshold)
  # plot 
  plot(NA, xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25), type = "n", axes = FALSE,
       xaxs = "i", yaxs = "i", xlab = "", ylab = "", asp = 1, ...)
  # plot axis 
  for (i in 1:ncol(chrstartend)) {
    chr.x <- seq(chrstartend[1L, i], chrstartend[2L, i], length.out = 100)
    chr.x <- c(chr.x[1], chr.x, chr.x[length(chr.x)])
    chr.y <- c(1.03, rep(1, 100), 1.03)
    chr.x.circ <- chr.y * sin(a * (chr.x - 1))
    chr.y.circ <- chr.y * cos(a * (chr.x - 1))
    lines(x = chr.x.circ, chr.y.circ, col = col.axis)
    # plot axis lable
    lchr.x.circ <- 1.15 * sin(a * (chr.x[50] - 1))
    lchr.y.circ <- 1.15 * cos(a * (chr.x[50] - 1))
    text(lchr.x.circ, lchr.y.circ, paste("Chr. " , i), col = col.lab, xpd = TRUE)
  }
  # plot lines
  if (max(xthr[, 5L]) > thr) {
    colinter <- seq(thr, max(xthr[, 5L]), length.out = 20)
    lcol <- length(colinter)
    colgrad <- rev(rainbow(n = lcol, start = 0.05, end = .16))
    for (l in order(xthr[, 5L])) {
      if (xthr[l, 5L] > thr) {
        x1.circ <- 0.975 * sin(a * (xthr[l, 2L] - 1))
        y1.circ <- 0.975 * cos(a * (xthr[l, 2L] - 1))
        x2.circ <- 0.975 * sin(a * (xthr[l, 4L] - 1))
        y2.circ <- 0.975 * cos(a * (xthr[l, 4L] - 1))
        x1x2 <- c(x1.circ, 0, x2.circ)
        y1y2 <- c(y1.circ, 0, y2.circ)
        x1x2.poly <- list(x1.circ)
        y1y2.poly <- list(y1.circ)
        # shmooth y transition from x1 to x2
        for (i in 2:99) {
          x1x2i.poly <- y1y2i.poly <- 0
          ifact <- i / 100
          const <- (1 - ifact)^2
          for (j in 1:3) {
            x1x2i.poly <- x1x2i.poly + const * x1x2[j]
            y1y2i.poly <- y1y2i.poly + const * y1y2[j]
            jfact <- j / 3
            const <- const * (1 - jfact) / jfact * ifact / (1 - ifact)
          }
          x1x2.poly[[i]] <- x1x2i.poly
          y1y2.poly[[i]] <- y1y2i.poly
        }
        x1x2.poly <- unlist(c(x1x2.poly, list(x2.circ)))
        y1y2.poly <- unlist(c(y1y2.poly, list(y2.circ)))
        lines(x1x2.poly, y1y2.poly, col = colgrad[findInterval(xthr[l, 5L], colinter)])
      }
    }
  }
  # Color key
  text(.825, -1.19, pos = 3, labels = "p-values", col = col.lab)
  rect(seq(.55, 1.1, length.out = lcol)[-1], rep(-1.25, lcol - 1),
       seq(.55, 1.1, length.out = lcol)[-lcol], rep(-1.2, lcol - 1),
       col = colgrad, border = NA, xpd = TRUE)
  lines(c(.555, 1.095), c(-1.27, -1.27), col = col.axis, xpd = TRUE)
  t.pos <- seq(.555, 1.095, length.out = 4)
  t.lab <- format(10^-seq(thr, colinter[length(colinter)], length.out = 4), digits = 3)
  for (i in 1:5) {
    lines(c(t.pos[i], t.pos[i]), c(-1.27, -1.29), col = col.axis, xpd = TRUE)
    text(t.pos[i], -1.35, pos = 1, labels = t.lab[i], col = col.lab, srt = 45, xpd = TRUE)
  }
}

