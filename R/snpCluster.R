#' @title Correlation distance interactions
#' @description Correlation distance interactions
#' 
#' @param snp A object of class \linkS4class{snpMatrix}.
#' @param iInx ...
#' 
#' @importFrom methods is
#' @export
iDistCor <- function(snp, iInx) {
  stopifnot(is(snp, "snpMatrix"))
  if (missing(iInx))
    stop("'iInx' must be specifid")
  if (any(colnames(iInx) != c("inx1", "inx2")))
    stop("Column of 'iInx' have to be named 'inx1' and 'inx2'")
  out <- corDists_x(snp@snpData, iInx[, 1L], iInx[, 2L])
  names.snp <- colnames(snp)
  namesX <- names.snp[iInx[, 1L]]
  inx.inxX <- which(iInx[, 1L] != iInx[, 2L])
  namesX[inx.inxX] <- paste(names.snp[iInx[inx.inxX, 1L]], 
                            names.snp[iInx[inx.inxX, 2L]], sep = ":")
  attr(out,"Labels") <- namesX
  attr(out,"Size") <- nrow(iInx)
  attr(out,"Diag") <- FALSE
  attr(out,"Upper") <- FALSE
  attr(out,"method") <- "1-abs(cor(snp))"
  attr(out,"call") <- match.call()
  class(out) <- "dist"
  out
}


#' @title Internal function
#' @description cluster at distance zero in hybclustX.cor interactions
#' 
#' @param snp A object of class \linkS4class{snpMatrix}.
#' @param iInx ...
#' @param mc.cores number of cores for parallel computing. See \code{mclapply}
#' 
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @importFrom stats optimise
#' @export
iIdenticals <- function(snp, iInx, mc.cores = 1L) {
  stopifnot(is(snp, "snpMatrix"))
  if (missing(iInx))
    stop("'iInx' must be specifid")
  if (any(colnames(iInx) != c("inx1", "inx2")))
    stop("Column of 'iInx' have to be named 'inx1' and 'inx2'")
  p <- nrow(iInx)
  s <- optimise(function(s, p, m) p * s + p / s * (p / s - 1) / 2 * s / m,
                interval = c(2, p - 1), p = p, m = mc.cores)$minimum
  step <- as.integer(p / (s + 1))
  preclust <- corPreIdenticals_x(snp@snpData, iInx[, 1L], iInx[, 2L], step)
  kidenticals <- mclapply(preclust, 
                         function(i, snp, a, b) corIdenticals_x(snp, a, b, i), 
                         snp = snp@snpData,  a = iInx[, 1L], b = iInx[, 2L], 
                         mc.cores = mc.cores)
  identclust <- joinCorIdenticals_x(p, preclust, kidenticals)
  names.snp <- colnames(snp)
  clusters <- data.frame(names1 = names.snp[iInx[, 1L]],
                         names2 = names.snp[iInx[, 2L]],
                         clusters = identclust[[1L]], 
                         stringsAsFactors = FALSE)
  medoids <- data.frame(names1 = names.snp[identclust[[2L]]$inx1],
                        names2 = names.snp[identclust[[2L]]$inx2], 
                        stringsAsFactors = FALSE)
  out <- list(clusters = clusters,
              medoids = medoids,
              iInx = identclust[[2L]])
  class(out) <- "epi_identicals"
  out
}


#' @title A Clustering Algorithm based on Randomized Search with linkage 
#' correlation distance
#' @description Partitioning (clustering) into k clusters "around medoids" by 
#' randomized search. 1-abs(cor) is used as distance.
#' 
#' @param snp A object of class \linkS4class{snpMatrix}.
#' @param iInx ...
#' @param k positive integer specifying the number of clusters, less than 
#' the number of observations.
#' @param maxNeigbours positive integer specifying the maximum number of 
#' randomized searches.
#' @param nLocal positive integer specifying the number of optimisation runs.
#' Columns have to be similar to \code{snp}.
#' @param mc.cores number of cores for parallel computing. See \code{mclapply} 
#' in package parallel for details.
#' 
#' @importFrom parallel mclapply
#' @export
iClarans <- function(snp, iInx, k, maxNeigbours = 100, nLocal = 10, 
                         mc.cores = 1) {
  stopifnot(is(snp, "snpMatrix"))
  if (missing(iInx))
    stop("'iInx' must be specifid")
  if (any(colnames(iInx) != c("inx1", "inx2")))
    stop("Column of 'iInx' have to be named 'inx1' and 'inx2'")
  if (missing(k))
    stop("k must be specifid")
  if (k < 2L)
    stop("k must be at least two")
  # cluster optimisation by clarans in parallel
  clarans <- function(i, snp, iInx, k, maxNeigbours) {
    # cluster optimisation by clarans
    out <- corClarans_x(snp@snpData, iInx[, 1L], iInx[, 2L], k, maxNeigbours)
    out
  } # clarans
  out.nLocal <- mclapply(1L:nLocal, 
                         clarans,
                         snp, iInx, k, maxNeigbours,
                         mc.cores = mc.cores)
  all.objectives <- sapply(1:nLocal, function(i, x) {x[[i]][[3L]]}, out.nLocal)
  out.opt <- out.nLocal[[which.min(all.objectives)]]
  names.snp <- colnames(snp)
  clusters <- data.frame(names1 = names.snp[iInx[, 1L]],
                         names2 = names.snp[iInx[, 2L]],
                         clusters = out.opt[[1L]], 
                         stringsAsFactors = FALSE)
  medoids <- data.frame(names1 = names.snp[out.opt[[2L]]$inx1],
                        names2 = names.snp[out.opt[[2L]]$inx2], 
                        stringsAsFactors = FALSE)
  # output
  out <- list(clusters = clusters,
              medoids = medoids,
              objective = out.opt[[3L]],
              all.objectives = all.objectives)
  class(out) <- "epi_k-medoids"
  out
}


#' @title A three step Clustering Algorithm 
#' correlation distance
#' @description Dendrogram for big data.
#' 
#' @param snp A object of class \linkS4class{snpMatrix}.
#' @param iInx ...
#' @param k positive integer specifying the number of clusters, less than 
#' the number of observations.
#' @param identicals lgical, if zero clustering ...
#' @param maxNeigbours positive integer specifying the maximum number of 
#' randomized searches.
#' @param nLocal positive integer specifying the number of optimisation runs.
#' Columns have to be similar to \code{snp}.
#' @param method see hclust
#' @param mc.cores number of cores for parallel computing. See \code{mclapply} 
#' in package parallel for details.
#' @param trace ...
#' @param ... ...
#' 
#' @importFrom parallel mclapply
#' @importFrom utils installed.packages
#' @importFrom stats hclust as.dendrogram
#' @export
iqtcatClust <- function(snp, iInx, k, identicals = TRUE, maxNeigbours = 100, 
                        nLocal = 10, method = "complete", mc.cores = 1, 
                        trace = FALSE, ...) {
  if (is.element("fastcluster", rownames(installed.packages()))) {
    hclust <- fastcluster::hclust
  }
  stopifnot(is(snp, "snpMatrix"))
  if (any(snp@snpData == 0x00))
    stop("missing data are not allowed, please impute first")
  if (missing(iInx))
    stop("'iInx' must be specifid")
  if (any(colnames(iInx) != c("inx1", "inx2")))
    stop("Column of 'iInx' have to be named 'inx1' and 'inx2'")
  names.snp <- colnames(snp)
  if (identicals) {
    # identicals
    if (trace)
      cat("Step 1: Search for identicals is running\n")
    identicalFit <- iIdenticals(snp, iInx)
    iInx <- identicalFit$iInx
  } else if (trace) {
    cat("Step 1: Search for identicals is switch off\n")
  }
  # CLARANS
  nvar  <- nrow(iInx)
  if (missing(k))
    k <- as.integer(nvar / 10000L)
  if (k >= 2L) {
    if (identicals && nrow(identicalFit$medoids) <= k * 2)
      stop("Number of medoids from pefect correlated clustering is < k * 2")
    if (trace)
      cat("Step 2: CLARANS is running, 'k' is:", k, "\n")
    clarFit <- iClarans(snp, iInx, k, maxNeigbours, nLocal, mc.cores)
    if (trace)
      cat("   objectives:",
          format(clarFit$all.objectives, sort = TRUE, digits = 4), "\n")
    # if cluster < 2 add to a bigger cluster
    clust.inx <- seq_len(k)
    cluster.size <- rep(NA, k)
    for (i in clust.inx)
      cluster.size[i] <- sum(clarFit$clusters$clusters == i)
    smallclust  <- which(cluster.size < 2)
    if (length(smallclust)) {
      clust.inx <- clust.inx[-smallclust]
      min.bigclust <- clust.inx[which.min(min(cluster.size[clust.inx]))]
      clarFit$clusters[clarFit$clusters$clusters %in% smallclust, ] <- min.bigclust
    }
    if (max(cluster.size) > 65536L)
      stop("Clusters from CLARANS are to big for hclust, choose larger 'k'")
    # HClust
    if (trace)
      cat("Step 3: HClust is running\n")
    hclust.sub <- function(i, snp, iInx, clarFit, method, ...) {
      inx.i <- which(clarFit$clusters$clusters == i)
      iInx <- iInx[inx.i, ]
      out <- as.dendrogram(hclust(iDistCor(snp, iInx), method, ...))
      out
    } # hclust.sub
    hclustFit <- mclapply(clust.inx, hclust.sub, 
                          snp, iInx, clarFit, method, ..., 
                          mc.cores = mc.cores)
    dendro <- do.call(merge, c(hclustFit, height = 1, adjust = "add.max"))
  } else {
    if (ncol(snp) > 65536L)
      stop("Data size is to big for hclust, choose larger 'k'")
    # HClust
    if (trace)
      cat("Step 2: CLARANS is switch off\nStep 3: HClust is running\n")
    dendro <- as.dendrogram(hclust(iDistCor(snp, iInx), method, ...))
  }
  if (identicals) {
    out <- list(dendrogram = dendro,
                clusters = identicalFit$clusters,
                medoids = identicalFit$medoids)
  } else {
    medos <- labels(dendro)
    clust <- 1:length(medos)
    names(clust) <- medos
    out <- list(dendrogram = dendro,
                clusters = clust,
                medoids = medos)
  }
  class(out) <- "iqtcatClust"
  out
}
