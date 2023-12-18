#' Detailed PCA
#'
#' @description
#' Detailed statistics on principal components analysis (PCA).
#'
#' @param X A complex matrix (or data frame) that serves as the data for principal
#' components analysis, featuring variables arranged in columns and observations in rows.
#'
#' @param scale A logical value (default = TRUE) indicating whether the variables should
#' be scaled to have unit variance before the analysis takes place.
#'
#' @param center A logical value (default = TRUE) indicating whether the variables should
#' be shifted to be zero centered.
#'
#' @return A list containing the following components: \describe{
#' \item{sdev}{The standard deviations of the principal components, calculated as the
#' square roots of the eigenvalues of the covariance/correlation matrix. The calculation
#' is actually done with the singular values of the data matrix.}
#' \item{eigval}{The eigenvalues of the covariance/correlation matrix, representing the
#' variance of the principal components. The calculation is actually done with the
#' singular values of the data matrix.}
#' \item{pctvar}{The percentage of variance explained by each principal component.}
#' \item{cumvar}{The cumulative percentage of variance explained by each principal component.}
#' \item{rotation}{The matrix of variable loadings (i.e., a matrix whose columns contain
#' the eigenvectors).}
#' \item{score}{The value of the rotated data (the centred (and scaled if requested) data
#' multiplied by the rotation matrix), representing the scores of the supplied data on all
#' principal components.}
#' \item{Psi}{The Psi index, which depends on the magnitude of the eigenvalues taken from
#' the correlation matrix of the data set.
#' \deqn{\Psi = \sum(\lambda_i-1)^2,}
#' where \eqn{\lambda} is the eigenvalue.}
#' \item{Phi}{The Phi statistic, which measures the average level of correlation among
#' the variables.
#' \deqn{\Phi = \sqrt{\frac{\sum\lambda_i^2-p}{p(p-1)}},}
#' where \eqn{\lambda} is the eigenvalue, and \eqn{p} is the number of variables.}
#' \item{correlation}{The correlations of the principal components with the variables
#' (Jackson, 1991).}
#' \item{indexload}{The index of the loadings (Vieira, 2012).}
#' \item{center}{The centering used, or \code{FALSE}.}
#' \item{scale}{The scaling used, or \code{FALSE}.}
#' }
#'
#' @references \itemize{
#' \item Gleason, T. C. and Staelin, R. (1975).
#' A Proposal for Handling Missing Data.
#' \emph{Psychometrika}, 40(2), 229--252.
#' \item Mardia, K. V., Kent, J. T. and Bibby, J. M. (1979).
#' Multivariate Analysis. London, UK: Academic Press.
#' \item Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988).
#' The New S Language. Wadsworth & Brooks/Cole.
#' \emph{Computer Science Series, Pacific Grove, CA}.
#' \item Jackson, J. Edward. (1991).
#' A User's Guide to Principal Components.
#' John Wiley & Sons, New York, USA.
#' \item Vieira, Vasco M. N. C. S. (2012).
#' Permutation Tests to Estimate Significances on Principal Components Analysis.
#' \emph{Computational Ecology and Software}, 2(2), 103--123.
#' \item Venables, W. N. and Ripley, B. D. (2013).
#' Modern Applied Statistics with S-PLUS.
#' Springer Science & Business Media.
#' }
#'
#' @export

detailedPCA <- function(X, center = TRUE, scale = TRUE) {

  X <- as.matrix(X)

  ## number of observations
  n <- nrow(X)

  ## number of variables
  p <- ncol(X)

  ## centering and scaling of matrix X
  X <- scale(X, center = center, scale = scale)

  ## the centering used
  centering.used <- attr(X, "scaled:center")
  if (is.null(centering.used)) {
    centering.used <- FALSE
  }

  ## the scaling used
  scaling.used <- attr(X, "scaled:scale")
  ## check if any variable has zero variance after scaling
  if(any(scaling.used == 0)) {
    stop("Cannot rescale a constant/zero column to unit variance!")
  }
  if (is.null(scaling.used)) {
    scaling.used <- FALSE
  }

  ## singular value decomposition of matrix X
  S <- svd(X, nu = 0)

  ## standard deviations
  sdev <- S$d / sqrt(max(1, n-1))

  ## eigenvalues
  eigval <- sdev^2

  ## If the number of observations is less than the number of variables, the covariance
  ## matrix becomes singular (non-invertible). In a singular matrix, at least one of the
  ## eigenvalues will be zero. Removing the last (smallest) eigenvalue helps address this
  ## issue and improves the stability of the eigenvalue decomposition.
  if (n < p) {
    eigval <- eigval[-length(eigval)]
  }

  ## the percentage of variance explained by each eigenvalue
  pctvar <- eigval / sum(eigval) * 100

  ## the cumulative percentage of variance explained by each eigenvalue
  cumvar <- cumsum(pctvar)

  ## the matrix of variable loadings
  dimnames(S$v) <- list(colnames(X), paste0("PC", seq_len(ncol(S$v))))
  rotation <- S$v

  ## the scores of the supplied data on the PCs
  score <- X %*% rotation

  ## Psi index
  Psi <- sum((eigval - 1)^2)

  ## Phi statistic
  Phi <- sqrt((sum(eigval^2) - p) / (p*(p-1)))

  ## the correlations of the PCs with the variables
  correlation <- t(rotation) * sqrt(eigval)

  ## the index of the loadings
  indexload <- t(rotation)^2 * eigval^2

  return(list(sdev = sdev, eigval = eigval,
              pctvar = pctvar, cumvar = cumvar,
              rotation = rotation, score = score,
              Psi = Psi, Phi = Phi,
              correlation = correlation, indexload = indexload,
              center = centering.used, scale = scaling.used))
}

