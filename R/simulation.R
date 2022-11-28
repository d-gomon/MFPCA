#' Simulate multivariate functional datav2
#' 
#' This function provides a unified simulation structure for multivariate 
#' functional data \eqn{f_1, \ldots, f_N} on one- or two-dimensional domains, 
#' based on a truncated multivariate Karhunen-Loeve representation: \deqn{f_i(t)
#' = \sum_{m = 1}^M \rho_{i,m} \psi_m(t).} The multivariate eigenfunctions 
#' (basis functions) \eqn{\psi_m}  are constructed from univariate orthonormal 
#' bases. There are two different concepts for the construction, that can be 
#' chosen by the parameter \code{type}: A split orthonormal basis (\code{split},
#' only one-dimensional domains) and weighted univariate orthonormal bases 
#' (\code{weighted}, one- and two-dimensional domains). The scores 
#' \eqn{\rho_{i,m}} in the Karhunen-Loeve representation are simulated 
#' independently from a normal distribution with zero mean and decreasing 
#' variance. See Details.
#' 
#' The parameter \code{type} defines how the eigenfunction basis for the 
#' multivariate Karhunen-Loeve representation is constructed: \itemize{ \item 
#' \code{type = "split"}: The basis functions of an underlying 'big' orthonormal
#' basis are split in \code{M} parts, translated and possibly reflected. This 
#' yields an orthonormal basis of multivariate functions with \code{M} 
#' elements. This option is implemented only for one-dimensional domains. \item 
#' \code{type = "weighted":} The multivariate eigenfunction basis consists of 
#' weighted univariate orthonormal bases.  This yields an orthonormal basis of 
#' multivariate functions with \code{M} elements. For data on two-dimensional 
#' domains (images), the univariate basis is constructed as a tensor product of 
#' univariate bases in each direction (x- and y-direction). }
#' 
#' Depending on \code{type}, the other parameters have to be specified as 
#' follows: \subsection{Split 'big' orthonormal basis}{ The parameters \code{M} 
#' (integer), \code{eFunType} (character string) and \code{ignoreDeg} (integer 
#' vector or \code{NULL}) are passed to the function \code{\link{eFun}} to 
#' generate a univariate orthonormal basis on a 'big' interval. Subsequently, 
#' the basis functions are split and translated, such that the \eqn{j}-th part 
#' of the split function is defined on the interval corresponding to 
#' \code{argvals[[j]]}. The elements of the multivariate basis functions are 
#' given by these split parts of the original basis functions multiplied by a 
#' random sign \eqn{\sigma_j \in \{-1,1\}, j = 1, \ldots, p}{\sigma_j in {-1,1},
#' j = 1, \ldots, p}.}
#' 
#' \subsection{Weighted orthonormal bases}{ The parameters \code{argvals, M, 
#' eFunType} and \code{ignoreDeg} are all lists of a similar structure. They are
#' passed element-wise to the function \code{\link{eFun}} to generate 
#' orthonormal basis functions for each element of the multivariate functional 
#' data to be simulated. In case of bivariate elements (images), the 
#' corresponding basis functions are constructed as tensor products of 
#' orthonormal basis functions in each direction (x- and y-direction).
#' 
#' If the \eqn{j}-th element of the simulated data should be defined on a 
#' one-dimensional domain, then \itemize{ \item \code{argvals[[j]]} is a list, 
#' containing one vector of observation points. \item \code{M[[j]]} is an 
#' integer, specifying the number of basis functions to use for this entry. 
#' \item  \code{eFunType[[j]]} is a character string, specifying the type of 
#' orthonormal basis functions to use for this entry (see \code{\link{eFun}} for
#' possible options). \item \code{ignoreDeg[[j]]} is a vector of integers, 
#' specifying the degrees to ignore when constructing the orthonormal basis 
#' functions. The default value is \code{NULL}. }
#' 
#' If the \eqn{j}-th element of the simulated data should be defined on a 
#' two-dimensional domain, then \itemize{ \item \code{argvals[[j]]} is a list, 
#' containing two vectors of observation points, one for each direction 
#' (observation points in x-direction and in y-direction). \item \code{M[[j]]} 
#' is a vector of two integers, giving the number of basis functions for each 
#' direction (x- and y-direction). \item \code{eFunType[[j]]} is a vector of two
#' character strings, giving the type of orthonormal basis functions for each 
#' direction (x- and y-direction, see \code{\link{eFun}} for possible options).
#' The corresponding basis functions are constructed as tensor products of 
#' orthonormal basis functions in each direction. \item \code{ignoreDeg[[j]]} is
#' a list, containing two integer vectors that specify the degrees to ignore 
#' when constructing the orthonormal basis functions in each direction. The 
#' default value is \code{NULL}. } The total number of basis functions (i.e. the
#' product of \code{M[[j]]} for all \code{j}) must be equal!}
#' 
#' @param type A character string, specifying the construction method for the 
#'   multivariate eigenfunctions (either \code{"split"} or \code{"weighted"}). 
#'   See Details.
#' @param argvals A list, containing the observation points for each element of 
#'   the multivariate functional data that is to be simulated. The length of 
#'   \code{argvals} determines the number of elements in the resulting simulated
#'   multivariate functional data. See Details.
#' @param M An integer (\code{type = "split"}) or a list of integers (\code{type
#'   = "weighted"}), giving the number of univariate basis functions to use. See 
#'   Details.
#' @param eFunType A character string (\code{type = "split"})   or a list of 
#'   character strings (\code{type = "weighted"}), specifying the type of 
#'   univariate orthonormal basis functions to use. See Details.
#' @param ignoreDeg A vector of integers (\code{type = "split"})   or a list of 
#'   integer vectors (\code{type = "weighted"}), specifying the degrees to 
#'   ignore when generating the univariate orthonormal bases. Defaults to 
#'   \code{NULL}. See Details.
#' @param eValType A character string, specifying the type of 
#'   eigenvalues/variances used for the simulation of the multivariate functions
#'   based on the truncated Karhunen-Loeve representation. See 
#'   \code{\link{eVal}} for details.
#' @param N An integer, specifying the number of multivariate functions to be 
#'   generated.
#'   
#' @return \item{simData}{A \code{\linkS4class{multiFunData}} object with 
#'   \code{N} observations, representing the simulated multivariate functional 
#'   data.} \item{trueFuns}{A \code{\linkS4class{multiFunData}} object with 
#'   \code{M} observations, representing the multivariate eigenfunction basis 
#'   used for simulating the data.} \item{trueVals}{A vector of numerics, 
#'   representing the eigenvalues used for simulating the data.}
#'   
#' @seealso \code{\linkS4class{multiFunData}}, \code{\link{eFun}}, 
#'   \code{\link{eVal}}, \code{\link{simFunData}}, \code{\link{addError}}, 
#'   \code{\link{sparsify}}.
#'   
#' @references C. Happ, S. Greven (2018): Multivariate Functional Principal 
#'   Component Analysis for Data Observed on Different (Dimensional) Domains. 
#'   Journal of the American Statistical Association, 113(522): 649-659. 
#'   
#' @importFrom stats rnorm
#'   
#' @export simMultiFunDatav2
#'   
#' @examples
#' oldPar <- par(no.readonly = TRUE)
#' 
#' # split
#' split <- simMultiFunData(type = "split", argvals = list(seq(0,1,0.01), seq(-0.5,0.5,0.02)),
#'                  M = 5, eFunType = "Poly", eValType = "linear", N = 7)
#' 
#' par(mfrow = c(1,2))
#' plot(split$trueFuns, main = "Split: True Eigenfunctions", ylim = c(-2,2))
#' plot(split$simData, main = "Split: Simulated Data")
#' 
#' # weighted (one-dimensional domains)
#' weighted1D <- simMultiFunData(type = "weighted",
#'                  argvals = list(list(seq(0,1,0.01)), list(seq(-0.5,0.5,0.02))),
#'                  M = c(5,5), eFunType = c("Poly", "Fourier"), eValType = "linear", N = 7)
#' 
#' plot(weighted1D$trueFuns, main = "Weighted (1D): True Eigenfunctions", ylim = c(-2,2))
#' plot(weighted1D$simData, main = "Weighted (1D): Simulated Data")
#' 
#' # weighted (one- and two-dimensional domains)
#' weighted <- simMultiFunData(type = "weighted",
#'                argvals = list(list(seq(0,1,0.01), seq(0,10,0.1)), list(seq(-0.5,0.5,0.01))),
#'                M = list(c(5,4), 20), eFunType = list(c("Poly", "Fourier"), "Wiener"),
#'                eValType = "linear", N = 7)
#' 
#' plot(weighted$trueFuns, main = "Weighted: True Eigenfunctions (m = 2)", obs = 2)
#' plot(weighted$trueFuns, main = "Weighted: True Eigenfunctions (m = 15)", obs = 15)
#' plot(weighted$simData, main = "Weighted: Simulated Data (1st observation)", obs = 1)
#' plot(weighted$simData, main = "Weighted: Simulated Data (2nd observation)", obs = 2)
#' 
#' par(oldPar)
simMultiFunDatav2 <- function(type, argvals, M, eFunType, ignoreDeg = NULL, eValType, N, custom_values)
{
  if(! all(is.character(type), length(type) == 1))
    stop("Parameter 'type' must be passed as a string.")
  
  if(! (is.list(argvals) & all(is.numeric(unlist(argvals)))) )
    stop("Parameter 'argvals' must be passed as a list of numerics.")
  
  if(! all(is.numeric(unlist(M))))
    stop("Parameter 'M' must contain only numerics.") 
  
  if(! all(is.character(unlist(eFunType))))
    stop("Parameter 'eFunType' must contain only strings.")
  
  if(!(is.null(ignoreDeg ) | all(is.numeric(ignoreDeg), ignoreDeg > 0)))
    stop("Parameter 'ignoreDeg' must be either NULL or a vector of positive numbers.") 
  
  if(! all(is.character(eValType), length(eValType) == 1))
    stop("Parameter 'eValType' must be passed as a string.")
  
  if(! all(is.numeric(N), length(N) == 1, N > 0))
    stop("Parameter 'N' must be passed as a positive number.") 
  
  # generate eigenfunctions
  trueFuns <- switch(type,
                     split = simMultiSplitv2(argvals, M, eFunType, ignoreDeg, eValType, N),
                     weighted = simMultiWeightv2(argvals, M, eFunType, ignoreDeg, eValType, N),
                     stop("Choose either 'split' or 'weighted' for the simulation of multivariate functional data.")
  )
  
  # number of eigenfunctions generated
  Mtotal <- nObs(trueFuns)
  
  # number of elements in multivariate functional basis
  p <- length(trueFuns)
  
  # generate eigenvalues and scores
  trueVals <- eValv2(Mtotal, eValType , custom_values)
  scores <- t(replicate(N, stats::rnorm(Mtotal, sd = sqrt(trueVals))))
  
  # generate individual observations
  simData  <- vector("list", p)
  
  for(j in seq_len(p))
  {
    X <- apply(trueFuns[[j]]@X, -1, function(v){scores %*% v})
    
    if(N == 1)
      dim(X) <- c(1, nObsPoints(trueFuns[[j]]))
    
    simData[[j]] <- funData(trueFuns[[j]]@argvals, X)
  } 
  
  return(list(simData = multiFunData(simData),
              trueFuns = trueFuns,
              trueVals = trueVals,
              scores = scores))
}


#' Simulate multivariate eigenfunctions based on a split 'big' ONB
#'
#' @keywords internal
simMultiSplitv2 <- function(argvals, M, eFunType, ignoreDeg = NULL, eValType, N)
{
  # consistency check
  if(any( c(length(M), length(eFunType), length(eValType)) != 1) )
    stop("argvals, M, eFunType, eValType must all be of length 1!")
  
  # number of elements
  p <- length(argvals)
  
  # "rearrange" argvalss
  x <- vector("list", length = length(argvals))
  splitVals <- rep(NA, length(argvals) + 1)
  
  x[[1]] <- unlist(argvals[[1]]) # convert to vector, if argvals[[1]] is a list
  splitVals[1:2] <- c(0, length(x[[1]]))
  
  for(i in 2:p)
  {
    x[[i]] <- unlist(argvals[[i]]) # convert to vector, if argvals[[i]] is a list
    x[[i]] <- argvals[[i]] - min(argvals[[i]]) + max(x[[i-1]])
    splitVals[i+1] <- splitVals[i]+length(x[[i]])
  }
  
  # generate "big" orthonormal system
  f <-  eFun(unlist(x), M, ignoreDeg = ignoreDeg, type = eFunType)
  
  # sample sign randomly
  s <- sample(c(-1,1), p, 0.5)
  
  # result object
  trueFuns  <- vector("list", p)
  
  for(j in seq_len(p))
    trueFuns[[j]] <- funData(argvals[[j]],  s[j] * f@X[,(1 + splitVals[j]):splitVals[j+1]])
  
  return(multiFunData(trueFuns))
}




#' Simulate multivariate eigenfunctions based on weighted orthonormal bases
#' 
#' @importFrom foreach "%do%"
#' @importFrom stats runif
#'
#' @keywords internal
simMultiWeightv2 <- function(argvals, M, eFunType, ignoreDeg = NULL, eValType, N)
{
  p <- length(argvals)
  
  # dimension for each component
  dimsSupp <- foreach::foreach(j = seq_len(p), .combine = "c")%do%{length(argvals[[j]])}
  
  if(any(dimsSupp > 2))
    stop("Function simMultiWeight: method is not implemented for objects of dimension > 2!")
  
  if(p > 1)
  {
    if(isTRUE(do.call(all.equal, lapply(M, prod))))
    {
      Mtotal <- prod(M[[1]])
    }
    else
      stop("Function simMultiWeight: basis dimensions must be equal!")
  }
  else
  {
    Mtotal <- prod(M[[1]])
  }
  
  # mixing parameters
  alpha <- stats::runif(p, 0.2, 0.8)
  weight <- sqrt(alpha / sum(alpha))
  
  # generate basis
  basis <- vector("list", p)
  
  for(j in seq_len(p))
  {
    if(dimsSupp[j] == 1) # one-dimensional
      basis[[j]] <- weight[j] * eFun(argvals[[j]][[1]], M = M[[j]], ignoreDeg = ignoreDeg[[j]], type = eFunType[[j]])
    else # dimsSupp[j] == 2, i.e. two-dimensional
      basis[[j]]  <- weight[j] * tensorProduct(eFun(argvals[[j]][[1]], M = M[[j]][1], ignoreDeg = ignoreDeg[[j]][[1]], type = eFunType[[j]][1]),
                                               eFun(argvals[[j]][[2]], M = M[[j]][2], ignoreDeg = ignoreDeg[[j]][[2]], type = eFunType[[j]][2]))
  }
  
  return(multiFunData(basis))
}


#' Generate a sequence of simulated eigenvalues
#' 
#' This function generates \eqn{M} decreasing eigenvalues.
#' 
#' The function implements three types of eigenvalues: \itemize{\item 
#' \code{"linear":} The eigenvalues start at  \eqn{1} and decrease linearly 
#' towards \eqn{0}: \deqn{\nu_m = \frac{M+1-m}{m}.}{\nu_m = (M+1-m)/m.} \item
#' \code{"exponential":} The eigenvalues start at \eqn{1}  and decrease
#' exponentially towards \eqn{0}: \deqn{\nu_m =
#' \exp\left(-\frac{m-1}{2}\right).}{\nu_m = exp(-(m-1)/2).}\item
#' \code{"wiener":} The eigenvalues correspond to the eigenvalues of the Wiener
#' process: \deqn{\nu_m = \frac{1}{(\pi/2 \cdot (2m-1))^2}.}{\nu_m = (pi/2 *
#' (2m-1))^(-2)} }
#' 
#' @param M An integer, the number of eigenvalues to be generated.
#' @param type A character string specifying the type of eigenvalues that should
#'   be calculated. See Details.
#'   
#' @return A vector containing the \code{M} decreasing eigenvalues.
#'   
#' @importFrom graphics points
#' 
#' @keywords internal

eValv2 <- function(M, type, custom_values)
{
  if(! all(is.numeric(M), length(M) == 1, M > 0))
    stop("Parameter 'M' must be passed as a positive number.") 
  
  if(! all(is.character(type), length(type) == 1))
    stop("Parameter 'type' must be passed as a string.")
  
  ret <- switch(type,
                linear = ((M+1) - (seq_len(M))) / M,
                exponential = exp(-((seq_len(M)-1) / 2)),
                wiener = 1/(pi/2 * (2 * (seq_len(M)) - 1))^2,
                custom = custom_values,
                stop("Choose either linear, exponential or wiener"))
  return(ret)
}