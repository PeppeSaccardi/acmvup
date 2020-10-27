#' @title Log-likelihood Function
#'
#' @description This function compute the log-likelihood function value
#' given the data and the covariates
#'
#' @param par The point in which the log-likelihood has to be computed
#' @param data The matrix of the observed data
#' @param lista_phi The list containing all the matrices of covariates to model each \code{phi} element
#' @param lista_d  The list containing all the matrices of covariates to model each \code{d} element
#' @param fnscale Scale coefficient: default value equal to \code{1}
#'
#' @return Log-Likelihood value
#' @export
#'
#' @examples
#' data <- matrix(rnorm(300), ncol=3)
#' lista_d <- list()
#' lista_phi <- list()

#' lista_d[[1]] <- matrix(c(rep(1,100),rnorm(100)), byrow = FALSE, ncol=2)
#' lista_d[[2]] <- matrix(c(rep(1,100),rnorm(200)), byrow = FALSE, ncol=3)
#' lista_d[[3]] <- matrix(rep(1,100), byrow = FALSE, ncol=1)

#' lista_phi[[1]] <- matrix(c(rep(1,100),rnorm(200)),byrow = FALSE, ncol=3)
#' lista_phi[[2]] <- matrix(rep(1,100),ncol=1)
#' lista_phi[[3]] <- matrix(c(rep(1,100),rnorm(100)),byrow = FALSE, ncol=2)

#' par <- rnorm(12)

#' optimal_loglik(par,data,lista_phi,lista_d)
#'
#'
optimal_loglik <- function(par,data,lista_phi,lista_d,fnscale = 1){
  # Information from the input
  n <- nrow(data)
  p <- ncol(data)
  q <- p*(p-1)/2
  l <- sum(unlist(lapply(lista_phi, ncol)))
  beta <- par[1:l]
  lambda <- par[-c(1:l)]

  # Definition of some useful lists and vectors that implement the idea of all the
  # indicator functions
  index_d <- list()
  index_phi <- list()
  count_d <- 1
  count_phi <- 1
  index_d[[1]] <- count_d:ncol(lista_d[[1]])
  index_phi[[1]] <- count_phi:ncol(lista_phi[[1]])
  for(j in 2:q){
    if( j <= p){
      count_d <- count_d + ncol(lista_d[[j-1]])
      index_d[[j]] <- (count_d):(count_d + ncol(lista_d[[j]])-1)
      count_phi <- count_phi + ncol(lista_phi[[j-1]])
      index_phi[[j]] <- count_phi:(count_phi + ncol(lista_phi[[j]])-1)
    }else{
      count_phi <- count_phi + ncol(lista_phi[[j-1]])
      index_phi[[j]] <- count_phi:(count_phi + ncol(lista_phi[[j]])-1)
    }
  }

  # Log-likelihood components initialization
  cs <- NULL                                        # Vector for the linear form piece in lambda
  for( j in 1:p) cs <- c(cs,colSums(lista_d[[j]]))
  qf <- 0                                           # Quadratic form initialization

  for(i in 1:n){
    # For all i we compute the vectors v, d^-1 and omega (w), that will be used to
    # compute iteratively the quadratic form piece of the Log-Likelihood
    v <- NULL
    phi <- NULL
    d <- NULL
    for(j in 1:q){
      if( j <= p){
        phi <- c(phi,sum(lista_phi[[j]][i,]*beta[index_phi[[j]]]))
        d <- c(d,sum(lista_d[[j]][i,]*lambda[index_d[[j]]]))

      }else{
        phi <- c(phi,sum(lista_phi[[j]][i,]*beta[index_phi[[j]]]))
      }
    }
    d_m1 <- 1/sqrt(exp(d))
    v <- data[i,1]
    for (j in 2:p) {
      ji <- 1+(j-1)*(j-2)/2
      jf <- (j-1)+(j-1)*(j-2)/2
      v <- c(v, data[i,j]+sum(data[i,1:(j-1)]*phi[ji:jf]))
    }
    qf <- qf + sum((d_m1*v)^2)

  }

  ll <- -0.5*(sum(lambda*cs)+qf)    # Final log-likelihood value

  return(fnscale*ll)                # Using fnscale = -1, can be useful for optimization

}
