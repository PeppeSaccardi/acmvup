#' @title Fitted Covariance Matrix
#'
#' @description This function compute the covariancce matrix using the optimal parameters
#' found by fitting the model, given two lists corresponding to the covariates used to
#' model each element of the unconstrained parametrization
#'
#' @param optimal_par The optimal parameters (MLE) fitted using the available data
#' @param lista_phi The list containing the row of the matrices of covariates to model each \code{phi}
#' element we are interested in retrive the Covariance matrix
#' @param lista_d  The list containing the row of the matrices of covariates to model each \code{d}
#' element we are interested in retrive the Covariance matrix
#'
#'
#' @return Log-Likelihood value
#' @export
#'
#' @examples
#'
#' lista_d <- list()
#' lista_phi <- list()
#' lista_d[[1]] <- 1
#' lista_d[[2]] <- 1
#' lista_d[[3]] <- c(1,rnorm(1))
#' lista_d[[4]] <- c(1,rnorm(2))
#'
#' lista_phi[[1]] <- c(1,rnorm(2))
#' lista_phi[[2]] <- 1
#' lista_phi[[3]] <- 1
#' lista_phi[[4]] <- c(1,rnorm(1))
#' lista_phi[[5]] <- 1
#' lista_phi[[6]] <- c(1,rnorm(3))
#'
#' optimal_par <- rnorm(19)
#'
#' S <- fitted_cov_matrix(optimal_par,lista_phi,lista_d)
#' S

fitted_cov_matrix <- function(optimal_par,lista_phi,lista_d){

  # Input information
  l <- sum(unlist(lapply(lista_phi, length)))
  r <- sum(unlist(lapply(lista_d, length)))
  beta <- optimal_par[1:l]
  lambda <- optimal_par[-c(1:l)]
  p <- length(lista_d)
  q <- length(lista_phi)

  # Useful lists for the computation of the decomposition vector phi and d
  index_d <- list()
  index_phi <- list()
  count_d <- 1
  count_phi <- 1
  index_d[[1]] <- count_d:length(lista_d[[1]])
  index_phi[[1]] <- count_phi:length(lista_phi[[1]])
  for(j in 2:q){
    if( j <= p){
      count_d <- count_d + length(lista_d[[j-1]])
      index_d[[j]] <- (count_d):(count_d + length(lista_d[[j]])-1)
      count_phi <- count_phi + length(lista_phi[[j-1]])
      index_phi[[j]] <- count_phi:(count_phi + length(lista_phi[[j]])-1)
    }else{
      count_phi <- count_phi + length(lista_phi[[j-1]])
      index_phi[[j]] <- count_phi:(count_phi + length(lista_phi[[j]])-1)
    }
  }

  # Computation of the decomposition parameters
  phi <- NULL                           # Vector of phi parameters
  d <- NULL                             # Vector of d^-2 parameters
  for(j in 1:q){
    if( j <= p){
      phi <- c(phi,sum(lista_phi[[j]]*beta[index_phi[[j]]]))

      d <- c(d,sum(lista_d[[j]]*lambda[index_d[[j]]]))

    }else{
      phi <- c(phi,sum(lista_phi[[j]]*beta[index_phi[[j]]]))
    }
  }
  d_m2 <- 1/exp(d)
  # Matrices Definition
  D_m2 <- diag(d_m2)                    # First, the diagonal matrix D^-2

  T <- matrix(0,p,p)
  index_lwr <- which(lower.tri(T, diag = FALSE), arr.ind=TRUE)
  T[index_lwr] <- phi                   # Second, the lower-unit triangular
  T <- T + diag(rep(1,p))                 # matrix T


  Sigma <- solve(t(T)%*%D_m2%*%T)       # We return the covariance matrix
                                        # as the inverse of the precision matrix
  return(Sigma)
}

# optimal_par <- rnorm(19)
# S <- fitted_cov_matrix(optimal_par,lista_phi,lista_d)
# S
# eigen(S,symmetric = TRUE)$values
# Fitted matrix
# lista_d <- list()
# lista_phi <- list()

# lista_d[[1]] <- 1
# lista_d[[2]] <- 1
# lista_d[[3]] <- c(1,rnorm(1))
# lista_d[[4]] <- c(1,rnorm(2))
#
# lista_phi[[1]] <- c(1,rnorm(2))
# lista_phi[[2]] <- 1
# lista_phi[[3]] <- 1
# lista_phi[[4]] <- c(1,rnorm(1))
# lista_phi[[5]] <- 1
# lista_phi[[6]] <- c(1,rnorm(3))
