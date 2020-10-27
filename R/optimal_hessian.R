#' @title Log-likelihood Hessian
#'
#' @description This function compute the hessian matrix of the log-likelihood
#' function given the data and the covariates
#'
#' @param par The point in which the hessian matrix has to be computed
#' @param data The matrix of the observed data
#' @param lista_phi The list containing all the matrices of covariates to model each \code{phi} element
#' @param lista_d  The list containing all the matrices of covariates to model each \code{d} element
#' @param fnscale Scale coefficient: default value equal to \code{1}
#'
#' @return Hessian matrix of the log-likelihood function
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

#' optimal_hessian(par,data,lista_phi,lista_d)
#'
#'
optimal_hessian <- function(par,data,lista_phi,lista_d, fnscale = 1){
  # Information from the input
  n <- nrow(data)
  p <- ncol(data)
  q <- p*(p-1)/2
  l <- sum(unlist(lapply(lista_phi, ncol)))
  r <- sum(unlist(lapply(lista_d, ncol)))
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

  idxn <- NULL
  for(j in 2:p) idxn <- c(idxn,j*(j-1)/2)

  phi_index_in_T <- list()
  phi_index_in_T[[1]] <- 1
  for(j in 2:(p-1)) phi_index_in_T[[j]] <- (idxn[j-1]+1):(idxn[j])

  iv <- NULL
  iy <- NULL
  for(j in 1:length(phi_index_in_T)){
    iv <- c(iv,rep(j+1, length(phi_index_in_T[[j]]) ) )
    iy <- c(iy,1:length(phi_index_in_T[[j]]))
  }

  lista_app <- list()                     # Only used for the Beta-Beta Block
  for (j in 1:length(phi_index_in_T)) {
    app <- NULL
    for (i in phi_index_in_T[[j]]) {
      app <- c(app,index_phi[[i]])
    }
    lista_app[[j]] <- app
  }

  regressors <- NULL                      # We need this matrix for the Block Beta-Beta
  for(i in 1:(length(lista_phi))) regressors <- cbind(regressors,lista_phi[[i]])

  iy_update <- NULL                       # Only used for the Beta-Beta Block
  iv_update <- NULL                       # Only used for the Beta-Beta Block
  for (i in 1:length(phi_index_in_T)) {
    for (j in 1:length(phi_index_in_T[[i]])) {
      iy_update <- c(iy_update,rep(j,length(index_phi[[phi_index_in_T[[i]][j]]])))
      iv_update <- c(iv_update,rep(i+1,length(index_phi[[phi_index_in_T[[i]][j]]])))
    }
  }
  # Hessian Blocks Initialization
  B_L_L <- matrix(rep(0,r^2),ncol = r)    # Block with the partial derivative Lambda-Lambda
  B_L_B <- matrix(rep(0,r*l),ncol = l)    # Block with the partial derivative Lambda-Beta
  B_B_B <- matrix(rep(0,l^2), ncol = l)   # Block with the partial derivative Beta-Beta

  for(i in 1:n){
    # For all i we compute the vectors v, d^-1 and omega (w), that will be used to
    # compute each hessian matrix block
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
    w <- d_m1*v

    # Block Lambda-Lambda iterative computation
    for (j in 1:p) {
      for (s in 1:ncol(lista_d[[j]])) {
        for (h in 1:ncol(lista_d[[j]])) {
          B_L_L[index_d[[j]][s],index_d[[j]][h]] <- B_L_L[index_d[[j]][s],index_d[[j]][h]] -
            0.5* (d_m1[j]^2)*(v[j]^2)*lista_d[[j]][i,s]*lista_d[[j]][i,h]
        }
      }
    }

    # Block Lambda-Beta iterative computation
    for (j in 1:q) {
      for(s in 1:ncol(lista_d[[iv[j]]])){
        for (t in 1:ncol(lista_phi[[j]])) {
          B_L_B[index_d[[iv[j]]][s],index_phi[[j]][t]] <- B_L_B[index_d[[iv[j]]][s],index_phi[[j]]][t]+
            w[iv[j]]*d_m1[iv[j]]*data[i,iy[j]]*lista_phi[[j]][i,t]*lista_d[[iv[j]]][i,s]
        }
      }
    }

    # Block Beta-Beta iterative computation
    for (j in 1:(p-1)) {
      for (u in lista_app[[j]]) {
        for (z in lista_app[[j]]) {
          B_B_B[u,z] <- B_B_B[u,z]-
            regressors[i,u]*regressors[i,z]*data[i,iy_update[u]]*data[i,iy_update[z]]*(d_m1[iv_update[u]])^2
        }
      }
    }

  }

  B_B_L <- t(B_L_B)      # Block with the partial derivative Beta-Lambda (Schwarz Thm)

  # Hessian matrix Assemblage
  H <- rbind(cbind(B_B_B,B_B_L),cbind(B_L_B,B_L_L))

  return(fnscale*H)      # Using fnscale = -1, can be useful for optimization
}



