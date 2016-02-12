
#' @title Create vector of non-zero growth rates for basal species
#'
#' @param amat Adjacency matrix of food web
#'
#' @return numeric vector of 1 if species is basal, and 0 otherwise
#' @export
#'

getR <- function(amat){
  r.i <- c()
  r.i[colSums(amat) == 0] <- 1
  r.i[colSums(amat) != 0] <- 0

  if(sum(r.i) == 0){r.i[sample(1:length(r.i), 1)] <- 1}
  return(r.i)
}



#' @title Basal growth function
#'
#' @description \code{G.i} returns population growth of basal species in the absence of predation
#'
#' @param r Growth rate in the absence of predators
#' @param B Biomass
#' @param K Carrying capacity
#'
#' @details testing a new tag
#' @return numeric vector of population sizes
#' @export

Gi <- function(r, B, K){return(r * B * (1 - (B/K)))}


#' @title Lotka-Volterra functional response
#'
#' @param B vector of biomasses
#' @param A binary adjacency matrix of the food web
#' @param B.0 half saturation constant
#' @param xpar control parameter that alters the form of the functional response
#'
#' @export
#'

Fij <- function(B, A, B.0, xpar){
  sum.bk <- rowSums(sapply(1:nrow(A), function(x){B[x] * A[x,]}))^(1+xpar)
  denom <- sum.bk + B.0^(1+xpar)

  F1 <- sapply(1:nrow(A), function(x){(B[x] * A[x,])^(1+xpar)})/denom

  return(F1)
}


#' @title Consumer-Interference functional response
#'
#' @param B vector of biomasses
#' @param A binary adjacency matrix of the food web
#' @param B.0 half saturation constant
#' @param xpar control parameter that alters the strength of interference
#'
#' @export
#'

Fbd <- function(B, A, B.0, xpar){
  sum.bk <- rowSums(sapply(1:nrow(A), function(x){B[x] * A[x,]}))
  denom <- sum.bk + (1 + (xpar * B)) * B.0

  F1 <- sapply(1:nrow(A), function(x){(B[x] * A[x,])})/denom

  return(F1)
}


#' @title Dynamic function for input into ode solver
#'
#' @param t time sequence for which output is wanted
#' @param states initial biomass values for the species in the system
#' @param par list of named parameters for input into the model
#'
#' @export

CRmod <- function(t,states,par){

  with(as.list(c(states, par)), {
    dB <- G.i(r = r.i, B = states, K = K) - x.i*states + rowSums((x.i * yij * FR(states, A, B.o, xpar = xpar) * states)) - rowSums((x.i * yij * t(FR(states, A, B.o, xpar = xpar)* states))/eij)

    list(c(dB))
  })

}


#' @title Extinction event function for the ODE system
#'
#' @param times time sequence for which output is wanted
#' @param states biomass values for the species in the system
#' @param parms additional required parameters for the event
#'
#' @export
#'
goExtinct <- function(times, states, parms){
  with(as.list(states), {
    for(i in 1:length(states)){
      if(states[i] < 10^-10){states[i] <- 0}else{states[i]}
    }
    return(c(states))
  })
}


#' @title Consumer-resource biomass dynamics model
#'
#' @param Adj binary adjacency matrix of the food web
#' @param t time sequence for which output is wanted
#' @param G function describing biomass growth in absence of predation
#' @param method function to be used in the ode solver
#' @param FuncRes name of the functional response to be used in the model
#' @param K carrying capacity
#' @param x.i mass specific metabolic rate
#' @param yij maximum rate at which species i assimilates species j
#' @param eij conversion efficiency
#' @param xpar control parameter for the functional response
#' @param B.o half saturation constant
#' @param plot logical determining whether or not a plot of biomass over time should be printed
#'
#' @export

CRsimulator <- function(Adj, t = 1:200, G = Gi, method = CRmod, FuncRes = Fij, K = 1, x.i = .5, yij = 6, eij = 1, xpar = .2, B.o =.5, ext = goExtinct, plot = FALSE){
  require(deSolve)

  grow <- get.r(Adj)

  par <- list(
    K = K,
    x.i = x.i,
    yij = yij,
    eij = 1,
    xpar = xpar,
    B.o = B.o,
    r.i = grow,
    A = Adj,
    G.i = G,
    FR = FuncRes
  )

  states <- runif(nrow(Adj), .5, 1)

  out <- ode(y=states, times=t, func=method, parms=par, events = list(func = ext, time = t))

  if(plot) print(matplot(out[,-1], typ = "l", lwd = 2, xlab = "Time", ylab = "Biomass"))

  return(out)
}





