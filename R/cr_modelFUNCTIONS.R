
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

  if(sum(r.i) == 0){
    warning("No basal species in simulation")
  }
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
    dB <- G.i(r = r.i, B = states, K = K) -                                   # growth
      x.i*states +                                                            # death
      rowSums((x.i * yij * FR(states, A, B.o, xpar = xpar) * states)) -       # consumption
      rowSums((x.i * yij * t(FR(states, A, B.o, xpar = xpar)* states))/eij)   # death by predation

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

CRsimulator <- function(Adj, states = NULL, t = 1:200, G = Gi, method = CRmod, FuncRes = Fij, K = 1, x.i = .5, yij = 6, eij = 1, xpar = .2, B.o =.5, ext = goExtinct, plot = FALSE){

  grow <- getR(Adj)

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

  if(is.null(states)){states <- runif(nrow(Adj), .5, 1)}

  out <- deSolve::ode(y=states, times=t, func=method, parms=par, events = list(func = ext, time = t))

  if(plot) print(matplot(out[,-1], typ = "l", lwd = 2, xlab = "Time", ylab = "Biomass"))

  return(out)
}


#' @title Dynamic network visualization
#'
#' @param mat initial adjacency matrix of the food web
#' @param dyn output of \code{CRsimulator}
#' @param path1 where to save the html plot
#'
#' @return An html file of the changes to the food web during the simulation.
#' @export
#'

netHTML <- function(mat, dyn, path1 = getwd()){
  if(!requireNamespace("animation", quietly = TRUE)){stop("This function requires the 'animation' package to be installed and loaded", call. = FALSE)}

  lay <- matrix(c(layout.sphere(graph.adjacency(mat))[,1], TrophInd(mat)$TL), ncol = 2)
  s <- matrix(0, nrow = nrow(dyn), ncol = ncol(mat))

  ani.options(interval = .25)
  saveHTML(
    {
      for(i in 1:50){
        fr <- Fij(dyn[i,-1], mat, .5, .2)
        strength <- melt(fr)[,3][melt(fr)[,3] > 0]
        fr[fr > 0 ] <- 1

        g.new <- graph.adjacency(t(fr))
        E(g.new)$weight <- strength/max(strength)*10
        s[i,c(which(dyn[i,-1] > 0))] <-log(dyn[i, c(which(dyn[i,] > 0)[-1])])+abs(min(log(dyn[i, c(which(dyn[i,] > 0)[-1])])))

        plot.igraph(g.new, vertex.size = s[i,], edge.width = E(g.new)$weight, layout = lay)
      }
    },
    img.name = paste(path1, "fwdyn", sep = ""), htmlfile = paste(path1, "fwdyn.html", sep = ""),
    interval = .25, nmax =500, ani.width = 500, ani.height = 500, outdir = path1
  )
}


#' @title Converts directed matrix to undirected
#'
#' @param tm Food web adjacency matrix
#'
#' @return A new adjacency matrix where if a_ij  = 1 so does a_ji
#'


conversion <- function(tm){
  for(i in 1:nrow(tm)){
    for(j in 1:ncol(tm)){
      if(tm[i,j] == 1 & tm[j,i] == 0){tm[j,i] <- 1}
    }
  }
  return(tm)
}


#' @title Get common food web structural indices
#'
#' @param dyn Matrix of biomass dynamics from \code{CRsimulator}
#' @param web Initial food web adjacency matrix
#'
#' @return A matrix of food web indices where each row is a time step in the dynamics simulation
#' @export
#'

WEBind <- function(dyn, web){
  if(!requireNamespace("rnetcarto", quietly = TRUE)){
    stop("This function requires the 'rnetcarto' package to be installed and loaded", call. = FALSE)
    }

  adj.list <- lapply(1:nrow(dyn), function(x){web[dyn[x, -1] > 0, dyn[x, -1] > 0]})
  g.list <- lapply(adj.list, graph.adjacency)

  # Number of species
  N <- sapply(adj.list, nrow)
  # Number of links
  Ltot <- sapply(adj.list, sum)
  # Link Density
  LD <- sapply(adj.list, function(x) sum(x)/nrow(x))
  # Connectance: Links / (N * (N - 1))
  C <- sapply(adj.list, function(x){sum(x)/(nrow(x) * (nrow(x) - 1))})
  # Web diameter
  D <- sapply(g.list, diameter)
  # Average path length
  APL <- sapply(g.list, average.path.length)
  # Clustering coefficient
  CC <- sapply(g.list, transitivity)
  # Modularity
  mod <- lapply(lapply(adj.list, conversion), function(x){
    if(nrow(x) > 2 && sum(x) > 2){rnetcarto::netcarto(x)}else{list(data.frame(module = c(0,0)),0)}
  })
  M <- sapply(mod, "[[", 2)
  nMod <- sapply(lapply(mod, "[[", 1), function(x) max(x$module) + 1)

  indices <- matrix(c(N, Ltot, LD, C, D, APL, CC, M, nMod), nrow = nrow(dyn))
  colnames(indices) <- c("N", "Ltot", "LD", "C", "D", "APL", "CC", "M", "nMod")
  return(indices)
}


#' @title Find the counts of three species configurations through time from \code{CRsimulator}
#'
#' @param dyn Matrix of biomass dynamics from \code{CRsimulator}
#' @param web Initial food web adjacency matrix
#'
#' @return A dataframe of counts of three species configurations where each row is a time step in the model
#' @export
#'

motifCounter3 <- function(dyn, web){
  adj.list <- lapply(1:nrow(dyn), function(x){web[dyn[x, -1] > 0, dyn[x, -1] > 0]})
  g.list <- lapply(adj.list, graph.adjacency)

  triad.count <- lapply(g.list, triad.census)
  triad.matrix <- matrix(unlist(triad.count), nrow = length(g.list), ncol = 16, byrow = T)
  colnames(triad.matrix) <- c("empty", "single", "mutual", "s5", "s4", "s1", "d4",
                              "d3", "s2", "s3","d8", "d2", "d1", "d5", "d7", "d6")

  triad.df <- as.data.frame(triad.matrix)

  motif.data.frame <- data.frame(s1 = triad.df$s1, s2 = triad.df$s2, s3 = triad.df$s3, s4 = triad.df$s4,
                                 s5 = triad.df$s5, d1 = triad.df$d1, d2 = triad.df$d2, d3 = triad.df$d3, d4 = triad.df$d4,
                                 d5 = triad.df$d5, d6 = triad.df$d6, d7 = triad.df$d7, d8 = triad.df$d8)

  return(motif.data.frame)
}


#' Change in trophic level through time
#'
#' @param dyn Matrix of biomass dynamics from \code{CRsimulator}
#' @param web Initial food web adjacency matrix
#'
#' @return Matrix of trophic position for each species in the food web
#' @export
#'


trophicChange <- function(dyn, web){
  if(!requireNamespace("NetIndices", quietly = TRUE)){
    stop("This function requires the 'NetIndices' package to be installed and loaded", call. = FALSE)
    }

  adj.list <- lapply(1:nrow(dyn), function(x){web[dyn[x, -1] > 0, dyn[x, -1] > 0]})
  til <- lapply(adj.list, NetIndices::TrophInd)

  m <- matrix(0, ncol = (ncol(dyn) - 1), nrow = nrow(dyn))
  for(x in 1:nrow(dyn)){
    m[x, dyn[x, -1] > 0] <- til[[x]]$TL
  }
  return(m)
}
