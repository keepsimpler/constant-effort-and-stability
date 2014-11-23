library(vegan)
library(bipartite)
library(igraph)
library(plyr)

library(doMC)  # 
registerDoMC()  # register Multi Cores
getDoParWorkers()  # get available Cores


###############################################################################
#' @title Generate a connected graph using package [igraph]
#'
#' @param s, size of network. 
#' if graph type is bipartite, s[1], s[2] represent size of two groups; else s is size of network
#' @param k, average degree for the network.
#' 1 < k < s for unipartite network, 1 < k < s[1]*s[2]/(s[1]+s[2]) for bipartite network.
#' @param gtype, Graph type generated: 'bipartite', 'sf', 'er', 'regular'.
#' @param maxtried, the maximum number of tried times. 
#' If have tried [maxtried] times, the function will return no matter whether the connected graph is generated.
#' @param ... the parms optional adhere to the implementation functions of [igraph]
#' @return the connected graph or NULL
#' @details .  
#' @import igraph
get.graph.connected <- function(s, k, gtype, maxtried = 100, expower = 2.5, ...) {
  if (gtype == 'bipartite' && is.na(s[2])) {  # the bipartite graph need size of two groups of nodes
    warning('sizes of TWO groups of nodes should be designated. 
            we have assumed the size of second group equal to the size of first group.')
    s[2] = s[1]  # if missed second size, we assume it equal to the first size.
  }
  count = 0
  repeat {  # generate a connected graph
    if (gtype == 'bipartite') {
      G = bipartite.random.game(s[1], s[2], type = 'gnm', m = ceiling(k * (s[1] + s[2])))
    } else if (gtype == 'sf') {
      G = static.power.law.game(s, k * s, exponent.out = expower)
    }
    else if (gtype == 'er') {
      G = erdos.renyi.game(s, p.or.m = k * s, type = 'gnm')
    }
    else if (gtype == 'regular') {
      G = k.regular.game(s, k)
    }
    else if (gtype == 'complete') {
      G = graph.full(s)
    }
    if (igraph::is.connected(G)) break  # until a connected graph is generated
    count = count + 1
    if (count == maxtried) {
      warning(paste('Tried', maxtried, 'times, But connected graph still cannot be generated.'))
      break
    }
  }
  G
}

#plot(G, layout = layout.bipartite)


#' @title assign interaction strength and interaction types for community matrix
#' @param graph, generated networks
#' @param beta0, intraspecific competition, used to decide the upper bound of local stability criteria
#' @param beta1.min, beta1.max, beta1.random.type, parameters of distribution of interspecific interaction strengths
#'        if random distribution type is 'norm', then [beta1.max] is SD of normal distribution
#'        if random distribution type is 'unif', then [beta1.max] is the upper bound of uniform distribution
#'        [beta1.min] which is the lower bound of uniform distribution, is conserved to equal 0
#' @param inter.type, interaction types among species:
#'        arbitrary, competition, mutualism, predatorprey, mixture
get.community.matrix <- function(graph, beta0 = 0, beta1.min = 0, beta1.max = 1, beta1.random.type = 'norm', inter.type) {
  if (class(graph) == 'igraph') {  # if is a igraph object of package <igraph> 
    graph = as.matrix(get.adjacency(graph))
  }
  s = dim(graph)[1]  # number of nodes
  edges = sum(graph > 0)  # number of edges
  M = graph
  if (inter.type == 'arbitrary') {
    if (beta1.random.type == 'norm')
      M[M > 0] = rnorm(edges, mean = 0, sd = beta1.max)
    else if (beta1.random.type == 'unif')
      M[M > 0] = runif(edges, min = - beta1.max, max = beta1.max)
  }
  else if (inter.type == 'competition') {
    if (beta1.random.type == 'norm')
      M[M > 0] = - abs(rnorm(edges, mean = 0, sd = beta1.max))
    else if (beta1.random.type == 'unif')
      M[M > 0] = runif(edges, min = - beta1.max, max = 0)
  }
  else if (inter.type == 'mutualism') {
    if (beta1.random.type == 'norm')
      M[M > 0] = abs(rnorm(edges, mean = 0, sd = beta1.max))
    else if (beta1.random.type == 'unif')
      M[M > 0] = runif(edges, min = 0, max = beta1.max)
  }
  else if (inter.type == 'predatorprey') {
    for (i in 1:s) {
      for (j in (i + 1) : s) {
        if (M[i, j] > 0) {  # if interaction exist
          if (sample(c(0, 1), 1) == 0) { # half probability
            if (beta1.random.type == 'norm') {
              M[i, j] = abs(rnorm(1, mean = 0, sd = beta1.max))
              M[j, i] = - abs(rnorm(1, mean = 0, sd = beta1.max))              
            }
            else if (beta1.random.type == 'unif') {
              M[i, j] = runif(edges, min = 0, max = beta1.max)
              M[j, i] = runif(edges, min = - beta1.max, max = 0)
            }
          }
          else {
            if (beta1.random.type == 'norm') {
              M[i, j] = - abs(rnorm(1, mean = 0, sd = beta1.max))
              M[j, i] = abs(rnorm(1, mean = 0, sd = beta1.max))                          
            }
            else if (beta1.random.type == 'unif') {
              M[i, j] = runif(edges, min = - beta1.max, max = 0)
              M[j, i] = runif(edges, min = 0, max = beta1.max)
            }
          }
        }
      }
    }
  }  # end of predatorprey
  diag(M) = rep(beta0, s)  # assign the diagonal elements
  M
}



#' @title get covariance matrix of multivariate OU process
#' @param phi, community matrix
#' @param C, the covariance matrix of environmental fluctuating
mou.vars <- function(phi, C) {
  s = dim(phi)[1]
  I = diag(1, s)
  - matrix(solve(kronecker(I, phi) + kronecker(phi, I)) %*% as.vector(C), nrow = s, ncol = s)
}

get.stability <- function(phi, C) {
  n = dim(phi)[1]
  eigs = eigen(phi)$values
  lev = max(Re(eigs))  # largest eigenvalue
  sev = max(Re(eigs))  # smallest eigenvalue
  sndev = sort(Re(eigs))[n - 1]  # second largest eigenvalue
  syn.eig = sev / sndev  # spectral gap between the second and smallest eigenvalue as sign of synchronization
  
  vars = mou.vars(phi, C)
  vars.sum = sum(vars)
  vars.max = sum(sqrt(diag(vars)))^2
  syn = vars.sum / vars.max
  
  sensitivity = - solve(phi)  # sensitivity matrix is the negative inverse of community matrix
  sens.sum = sum(sensitivity)
  sens.self = diag(sensitivity)
  sens.self.sum = sum(sens.self)
  sens.other = diag(1/sens.self) %*% sensitivity
  sens.other.sum = sum(sens.other)
  
  c(lev = lev, sev = sev, sndev = sndev, syn.eig = syn.eig, vars.sum = vars.sum, vars.max = vars.max, syn = syn,
    sens.self.sum = sens.self.sum, sens.other.sum = sens.other.sum, sens.sum = sens.sum)
}









#' @title generate Correlated Random UNIFormly variables
crunif <- function(rows, cols, cor) {
  # generate normals, check correlations
  X = array(rnorm(rows * cols), dim = c(rows, cols))
  #cor(X)
  
  # desired correlation
  #M = c(1.0, 0.7, 0.6, 0.6, 0.7, 1.0, 0.6, 0.6, 0.6, 0.6, 1.0, 0.8, 0.6, 0.6, 0.8, 1.0)
  #dim(M) = c(cols, cols)
  cor = 2 * sin(pi * cor / 6)
  
  # adjust correlations for uniforms
#   for (i in 1:(cols-1)) {
#     for (j in (i + 1) : cols) {
#       if (i != j) {
#         M[i, j] = 2 * sin(pi * M[i, j] / 6)
#         M[j, i] = 2 * sin(pi * M[j, i] / 6)
#       }
#     }
#   }
  
  # induce correlation, check correlations
  C = chol(cor)
  Y = X %*% C
  #cor(Y)
  
  # create uniforms, check correlations
  Y[, 1:cols] = pnorm(Y[, 1:cols])
  #cor(Y)
  Y
  
  # plot results (marginals)
#   par(mfrow = c(2, 2))
#   for (i in 1:cols) {
#     hist(Y[, i], main = paste('Y', i), xlab = '')
#   }
}
