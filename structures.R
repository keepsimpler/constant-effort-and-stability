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
graph.connected <- function(s, k, gtype, maxtried = 100, expower = 2.5, ...) {
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


get.community.matrix <- function(graph, beta0 = 0, beta1.min = 0, beta1.max = 1) {
  if (class(graph) == 'igraph') {  # if is a igraph object of package <igraph> 
    graph = as.matrix(get.adjacency(graph))
  }
  s = dim(graph)[1]  # number of nodes
  edges = sum(graph > 0)  # number of edges
  B = graph
  B[B > 0] = runif(edges, min = -beta1.max, max = beta1.max)
  diag(B) = rep(beta0, s)  # assign the diagonal elements
  B
}



#' @title get covariance matrix of multivariate OU process
mou.vars <- function(phi, C) {
  s = dim(phi)[1]
  I = diag(1, s)
  - matrix(solve(kronecker(I, phi) + kronecker(phi, I)) %*% as.vector(C), nrow = s, ncol = s)
}

get.vars.sum.and.syn <- function(vars) {
  c(vars.sum = sum(vars), vars.max = sum(sqrt(diag(vars)))^2, syn = sum(vars) / sum(sqrt(diag(vars)))^2)
}







gen.random.matrix.constant.colsums <- function(n) {
  M = matrix(runif(n * n, min = 0, max = 1), nrow = n, ncol = n)
  M = M / rowSums(M)
  M = t(M)
}

gen.random.matrix.constant.colsums.2 <- function(n) {
  M = matrix(runif(n * n), nrow = n, ncol = n)
  diag(M) = 0
  M = M / rowSums(M)
  M = t(M)
}

gen.random.matrix <- function(n) {
  M = matrix(runif(n * n), nrow = n, ncol = n)
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
