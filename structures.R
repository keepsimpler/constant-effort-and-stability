library(vegan)
library(bipartite)
library(igraph)
library(plyr)

#' @title random rewiring algorithm that change degree distribution
#'        let node with richer neighbors get richer and increase degree heterogeneity
#' @param degrees, the old degree distribution
#' @param ntried, maximum tried times
#' @return the new degree distribution which is more heterogeneous,
#'         or, NULL if maximum tried times approach which new distribution still not be found
degree.richer.richer <- function(degrees, ntried = 1000) {
  n = length(degrees)  # node number
  flag = FALSE  # is link rewiring success?
  # randomly sample two different nodes whose neighbors are larger than 1 and less than n-1
  for (i in 1:ntried) {
    two.rand = sample(1:n, size = 2, replace = FALSE)
    if ( all(degrees[two.rand] > 1) & all(degrees[two.rand] < n - 1) ) {
      flag = TRUE
      break
    }
  }
  if (flag == TRUE) {  # if two different nodes are sampled
    a = two.rand[1]
    b = two.rand[2]
    # add one neighbor to the node which has more neighbors, and
    # substract one neighbor from the node which has less neighbors
    if (degrees[a] >= degrees[b]) { 
      degrees[a] = degrees[a] + 1
      degrees[b] = degrees[b] - 1
    }
    else if (degrees[b] > degrees[a]) {
      degrees[b] = degrees[b] + 1
      degrees[a] = degrees[a] - 1
    }
    return(degrees)
  }
  else {
    warning('Max tried reached, but still can not rewire link to nodes with more neighbors!')
    return(NULL)
  }
}

#' @title generate graphs with different heterogeneity
#' @param n, number of nodes
#' @param k, average degree
#' @param ntried, maximum tried times for link rewiring algorithm, which be transfered to function [degree.richer.richer]
#' @param rep1, num of rounds, one round is the progress from the least heterogeneity (regular degree) to the max-possible hetero.  
#' @param rep2, number of random graphs generated for the same degree sequence
get.graphs.hetero <- function(n, k, ntried = 1000, rep1 = 10, rep2 = 10) {
  ret = list()
  for (i in 1:rep1) {
    degrees = rep(k, n)  # initialize a regular network
    repeat {
      degrees = degree.richer.richer(degrees, ntried)  # get new degree distribution by link rewiring
      if (is.null(degrees)) break  # if new degree distribution cannot be found, then exit loops
      if (is.graphical.degree.sequence(degrees)) {  # if new degree sequence can be realized in a simple graph
        for (i in 1:rep2) { # generate [rep2] random graphs according to the new degree sequence
          G = degree.sequence.game(out.deg = degrees, method = 'vl')    
          ret[[length(ret) + 1]] = G        
        }
      }
    }    
  }
  ret
}

#' @title get degree heterogeneity of a simple undirected graph
#' @param graph, can be a igraph object or an adjacency matrix
get.degree.hetero <- function(graph) {
  if (class(graph) == 'igraph') {  # if is a igraph object of package <igraph> 
    degrees = degree(graph)
  }
  else {
    degrees = rowSums(graph)    
  }
  degrees = sort(degrees, decreasing = TRUE)
  edges = sum(degrees)
  degrees = degrees / edges
  hetero = - sum( degrees * log(degrees) )
  hetero
}


#' @title generate random graphs with a level of modularity 
#' @param n, number of vertices
#' @param g, number of groups
#' @param k, average degree
#' @param q, modularity level
get.graph.modularity <- function(n, g, k, q) {
  k1 = g * g * k / (g + g * (g - 1) * q)
  n1 = n / g
  M = NULL
  for (i in 1:g) {
    M1 = NULL
    for (j in 1:g) {
      if (i == j) 
        M2 = as.matrix(get.adjacency(erdos.renyi.game(n1, n1 * k1 / 2, type = 'gnm')))
      else
        M2 = as.matrix(get.adjacency(erdos.renyi.game(n1, n1 * k1 * q / 2, type = 'gnm')))
      M1 = cbind(M1, M2)
    }
    M = rbind(M, M1)
  }
  M
}




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
#' @param isfixed.colsums, need the column sums are identical with each other?
get.community.matrix <- function(graph, beta0 = 0, beta1.min = 0, beta1.max = 1, beta1.random.type = 'norm',
                                 inter.type, isfixed.colsums = FALSE) {
  if (class(graph) == 'igraph') {  # if is a igraph object of package <igraph> 
    graph = as.matrix(get.adjacency(graph))
  }
  s = dim(graph)[1]  # number of nodes
  edges = sum(graph > 0)  # number of edges
  M = graph
  if (inter.type == 'arbitrary') {  # random matrix whose elements with normal or uniform distribution with zero mean
    if (beta1.random.type == 'norm')
      M[M > 0] = rnorm(edges, mean = 0, sd = beta1.max)
    else if (beta1.random.type == 'unif')
      M[M > 0] = runif(edges, min = - beta1.max, max = beta1.max)
  }
  else if (inter.type == 'competition') {  # random matrix whose elements with negative half-normal or uniform distribution
    if (beta1.random.type == 'norm')
      M[M > 0] = - abs(rnorm(edges, mean = 0, sd = beta1.max))
    else if (beta1.random.type == 'unif')
      M[M > 0] = runif(edges, min = - beta1.max, max = - beta1.min)
  }
  else if (inter.type == 'mutualism') {  # random matrix whose elements with positive half-normal or uniform distribution
    if (beta1.random.type == 'norm')
      M[M > 0] = abs(rnorm(edges, mean = 0, sd = beta1.max))
    else if (beta1.random.type == 'unif')
      M[M > 0] = runif(edges, min = beta1.min, max = beta1.max)
  }
  else if (inter.type == 'predatorprey') {
    for (i in 1 : (s - 1)) {
      for (j in (i + 1) : s) {
        if (M[i, j] > 0) {  # if interaction exist
          if (sample(c(0, 1), 1) == 0) { # half probability
            if (beta1.random.type == 'norm') {
              M[i, j] = abs(rnorm(1, mean = 0, sd = beta1.max))
              M[j, i] = - abs(rnorm(1, mean = 0, sd = beta1.max))              
            }
            else if (beta1.random.type == 'unif') {
              M[i, j] = runif(1, min = beta1.min, max = beta1.max)
              M[j, i] = runif(1, min = - beta1.max, max = - beta1.min)
            }
          }
          else {
            if (beta1.random.type == 'norm') {
              M[i, j] = - abs(rnorm(1, mean = 0, sd = beta1.max))
              M[j, i] = abs(rnorm(1, mean = 0, sd = beta1.max))                          
            }
            else if (beta1.random.type == 'unif') {
              M[i, j] = runif(1, min = - beta1.max, max = - beta1.min)
              M[j, i] = runif(1, min = beta1.min, max = beta1.max)
            }
          }
        }
      }
    }
  }  # end of predatorprey
  
  if (isfixed.colsums == TRUE) {  # if the hypothesis of constant interaction effort is assumed
    M = M / rowSums(M)
    M = t(M)
  }
  
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

#' @title get different stability measurements based on the community matrix and the environmental fluctuating covariance matrix
#' @param phi, community matrix
#' @param C, the covariance matrix of environmental fluctuating
get.stability <- function(phi, C, r) {
  n = dim(phi)[1]
  eigs = eigen(phi)$values
  lev = max(Re(eigs))  # largest eigenvalue
  sev = min(Re(eigs))  # smallest eigenvalue
  sndev = sort(Re(eigs))[n - 1]  # second largest eigenvalue
  syn.eig = sev / sndev  # spectral gap between the second and smallest eigenvalue as sign of synchronization
  
  vars = mou.vars(phi, C)  # the variance-covariance matrix
  vars.sum = sum(vars)  # temporal variability of the total community
  vars.self.sum = sum(diag(vars))
  vars.max = sum(sqrt(diag(vars)))^2  # the maximum variance of the total community when all species are perfectly synchronized with each other 
  syn = vars.sum / vars.max  # measure of synchronization among all species 
  
  sensitivity = - solve(phi)  # sensitivity matrix is the negative inverse of community matrix
  sens.sum = sum(sensitivity)  
  sens.self = diag(sensitivity)
  sens.self.sum = sum(sens.self)
  sens.other = diag(1/sens.self) %*% sensitivity
  sens.other.sum = sum(sens.other)
  sens.det = det(sensitivity)
  
  nstar = get.nstar.from.community.matrix(phi, r)
  
  c(lev = lev, sev = sev, sndev = sndev, syn.eig = syn.eig, vars.sum = vars.sum, vars.self.sum = vars.self.sum, vars.max = vars.max, syn = syn,
    sens.self.sum = sens.self.sum, sens.other.sum = sens.other.sum, sens.sum = sens.sum, sens.det = sens.det, nstar = nstar)
}


## get stabilities by graphs which have different heteros and by variance of interaction strengths
#' @param C, error variance-covariance matrix reflecting environmental fluctuations
get.stabilities.by.graphs <- function(graphs, beta0 = -2, beta1.mean = 1, beta1.sds = NULL,
                                      inter.type = 'competition', C, isfixed.colsums, r) {
  if (is.null(beta1.sds)) beta1.sds = seq(0, beta1.mean, beta1.mean / 10)
  ldply(beta1.sds, function(beta1.sd) {
    ldply(graphs, .parallel = TRUE, function(graph) {
      hetero = get.degree.hetero(graph)
      print(paste(hetero, ',', beta1.sd))
      beta1.min = beta1.mean - beta1.sd
      beta1.max = beta1.mean + beta1.sd
      M = get.community.matrix(graph, beta0 = beta0, beta1.min = beta1.min, beta1.max = beta1.max, beta1.random.type = 'unif',
                               inter.type = inter.type, isfixed.colsums = isfixed.colsums)
      if (max(Re(eigen(M)$values)) >= 0) warning('Not stable!')
      c(get.stability(M, C, r), hetero = hetero, beta1.sd = beta1.sd)
    })
  })
}

#' @title get variance-covariance matrix that reflects environmental fluctuations
#' @param n, number of species
#' @param rho, correlations among sensitivities of species to environmental fluctuations
get.env.flucs <- function(n, rho) {
  C = matrix(rho, nrow = n, ncol = n)
  diag(C) = 1
  C
}

#' @title For the LV1 model, the species densities in equilibrium can be derived from the community matrix 
#'        if the intrinsic growth rates are known
get.nstar.from.community.matrix <- function(phi, r) {
  n = length(r)
  I = diag(1, n)
  M = - solve(phi) %*% diag(r)
  M = M - I
  nstar = Null(t(M))
  if (dim(nstar)[2] > 0)
    nstar = c(nstar[, 1])
  else
    nstar = rep(0, n)
  nstar
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
