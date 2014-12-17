library(doMC)  # 
registerDoMC(20)  # register Multi Cores
getDoParWorkers()  # get available Cores


## Generate random community matrix with constant column sums, that reflect the constant effort hypothesis
n = 30
k = 5
r = rep(1, n)

## generate random networks with different degree heterogeneitise
graphs = get.graphs.hetero(n = n, k = k, rep1 = 10, rep2 = 1)

## get degree heterogeneities of graphs
heteros = ldply(graphs, function(graph) c(hetero = get.degree.hetero(graph)))
hist(heteros$hetero, breaks = 50)

## generate random networks with different modularity
#g = 2  # number of groups
#qs = seq(0, 2, 0.1)  # modularity level, measured by the propotion between intra-group edges and inter-group edges
#qs = rep(qs, 5)  # repeat several times?

#graphs.modularity = llply(qs, function(q) get.graph.modularity(n, g, k, q))

rhos = seq(0, 0.9, 0.3)
beta0 = - 1
beta1.max = 1 / sqrt(n * k)
beta1.means = seq(0.05, 0.2, 0.05)
isfixed.colsums = FALSE

stab = ldply(rhos, function(rho) {
  C = get.env.flucs(n, rho)
  ldply(beta1.means, function(beta1.mean) {
    beta1.sds = seq(0, beta1.mean, length.out = 2)
    get.stabilities.by.graphs(graphs, beta0 = beta0,  beta1.mean = beta1.mean, beta1.sds = beta1.sds,
                              inter.type = 'competition', C, isfixed.colsums = isfixed.colsums, r)
  })
})
graphs.length = length(graphs)
stab$beta1.means = rep(beta1.means, times = length(rhos), each = 2 * graphs.length)
stab$rhos = rep(rhos, each = 2 * graphs.length * length(beta1.means))

ldply(graphs, function(graph) {
  M = get.community.matrix(graph, beta0 = beta0, beta1.min = beta1.min, beta1.max = beta1.max, beta1.random.type = 'unif',
                           inter.type = inter.type, isfixed.colsums = isfixed.colsums)
  c(rankMatrix(M), c(Null(-solve(M) - diag(1, 30))) )
})

#stab = get.stabilities.by.graphs(graphs, beta0 = beta0,  beta1.mean = beta1.mean, beta1.sds = beta1.sds,
#                                           inter.type = 'competition', C, isfixed.colsums = isfixed.colsums, r)
stab = stab.notfixed.rho10
stab = stab %.% filter(beta1.sd == 0. & lev < -1e-2) %.% select(vars.sum, vars.self.sum, vars.max, syn, hetero)
title = 'mutualism beta1.sd = 0.1 and notfixed and rho = 0.'
pairs(stab, upper.panel = panel.smooth2, lower.panel = panel.cor, diag.panel = panel.hist, main = title)

## for modularity
stab = stab.modular
stab$qs = rep(qs, 2)
stab = stab %.% filter(beta1.sd == 0. & lev < -1e-2) %.% select(vars.sum, vars.self.sum, vars.max, syn, hetero, qs, beta1.sd)

graph = graph.connected(n, k, gtype = 'er')
M = get.community.matrix(graph)
plot(eigen(M)$values)

tmp = ldply(1:48, .parallel = TRUE, function(i) {
  graph = get.graph.connected(n, k, gtype = 'regular')
  M = get.community.matrix(graph, beta0 = -2, beta1.max = 1, beta1.random.type = 'unif',
                           inter.type = 'predatorprey', isfixed.colsums = FALSE)
  get.stability(M, C)
})






rows = 1000
cols = 1000
M = array(rep(1, rows * cols), dim = c(rows, cols))
M = diag(1, rows)

Y = crunif(rows, cols, M)
plot(eigen(Y)$values)
hist(colSums(Y))
