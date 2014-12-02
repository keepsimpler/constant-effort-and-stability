library(doMC)  # 
registerDoMC(20)  # register Multi Cores
getDoParWorkers()  # get available Cores


## Generate random community matrix with constant column sums, that reflect the constant effort hypothesis
n = 30
k = 3
C = diag(1, n)

## generate random networks with different degree heterogeneitise
graphs = get.graphs.hetero(n = n, k = k, rep1 = 10, rep2 = 1)

## get degree heterogeneities of graphs
heteros = ldply(graphs, function(graph) c(hetero = get.degree.hetero(graph)))
hist(heteros$hetero, breaks = 50)

stabilities = ldply(graphs, .parallel = TRUE, function(graph) {
  print(get.degree.hetero(graph))
  M = get.community.matrix(graph, beta0 = -2, beta1.max = 1, beta1.random.type = 'unif',
                           inter.type = 'predatorprey', isfixed.colsums = TRUE)
  get.stability(M, C)
})

stabilities = cbind(stabilities, heteros)


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
