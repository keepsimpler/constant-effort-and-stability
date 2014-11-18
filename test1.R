## Generate random community matrix with constant column sums, that reflect the constant effort hypothesis
n = 50
k = 3
C = diag(1, n)
graph = graph.connected(n, k, gtype = 'er')
M = get.community.matrix(graph)
plot(eigen(M)$values)


tmp = ldply(1:48, .parallel = TRUE, function(i) {
  graph = graph.connected(n, k, gtype = 'regular')
  M = get.community.matrix(graph, beta0 = -2)
  lev = max(Re(eigen(M)$values))
  vars = mou.vars(M, C)
  c(lev = lev, vars.sum = sum(vars), vars.max = sum(sqrt(diag(vars)))^2, syn = sum(vars) / sum(sqrt(diag(vars)))^2)
})










rows = 1000
cols = 1000
M = array(rep(1, rows * cols), dim = c(rows, cols))
M = diag(1, rows)

Y = crunif(rows, cols, M)
plot(eigen(Y)$values)
hist(colSums(Y))