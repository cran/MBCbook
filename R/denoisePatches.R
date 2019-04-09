denoisePatches <- function(Y,out,P,sigma=10){
  Y = as.matrix(Y)
  p = ncol(Y); G = out@bestResult@nbCluster
  X = matrix(0,nrow(Y),ncol(Y))
  mu = out@bestResult@parameters@mean
  for (g in 1:G){
    Sg = out@bestResult@parameters@variance[[g]]
    SSg = Sg %*% solve(Sg + sigma^2*diag(p))
    X = X + P[,g] * (matrix(1,nrow(Y),1) %*% mu[g,] + (Y - matrix(1,nrow(Y),1) %*% mu[g,]) %*% SSg)
  }
  X
}