constrEM <-
function(X,K,C,maxit=30){
  
  f <- function(X,mu,S){
    if (is.vector(X)) X = t(X) 
    n = nrow(X); p = ncol(X)
    if (n> 1) Xbar = X - matrix(1,n,1) %*% mu
    else Xbar = X - mu
    D = 1/(2*pi)^(p/2) * 1 / det(S)^(1/2) *
      exp(-1/2*diag(Xbar%*%solve(S)%*%t(Xbar)))  
  }
  
  # Initialization
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  C = as.numeric(as.factor(C))
  L = length(unique(C))
  ll = c(-Inf)
  T = matrix(NA,L,K)
  prop = rep(1/K,K)
  mu = matrix(rnorm(K*p,0,2),ncol=p)
  S = array(NA,c(K,p,p))
  for (k in 1:K) S[k,,] = diag(p)
  
  # Global statistics
  Cmeans = matrix(NA,L,p)
  Cn = rep(NA,L)
  for (l in 1:L){
    Cn[l] = sum(C==l)
    if (Cn[l] > 1) Cmeans[l,] = colMeans(X[C==l,]) 
    else Cmeans[l,] = X[C==l,]
  }
  
  # EM loop
  for (it in 1:maxit){
    cat('.')
    # E step
    for (l in 1:L){
      for (k in 1:K){
        T[l,k] = prop[k] * prod(f(X[C=l,],mu[k,],S[k,,]))
      }
    }
    T = T / rowSums(T) %*% matrix(1,1,K)
    
    # M step
    for (k in 1:K){
      prop[k] = sum(T[,k]) / L
      mu[k,] = colSums(T[,k]*Cmeans*Cn) / sum(T[,k]*Cn)
      Sk = matrix(0,p,p)
      for (l in 1:L){
        XX = X[C==l,] - matrix(1,Cn[l],1) %*% mu[k,]
        Sk = Sk + T[l,k] * t(XX)%*%XX
      }
      S[k,,] = Sk / sum(T[,k]*Cn)
    }
  }
  cat('\n')
  
  # Results
  cls = rep(NA,n)
  for (l in 1:L) cls[C==l] = which.max(T[l,])
  res = list(cls=cls,T=T,prop=prop,mu=mu,S=S,ll=ll)
  return(res)
}
