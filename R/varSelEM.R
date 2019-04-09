varSelEM <- function(X,G,maxit=100,eps=1e-6){
  # Initialization
  n = nrow(X); p = ncol(X)
  w = matrix(NA,n,G); a = b = c = u = v = array(NA,c(n,G,p))
  mu = mvrnorm(G,mu=colMeans(X),Sigma=diag(p)) # Mean and sd for relevant var.
  #mu = kmeans(X,G)$centers
  sigma = matrix(1,G,p)                     
  lambda = rbind(rep(0,p),rep(1,p))           # Mean and sd for irrelevant var.
  rho = rep(0.5,p); alpha = rep(1/G,G)
  
  cond = TRUE; i = 1; ll = c(-Inf)
  par(mfrow=c(1,5))
  while(cond){ cat('.'); i = i+1
    # E step
    for (g in 1:G){
      for (j in 1:p){
        a[,g,j] = rho[j] * dnorm(X[,j],mu[g,j],sqrt(sigma[g,j]))
        b[,g,j] = (1-rho[j]) * dnorm(X[,j],lambda[1,j],sqrt(lambda[2,j]))
        a[(a<.Machine$double.eps)] = .Machine$double.eps
        b[(b<.Machine$double.eps)] = .Machine$double.eps
        c[,g,j] = a[,g,j] + b[,g,j]
      }
      w[,g] = alpha[g] * apply(c[,g,],1,'prod')
    }
    w = w / rowSums(w)%*%matrix(1,1,G)
    for (g in 1:G){
      for (j in 1:p){
        u[,g,j] = a[,g,j] / c[,g,j] * w[,g]
        v[,g,j] = w[,g] - u[,g,j]
      }
    }
    #print(apply(b,3,'sum'))
    
    # M step
    alpha = colSums(w) / n
    for (g in 1:G){
      mu[g,] = colSums(u[,g,] * X) / colSums(u[,g,])
      sigma[g,] = colSums(u[,g,] * (X-matrix(1,n,1)%*%mu[g,])^2) /  colSums(u[,g,])
    }
    lambda[1,] = apply(apply(v,c(1,3),'sum') * X,2,'sum') / apply(apply(v,c(1,3),'sum'),2,'sum')
    lambda[2,] = apply(apply(v,c(1,3),'sum') * (X-matrix(1,n,1)%*%lambda[1,])^2,2,'sum') / apply(apply(v,c(1,3),'sum'),2,'sum')
    rho = apply(apply(u,c(1,3),'sum'),2,'sum') / n
    sigma[sigma<1e-3] = 1e-3; lambda[2,][lambda[2,]<1e-3] = 1e-3
    if (i %% 2 == 0) barplot(rho,ylim=c(0,1),xlab="Dimensions",ylab="Probability",main=paste('Iteration',i))
    #print(rho)
    #plot(X,col=max.col(w),pch=(16:19)[max.col(w)]); Sys.sleep(0.25)
    
    # Likelihood
    lli = matrix(1,n,G)
    for (g in 1:G){
      for (j in 1:p)
      lli[,g] = lli[,g] * ((rho[j] * dnorm(X[,j],mu[g,j],sigma[g,j]) + (1-rho[j]) * dnorm(X[,j],lambda[1,j],lambda[2,j]))) 
    }
    ll[i] = log(sum(matrix(1,n,1)%*%alpha * lli))
   
    
    # Break condition
    cond = ((ll[i-1] - ll[i] < 0) & abs((ll[i-1] - ll[i]) / ll[i]) > eps) & (i < maxit)
  }
  cat('\n')
  plot(ll,type='b')
  list(mu=mu,sigma=sigma,lambda=lambda,alpha=alpha,rho=rho,P=w,cls=max.col(w),ll=ll)
}