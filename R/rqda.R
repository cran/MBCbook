rqda <- function(X,lbl,Y,maxit=50,disp=FALSE,...){
  n = nrow(X)
  p = ncol(X)
  Ti = matrix(NA,n,2)
  S = array(NA,c(2,p,p))
  S[1,,] = S[2,,] = diag(p)
  m = rbind(colMeans(X[lbl==1,]),colMeans(X[lbl==2,]))
  #m = mvrnorm(2,mu=colMeans(X),Sigma=diag(p))
  nu = rep(n/2,2)
  gamma = rep(0.5,2)
  ll = c()
  
  # Learning
  for (i in 1:maxit){
    cat('.')
    #browser()
    # E step
    Ti[,1] = dmvnorm(X,m[1,],S[1,,]) * ((lbl==1)*gamma[1] + (lbl==2)*(1-gamma[1]))
    Ti[,2] = dmvnorm(X,m[2,],S[2,,]) * ((lbl==2)*gamma[2] + (lbl==1)*(1-gamma[2]))
    Ti = Ti / t(matrix(1,2,1) %*% rowSums(Ti))
    if (disp & i <= 30 & (i<=5 | i%%5==0)){
      par(mfrow=c(1,2))
      #barplot(Ti[,1],col='lightblue',ylim=c(0,1),ylab=expression(paste('P(y=1|x,',theta,')')))
      plot(Ti[,1],type='h',col='lightblue',ylim=c(0,1),ylab=expression(paste('P(y=1|x,',theta,')')))
      # plot(X,col=max.col(Ti),pch=(18:19)[max.col(Ti)],main="Estimated labels"); Sys.sleep(0.25)
      prec <- 150; pastel <- .9
      x1 <- seq(-2.1,3.1,length=prec); x2 <- seq(-3.1,2.1,length=prec)
      s <- expand.grid(x1,x2); s <- as.data.frame(s)
      P = matrix(NA,nrow(s),2)
      P[,1] = nu[1] * dmvnorm(s,m[1,],S[1,,])
      P[,2] = nu[2] * dmvnorm(s,m[2,],S[2,,])
      P = P / t(matrix(1,2,1) %*% rowSums(P))
      plot(X,type='n',xlim=c(-2,3),ylim=c(-3,2),xlab='',ylab='')
      points(s,type='p',pch=16,col=c(rgb(1,pastel,pastel),rgb(pastel,pastel,1))[max.col(P)])
      points(X,col=lbl+1,pch=19,xlim=c(-2,3),ylim=c(-3,2))
      for (k in 1:2){
        PP <- P[,k] - apply(P[,-k,drop = FALSE],1, max)
        contour(x1,x2,matrix(PP,prec,prec),level=0,add=1,col="black",lwd=3,drawlabels=0);
      }
    }
    
    # M step
    for (k in 1:2){
      nu[k] = sum(Ti[,k])
      m[k,] = 1/nu[k] * colSums(Ti[,k]%*%matrix(1,1,p) * X)
      Xbark = X - matrix(1,n,1)%*%m[k,]
      S[k,,] = 1/nu[k] * t(Ti[,k]*Xbark) %*% as.matrix(Xbark)
      gamma[k] = 1/nu[k] * sum(Ti[lbl==k,k])
    }
    #print(m)
    
    # Loglikelihood
    ll = c(ll,sum(log(dmvnorm(X[max.col(Ti)==1,],m[1,],S[1,,])))+sum(log(dmvnorm(X[max.col(Ti)==2,],m[2,],S[2,,]))))
    if (i>1) if(abs(ll[i-1]-ll[i])<1e-6) break
  }
  #cat('\n')
  
  # Classification
  Pi = matrix(NA,nrow(Y),2)
  Pi[,1] = nu[1] * dmvnorm(Y,m[1,],S[1,,])
  Pi[,2] = nu[2] * dmvnorm(Y,m[2,],S[2,,])
  #browser()
  Pi = Pi / t(matrix(1,2,1) %*% rowSums(Pi))
  
  # Result
  list(nu=nu,mu=m,S=S,gamma=gamma,Ti=Ti,Pi=Pi,cls=max.col(Pi),ll=ll)
}