reconstructImage <- function(X,nl,nc){
  X = t(X)
  f = sqrt(dim(X)[1])
  TT = array(NA,c(nl,nc,f^2))
  for (x in 1:f){
    for (y in 1:f){
      i = f*(y-1)+x
      TT[x:(x+nl-f),y:(y+nc-f),i] = matrix(X[i,],nrow=nl-f+1)
    }
  }
  Im = apply(TT,c(1,2),mean,na.rm=TRUE)
}