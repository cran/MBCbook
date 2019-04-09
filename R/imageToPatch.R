imageToPatch <- function(Im,f){
  nl = dim(Im)[1]; nc = dim(Im)[2]
  X = matrix(0,(nl-f+1)*(nc-f+1),f^2)
  for (i in 1:f){
    for (j in 1:f){
      v = Im[i:(i+nl-f),j:(j+nc-f)]
      X[,i+f*(j-1)] = as.vector(v) 
    }
  }
  as.data.frame(X)
}