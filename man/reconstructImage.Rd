\name{reconstructImage}
\alias{reconstructImage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reconstructing an image from a patch decomposition
}
\description{
A simple way of reconstructing an image from a patch decomposition.
}
\usage{
reconstructImage(X,nl,nc)
}
\arguments{
  \item{X}{
the matrix of patches to be used for reconstructing the image.
}
  \item{nl}{
the number of rows of the image.
}
  \item{nc}{
the number of columns of the image.
}
}
\value{
an image is returned as a matrix object, that can be display with the \code{\link{imshow}} function.
}
\author{
C. Bouveyron & J. Delon
}
\examples{
Im = diag(16) 
ImNoise = Im + rnorm(256,0,0.1)
X = imageToPatch(ImNoise,4)
out = mixmodCluster(X,10,model=mixmodGaussianModel(family=c("spherical")))
res = mixmodPredict(X,out@bestResult)
Xdenoised = denoisePatches(X,out,P = res@proba,sigma = 0.1) 
ImRec = reconstructImage(Xdenoised,16,16)
par(mfrow=c(1,3)); imshow(Im); imshow(ImNoise); imshow(ImRec)
}
