\name{imageToPatch}
\alias{imageToPatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform an image into a collection of patches
}
\description{
Transform an image into a collection of small images (patches) that cover the original image.
}
\usage{
imageToPatch(Im,f)
}
\arguments{
  \item{Im}{
the image for which one wants to extract local patches.
}
  \item{f}{
the size of the desired patches (fxf).
}
}
\value{
A data frame of all extracted patches is returned.
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
