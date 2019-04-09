\name{denoisePatches}
\alias{denoisePatches}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Denoising of image patches
}
\description{
Denoising of image patches based on the clustering of patches.
}
\usage{
denoisePatches(Y,out,P,sigma=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
a data frame containing as rows the image patches to denoise
}
  \item{out}{
the mixmodCluster object that contains mixture parameters
}
  \item{P}{
the posterior probabilities that patches belong to the clusters
}
  \item{sigma}{
the noise standard deviation
}
}
\value{
A data fame of the denoised patches is returned.
}
\note{
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
