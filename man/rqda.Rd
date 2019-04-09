\name{rqda}
\alias{rqda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust (quadratic) discriminant analysis
}
\description{
Robust (quadratic) discriminant analysis implements a discriminant analysis method which is robust to label noise. This function implements the method described in Lawrence and Scholkopf (2003, ISBN:1-55860-778-1).
}
\usage{
rqda(X,lbl,Y,maxit=50,disp=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a data frame containing the learning observations.
}
  \item{lbl}{
the class labels of the learning observations.
}
  \item{Y}{
a data frame containing the new observations to classify.
}
  \item{maxit}{
the maximum number of iterations.
}
  \item{disp}{
logical, if \code{TRUE}, several plots are displayed.
}
  \item{...}{
additional arguments to provide to subfunctions.
}
}
\value{
A list is returned with the following elements:
\item{nu}{the estimated class proportions.}
\item{mu}{the estimated class means.}
\item{S}{the estimated covariance matrices.}
\item{gamma}{the estimated purity level of the labels.}
\item{Ti}{the posterior probabilties of the labels knowing the observed labels for the learning observations.}
\item{Pi}{the class posterior probabilities of the observations to classify.}
\item{cls}{the class assignments of the observations to classify.}
\item{ll}{the log-likelihood value.}
}
\references{
Lawrence, N., and Scholkopf, B., Estimating a kernel Fisher discriminant in the presence of label noise, Pages 306–313 of: Proceedings of the Eighteenth International Conference on Machine Learning. ICML’01. San Francisco, CA, USA, 2001 (ISBN:1-55860-778-1).
}
\author{
C. Bouveyron
}
\examples{
n = 50
m1 = c(0,0); m2 = 1.5*c(1,-1)
S1 = 0.1*diag(2); S2 = 0.25 * diag(2)
X = rbind(mvrnorm(n,m1,S1),mvrnorm(2*n,m2,S2))
cls = rep(1:2,c(n,2*n))

# Label perturbation
ind = rbinom(3*n,1,0.4); lb = cls
lb[ind==1 & cls==1] = 2
lb[ind==1 & cls==2] = 1

# Classification with RQDA
res = rqda(X,lb,X)
table(cls,res$cls)
}