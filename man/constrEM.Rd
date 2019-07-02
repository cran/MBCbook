\name{constrEM}
\alias{constrEM}

%%%
\title{
Semi-supervised clustering with must-link constraints
}
\description{
Semi-supervised clustering with must-link constraints allows to cluster data for which must-link constraints are available. This function implements the method described in Shental et al. (2003, ISBN:9781615679119).
}

%%%
\usage{
constrEM(X, K, C, maxit = 30)
}
\arguments{
  \item{X}{
a data frame of observations, assuming the rows are the observations and the columns the variables. Note that NAs are not allowed.
}
  \item{K}{
the number of desired groups.
}
  \item{C}{
a vector encoding the must-link constraints through chuncklets. This vector has to be of the length of the number of observations. Two observations that have to be in the same group must be in the same chuncklet. For instance, the chuncklet vector (1,2,3,4,3,5) indicate that 3rd and the 5th observations have a must-link constraint. If there is no must-link constraints, this vector should be simply 1:nrow(X).
}
  \item{maxit}{
the maximum number of iterations.
}
}

%%%
\value{
A list is returned with the following fields:
\item{cls}{a vector containg the group memberships of the observations.}
\item{T}{the posterior probabilities that the observations belong to the K groups.}
\item{prop}{the estimated mixture proportions.}
\item{mu}{the estimated mixture means.}
\item{S}{the estimated mixture covariance matrices.}
\item{ll}{the log-likelihood value at convergence.}
}

%%%
\references{
This function implements the method described in Shental, N., Bar-Hillel, A., Hertz, T., and Weinshall, D., Computing Gaussian mixture models with EM using equivalence constraints, Proceedings of the 16th International Conference on Neural Information Processing Systems, pages 465–472, 2003 (ISBN:9781615679119).
}
\author{
C. Bouveyron
}
\examples{
# Simulation of some data
set.seed(123)
n = 200
m1 = c(0,0); m2 = 4*c(1,1); m3 = 4*c(1,1)
S1 = diag(2); S2 = rbind(c(1,0),c(0,0.05))
S3 = rbind(c(0.05,0),c(0,1))
X = rbind(mvrnorm(n,m1,S1),mvrnorm(n,m2,S2),mvrnorm(n,m3,S3))
cls = rep(1:3,c(n,n,n))

# Encoding the constraints through chunklets
# Observations 397 and 408 are in the same chunklet
a = 398
b = 430
C = c(1:(b-1),a,b:(nrow(X)-1))

# Clustering with constrEM
res = constrEM(X,K=3,C,maxit=20)
}
