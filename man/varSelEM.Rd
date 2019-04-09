\name{varSelEM}
\alias{varSelEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A variable selection algorithm for clustering
}
\description{
A variable selection algorithm for clustering which implements the method described in Law et al. (2004) <doi:10.1109/TPAMI.2004.71>.
}
\usage{
varSelEM(X,G,maxit=100,eps=1e-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a data frame containing the observations to cluster.
}
  \item{G}{
the expected number of groups (integer).
}
  \item{maxit}{
the maximum number of iterations (integer). The default value is 100.
}
  \item{eps}{
the convergence threshold. The default value is 1e-6.
}
}
\value{
A list is returned with the following elements:
\item{mu}{the group means for relevant variables.}
\item{sigma}{the group variances for relevant variables.}
\item{lambda}{the group means for irrelevant variables}
\item{alpha}{the group variances for irrelevant variables.}
\item{rho}{the feature saliency.}
\item{P}{the group posterior probabilities.}
\item{cls}{the group memberships.}
\item{ll}{the log-likelihood value.}
}
\references{
Law, M. H., Figueiredo, M. A. T., and Jain, A. K., Simultaneous feature selection and clustering using mixture models, IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 26, pp. 1154–1166, 2004 <doi:10.1109/TPAMI.2004.71>.
}
\author{
C. Bouveyron
}
\examples{
data(wine27)
X = scale(wine27[,1:27]) 
cls = wine27$Type

# Clustering and variable selection with VarSelEM
res = varSelEM(X,G=3)

# Clustering table
table(cls,res$cls)
}
