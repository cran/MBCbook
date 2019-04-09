\name{PoliticalBlogs}
\alias{PoliticalBlogs}
\docType{data}
\title{
The political blog data set
}
\description{
The political blog data set shows the linking structure in online blogs which commentate on French political issues; the data were collected by Observatoire Presidentielle in October 2006. The data were first used by Latouche et al. (2011) <doi:10.1214/10-AOAS382>.
}
\usage{data("PoliticalBlogs")}
\format{
A large network object, which can be managed with the network library, with 196 nodes.
}
\references{
P. Latouche, E. Birmelé, and C. Ambroise. "Overlapping stochastic block models with application to the French political blogosphere". In : Annals of Applied Statistics 5.1, p. 309-336, 2011 <doi:10.1214/10-AOAS382>.
}
\examples{
data(PoliticalBlogs)
\donttest{
# Visualization with the network library
library(network)
plot(PoliticalBlogs)
}
}
\keyword{datasets}
