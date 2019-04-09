\name{NIR}
\alias{NIR}
\docType{data}
\title{
The chemometrics near-infrared (NIR) data set
}
\description{
The chemometrics near-infrared (NIR) data set has 202 observations and 2801 variables: 2800 near-infrared wavelength measures and 1 class variable. The data were obtained from the analysis of three types of textiles. The data set was first introduce in Devos et al. (2009) <doi:10.1016/j.chemolab.2008.11.005>.
}
\usage{data("velibCount")}
\format{
  A data frame with 202 observations and 2801 variables. The first variable indicates the class-memberships of the observations.
}
\references{
Devos, O., Ruckebusch, C., Durand, A., Duponchel, L., and Huvenne, J.-P., Support vector machines (SVM) in near infrared (NIR) spectroscopy: Focus on parameters optimization and model interpretation, Chemometrics and Intelligent Laboratory Systems, 96, 27–33, 2009 <doi:10.1016/j.chemolab.2008.11.005>.
}
\examples{
data(NIR)
matplot(t(NIR[,-1]),type='l',col=NIR[,1])
}
\keyword{datasets}
