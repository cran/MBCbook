\name{velib2D}
\alias{velib2D}
\docType{data}
\title{
The bivariate Vélib data set
}
\description{
The bivariate Vélib data set contains data from the bike sharing system of Paris, called Vélib. The data are loading profiles and percentage of broken docks of the bike stations over one week. The data were collected every hour during the period Sunday 1st Sept. - Sunday 7th Sept., 2014. The data were first used in Bouveyron et al. (2015) <doi:10.1214/15-AOAS861>.
}
\usage{data("velib2D")}
\format{
The format is:

- availableBikes: the loading profiles (nb of available bikes / nb of bike docks) of the 1189 stations at 181 time points.

- brokenDockss: the percentage of broken docks of the 1189 stations at 181 time points.

- position: the longitude and latitude of the 1189 bike stations.

- dates: the download dates.

- bonus: indicates if the station is on a hill (bonus = 1).

- names: the names of the stations.
}
\source{
The real time data are available at https://developer.jcdecaux.com/ (with an api key).
}
\references{
The data were first used in C. Bouveyron, E. Côme and J. Jacques, The discriminative functional mixture model for the analysis of bike sharing systems, The Annals of Applied Statistics, vol. 9 (4), pp. 1726-1760, 2015 <doi:10.1214/15-AOAS861>.
}
\examples{
data(velib2D)
}
\keyword{datasets}
