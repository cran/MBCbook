\name{wine27}
\alias{wine27}
\docType{data}
\title{
The (27-dimensional) Italian Wine data set
}
\description{
The (27-dimensional) Italian Wine data set is the result of a chemical analysis of 178 wines grown in the same region in Italy but derived from three different cultivars. The analysis determined the quantities of 27 constituents found in each of the three types of wines.
}
\usage{data("wine27")}
\format{
A data frame with 178 observations on the following 29 variables.
  \describe{
    \item{\code{Alcohol}}{a numeric vector}
    \item{\code{Sugar.free_extract}}{a numeric vector}
    \item{\code{Fixed_acidity}}{a numeric vector}
    \item{\code{Tartaric_acid}}{a numeric vector}
    \item{\code{Malic_acid}}{a numeric vector}
    \item{\code{Uronic_acids}}{a numeric vector}
    \item{\code{pH}}{a numeric vector}
    \item{\code{Ash}}{a numeric vector}
    \item{\code{Alcalinity_of_ash}}{a numeric vector}
    \item{\code{Potassium}}{a numeric vector}
    \item{\code{Calcium}}{a numeric vector}
    \item{\code{Magnesium}}{a numeric vector}
    \item{\code{Phosphate}}{a numeric vector}
    \item{\code{Chloride}}{a numeric vector}
    \item{\code{Total_phenols}}{a numeric vector}
    \item{\code{Flavanoids}}{a numeric vector}
    \item{\code{Nonflavanoid_phenols}}{a numeric vector}
    \item{\code{Proanthocyanins}}{a numeric vector}
    \item{\code{Color_Intensity}}{a numeric vector}
    \item{\code{Hue}}{a numeric vector}
    \item{\code{OD280.OD315_of_diluted_wines}}{a numeric vector}
    \item{\code{OD280.OD315_of_flavanoids}}{a numeric vector}
    \item{\code{Glycerol}}{a numeric vector}
    \item{\code{X2.3.butanediol}}{a numeric vector}
    \item{\code{Total_nitrogen}}{a numeric vector}
    \item{\code{Proline}}{a numeric vector}
    \item{\code{Methanol}}{a numeric vector}
    \item{\code{Type}}{a factor with levels \code{Barbera}, \code{Barolo}, \code{Grignolino}}
    \item{\code{Year}}{a numeric vector}
  }
}
\details{
This data set is an expended version of the popular one from the UCI machine learning repository (http://archive.ics.uci.edu/ml/datasets/Wine).
}
\examples{
data(wine27)
}
\keyword{datasets}
