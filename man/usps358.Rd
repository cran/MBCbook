\name{usps358}
\alias{usps358}
\docType{data}
\title{
The handwritten digits usps358 data set
}
\description{
The handwritten digits usps358 data set is a subset of the famous USPS data from UCI, which contains only the 1 756 images of the digits 3, 5 and 8.
}
\usage{data("usps358")}
\format{
  A data frame with 1756 observations on the following 257 variables: \code{cls} is a numeric vector encoding the class of the digits, \code{V1} to \code{V256} are numeric vectors corresponding to the pixels ot the 8x8 images.
}
\source{
The data set is a subset of the famous USPS data from UCI (https://archive.ics.uci.edu/ml/index.php). The usps358 data set contains only the 1 756 images of the digits 3, 5 and 8 which are the most difficult digits to discriminate.
}
\examples{
data(usps358)
}
\keyword{datasets}
