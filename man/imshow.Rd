\name{imshow}
\alias{imshow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Display an image
}
\description{
A simple way of displaying an image, using the \code{image} function.
}
\usage{
imshow(x,col=palette(gray(0:255/255)),useRaster = TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the image to display as a matrix.
}
  \item{col}{
the color palette to use when displaying the image.
}
  \item{useRaster}{
logical; if TRUE a bitmap raster is used to plot the image instead of polygons. The grid must be regular in that case, otherwise an error is raised. For the behaviour when this is not specified, see the ‘Details’ section of the \code{image} function.
}
  \item{...}{
additionial arguments to provide to subfunctions.}
}
\seealso{
\code{\link{image}}
}
\examples{
Im = diag(16)
imshow(Im)
}