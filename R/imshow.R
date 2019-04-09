imshow <- function(x,col=palette(gray(0:255/255)),useRaster = TRUE,...){
  image(t(x)[,nrow(x):1],col=col,useRaster=useRaster,xaxt='n',yaxt='n',...)
}
