# EXAMPLE 1
  cf <- function(t) exp(-t^2/2)
  options<-list()
  options$isPlot <- FALSE
  options$isInterp <- TRUE
  result   <- cf2DistGP(cf,c(),c(),options)
  xGiven   <- result$x
  pdfGiven <- result$pdf
  pdf <- function(x) InterpPDF(x,xGiven,pdfGiven)
  x   <- seq(-10,10,length.out=1001)
  plot (x,pdf(x),"l")

# EXAMPLE 2
  cf <- function(t) exp(-t^2/2)
  options <- list()
  options$isPlot <- FALSE
  options$isInterp <- TRUE
  result <- cf2DistGP(cf = cf, options = options)
  pdf    <- function(x) InterpPDF(x,result)
  x      <- seq(-10,10,length.out=1001)
  plot(x,pdf(x),"l")
