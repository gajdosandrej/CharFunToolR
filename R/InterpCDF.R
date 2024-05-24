InterpCDF<-function (x,xGiven,cdfGiven){
        if(is.list(xGiven)){

                result <- xGiven
                xGiven <- sort(result$x)
                cdfGiven <- sort(result$cdf)
        }
        else{

                errorCondition("Missing Inputs")
        }

        szx <- dim(x)
        x <- x[]
        id0 <- x < min(xGiven)
        id1 <- x > max(xGiven)
        id <- (x >= min(xGiven) & x <= max(xGiven))
        cdf <- rep(0,length(x))

        if(any(id0)){
                cdf[id0] <- 0
        }

        if(any(id1)){
                cdf[id1] <- 1
        }

        if(any(id)){
                cdf[id] <- interpBarycentric(xGiven, cdfGiven, x[id])[[2]]
        }

        cdf <- pmax(0,pmin(1,cdf))
        dim(cdf) <- szx

        return(cdf)
}

