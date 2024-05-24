InterpPDF<-function (xNew,xGiven,pdfGiven){
        if (is.list(xGiven)) {
                result <- xGiven
                xGiven <- result$x
                pdfGiven <- result$pdf
        }
        else{
                errorCondition('Missing Inputs')
        }
        szx <- dim(xNew)
        xNew <- xNew[]
        id <- (xNew >= min(xGiven) & xNew <= max(xGiven))
        pdf<-rep(0,length(xNew))
       #pdf1<- interpBarycentric(xGiven,pdfGiven,xNew[id])


        if (any(id)) {
                pdf[id] <- interpBarycentric(xGiven,pdfGiven,xNew[id])[[2]]

        }
        pdf <- pmax(0,pdf)
        dim(pdf)<-szx

        return(pdf)
}
