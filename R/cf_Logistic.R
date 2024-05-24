cf_Logistic <- function(t, mu, beta, coef, niid){

        ## CHECK THE INPUT PARAMETERS

        if(missing(niid)){
                niid <- vector()
        }

        if(missing(coef)){
                coef <- vector()
        }

        if(missing(beta)){
                beta <- vector()
        }

        if(missing(mu)){
                mu   <- vector()
        }

        if(length(mu)==0){
                mu   <- 0
        }

        if(length(beta)==0){
                beta <- 1
        }

        if(length(coef)==0){
                coef <- 1
        }

        ## Equal size of the parameters

        l_max <- max(c(length(coef), length(mu), length(beta)))
        if (l_max > 1) {
                if (length(coef) == 1) {
                        coef <- rep(coef, l_max)
                }
                if (length(mu) == 1) {
                        mu   <- rep(mu, l_max)
                }
                if(length(beta) == 1){
                        beta <- rep(beta, l_max)
                }
                if ((any(lengths(list(coef, mu, beta)) < l_max))) {
                        stop("Input size mismatch.")
                }
        }

        ## Characteristic function

        szt <- dim(t)
        t   <- c(t)
        cf  <- matrix()
        cf  <- pi * t %*% t(coef * beta) * exp(1i * t %*% t(coef * mu)) * pracma::csch(pi * t %*% t(coef * beta))
        cf  <- apply(cf, 1, prod)
        dim(cf)  <- szt
        cf[t==0] <- 1
        if(length(niid)!=0){
                if(length(niid)==1){
                        cf <- cf ^ niid
                }
                else{
                        errorCondition("niid should be a scalar (positive integer) value")
                }
        }

        return(cf)
 }
