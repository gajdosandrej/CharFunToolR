# EXAMPLE 1:
## CF of Bartlett distribution with the specified degrees of freedom
df    <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3)
t     <- seq(-1,1,length.out=201)
plotReIm(function(t) cfTest_Bartlett(t,df),t,
title='CF of the Bartlett distribution with specified dfs')

# EXAMPLE 2:
#PDF/CDF of Bartlett distribution with the specified degrees of freedom
df    <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3)
cf    <- function(t) cfTest_Bartlett(t,df)
x     <- seq(0,40,length.out=200)
prob  <- c(0.9, 0.95, 0.99)
options<-list()
options$xMin <- 0
result <- cf2DistGP(cf,x,prob,options)
result


# EXAMPLE 3:
#Quantiles of Bartlett distribution computed by the algorithm cf2DistGP
k     <- 13
df    <-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 19, 29, 49, 99)
prob  <- c(0.9, 0.95, 0.99)
Table<-matrix(data = 0,nrow = 13,ncol = 3)
options<-list()
options$isPlot <- false
for (i in 1:13){
print(df[i])
nu <- df[i]*rep(1,k)
cf <- function(t) cfTest_Bartlett(t,nu)
result <-cf2DistGP(cf,c(),prob,options)
print(result$qf)
Table[i,] <- result$qf
}
print(Table)

# EXAMPLE 4:
#Quantiles of Bartlett distribution computed by the algorithm cf2QF
#k     <- 5
#df    <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 19, 29, 49, 99)
#prob  <- c(0.9, 0.95, 0.99)
#Table <- matrix(data = 0,nrow = 13,ncol = 3)
#options<-list()
#options$crit = 1e-10
#options$maxiter = 10
#for (i in 1:13){
#print(df[i])
#nu <- df[i]*rep(1,k)
#cf <- function(t) cfTest_Bartlett(t,nu)
#result<- cf2DistGP(cf,c(),prob,options)
#print(result$qf)
#Table[i,] <- result$qf
#}
#print(Table)


# Computing the exact distribution of the Bartlett's test statistic
# EXAMPLE 5:
k         <- 15
nu_1      <- c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3)
nu        <- sum(nu_1)
alpha     <- list()
alpha[[1]]  <- nu_1/2
weight    <- list()
weight[[1]] <- alpha[[1]]/sum(alpha[[1]])
c         <- nu * log(k * prod(weight[[1]] ^ weight[[1]]))
b         <- 1 + 1 / (3 * (k - 1)) * (sum(1 / nu_1) - 1 / nu)
shift     <- c / b
coef      <- -nu / b

# Characteristic function of the Bartlett's test statistic
cf_logR <- function(t) cf_LogRV_MeansRatioW(t,k,alpha,weight,coef)
cf      <- function(t) exp(1i * t * shift) * cf_logR(t)

# Evaluate the distribution function by using the algorithm cf2DistGP
t     <- seq(-1,1,length.out=201)
x       <- seq(0,40)
prob    <- c(0.9, 0.95, 0.99)
options<-list()
options$xMin <- 0
plotReIm(cf,t,title='CF of the Bartlett distribution with specified dfs')
result  <- cf2DistGP(cf,x,prob,options)
result$qf


# EXAMPLE 6
k<-3
nu_1      <- c(9,9,9)
nu        <- sum(nu_1)
alpha     <- list()
alpha[[1]]  <- nu_1/2
weight    <- list()
weight[[1]] <- alpha[[1]]/sum(alpha[[1]])
c         <- nu * log(k * prod(weight[[1]] ^ weight[[1]]))
b         <- 1 + 1 / (3 * (k - 1)) * (sum(1 / nu_1) - 1 / nu)
shift     <- c / b
coef      <- -nu / b
cf_logR <- function(t) cf_LogRV_MeansRatioW(t,k,alpha,weight,coef)
cf      <- function(t) exp(1i * t * shift) * cf_logR(t)

# Evaluate the distribution function by using the algorithm cf2DistGP
t     <- seq(-1,1,length.out=201)
x       <- seq(0,40)
prob    <- c(0.9, 0.95, 0.99)
options<-list()
options$xMin <- 0
plotReIm(cf,t,title='CF of the Bartlett distribution with specified dfs')
result  <- cf2DistGP(cf,x,prob,options)
result$qf

#example 7
k     <- 3
df    <- c(9,9,9)
prob  <- 0.95
#Table <- matrix(data = 0,nrow = 3,ncol = 1)
options<-list()
options$crit = 1e-10
options$maxiter = 10

        #print(df[i])
        cf <- function(t) cfTest_Bartlett(t,df)
        result<- cf2DistGP(cf,c(),prob,options)

        print(result$qf)


#Example 8
group1 <- c(85, 86, 88, 75, 78, 94, 98, 79, 71, 80)
group2 <- c(91, 92, 93, 85, 87, 84, 82, 88, 95, 96)
group3 <- c(79, 78, 88, 94, 92, 85, 83, 85, 82, 81)

data_list <- list(Group1 = group1, Group2 = group2, Group3 = group3)
bartlett_result <- bartlett.test(data_list)

print(bartlett_result)

#Example 9
k<-3
nu_1<-c(9,9,9)
nu<-sum(nu_1)
n<-30
group1 <- c(85, 86, 88, 75, 78, 94, 98, 79, 71, 80)
group2 <- c(91, 92, 93, 85, 87, 84, 82, 88, 95, 96)
group3 <- c(79, 78, 88, 94, 92, 85, 83, 85, 82, 81)

s_1<-c(var(group1),var(group2),var(group2))
s<-sum(s_1)
c         <- 1 + (1 /( 3 * (k - 1))) * (sum(1 / nu_1) - (1 / nu))
c
T<-(n-k)*log(s^2)-sum(nu_1*log(s_1^2))
B<-3.30244
T
T/c
#https://www.statology.org/bartletts-test/
#https://en.wikipedia.org/wiki/Bartlett%27s_test


#skuska
k     <- 5
df    <- c(7, 9, 4, 11, 12)
prob  <- c(0.9, 0.95, 0.99)
options<-list()
options$isPlot <- FALSE
options$xMin <- 0
result <- cf2DistGP(cf,x,prob,options)
result$qf

