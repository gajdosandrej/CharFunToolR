## EXAMPLE 1
# CF of a linear combination of K=50 independent log-Beta RVs
coef <- 1 / ((((1:50) - 0.5) * pi) ^ 2)
plot(
        1:50,
        coef,
        xlab = "",
        ylab = "",
        type = "p",
        pch = 20,
        col = "blue",
        cex = 1,
        main = expression('Coefficients of the linear combination of' ~ chi ^ 2 ~ 'RVs with DF=1')
)
lines(1:50, coef, col = "blue")
alpha <- 5 / 2
beta <- 3 / 2
t <- seq(-100, 100, (100 - (-100)) / 200)
plotGraf(function(t)
        cf_LogRV_Beta(t, alpha, beta, coef),
        t,
        title = 'Characteristic function of a linear combination of log-Beta RVs')

#library(ggplot2)
#df <- data.frame(x_axis = t, y1_axis = Re(cf), y2_axis = Im(cf))
#ggplot(df, aes(x = x_axis)) +
#  geom_line(aes(y = y1_axis), colour="blue") +
#  geom_line(aes(y = y2_axis), colour = "red") +
#  ylab(label="") +
#  xlab("") +
#  ggtitle("Characteristic function of a linear combination of log-Beta RVs") +
#  theme(plot.title = element_text(hjust = 0.5))
#   figure; plot(t,real(cf),t,imag(cf)); grid on;
#   title('Characteristic function of a linear combination of log-Beta RVs')
#

##  EXAMPLE 2:
# PDF/CDF from the CF by cf2DistGP
alpha <- 5 / 2
beta  <- 3 / 2
coef <- 1 / ((((1:50) - 0.5) * pi) ^ 2)
cf <- function(t) {
        cf_LogRV_Beta(t, alpha, beta, coef)
}
options <- list()
options$xMax <- 0
result <- cf2DistGP(cf = cf, options = options)

## EXAMPLE 3:
# Distribution of log(R), where R=geometric/arithmetic mean of Gamma RVs
# Let X_1,...,X_n are iid RVs, X_j ~ Gamma(A,B), where A > 0 is the known
# shape parameter and B > 0 is the (unknown, common) rate parameter.
# Let R = geometricmean(X)/mean(X). According to Glaser (JASA 1976)
# log(R) ~ (1/n) * sum_{j=1}^{n-1} log(Y_j), Y_j ~ Beta(alpha,beta_j),
# where alpha = A and beta_j = j/n for j = 1,...,n-1. That is, log(R) is
# distributed as linear combination of independent logBeta random
# variables log(Y_j).
# Here we evaluate the PDF/CDF of W = -log(R) (i.e. minus of log(R))
n <- 10
alpha <- 1 #A = 1, i.e. X_j are from exponential distribution
beta <- (1:n - 1) / n
coef <- -1 / n
t <- seq(from = -25,
         to = 25,
         length.out = 201)
plotGraf(function(t)
        cf_LogRV_Beta(t, alpha, beta, coef),
        t,
        title = 'Characteristic function of a linear combination of log-Beta RVs')
prob <- c(0.9, 0.95, 0.99)
options <- list()
options.xMin = 0
result <- cf2DistGP(cf = cf, prob = prob, options = options)
str(result)

## EXAMPLE 4 (Distribution of Wilks Lambda statistic)
# If E ~ Wp(m,Sigma) and H ~ Wp (n,Sigma) with m >= p, then the
# distribution of Lambda = L = det(E)/det(E + H) is Wilks Lambda
# distribution, denoted by L ~ Lambda(p,m,n), with L in (0,1).
# It holds that L ~ Prod_{i=1}^p B_i, where B_i follow independent
# Beta distributions with Bi ~ B{(m + 1 - i)/2, n/2)}, for i = 1,...,p.
# Let W = -log(L) and cdf_W(x) = Prob(W <= x) = Prob(-sum(log(Bi)) <= x).
# Then, cdf_L(u) = Prob(L <= u) = Prob(W > -log(u)) =
# 1 - Prob(W <= -log(u)) = 1 - cdf_W(x), where x = -log(u) and u =
# exp(-x). Moreover, pdf_Lambda(u) = pdf_W(x)/exp(-x) = pdf_W(-log(u))/u.
# The Lambda statistic is used to test null (significance) hypothesis
# expressed by the matrix H. The null hypothesis is rejected for small
# values of the observed statistic Lambda, or large value -log(Lambda).
p <- 10
m <- 20
n <- 7
i  <- 1:p
alpha <- (m + 1 - i) / 2
beta  <- n / 2
coef  <- -1
t <- seq(from = -5,
         to = 5,
         length.out = 201)
plotGraf(function(t)
        cf_LogRV_Beta(t, alpha, beta, coef),
        t,
        title = 'Characteristic function of a linear combination of log-Beta RVs')
prob <- c(0.99, 0.95, 0.9)
x <- seq(from = 2,
         to = 8,
         length.out = 100)
options <- list()
options.xMin = 0

# Distribution of -log(Lambda)
result <- cf2DistGP(
        cf = cf,
        x = x,
        prob = prob,
        options = options
)
str(result)

# PDF of Lambda
# If W = -log(Lambda), then pdf_Lambda(u) = pdf_W(x)/exp(-x),
# where x = -log(u) and u = exp(-x)
plot.new()
plot.window(xlim = c(0, 0.15), ylim = c(0, 40))
axis(1)
axis(2)
#title(main=paste("PDF of Wilks", expression(Lambda), "distribution with p=10, m = 20, n=7"))
title(main = expression("PDF of Wilks" ~ Lambda ~ "distribution with p=10, m = 20, n=7"))
title(xlab = expression(Lambda))
title(ylab = "PDF")
grid(lty = "solid")
lines(exp(-result$x),
      result$pdf / exp(-result$x),
      col = "blue",
      lwd = 2)

# CDF of Lambda
# If W = -log(Lambda), then cdf_Lambda(u) = 1-cdf_W(x),
# where x = -log(u) and u = exp(-x)
plot.new()
plot.window(xlim = c(0, 0.15), ylim = c(0, 1))
axis(1)
axis(2)
#title(main=paste("PDF of Wilks", expression(Lambda), "distribution with p=10, m = 20, n=7"))
title(main = expression("CDF of Wilks" ~ Lambda ~ "distribution with p=10, m = 20, n=7"))
title(xlab = expression(Lambda))
title(ylab = "PDF")
grid(lty = "solid")
lines(exp(-result$x), 1 - result$cdf, col = "blue", lwd = 2)
print("Quantiles of Wilks Lambda distribution")
print(paste("alpha =", prob))
print(paste("quantile =", exp(-result$qf)))
