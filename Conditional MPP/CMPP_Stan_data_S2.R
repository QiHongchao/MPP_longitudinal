N <- nrow(current)  ##number of observations in current dataset
N0 <- nrow(historical)  ##number of observations in historical dataset
K <- 3  ##number of fixed effects in current data
K0 <- 2  ##number of fixed effects in historical data
J <- nlevels(as.factor(current[, ncol(current)]))  ##id is in the last column, number of subject, current
J0 <- nlevels(as.factor(historical[, ncol(historical)]))  ##id is in the last column, number of subject, historical
L <- 2  ##number of random effects
x <- as.matrix(current[, 1:K])  ##design matrix of fixed effects, current
z <- as.matrix(current[, 1:L])  ##design matrix of random effects, current
y <- current[, ncol(current)-1]  ##response was put in the second last column, current
x0 <- as.matrix(historical[, 1:K0])  ##design matrix of fixed effects, historical
z0 <- as.matrix(historical[, 1:L])  ##design matrix of random effects, historical
y0 <- historical[, ncol(historical)-1]  ##response was put in the second last column, historical
num_power <- num_power
logsc <- logsc[,c(1,4)]  ##log scaling constant
sub_index0 <- sub_index_specification(data = historical)  ##subject index for historical data
sub_index <- sub_index_specification(data = current)  ##subject index for current data
Nobs <- 6
var_betamutual <- N * solve(t(as.matrix(current[,c("int", "time", "interaction")]))%*%
                              as.matrix(current[,c("int", "time", "interaction")]))[1:K0, 1:K0]
var_betatrt <- N * solve(t(as.matrix(current[,c("int", "time", "interaction")]))%*%
                           as.matrix(current[,c("int", "time", "interaction")]))[K, K]