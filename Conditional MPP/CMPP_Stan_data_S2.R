N <- nrow(current)  ##number of observations in current dataset
N0 <- nrow(historical)  ##number of observations in historical dataset
K <- 3  ##number of fixed effects in current data
K0 <- 2  ##number of fixed effects in historical data
J <- nlevels(as.factor(current$id))  ##number of subject, current
J0 <- nlevels(as.factor(historical$id))  ##number of subject, historical
L <- 2  ##number of random effects
x <- model.matrix(~ time + time:group, data = current)  ##design matrix of fixed effects, current
z <- model.matrix(~ time, data = current)  ##design matrix of random effects, current
y <- current$response  ##response, current
x0 <- model.matrix(~ time, data = historical)  ##design matrix of fixed effects, historical
z0 <- model.matrix(~ time, data = historical)  ##design matrix of random effects, historical
y0 <- historical$response  ##response, historical
num_power <- num_power
logsc <- logsc[,c(1,4)]  ##log scaling constant
sub_index0 <- sub_index_specification(data = historical)  ##subject index for historical data
sub_index <- sub_index_specification(data = current)  ##subject index for current data
Nobs <- 6
var_betamutual <- N * solve(t(x)%*%x)[1:K0, 1:K0]
var_betatrt <- N * solve(t(x)%*%x)[K, K]