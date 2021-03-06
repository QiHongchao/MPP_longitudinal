N0 <- nrow(historical)  ##number of observations, historical
N <- nrow(current)  ##number of observations, current
K0 <- 2 ##number of fixed effects in historical study, there is no interaction
K <- 3  ##number of fixed effects in current study
J0 <- nlevels(as.factor(historical$id))  ##number of subject, historical
J <- nlevels(as.factor(current$id))  ##number of subject, current
L <- 2  ##number of random effects, we assume same random effects in historical and current data
x0 <- model.matrix(~ time, data = historical)  ##design matrix of fixed effects, historical
z0 <- model.matrix(~ time, data = historical)  ##design matrix of random effects, historical
y0 <- historical$response  ##response, historical
x <- model.matrix(~ time + time:group, data = current)  ##design matrix of fixed effects, current
z <- model.matrix(~ time, data = current)  ##design matrix of random effects, current
y <- current$response  ##response, current
num_power <- num_power  ##number of power values
logsc <- logsc[,c(1,4)]  ##Log scaling constant
sub_index0 <- sub_index_specification(data = historical)  ##subject index for historical data
sub_index <- sub_index_specification(data = current)  ##subject index for current data
var_betamutual <- N * solve(t(x)%*%x)[1:K0, 1:K0]
var_betatrt <- N * solve(t(x)%*%x)[K, K]