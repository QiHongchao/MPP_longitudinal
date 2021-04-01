N <- nrow(current)  ##number of observations in current dataset
N0 <- nrow(historical)  ##number of observations in historical dataset
K <- 3  ##number of fixed effects in current data
K0 <- 2  ##number of fixed effects in historical data
J <- nlevels(as.factor(current$id)) ##number of current subjects
J0 <- nlevels(as.factor(historical$id))  ##number of historical subjects
L <- 2  ##number of random effects
x <- model.matrix(~ time + time:group, data = current)  ##design matrix of fixed effects, current
z <- model.matrix(~ time, data = current)  ##design matrix of random effects, current
y <- current$reponse  ##response, current
x0 <- model.matrix(~ time, data = historical)  ##design matrix of fixed effects, historical
z0 <- model.matrix(~ time, data = historical)  ##design matrix of random effects, historical
y0 <- historical$response  ##response, historical
sub_index0 <- sub_index_specification(data = historical)  ##subject index for historical data
sub_index <- sub_index_specification(data = current)  ##subject index for current data
var_beta01 <- N0 * solve(t(x0)%*%x0)
var_betatrt <- N * solve(t(x)%*%x)[K,K]
