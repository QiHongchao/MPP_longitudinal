N <- nrow(historical)  ##number of observations
K <- 2  ##number of fixed effects, there is no interaction term in step 1
J <- nlevels(as.factor(historical$id))  ##number of subject
L <- 2  ##number of random effects
x <- model.matrix(~ time, data = historical)  ##design matrix of fixed effects
z <- model.matrix(~ time, data = historical)  ##design matrix of random effects
y <- historical$response  ##response
sub_index <- sub_index_specification(data = historical)  ##subject index for historical data 
var_beta <- N * solve(t(x)%*%x)