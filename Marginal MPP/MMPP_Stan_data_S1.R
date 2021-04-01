N <- nrow(historical)  ##number of observations
K <- 2  ##number of fixed effects, in historical control arm, there is only intercept and time effect
J <- nlevels(as.factor(historical$id))  ##number of subject
L <- 2  ##number of random effects
x <- model.matrix(~ time, data = historical)  ##design matrix of fixed effects
z <- model.matrix(~ time, data = historical)  ##design matrix of random effects
y <- historical$response  ##response
sub_index <- sub_index_specification(data = historical)  
var_beta <- N * solve(t(x)%*%x)