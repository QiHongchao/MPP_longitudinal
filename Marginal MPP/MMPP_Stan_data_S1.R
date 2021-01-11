N <- nrow(historical)  ##number of observations
K <- 2  ##number of fixed effects, in historical control arm, there is only intercept and time effect
J <- nlevels(as.factor(historical[, ncol(historical)]))  ##id is put in the last column, number of subject
L <- 2  ##number of random effects
x <- as.matrix(historical[, 1:K])  ##design matrix of fixed effects
z <- as.matrix(historical[, 1:L])  ##design matrix of random effects
y <- historical[, ncol(historical)-1]  ##response was put in the second last column
sub_index <- sub_index_specification(data = historical)  
var_beta <- N * solve(t(x)%*%x)