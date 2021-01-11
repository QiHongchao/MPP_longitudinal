N <- nrow(historical)  ##number of observations
K <- 2  ##number of fixed effects, there is no interaction term in step 1
J <- nlevels(as.factor(historical[, ncol(historical)]))  ##id is in the last column, number of subject
L <- 2  ##number of random effects
x <- as.matrix(historical[, 1:K])  ##design matrix of fixed effects
z <- as.matrix(historical[, 1:L])  ##design matrix of random effects
y <- historical[, ncol(historical)-1]  ##response was put in the second last column
sub_index <- sub_index_specification(data = historical)  
var_beta <- N * solve(t(as.matrix(historical[,c("int", "time")]))%*%
                        as.matrix(historical[,c("int", "time")]))