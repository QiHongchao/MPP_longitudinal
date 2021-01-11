N <- nrow(current)  ##number of observations of current data
K <- 3  ##number of fixed effects
J <- nlevels(as.factor(current[, ncol(current)]))  ##id is in the last column, number of subject
L <- 2  ##number of random effects
x <- as.matrix(current[, 1:K])  ##design matrix of fixed effects, current data
z <- as.matrix(current[, 1:L])  ##design matrix of random effects, current data
y <- current[, ncol(current)-1] 
sub_index <- sub_index_specification(data = current)  ##subject index for current data
var_beta <- N * solve(t(as.matrix(current[,c("int", "time", "interaction")]))%*%
                        as.matrix(current[,c("int", "time", "interaction")]))