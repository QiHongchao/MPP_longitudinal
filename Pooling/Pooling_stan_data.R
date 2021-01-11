N <- nrow(pooling_data)  ##number of observations of pooled data
K <- 3  ##number of fixed effects
J <- nlevels(as.factor(pooling_data[, ncol(pooling_data)]))  ##id is in the last column, number of subject
L <- 2  ##number of random effects
x <- as.matrix(pooling_data[, 1:K])  ##design matrix of fixed effects, pooled data
z <- as.matrix(pooling_data[, 1:L])  ##design matrix of random effects, pooled data
y <- pooling_data[, ncol(pooling_data)-1] 
sub_index <- sub_index_specification(data = pooling_data)  ##subject index for pooled data
var_beta <- N * solve(t(as.matrix(pooling_data[,c("int", "time", "interaction")]))%*%
                        as.matrix(pooling_data[,c("int", "time", "interaction")]))
