N <- nrow(pooling)  ##number of observations of pooling data
J <- nlevels(as.factor(pooling$id))
x <- as.matrix(pooling[, c("int", "time", "interaction")])  ##design matrix of fixed effects, pooling data
z <- as.matrix(pooling[, c("int", "time")])  ##design matrix of random effects, pooling data
y <- pooling$response
K <- ncol(x) ##number of fixed effects
L <- ncol(z)  ##number of random effects
sub_index <- sub_index_specification(data = pooling)  ##subject index for pooling data