N <- nrow(pooling_data)  ##number of observations of pooled data
K <- 3  ##number of fixed effects
J <- nlevels(as.factor(pooling_data$id))  ##number of subject
L <- 2  ##number of random effects
x <- model.matrix(~ time + time:group, data = pooling_data)  ##design matrix of fixed effects, pooled data
z <- model.matrix(~ time, data = pooling_data)   ##design matrix of random effects, pooled data
y <- pooling_data$response  ##response
sub_index <- sub_index_specification(data = pooling_data)  ##subject index for pooled data
var_beta <- N * solve(t(x)%*%x)