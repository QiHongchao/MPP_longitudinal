N <- nrow(current)  ##number of observations of current data
K <- 3  ##number of fixed effects
J <- nlevels(as.factor(current$id))  ##number of subjects
L <- 2  ##number of random effects
x <- model.matrix(~ time + time:group, data = current)  ##design matrix of fixed effects, current data
z <- model.matrix(~ time, data = current)  ##design matrix of random effects, current data
y <- current$response  ##response
sub_index <- sub_index_specification(data = current)  ##subject index for current data
var_beta <- N * solve(t(x)%*%x)
