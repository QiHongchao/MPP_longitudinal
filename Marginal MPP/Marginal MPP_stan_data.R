N0 <- nrow(historical)  ##number of observations, historical
N <- nrow(current)  ##number of observations, current
J0 <- nlevels(as.factor(historical$id))  ##id is put in the last column, number of subject, historical
J <- nlevels(as.factor(current$id))  ##id is put in the last column, number of subject, current
x0 <- as.matrix(historical[, c("int", "time")])  ##design matrix of fixed effects, historical
z0 <- as.matrix(historical[, c("int", "time")])  ##design matrix of random effects, historical
y0 <- historical$response  ##response was put in the second last column, historical
x <- as.matrix(current[, c("int", "time", "interaction")])  ##design matrix of fixed effects, current
z <- as.matrix(current[, c("int", "time")])  ##design matrix of random effects, current
y <- current$response  ##response was put in the second last column, current
K0 <- ncol(x0) ##number of fixed effects in historical study, there is no interaction
K <- ncol(x)  ##number of fixed effects in current study
L <- ncol(z0)  ##number of random effects, we assume same random effects in historical and current data
num_power <- num_power  ##number of power values
logsc <- logsc  ##Log scaling constant
sub_index0 <- sub_index_specification(data = historical)  ##subject index for historical data
sub_index <- sub_index_specification(data = current)  ##subject index for current data