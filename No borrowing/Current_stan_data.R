N <- nrow(current)  ##number of observations of current data
J <- nlevels(as.factor(current$id))  ##id is in the last column, number of subject
x <- as.matrix(current[, c("int", "time", "interaction")])  ##design matrix of fixed effects, current
z <- as.matrix(current[, c("int", "time")])  ##design matrix of random effects, current
y <- current$response 
K <- ncol(x)  ##number of fixed effects in current study
L <- ncol(z)  ##number of random effects, we assume same random effects in historical and current data
sub_index <- sub_index_specification(data = current)  ##subject index for current data