##Generate current data
current_data_gen<-function(heterogeneity,fe,trt,ntp,ngroup,varb,varw,num_control,num_treat){
  ##Between-study heterogeneity
  if (heterogeneity == 0) {cov_study <- diag(0, 2)}
  if (heterogeneity == 1) {cov_study <- diag(c(0.01, 0))}
  if (heterogeneity == 2) {cov_study <- diag(c(0.09, 0))}
  if (heterogeneity == 3) {cov_study <- diag(c(0.16, 0))}
  if (heterogeneity == 4) {cov_study <- diag(0.01, 2)}
  if (heterogeneity == 5) {cov_study <- diag(0.09, 2)}
  if (heterogeneity == 6) {cov_study <- diag(0.16, 2)}
  fe_current <- rmvnorm(1, fe, cov_study)
  ##Design matrix for fixed effect
  time <- rep(seq(0, 1, 1/(ntp-1)),ngroup)
  group <- rep(c(0, 1), each = ntp)
  dm_fe <- model.matrix( ~ time + time:group)
  ##Design matrix for random effects
  dm_re <- model.matrix( ~ time[1:ntp])
  ##Control arm
  data_control <- data.frame(time = rep(seq(0, 1, 1/(ntp-1)),num_control), group = 0)
  for (i in 1:num_control){
    data_control$response[(ntp*(i-1)+1):(ntp*i)]<- rmvnorm(1,dm_fe[1:ntp,]%*%c(fe_current, trt), 
                                                   dm_re%*%diag(varb)%*%t(dm_re)+varw*diag(1,ntp))
  }
  ##Treatment arm
  data_treat <- data.frame(time = rep(seq(0, 1, 1/(ntp-1)),num_treat), group = 1)
  for (i in 1:num_treat){
    data_treat$response[(ntp*(i-1)+1):(ntp*i)]<- rmvnorm(1,dm_fe[(ntp+1):(ntp*2),]%*%c(fe_current, trt), 
                                                         dm_re%*%diag(varb)%*%t(dm_re)+varw*diag(1,ntp))
  }
  ##Data
  data<-rbind(data_control,data_treat)
  data$id<-sort(rep(1:(num_control+num_treat),ntp))
  ##Return the dataset
  return(data)
}

##Genreate historical data
historical_data_gen<-function(heterogeneity,fe,ntp,varb,varw,num_control){
  ##Between-study heterogeneity
  if (heterogeneity == 0) {cov_study <- diag(0, 2)}
  if (heterogeneity == 1) {cov_study <- diag(c(0.01, 0))}
  if (heterogeneity == 2) {cov_study <- diag(c(0.09, 0))}
  if (heterogeneity == 3) {cov_study <- diag(c(0.16, 0))}
  if (heterogeneity == 4) {cov_study <- diag(0.01, 2)}
  if (heterogeneity == 5) {cov_study <- diag(0.09, 2)}
  if (heterogeneity == 6) {cov_study <- diag(0.16, 2)}
  
  ##Design matrix for fixed effect
  time <- seq(0, 1, 1/(ntp-1)) ##time point start from 0, equally spaced, time interval is 1/5
  dm_fe <- model.matrix( ~ time)
  ##Design matrix for random effects
  dm_re <- dm_fe
  ##Historical control arm
  data_control <- data.frame(time = rep(seq(0, 1, 1/(ntp-1)),num_control), group = 0)
  
  fe_hist <- fe + t(rmvnorm(1, rep(0,ncol(dm_fe)), cov_study))
  
  for (i in 1:num_control) {
    data_control$response[(ntp*(i-1)+1):(ntp*i)] <- rmvnorm(1,dm_fe%*%fe_hist, 
                                                          dm_re%*%diag(varb)%*%t(dm_re)+varw*diag(1,ntp))
  }
  data_control$id <- sort(rep(1:num_control,ntp))
  ##Return the dataset
  return(data_control)
}