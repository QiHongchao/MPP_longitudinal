##In oder to simplify the case, I simulate an RCT with equal number of patients in treatment and control arm, respectively.
##Design matrix: we assume current and historical study have the same design matrix.
##Design matrix for fixed effects and random effects
design_matrix_gen_current<-function(ntp,group=2,num_fe=3,num_re=2,fe=c(2,1),varb1,varb2,varw,num_control,num_treat){
  ##with treatment effect or not
  if (get("heterogeneity_trt", envir = .GlobalEnv)[2] == "No") {fe[3] <- 0}
  if (get("heterogeneity_trt", envir = .GlobalEnv)[2] == "Yes") {fe[3] <- interaction_time_trt}
  
  ##Design matrix for fixed effect
  dm_fe<-matrix(NA,ntp*group,num_fe)
  dm_fe[,1]<-1
  dm_fe[,2]<-rep(seq(0, 1, 1/(ntp-1)),group) ##time point start from 0, equally spaced, time interval is 1
  dm_fe[,3]<-dm_fe[,2]*c(rep(0,ntp),rep(1,ntp))
  ##Design matrix for random effects
  dm_re<-matrix(NA,ntp,num_re)
  dm_re[,1]<-1
  dm_re[,2]<-seq(0, 1, 1/(ntp-1))
  ##Control arm
  control<-rmvnorm(num_control,dm_fe[1:ntp,]%*%fe, dm_re%*%(c(varb1,varb2)*diag(1,num_re))%*%t(dm_re)+varw*diag(1,ntp))
  dm_fe_control<-matrix(NA,ntp*num_control,num_fe+2)
  dm_fe_control[,1]<-1
  dm_fe_control[,2]<-rep(seq(0, 1, 1/(ntp-1)),num_control)
  dm_fe_control[,3]<-0
  for (i in 1:num_control){
    dm_fe_control[(ntp*(i-1)+1):(ntp*i),num_fe+1]<-control[i,]
  }
  ##Treatment arm
  treat<-rmvnorm(num_treat,dm_fe[(ntp+1):(ntp*2),]%*%fe, dm_re%*%(c(varb1,varb2)*diag(1,num_re))%*%t(dm_re)+varw*diag(1,ntp))
  dm_fe_treat<-matrix(NA,ntp*num_control,num_fe+2)
  dm_fe_treat[,1]<-1
  dm_fe_treat[,2]<-rep(seq(0, 1, 1/(ntp-1)),num_treat)
  dm_fe_treat[,3]<-dm_fe_treat[,2]
  for (i in 1:num_treat){
    dm_fe_treat[(ntp*(i-1)+1):(ntp*i),num_fe+1]<-treat[i,]
  }
  ##Response
  response<-rbind(dm_fe_control,dm_fe_treat)
  response[,num_fe+2]<-sort(rep(1:(num_control+num_treat),ntp))
  colnames(response)<-c("int","time","interaction","response","id")
  ##Return the dataset
  return(as.data.frame(response))
}

design_matrix_gen_hc<-function(ntp,num_fe=2,num_re=2,fe=c(2,1),varb1,varb2,varw,num_control){
  ##standard deviation of time trend
  if (get("heterogeneity_trt", envir = .GlobalEnv)[1] == "No") {cov_study <- matrix(rep(0, 4), 2)}
  if (get("heterogeneity_trt", envir = .GlobalEnv)[1] == "Low+RI") {cov_study <- matrix(c(0.01, 0, 0, 0), 2)}
  if (get("heterogeneity_trt", envir = .GlobalEnv)[1] == "Moderate+RI") {cov_study <- matrix(c(0.09, 0, 0, 0), 2)}
  if (get("heterogeneity_trt", envir = .GlobalEnv)[1] == "High+RI") {cov_study <- matrix(c(0.16, 0, 0, 0), 2)}
  if (get("heterogeneity_trt", envir = .GlobalEnv)[1] == "Low+RIS") {cov_study <- diag(rep(0.01, 2))}
  if (get("heterogeneity_trt", envir = .GlobalEnv)[1] == "Moderate+RIS") {cov_study <- diag(rep(0.09, 2))}
  if (get("heterogeneity_trt", envir = .GlobalEnv)[1] == "High+RIS") {cov_study <- diag(rep(0.16, 2))}
  
  ##Design matrix for fixed effect
  dm_fe<-matrix(NA,ntp,num_fe)
  dm_fe[,1]<-1
  dm_fe[,2]<-seq(0, 1, 1/(ntp-1)) ##time point start from 0, equally spaced, time interval is 1/5
  ##Design matrix for random effects
  dm_re<-matrix(NA,ntp,num_re)
  dm_re[,1]<-1
  dm_re[,2]<-seq(0, 1, 1/(ntp-1))
  ##Control arm
  control<-rmvnorm(num_control,dm_fe[1:ntp,]%*%fe + dm_fe[1:ntp,]%*%as.numeric(rmvnorm(1, rep(0,ncol(dm_fe)), cov_study)), 
                   dm_re%*%(c(varb1,varb2)*diag(1,num_re))%*%t(dm_re)+varw*diag(1,ntp))
  dm_fe_control<-matrix(NA,ntp*num_control,num_fe+2)
  dm_fe_control[,1]<-1
  dm_fe_control[,2]<-rep(seq(0, 1, 1/(ntp-1)),num_control)
  for (i in 1:num_control){
    dm_fe_control[(ntp*(i-1)+1):(ntp*i),num_fe+1]<-control[i,]
  }
  ##Response
  response<-dm_fe_control
  response[,num_fe+2]<-sort(rep(1:num_control,ntp))
  colnames(response)<-c("int","time","response","id")
  ##Return the dataset
  return(as.data.frame(response))
}