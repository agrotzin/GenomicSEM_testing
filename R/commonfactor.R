commonfactor <-function(covstruc,estimation="DWLS"){
  time<-proc.time()
  
  #function to create lavaan syntax for a 1 factor model given k phenotypes
  write.Model1 <- function(k) {  
    Model1 <- ""
    for (i in 1) {
      linestart <- paste("F1"," =~ NA*",colnames(S_LD)[i], sep = "")  
      if (k-i > 0) {
        linemid <- ""
        for (j in (i+1):k) {
          linemid <- paste(linemid, " + ", colnames(S_LD)[j], sep = "")
        }
      } else {linemid <- ""}
      
    }
    Model1 <- paste(Model1, linestart, linemid, " \n ", "F1 ~~ 1*F1", " \n ", sep = "")
    
    return(Model1)
  } 
  

  ##read in the LD portion of the V (sampling covariance) matrix
  V_LD<-as.matrix(covstruc[[1]])
  
  ##read in the LD portion of the S (covariance) matrix
  S_LD<-as.matrix(covstruc[[2]])
  
  ##k = number of phenotypes in dataset (i.e., number of columns in LD portion of S matrix)
  k<-ncol(S_LD)
  
  ##return error if model is underidentified 
  if(k == 2){
    stop("Their are only 2 variables in the genetic covariance matrix so the common factor model will be under identified (df = -1). 
         You can either specify a common factor model with constrained factor loadings with the user model function or rerun ldsc with at least one additional variable.")
  }       
  
  ##create the 1 factor model with k # of indicators
  Model1 <- write.Model1(k)
  
  ##pull the column names specified in the munge function
  rownames(S_LD)<-colnames(S_LD)
  
  ##smooth to near positive definite if either V or S are non-positive definite
  ##smooth to near positive definite if either V or S are non-positive definite
  ks<-nrow(S_LD)
  S_LDb<-S_LD
  smooth1<-ifelse(eigen(S_LD)$values[ks] <= 0, S_LD<-as.matrix((nearPD(S_LD, corr = FALSE))$mat), S_LD<-S_LD)
  LD_sdiff<-max(abs(S_LD-S_LDb))
  
  kv<-nrow(V_LD)
  V_LDb<-V_LD
  smooth2<-ifelse(eigen(V_LD)$values[kv] <= 0, V_LD<-as.matrix((nearPD(V_LD, corr = FALSE))$mat), V_LD<-V_LD)
  LD_sdiff2<-max(abs(V_LD-V_LDb))
  
  SE_pre<-matrix(0, k, k)
  SE_pre[lower.tri(SE_pre,diag=TRUE)] <-sqrt(diag(V_LDb))
  
  SE_post<-matrix(0, k, k)
  SE_post[lower.tri(SE_post,diag=TRUE)] <-sqrt(diag(V_LD))
  
  Z_pre<-S_LDb/SE_pre
  Z_post<-S_LD/SE_post
  Z_diff<-(Z_pre-Z_post)
  Z_diff[which(!is.finite(Z_diff))]<-0
  Z_diff<-max(Z_diff)
  rm(V_LDb,S_LDb,Z_pre,Z_post)
  
  ##run model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
  #transform V_LD matrix into a diagonal weight matrix: 
  z<-(k*(k+1))/2
  
  ##save the ordering
  order <-(1:nrow(V_LD))
  
  ##reorder the weight (inverted V_LD) matrix
  V_Reorder<-V_LD[order,order]
  W_Reorder<-diag(z)
  diag(W_Reorder)<-diag(V_Reorder)
  W_Reorder<-solve(W_Reorder)

  print("Running Model")    
  if(estimation == "DWLS"){
    ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
    empty<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2, optim.dx.tol = .01))
  }
  
  if(estimation == "ML"){
    empty<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "ML", sample.nobs = 200, optim.dx.tol = .01,sample.cov.rescale=FALSE))
  }
  
  empty$warning$message[1]<-ifelse(is.null(empty$warning$message), empty$warning$message[1]<-0, empty$warning$message[1])
  
  if(class(empty$value)[1] == "simpleError" | grepl("solution has NOT",  as.character(empty$warning)) == TRUE){
    print("The common factor initially failed to converge. A lower bound of 0 on residual variances has been added to try and troubleshoot this.")
    
    #create unique combination of letters for residual variance parameter labels
    n<-combn(letters,4)[,sample(1:14000, k, replace=FALSE)]
    
    Model3<-""
    for (p in 1:k) {
      linestart3a <- paste(colnames(S_LD)[p], " ~~ ",  paste(n[,p],collapse=""), "*", colnames(S_LD)[p], sep = "")
      linestart3b <- paste(paste(n[,p],collapse=""), " > .0001", sep = "")
      Model3<-paste(Model3, linestart3a, " \n ", linestart3b, " \n ", sep = "")}
    
    Model1<-paste(Model1,Model3)
    
    if(estimation == "DWLS"){
      empty<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2, optim.dx.tol = .01))
    }
    
    if(estimation == "ML"){
      empty<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "ML", sample.nobs = 200, optim.dx.tol = .01,sample.cov.rescale=FALSE))
    }
    
  }else{}
  
  if(class(empty$value)[1] == "simpleError"){
    print("The common factor model failed to converge on a solution. Please try specifying an alternative model using the usermodel function.")
  }
  
  if(!(is.null(empty$warning))){
    if(grepl("solution has NOT",  as.character(empty$warning)) == TRUE){
      print("The common factor model failed to converge on a solution. Please try specifying an alternative model using the usermodel function.")
    }}
  
  #pull the delta matrix (this doesn't depend on N)
  S2.delt <- lavInspect(Model1_Results, "delta")
  
  ##weight matrix from stage 2
  S2.W <- lavInspect(Model1_Results, "WLS.V") 
  
  #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
  bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt) 
  
  #create the "lettuce" part of the sandwich
  lettuce <- S2.W%*%S2.delt
  
  #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
  Ohtt <- bread %*% t(lettuce)%*%V_Reorder%*%lettuce%*%bread  
  
  #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
  SE <- as.matrix(sqrt(diag(Ohtt)))
  
  #calculate model chi-square
  Eig<-as.matrix(eigen(V_LD)$values)
  Eig2<-diag(z)
  diag(Eig2)<-Eig
  
  #Pull P1 (the eigen vectors of V_eta)
  P1<-eigen(V_LD)$vectors
  
  implied<-as.matrix(fitted(Model1_Results))[1]
  implied_order<-colnames(S_LD)
  implied[[1]]<-implied[[1]][implied_order,implied_order]
  implied2<-S_LD-implied[[1]]
  eta<-as.vector(lowerTriangle(implied2,diag=TRUE))
  Q<-t(eta)%*%P1%*%solve(Eig2)%*%t(P1)%*%eta
  chisq<-Q
  
  ##df of user model
  df <- lavInspect(Model1_Results, "fit")["df"]
  
  #caculate CFI
  resid_CFI<-S_LD
  diag(resid_CFI)<-0
  
  #eta: the vector of unique elements of the residual matrix
  eta_CFI<-as.vector(lowerTriangle(resid_CFI,diag=TRUE))
  
  #matrix algebra weighting the vector of residuals by the precision of those residuals (i.e., P1 and Eig)
  CFI_chi<-t(eta_CFI)%*%P1%*%solve(Eig2)%*%t(P1)%*%eta_CFI
  
  ##df of independence Model
  dfCFI <- (((k * (k + 1))/2) - k)
  
  #calculate CFI
  CFI<-as.numeric(((CFI_chi-dfCFI)-(chisq-df))/(CFI_chi-dfCFI))
  CFI<-ifelse(CFI > 1, 1, CFI)
  
  ##transform the S covariance matrix to S correlation matrix
  D=sqrt(diag(diag(S_LD)))
  S_Stand=solve(D)%*%S_LD%*%solve(D)
  rownames(S_Stand)<-rownames(S_LD)
  colnames(S_Stand)<-colnames(S_Stand)
  
  #obtain diagonals of the original V matrix and take their sqrt to get SE's
  Dvcov<-sqrt(diag(V_LD))
  
  #calculate the ratio of the rescaled and original S matrices
  scaleO=as.vector(lowerTriangle((S_Stand/S_LD),diag=T))
  
  ## MAke sure that if ratio in NaN (devision by zero) we put the zero back in: ### TEMP STUPID MICHEL FIX!
  scaleO[is.nan(scaleO)] <- 0
  
  #rescale the SEs by the same multiples that the S matrix was rescaled by
  Dvcovl<-as.vector(Dvcov*t(scaleO))
  
  #obtain the sampling correlation matrix by standardizing the original V matrix
  Vcor<-cov2cor(V_LD)
  
  #rescale the sampling correlation matrix by the appropriate diagonals
  V_stand<-diag(Dvcovl)%*%Vcor%*%diag(Dvcovl)
  V_stand2<-diag(z)
  diag(V_stand2)<-diag(V_stand)
  
  ### make sure no value on the diagonal of V is 0 
  diag(V_stand2)[diag(V_stand2) == 0] <- 2e-9
  
  W_stand<-solve(V_stand2[order,order])
  
  if(estimation == "DWLS"){
    emptystand<-.tryCatch.W.E(Fit_stand <- sem(Model1, sample.cov = S_Stand, estimator = "DWLS", WLS.V = W_stand, sample.nobs = 2, optim.dx.tol = .01)) 
  }
  
  if(estimation == "ML"){
    emptystand<-.tryCatch.W.E(Fit_stand <- sem(Model1, sample.cov = S_Stand, estimator = "ML",  sample.nobs = 200, optim.dx.tol = .01,sample.cov.rescale=FALSE)) 
  }
  
  ##perform same procedures for sandwich correction as in the unstandardized case
  delt_stand <- lavInspect(Fit_stand, "delta") 
  W_stand <- lavInspect(Fit_stand, "WLS.V") 
  bread_stand <- solve(t(delt_stand)%*%W_stand %*%delt_stand)
  lettuce_stand <- W_stand%*%delt_stand
  Vcov_stand<-as.matrix(V_stand[order,order])
  Ohtt_stand <- bread_stand %*% t(lettuce_stand)%*%Vcov_stand%*%lettuce_stand%*%bread_stand
  SE_stand <- as.matrix(sqrt(diag(Ohtt_stand)))
  
  unstand<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
  unstand<-subset(unstand, unstand$free != 0)                    
  unstand$free<-NULL
  
  stand<-data.frame(inspect(Fit_stand,"list")[,c(8,14)])
  stand<-subset(stand, stand$free != 0)
  stand$free<-NULL

  AIC<-(Q + 2*lavInspect(Model1_Results, "fit")["npar"])
  SRMR<-lavInspect(Model1_Results, "fit")["srmr"]
  
  modelfit<-cbind(chisq,df,AIC,CFI,SRMR)
  results<-cbind(unstand,SE,stand,SE_stand)
  
  ##name the columns of the results file
  colnames(results)=c("lhs","op","rhs","Unstandardized_Estimate","Unstandardized_SE","Standardized_Est","Standardized_SE")
  
  ##name model fit columns
  colnames(modelfit)=c("chisq","df","AIC","CFI","SRMR")
  modelfit<-data.frame(modelfit)
  modelfit$p_chisq<-ifelse(modelfit$chisq != 'NA', modelfit$p_chisq<-pchisq(modelfit$chisq, modelfit$df,lower.tail=FALSE), modelfit$p_chisq<-NA)
  modelfit$chisq<-ifelse(modelfit$df == 0, modelfit$chisq == NA, modelfit$chisq)  
  modelfit$AIC<-ifelse(modelfit$df == 0, modelfit$AIC == NA, modelfit$AIC)  
  modelfit$p_chisq<-ifelse(modelfit$df == 0, modelfit$p_chisq == NA, modelfit$p_chisq)  
  
  order<-c(1,2,6,3,4,5)
  modelfit<-modelfit[,order]
  
  time_all<-proc.time()-time
  print(time_all[3])
  
  
  if(LD_sdiff > 0){
    print(paste("The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was ", LD_sdiff, "As a result of the smoothing, the largest Z-statistic change for the genetic covariances was ", Z_diff, ". We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS.", sep = " "))
  }
  
  if(LD_sdiff > .025){
    warning("A difference greater than .025 was observed pre- and post-smoothing in the genetic covariance matrix. This reflects a large difference and results should be interpreted with caution!! This can often result from including low powered traits, and you might consider removing those traits from the model. If you are going to run a multivariate GWAS we strongly recommend setting the smooth_check argument to true to check smoothing for each SNP.")
  }
  
  if(Z_diff > .025){
    warning("A difference greater than .025 was observed pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a large difference and results should be interpreted with caution!! This can often result from including low powered traits, and you might consider removing those traits from the model. If you are going to run a multivariate GWAS we strongly recommend setting the smooth_check argument to true to check smoothing for each SNP.")
  }
  
  if(LD_sdiff2 > 0){
    print(paste("The V matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was ", LD_sdiff2, "As a result of the smoothing, the largest Z-statistic change for the genetic covariances was ", Z_diff,  ". We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS.", sep = " "))
  }
  
  if(modelfit$CFI < 0){
    warning(paste("CFI estimates below 0 should not be trusted, and indicate that the other model fit estimates should be interpreted with caution. A negative CFI estimates typically appears due to negative residual variances."))
  }
  
  results$p_value<-2*pnorm(abs(as.numeric(results$Unstandardized_Estimate)/as.numeric(results$Unstandardized_SE)),lower.tail=FALSE)
  results$p_value<-ifelse(results$p_value == 0, "< 5e-300", results$p_value)
  
  return(list(modelfit=modelfit,results=results))
  
}
