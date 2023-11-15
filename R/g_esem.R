g_esem <-function(covstruc, model = "", rotation = "geomin", rotation.args=list(NULL),std.lv=TRUE, toler=NULL){ 
  time<-proc.time()
  ##determine if the model is likely being listed in quotes and print warning if so
  test<-c(str_detect(model, "~"),str_detect(model, "="),str_detect(model, "\\+"))
  if(any(test) != TRUE){
    warning("Your model name may be listed in quotes; please remove the quotes and try re-running if the function has returned an error about not locating the ReorderModel.")
  }
  
  ##read in the LD portion of the V (sampling covariance) matrix
  V_LD<-as.matrix(covstruc[[1]])
  
  ##read in the LD portion of the S (covariance) matrix
  S_LD<-as.matrix(covstruc[[2]])
  
  ##k = number of phenotypes in dataset (i.e., number of columns in LD portion of S matrix)
  k<-ncol(S_LD)
  
  ##size of V matrix used later in code to create diagonal V matrix
  z<-(k*(k+1))/2
  
  Model1<-model
  
  ##pull the column names specified in the munge function
  S_names<-colnames(S_LD)
  rownames(S_LD)<-colnames(S_LD)
  
  ##name columns of V to remove any variables not used in the current analysis
  y<-expand.grid(S_names,S_names)
  y<-y[!duplicated(apply(y,1,function(x) paste(sort(x),collapse=''))),]
  V_Names<-paste(y$Var1,y$Var2,sep=" ")
  colnames(V_LD)<-V_Names
  rownames(V_LD)<-V_Names
  
  ##determine whether all variables in S are in the model
  ##if not, remove them from S_LD and V_LD for this particular run
  remove2<-c()
  w<-1
  
  ##also for exact cases
  for(i in 1:length(S_names)){
    S_names[[i]]<-paste0("\\b", S_names[[i]],"\\b",sep="")
  }
  
  for(i in 1:length(S_names)){
    b<-grepl(S_names[i], model)
    if(b == FALSE){
      remove<-paste0("\\b", colnames(S_LD)[i],"\\b",sep="")
      remove2[w]<-i
      V_LD <- V_LD[-grep(pattern=remove[1],row.names(V_LD)),-grep(pattern=remove[1],colnames(V_LD))]
      w<-w+1
      if (!(is.matrix(V_LD))) {
        stop("None of the trait names in the LDSC output match names in the model")
      }
    }else{}
  }
  
  if(is.null(remove2) == FALSE){
    S_LD<-S_LD[-remove2,-remove2]
  }
  
  ##redefine k and z and model names after removing non-used variables
  k<-ncol(S_LD)
  z<-(k*(k+1))/2
  
  ##smooth to near positive definite if either V or S are non-positive definite
  S_LDb<-S_LD
  smooth1<-ifelse(eigen(S_LD)$values[nrow(S_LD)] <= 0, S_LD<-as.matrix((nearPD(S_LD, corr = FALSE))$mat), S_LD<-S_LD)
  LD_sdiff<-max(abs(S_LD-S_LDb))
  
  V_LDb<-V_LD
  smooth2<-ifelse(eigen(V_LD)$values[nrow(V_LD)] <= 0, V_LD<-as.matrix((nearPD(V_LD, corr = FALSE))$mat), V_LD<-V_LD)
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
  #transform V_LD matrix into a weight matrix: 
  W <- solve(V_LD)

  empty3<-.tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_LD, estimator = "WLS", WLS.V = W, sample.nobs = 2,warn=FALSE,std.lv=std.lv, optim.dx.tol = +Inf,optim.force.converged=TRUE,control=list(iter.max=1),rotation=rotation,rotation.args=rotation.args))
  
  if(class(empty3$value) != "lavaan"){
    warning(paste("The function has stopped due to convergence issues for your primary model. Please contact us with your specific model and variables used or try specifying an alternative model"))
  }

#save internal ordering of S and V
order <- .rearrange_esem(k = k, fit = ReorderModel, names = rownames(S_LD))

#reorder S                         
S_LD<-S_LD[order[[1]],order[[1]]]

#create weight matrix using inverted, reordered sampling covariance matrix
W_Reorder<-solve(V_LD[order[[2]],order[[2]]])

  print("Running model")
  
  ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
  empty4<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "WLS", std.lv=std.lv,WLS.V = W_Reorder, sample.nobs = 2,optim.dx.tol = +Inf,rotation=rotation,rotation.args=rotation.args))
  
  empty4$warning$message[1]<-ifelse(is.null(empty4$warning$message), empty4$warning$message[1]<-0, empty4$warning$message[1])

  if(class(empty4$value)[1] == "simpleError"){
    warning("The model failed to converge on a solution. Please try specifying an alternative model")}
  
      #calculate model chi-square: still need this becasue that part is contingent on number of observations
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

      ##df of user model
      df<-lavInspect(Model1_Results, "fit")["df"]
        
        if(!(is.character(Q))){
          chisq<-Q
          AIC<-(Q + 2*lavInspect(Model1_Results, "fit")["npar"])}else{chisq<-Q
          AIC<-NA}
        
        SRMR<-lavInspect(Model1_Results, "fit")["srmr"]
        
        modelfit<-cbind(chisq,df,AIC,SRMR)
        
    ##pull results and name the columns of the results file
    results<-data.frame(inspect(Model1_Results, "list")[,c("lhs", "op", "rhs", "est", "se", "est.std", "se.std","efa")])
 
    colnames(results)=c("lhs","op","rhs","Unstand_Est","Unstand_SE", "Stand_Est", "Stand_SE","efa")
 
    ##name model fit columns
    colnames(modelfit)=c("chisq","df","AIC","SRMR")
    
    modelfit<-data.frame(modelfit)
    
    if(!(is.character(modelfit$chisq)) & !(is.factor(modelfit$chisq))){
      modelfit$chisq<-as.numeric(as.character(modelfit$chisq))
      modelfit$df<-as.numeric(as.character(modelfit$df))
      modelfit$p_chisq<-ifelse(!(is.character(modelfit$chisq)), modelfit$p_chisq<-pchisq(modelfit$chisq, modelfit$df,lower.tail=FALSE), modelfit$p_chisq<-NA)
      modelfit$chisq<-ifelse(modelfit$df == 0, modelfit$chisq == NA, modelfit$chisq)  
      modelfit$AIC<-ifelse(modelfit$df == 0, modelfit$AIC == NA, modelfit$AIC)  
      modelfit$p_chisq<-ifelse(modelfit$df == 0, modelfit$p_chisq == NA, modelfit$p_chisq)
      modelfit$SRMR<-ifelse(modelfit$df == 0, modelfit$SRMR == NA, modelfit$SRMR)
      order<-c(1,2,5,3,4)
      modelfit<-modelfit[,order]
      }
    
    if(modelfit$df == 0){
      print("Model fit statistics are all printed as NA as you have specified a fully saturated model (i.e., df = 0)")
    }
    
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

    results$p_value<-2*pnorm(abs(results$Unstand_Est/results$Unstand_SE),lower.tail=FALSE)
    results$p_value<-ifelse(results$p_value == 0, "< 5e-300", results$p_value)
    results$Unstand_SE<-ifelse(results$Unstand_SE == 0, "",results$Unstand_SE)
    results$Stand_SE<-ifelse(results$Stand_SE == 0, "",results$Stand_SE)
    results$Unstand_SE<-as.numeric(results$Unstand_SE)
    results$Stand_SE<-as.numeric(results$Stand_SE)
    results$p_value<-ifelse(results$Stand_SE == "", "",results$p_value)
    
    #remove parameter constraints to avoid confusion
    results<-subset(results, results$op != ">" & results$op != "<")

    if(empty4$warning$message[1] != 0){
      warning(paste0("The model produced the following warning: ", empty4$warning$message[1],sep=""))
    }  
    
    time_all<-proc.time()-time
    print(time_all[3])
    
   return(list(modelfit=modelfit,results=results))
}
