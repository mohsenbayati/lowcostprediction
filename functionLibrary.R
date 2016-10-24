
library(AUC)
library(Matrix)
library(pROC)  # Only used for delong CI calculations
library(caret) # Only used for delong to create cross-validation folds
library(glmnet)
library(mice)
#---------------------------------------------------------------------
myauc = function(outcome, predictions){
  
  if (sum(outcome)==0){
    warning("AUC is called with all zeros outcome vector")
    return(0.5)
  }else{
    return(AUC:::auc(AUC:::roc(predictions, factor(outcome))))
  }
}
#---------------------------------------------------------------------
# Function fills in missing values with column means.
MeanImpute = function(x) {
  
  return(apply(x, 2, function(a) {a[which(is.na(a))] = sum(a[!is.na(a)])/length(which(!is.na(a)))
  return(a)}))
}
#--------------------------------------------------------------------
# This function subsamples controls so the ratio of controls to cases would be "ratio"
BalanceData = function(XY, ratio = 3){
  d = ncol(XY)
  indCase = XY[,d]==1
  XYCase = XY[indCase,]
  XYControl = XY[!indCase,]
  
  nCase = nrow(XYCase)
  nControl = ratio*nCase
  XYControl = XYControl[sample(nrow(XYControl),nControl),]
  XY_train = rbind(XYCase,XYControl)  
}
#---------------------------------------------------------------------
# x is not really used, it is put in place so the function would be called via apply() function
HelperFunction1 = function(x,XY,LASSO=FALSE){
  
  n = nrow(XY)
  p = ncol(XY) - 1
  boot = sample(n,replace=TRUE)
  XY_b = XY[boot,]
  
  #----------- balance the training set to have 3 to 1 control to case ratio
  XY_train = BalanceData(XY_b)
  #--------------------------------------------------------------------------
  
  if (!LASSO){
    model_b = glm(Response~.,data=XY_train,family=binomial)
    
    # Predicting on the pre-balancing but the same data
    pred_b = predict(model_b, newdata=XY_b,type="response")
    
    # Predicting on the original (pre-bootstrap) data
    pred_b0 = predict(model_b, newdata =XY,type="response" )
  }else{
    
    model_b = cv.glmnet(x=as.matrix(XY_train[,1:p]),y=XY_train[,p+1],family="binomial",type="auc",nfolds = 5)
    
    # Predicting on the pre-balancing but the same data
    pred_b = predict(model_b, newx=as.matrix(XY_b[,1:p]),type="response")    
    
    # Predicting on the original (pre-bootstrap) data
    pred_b0 = predict(model_b, newx=as.matrix(XY[,1:p]),type="response")    
    # Predicting on the original data (but pre-balancing)    
  }
  
  AUC_B = myauc(XY_b[,p+1],pred_b)
  AUC_b0 = myauc(XY[,p+1],pred_b0)
  
  return(AUC_B-AUC_b0)
}
#---------------------------------------------------------------------

RunInnerBootstrap = function(XY,B=100,LASSO=FALSE){
  n = nrow(XY)
  p = ncol(XY)-1
  colnames(XY)[p+1]="Response"

  optimism_all = apply(matrix(1:B,nrow=1,ncol=B),2,function(x){HelperFunction1(x,XY,LASSO)})
  
  optimism = mean(optimism_all)
  #-----------------------
  # Now build train/test on the original data
  
  #----------- balance the training set to have 3 to 1 control to case ratio
  XY_train = BalanceData(XY)
  #--------------------------------------------------------------------------
  
  if (!LASSO){  
    model = glm(Response~.,data=XY_train,family=binomial)
    
    # Predicting on the original data (but pre-balancing)
    pred = predict(model, newdata=XY,type="response")
  }else{
    model = cv.glmnet(x=as.matrix(XY_train[,1:p]),y=XY_train[,p+1],family="binomial",type="auc",nfolds = 5)
    
    # Predicting on the original data (but pre-balancing)
    pred = predict(model, newx=as.matrix(XY[,1:p]),type="response")    
    
  }
  
  outputAUC = myauc(XY[,p+1],pred) - optimism
  return(outputAUC)
}
#---------------------------------------------------------------------
BCa = function(XY,nboot,fun,alpha=0.025,quiet=T){
  
  st=proc.time()
  theta0 = fun(XY)
  n = nrow(XY)
  values = rep(0,nboot)
  for (b in 1:nboot){
    values[b] = fun(XY[sample(n,replace=TRUE),])
    if ((!quiet)&(b%%10==0)){
      et=proc.time()	
      cat('BCa b=',b,'took',et[3]-st[3],'secs\n')
      st=proc.time()	
    }
  }
  
  # Using formulae 14.9-14.15 in Efron Tibshirani 1993
  
  z_alpha = qnorm(alpha)
  z1_alpha = qnorm(1-alpha)
  
  z0h = qnorm(mean(values< theta0))
  
  theta_hat = rep(0,n)
  for (i in 1:nrow(XY)){
    theta_hat[i] = fun(XY[-i,])
    
    if ((!quiet)&(i%%10==0)){
      et=proc.time()
      cat('BCa jk=',i,'took',et[3]-st[3],'secs\n')
      st=proc.time()    
    }
    
  }
  theta_star = mean(theta_hat)
  jk_numerator = sum((theta_star-theta_hat)^3)
  jk_denominator = 6 * (sqrt( sum( (theta_hat-theta_star)^2)  ))^3
  ah = jk_numerator/jk_denominator
  
  a1 = pnorm(z0h+ (z0h+z_alpha)/(1-ah*(z0h+z_alpha)))
  a2 = pnorm(z0h+ (z0h+z1_alpha)/(1-ah*(z0h+z1_alpha)))
  
  output = quantile(values,probs=c(a1,a2))
  
  return(output)
}
#---------------------------------------------------------------------
BootstrapCIpercentile = function(XY,nboot,fun, alpha=0.025){
  n = nrow(XY)
  values = rep(0,nboot)
  for (b in 1:nboot){
    values[b] = fun(XY[sample(n,replace=TRUE),])
  }
  
  output = quantile(values,probs=c(alpha,1-alpha))
  
  return(output)
}
#---------------------------------------------------------------------
delongCI = function(XY, alpha=0.025, train_ratio = 0.5, LASSO = FALSE){
  n = nrow(XY)
  p = ncol(XY)-1
  colnames(XY)[p+1]="Response"
  
  train_indices = sample(n,round(n*train_ratio))
  XY_train = XY[train_indices,]
  XY_test = XY[-train_indices,]
  
  XY_train_bal = BalanceData(XY_train)
  #--------------------------------------------------------------------------
  if (!LASSO){
    model = glm(Response~.,data=XY_train_bal,family=binomial)
    pred = predict(model, newdata=XY_test,type="response")
  }else{
    model = cv.glmnet(x=as.matrix(XY_train[,1:p]),y=XY_train[,p+1],family="binomial",type="auc",nfolds = 5)
    pred = predict(model, newx=as.matrix(XY_test[,1:p]),type="response")    
  }
  
  rocObject = pROC:::roc(XY_test[,p+1],pred)
  
  output = pROC:::ci.auc(rocObject,method="delong",conf.level=1-2*alpha)
  
  output = as.numeric(output)[c(1,3)]
  return(output)
}
#---------------------------------------------------------------------
MultipleIMputationCI = function(XY, nMI=5, impMethod = "pmm", alpha=0.025, train_ratio = 0.5, LASSO = FALSE, miceMaxIt = 50){
# Methods: pmm or sample. Write methods(mice) to see a list of all options
    
  AUCs = rep(0,nMI)
  
  n = nrow(XY)
  p = ncol(XY)-1

  mice_output = mice(XY[,1:p],m=nMI,maxit=miceMaxIt,meth=impMethod,print=FALSE)
    
  for (i in 1:nMI){
    
    Xsc = as.data.frame(scale(as.matrix(complete(mice_output,i)),center = T, scale = T)) # Standardize
    XY_i = cbind(Xsc,XY[,p+1])
    colnames(XY_i)[p+1]="Response"
    
    train_indices = sample(n,round(n*train_ratio))
    XY_train = XY_i[train_indices,]
    XY_test = XY_i[-train_indices,]
    
    XY_train_bal = BalanceData(XY_train)
    #--------------------------------------------------------------------------
    if (!LASSO){
      model = glm(Response~.,data=XY_train_bal,family=binomial)
      pred = predict(model, newdata=XY_test,type="response")
    }else{
      model = cv.glmnet(x=as.matrix(XY_train[,1:p]),y=XY_train[,p+1],family="binomial",type="auc",nfolds = 5)
      pred = predict(model, newx=as.matrix(XY_test[,1:p]),type="response")    
    }
    
    AUCs[i] = pROC:::auc(XY_test[,p+1],pred)
  }
  
  output = quantile(AUCs,probs=c(alpha,1-alpha))
  
  return(output)
}
