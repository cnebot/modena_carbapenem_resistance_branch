
# This file contains Functions to estimate VAR order n=2,3 or 4
# arguments:
#           X is a matrix of n variables in columns  
#       ir_CUM= FALSE (by default) noncumulative Impulse response, if =TRUE then cumulative Impulse response
#       force.degree imposes the degree of the variables to run VAR. if NULL, it makes the variables to be stationary previously testing Unitroot  
#       force.p imposes order p in VAR estimation, if NULL it estimates the best lag being 6 the lag.max
#       if Amat or Bmat are not NULL then SVAR is also computed. 
#        Amat is the matrix of contemporaneous relationships between variables
#        Bmat is the matrix of dependencies between errors in equations



require(vars)
eval_var2<-function(X,ir_CUM=F, force.degree=NULL, force.p=NULL, Amat=NULL, Bmat=NULL){
  
  x<-X[,1]
  y<-X[,2]
  names<-colnames(X)
  
  SVAR.switch<-1 
  if(is.null(Amat) && is.null(Bmat)) SVAR.switch<-0  # If no A and B matrix Switch OFF, deactivating SVAR estimation
  

  # Transfom order of series force.degree or just making them stationary previously testing unitroot 
  
  if(is.null(force.degree)){
    
    dx<-0;dy<-0
     
  
  if(unitroot_ndiffs(ts(x), alpha=0.05)!=0){dx<-unitroot_ndiffs(ts(x), alpha=0.05);x<-diff(ts(x), differences=dx)}
  if(unitroot_ndiffs(ts(y), alpha=0.05)!=0){dy<-unitroot_ndiffs(ts(y), alpha=0.05);y<-diff(ts(y), differences=dy)}

    degree<-c(dx,dy)
    
    }else{ 
  
    degree<-force.degree
    dx=degree[1]
    dy=degree[2]
    }
  
  maxd<-max(dx,dy)
  
  if(dx<maxd){x<-head(x,-I(maxd-dx))}
  if(dy<maxd){y<-head(y,-I(maxd-dy))}    
  
  
  assign(names[1],x)
  assign(names[2],y)
  
  names(degree)<-names
  # Select order and estimate VAR
  
  matrix<-cbind(ts(x),ts(y))
  colnames(matrix)<-names
  
  
  var.fit.1 <- vars::VAR(matrix, lag.max = 6, p=force.p)
  
  
 if(SVAR.switch==1){
  svar.fit<-SVAR(x = var.fit.1, estmethod = "scoring", Amat = Amat, Bmat = Bmat, max.iter = 500, maxls = 5000, conv.crit = 1.0e-8)
  var.fit.1<-svar.fit$var
  # Some specific results from SVAR estimation
  SVAR.matrices<-list(A=  svar.fit$A, Ase=  svar.fit$Ase,  B=svar.fit$B, Bse=svar.fit$Bse, sigma.U=svar.fit$Sigma.U)
  
}
  var.fit<-var.fit.1
  
  # Some test from var estimation
  
  serial.test<-serial.test(var.fit, lags.pt=10, type="PT.asymptotic")
  
  normality.test <- normality.test(var.fit)
  
  var.stability <- stability(var.fit, type = "Rec-CUSUM")
  

  
  
  # Impulse response function
  
  impulse.response <- irf(var.fit , n.ahead = 6, ortho = TRUE,
                          cumulative = ir_CUM, boot = TRUE, ci = 0.95, runs = 100)
  
  # Variance decomposition of the forescats
  
  var.decomp<-fevd(var.fit)
  
  # causality analysis
  
  causality<- list(
    cause.x=causality(var.fit,  cause=names[1],boot=TRUE, boot.runs=1000)$Granger,
    cause.y=causality(var.fit,  cause=names[2],boot=TRUE, boot.runs=1000)$Granger
    )
  
  # lags of the model 
  lags<-var.fit$p
  if(SVAR.switch==0) SVAR.matrices<-c("Not a SVAR")
  
  return(
    list(
      degree=degree,
      var.fit=var.fit,
      causality=causality,
      impulse.response=impulse.response,
      serial.test=serial.test,
      normality.test=normality.test,
      var.stability=var.stability,
      var.decomp=var.decomp,
      lags=lags,
      SVAR.matrices=SVAR.matrices
    )
  )
}

eval_var3<-function(X,ir_CUM=F, force.degree=NULL, force.p=NULL, Amat=NULL, Bmat=NULL){
  
  x<-X[,1]
  y<-X[,2]
  z<-X[,3]

  names<-colnames(X)
  
  SVAR.switch<-1 
  if(is.null(Amat) && is.null(Bmat)) SVAR.switch<-0  # If no A and B matrix Switch OFF deactivating SVAR estimation
  
  
  # Transfomr order of series force.degree or just making them stationary previously testing unitroot 
  
  
  if(is.null(force.degree)){
    
    dx<-0;dy<-0;dz<-0
    
    if(unitroot_ndiffs(ts(x), alpha=0.05)!=0){dx<-unitroot_ndiffs(ts(x), alpha=0.05);x<-diff(ts(x), differences=dx)}
    if(unitroot_ndiffs(ts(y), alpha=0.05)!=0){dy<-unitroot_ndiffs(ts(y), alpha=0.05);y<-diff(ts(y), differences=dy)}
    if(unitroot_ndiffs(ts(z), alpha=0.05)!=0){dz<-unitroot_ndiffs(ts(z), alpha=0.05);z<-diff(ts(z), differences=dz)}
    
    degree<-c(dx,dy,dz)
  }else{ 
    
    degree<-force.degree
    dx=degree[1]
    dy=degree[2]
    dz=degree[3]
  }
  
  maxd<-max(da,dy,dz)
  
  if(dx<maxd){x<-head(x,-I(maxd-dx))}
  if(dy<maxd){y<-head(y,-I(maxd-dy))}    
  if(dz<maxd){z<-head(z,-I(maxd-dz))}    
  
  
   
  assign(names[1],x)
  assign(names[2],y)
  assign(names[3],z)
 
  names(degree)<-names
  # Select order and estimate VAR
  
  matrix<-cbind(ts(x),ts(y),ts(z))
  colnames(matrix)<-names
  
  var.fit.1 <- vars::VAR(matrix, lag.max=6,p=force.p)
  
  if(SVAR.switch==1){
    svar.fit<-vars::SVAR(x = var.fit.1, estmethod = "scoring", Amat = Amat, Bmat = Bmat, max.iter = 500, maxls = 5000, conv.crit = 1.0e-8)
    var.fit.1<-svar.fit$var
    # Some specific results from SVAR estimation
    SVAR.matrices<-list(A=  svar.fit$A, Ase=  svar.fit$Ase,  B=svar.fit$B, Bse=svar.fit$Bse, sigma.U=svar.fit$Sigma.U)
    
  }
  var.fit<-var.fit.1
  
  
  # Some test form var estimation
  
  serial.test<-serial.test(var.fit, lags.pt=10, type="PT.asymptotic")
  
  normality.test <- normality.test(var.fit)
  
  var.stability <- stability(var.fit, type = "Rec-CUSUM")
  
  # Impulse response function
  
  impulse.response <- irf(var.fit , n.ahead = 6, ortho = TRUE,
                          cumulative = ir_CUM, boot = TRUE, ci = 0.95, runs = 100)
  
  # Variance decomposition of the forescats
  
  var.decomp<-fevd(var.fit)
  
  # causality analisys
  
  causality<- list(
    cause.x=causality(var.fit,  cause=names[1],boot=TRUE, boot.runs=1000)$Granger,
    cause.y=causality(var.fit,  cause=names[2],boot=TRUE, boot.runs=1000)$Granger,
    cause.z=causality(var.fit,  cause=names[3],boot=TRUE, boot.runs=1000)$Granger
    )
  
  # lags of the model 
  lags<-var.fit$p
  if(SVAR.switch==0) SVAR.matrices<-c("Not a SVAR")
  
  return(
    list(
      degree=degree,
      var.fit=var.fit,
      causality=causality,
      impulse.response=impulse.response,
      serial.test=serial.test,
      normality.test=normality.test,
      var.stability=var.stability,
      var.decomp=var.decomp,
      lags=lags,
      SVAR.matrices=SVAR.matrices
    )
  )
}

eval_var4<-function(X,ir_CUM=F, force.degree=NULL, force.p=NULL, Amat=NULL, Bmat=NULL){
  
  x<-X[,1]
  y<-X[,2]
  z<-X[,3]
  zz<-X[,4]
  names<-colnames(X)
  
  SVAR.switch<-1 
  if(is.null(Amat) && is.null(Bmat)) SVAR.switch<-0  # If no A and B matrix Switch OFF deactivating SVAR estimation
  
  # Make series to be stationary
  
  if(is.null(force.degree)){
    
    dx<-0;dy<-0;dz<-0;dzz<-0
    
    if(unitroot_ndiffs(ts(x), alpha=0.05)!=0){dx<-unitroot_ndiffs(ts(x), alpha=0.05);x<-diff(x, differences=dx)}
    if(unitroot_ndiffs(ts(y), alpha=0.05)!=0){dy<-unitroot_ndiffs(ts(y), alpha=0.05);y<-diff(y, differences=dy)}
    if(unitroot_ndiffs(ts(z), alpha=0.05)!=0){dz<-unitroot_ndiffs(ts(z), alpha=0.05);z<-diff(z, differences=dz)}
    if(unitroot_ndiffs(ts(zz), alpha=0.05)!=0){dzz<-unitroot_ndiffs(ts(zz), alpha=0.05);zz<-diff(zz, differences=dzz)}
    
    degree<-c(dx,dy,dz,dzz)
  }else{ 
    
    degree<-force.degree
    dx=degree[1]
    dy=degree[2]
    dz=degree[3]
    dzz=degree[3]
  }
  maxd<-max(dx,dy,dz,dzz)


  if(dx<maxd){x<-head(x,-I(maxd-dx))}
  if(dy<maxd){y<-head(y,-I(maxd-dy))}    
  if(dz<maxd){z<-head(z,-I(maxd-dz))}    
  if(dzz<maxd){zz<-head(zz,-I(maxd-dzz))}    
  
  assign(names[1],x)
  assign(names[2],y)
  assign(names[3],z)
  assign(names[4],zz)
  
  names(degree)<-names
  # Select order and estimate VAR
  
  matrix<-cbind(ts(x),ts(y),ts(z),ts(zz))
  colnames(matrix)<-names
  
  
  var.fit.1 <- vars::VAR(matrix, lag.max = 6, p=force.p)
  
  if(SVAR.switch==1){
    svar.fit<-vars::SVAR(x = var.fit.1, estmethod = "scoring", Amat = Amat, Bmat = Bmat, max.iter = 500, maxls = 5000, conv.crit = 1.0e-8)
    var.fit.1<-svar.fit$var
    # Some specific results from SVAR estimation
    SVAR.matrices<-list(A=  svar.fit$A, Ase=  svar.fit$Ase,  B=svar.fit$B, Bse=svar.fit$Bse, sigma.U=svar.fit$Sigma.U)
    
  }
  var.fit<-var.fit.1
  
  
  # Some test form var estimation
  
  serial.test<-serial.test(var.fit, lags.pt=10, type="PT.asymptotic")
  
  normality.test <- normality.test(var.fit)
  
  var.stability <- stability(var.fit, type = "Rec-CUSUM")
  
  # Impulse response function
  
  impulse.response <- irf(var.fit , n.ahead = 6, ortho = TRUE,
                          cumulative = ir_CUM, boot = TRUE, ci = 0.95, runs = 100)
  
  # Variance decomposition of the forescats
  
  var.decomp<-fevd(var.fit)
  
  # causality analisys
  
  causality<- list(
    cause.x=causality(var.fit,  cause=names[1],boot=TRUE, boot.runs=1000)$Granger,
    cause.y=causality(var.fit,  cause=names[2],boot=TRUE, boot.runs=1000)$Granger,
    cause.z=causality(var.fit,  cause=names[3],boot=TRUE, boot.runs=1000)$Granger,  
    cause.zz=causality(var.fit,  cause=names[4],boot=TRUE, boot.runs=1000)$Granger)

  # lags of the model 
  lags<-var.fit$p
  
  if(SVAR.switch==0) SVAR.matrices<-c("Not a SVAR") 
  
  
  return(
    list(
      degree=degree,
      var.fit=var.fit,
      causality=causality,
      impulse.response=impulse.response,
      serial.test=serial.test,
      normality.test=normality.test,
      var.stability=var.stability,
      var.decomp=var.decomp,
      lags=lags,
      SVAR.matrices=SVAR.matrices
      )
  )
}

list.eval_var<-list(eval_var2,eval_var3, eval_var4)