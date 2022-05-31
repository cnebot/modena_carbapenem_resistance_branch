#Loading the necessary packages. These should be installed to R previously.
require(fUnitRoots)
require(vars)

# Variables
vars<- c("cskp_all_clinical", "gmneg_all_clinical", "all_monitor_abx",
         "ecoli_esbl_all_clinical")

varspick<- ts_p2 %>% dplyr::select(vars) %>% as_tibble() %>% dplyr::select(-date)
ts_matrix <- ts(varspick,start=c(2008,1),frequency=12)
#Creation of the data matrix
ts_matrix<-ts.intersect(varspick) #List of the time series to be modelled in VAR in parenthesis. All should be of the same length. Theoretically, any number of time-series may be included, though very high number of series leads to model instability.
cor(ts_matrix) #Printing the correlation matrix

end_year<-2013 #Defining the last year in the series. This cannot exceed the final year in the input time series.
end_month<-12 #Defining the last month in the series. By default, the rolling windows start from this month back to January.
IRF<-IRF1<-IRF2<-IRF3<-matrix(0,13,end_month+1) #Creation of the output matrices for the impulse response functions, their means and confidence intervals of all rolling windows. The number of output matrices can be increased as desired by copying the rows and creating further output matrices with different names. Note that R is case-sensitive and do not mistake the matrix names (all caps) for the impulse-response function (lower case irf) below.
IRF_mean<-IRF1_mean<-IRF2_mean<-IRF3_mean<-IRF4_mean<-IRF5_mean<-matrix(0,13,end_month+1)
IRF_upper<-IRF_lower<-matrix(0,13,end_month)
IRF1_upper<-IRF1_lower<-matrix(0,13,end_month)
IRF2_upper<-IRF2_lower<-matrix(0,13,end_month)
IRF3_upper<-IRF3_lower<-matrix(0,13,end_month)

varmodel_list<-list() #Defines a list in which all the models with the rolling windows are later saved.
Portmanteau<-NULL #Defines a vector in which the results of the Ljung-Box Portmanteau test will be collected.
#Defining the cycle for the rolling windows


for (j in end_month:1){
  X<-window(ts_matrix, start=c(2008,1), end=c(end_year, j))  #This truncates the data matrix to form the rolling windows.

  VARtype<-as.vector(0) #This defines a vector which monitors whether there is trend in any of the series.
  #A cycle to test stationarity for the series in the data matrix
  for (i in 1:ncol(X)){
    print(colnames(X)[i]) #Prints the name of the time-series tested for stationarity to the R console.
    options(warn=-1)
    unitroot_trend<-adfTest(X[,i], type="c")@test$p.value #Checks for trend-stationarity of all series in the rolling window separately; if significant, means lack of trend in the series.
    unitroot_const<-adfTest(X[,i], type="ct")@test$p.value #Checks for stationarity of all series in the rolling window separately; if significant, means stationarity of the series and lack of trend.
    #If both this and the former significant, the model type is with constant („const”); if this is significant, but the former is not, then the model type is with trend and constant together („both”); if both this and the former is not significant, then the series is not stationary and a VAR model cannot be estimated, in this case „non-stationary series message is printed. This is set automatically below.
    ifelse((unitroot_const>0.05), print("non-stationary series"),
           ifelse((unitroot_trend>0.05), VARtype<-VARtype+1, VARtype<-VARtype))
    print(unitroot_trend) #Prints the probability (p value) for trend-stationarity to the R console.
    print(unitroot_const) #Prints the probability (p value) for stationarity.
    options(warn=0)}
  
  if(VARtype>0){type1<-"both"} else {type1<-"const"} #Sets the type of VAR model, i.e. with constant or with constant AND trend.
  print(type1) #Prints the type calculated from the stationarity check to the R console.

  
  p<-VARselect(X, lag.max=12, season=12, type=type1)$selection #Calculates the optimal lag.

  print(p) #Prints the optimal lag of the VAR model recommended by the different criteria. Akaike Information Criterion (AIC) is the best estimator for most monthly time-series, and the model in this script uses the lag suggested by the AIC.
  ifelse(p[1]>9, p[1]<-2, p[1]<-p[1]) #A control point to stop building models with very high lags (more than 9), which are expected to be invalid models. If the optimal lag calculated from Akaike Information Criteria is higher than nine, it is artificially set at lag 2. By changing 2 to any other number the lag number can be set manually here.
  ifelse(p[1]==1, p[1]<-2, p[1]<-p[1]) #A control point suppressing lag 1
  var_model<-VAR(X, p=p[1], season=12, type=type1) #This estimates the VAR model with the lag and type set above. A yearly (12-month) seasonality is assumed.
  lags[j]<-var_model$p #Collects the lags of the different rolling windows.

  print(c("month: ", j), quote=F) #Prints the end month of the present rolling window.
  print(c("lag: ", lags[j]), quote=F) #Prints the lag used in the model.

  print(serial.test(var_model, lags.pt=12, type="PT.asymptotic")) #Tests the independence of residuals using an asymptotic Portmanteau test; if this is significant, the models is invalid. In this case another lag for the VAR model may be tried, to be set in the control point against very high lags. Note that this sets a uniform lag for all rolling windows.
  print(normality.test(var_model)$jb.mul) #Test the normality of residuals using a multivariate Jarque-Bera test, a multivariate skewness test and a multivariate kurtosis test; best if not significant, but not crucial for valid models. Be prepared that it will be significant in models with a high number of variables.
  print(arch.test(var_model, lags.multi = 5)$arch.mul) #Shows the autoregressive conditional heteroskedasticity with a Lagrange multiplier; best if not significant, but not crucial for valid models.
  # The results of the test may also be plotted by plot(arch.test(var_macicarb, lags.multi = 5)) , but this makes the run needing continuous user attention and generates a lot of extra plots. The same applied to stability tests plot(stability(var_macicarb), plot.type=c("s")) .
  #Interpreting the models using impulse-response functions.
}
  #variables Z, Z1, etc. is the impulse response function for impulse in quotation marks affecting the response variable in quotation marks. The impulse and the response variables should be the names from the model, i.e. those in the data matrix. Z numbers can be extended as desired by copying the command lines below. Alternatively, both impulse and response may be defined by column numbers of the data matrix (X), by replacing =”timeseries1” with =colnames(X)[a], where [a] is the number of the variable in the data matrix.
  Portmanteau[j]<-serial.test(var_model, lags.pt=12, type="PT.asymptotic")$serial$p.value #Collects the results of the Ljung-Box Portmanteau test for serial correlation
  varmodel_list[[j]]<-var_model
  Z<-irf(var_model, n.ahead=12, impulse=colnames(X)[1], response=colnames(X)[2], ortho=TRUE, cumulative=TRUE, boot=TRUE, runs=100) #This calculates the impulse-response function. By changing the index in colnames(X)[], we can define which variable is the impulse and which is the reponse; in this case the first variable in the matrix is the impulse and the second is the response. Not below in the next cycle (Z1) that they are exchanged to test for reciprocal effect. Alternatively, impulse="" is also possible when the name of the series for impulse (or for reponse) should be written between the quotation marks. Presently this script provides opportunity to test a three-variable model's all potential combinations.
  plot(Z) #This plots the impulse-response function. The impulse has a significant effect on the response if, and at the lags where, the zero is outside the confidence intervals (red dotted lines).
  #This cycle below stores the significance, mean and confidence intervals of the impulse response in the respective IRF variable; note that the first value corresponds to lag 0, not to lag 1, of the response horizon.
  for (k in 1:13) {
    ifelse(as.vector(Z$Upper[[1]][k])<0,
           ifelse(as.vector(Z$Lower[[1]][k])<0, IRF[k,j]<--1, IRF[k,j]<-0),
           ifelse(as.vector(Z$Lower[[1]][k])>0, IRF[k,j]<-1, IRF[k,j]<-0))
    IRF[k,(end_month+1)]<-sum(IRF[k,1:end_month])/end_month
    IRF_mean[k,j]<-Z$irf[[1]][k,]
    IRF_upper[k,j]<-Z$Upper[[1]][k,]
    IRF_lower[k,j]<-Z$Lower[[1]][k,]
    IRF_mean[k, (end_month+1)]<-max(IRF_mean[k,])

  #This is the end of the cycle for a single impulse-reponse pair.
  Z1<-irf(var_model, n.ahead=12, impulse=colnames(X)[2], response=colnames(X)[1], ortho=TRUE, cumulative=TRUE, boot=TRUE, runs=100)
  plot(Z1)
  for (k in 1:13) {
    ifelse(as.vector(Z1$Upper[[1]][k])<0,  
           ifelse(as.vector(Z1$Lower[[1]][k])<0, IRF1[k,j]<--1, IRF1[k,j]<-0),
           ifelse(as.vector(Z1$Lower[[1]][k])>0, IRF1[k,j]<-1, IRF1[k,j]<-0))
    IRF1[k, (end_month+1)]<-sum(IRF1[k,1:end_month])/end_month
    IRF1_mean[k,j]<-Z1$irf[[1]][k,]
    IRF1_upper[k,j]<-Z1$Upper[[1]][k,]
    IRF1_lower[k,j]<-Z1$Lower[[1]][k,]
    IRF1_mean[k, (end_month+1)]<-max(IRF1_mean[k,])
  }
  Z2<-irf(var_model, n.ahead=12, impulse= colnames(X)[1], response= colnames(X)[3], ortho=TRUE, cumulative=TRUE, boot=TRUE, runs=100)
  plot(Z2)
  for (k in 1:13) {
    ifelse(as.vector(Z2$Upper[[1]][k])<0,  
           ifelse(as.vector(Z2$Lower[[1]][k])<0, IRF2[k,j]<--1, IRF2[k,j]<-0),
           ifelse(as.vector(Z2$Lower[[1]][k])>0, IRF2[k,j]<-1, IRF2[k,j]<-0))
    IRF2[k, (end_month+1)]<-sum(IRF2[k,1:end_month])/end_month
    IRF2_mean[k,j]<-Z1$irf[[1]][k,]
    IRF2_upper[k,j]<-Z1$Upper[[1]][k,]
    IRF2_lower[k,j]<-Z1$Lower[[1]][k,]
    IRF2_mean[k, (end_month+1)]<-max(IRF2_mean[k,])
  }
  Z3<-irf(var_model, n.ahead=12, impulse= colnames(X)[3], response= colnames(X)[1], ortho=TRUE, cumulative=TRUE, boot=TRUE, runs=100)
  plot(Z3)
  for (k in 1:13) {
    ifelse(as.vector(Z3$Upper[[1]][k])<0,  
           ifelse(as.vector(Z3$Lower[[1]][k])<0, IRF3[k,j]<--1, IRF3[k,j]<-0),
           ifelse(as.vector(Z3$Lower[[1]][k])>0, IRF3[k,j]<-1, IRF3[k,j]<-0))
    IRF3[k, (end_month+1)]<-sum(IRF3[k,1:end_month])/end_month
    IRF3_mean[k,j]<-Z1$irf[[1]][k,]
    IRF3_upper[k,j]<-Z1$Upper[[1]][k,]
    IRF3_lower[k,j]<-Z1$Lower[[1]][k,]
    IRF3_mean[k, (end_month+1)]<-max(IRF3_mean[k,])
  }
  Z4<-irf(var_model, n.ahead=12, impulse= colnames(X)[2], response= colnames(X)[3], ortho=TRUE, cumulative=TRUE, boot=TRUE, runs=100)
  plot(Z4)
  for (k in 1:13) {
    ifelse(as.vector(Z4$Upper[[1]][k])<0,  
           ifelse(as.vector(Z4$Lower[[1]][k])<0, IRF4[k,j]<--1, IRF4[k,j]<-0),
           ifelse(as.vector(Z4$Lower[[1]][k])>0, IRF4[k,j]<-1, IRF4[k,j]<-0))
    IRF4[k, (end_month+1)]<-sum(IRF4[k,1:end_month])/end_month
    IRF4_mean[k,j]<-Z1$irf[[1]][k,]
    IRF4_upper[k,j]<-Z1$Upper[[1]][k,]
    IRF4_lower[k,j]<-Z1$Lower[[1]][k,]
    IRF4_mean[k, (end_month+1)]<-max(IRF4_mean[k,])
  }
  Z5<-irf(var_model, n.ahead=12, impulse= colnames(X)[3], response= colnames(X)[2], ortho=TRUE, cumulative=TRUE, boot=TRUE, runs=100)
  plot(Z5)
  for (k in 1:13) {
    ifelse(as.vector(Z5$Upper[[1]][k])<0,  
           ifelse(as.vector(Z5$Lower[[1]][k])<0, IRF5[k,j]<--1, IRF5[k,j]<-0),
           ifelse(as.vector(Z5$Lower[[1]][k])>0, IRF5[k,j]<-1, IRF5[k,j]<-0))
    IRF5[k, (end_month+1)]<-sum(IRF5[k,1:end_month])/end_month
    IRF5_mean[k,j]<-Z1$irf[[1]][k,]
    IRF5_upper[k,j]<-Z1$Upper[[1]][k,]
    IRF5_lower[k,j]<-Z1$Lower[[1]][k,]
    IRF5_mean[k, (end_month+1)]<-max(IRF5_mean[k,])
  }
}

print(c("lags: ", lags), quote=F) #Prints the lags for all rolling windows from the shortest to the longest.
print(Portmanteau) #prints the results of the Ljung-Box Portmanteau test for the rolling windows
print(c("effect of", Z$impulse, "on", Z$response), quote=F) #Prints the name of the impulse and the response variable.
print(IRF) #Prints the significance of reponse to all rolling windows (in columns) and all lags (in rows, note that the first row correspond to lag 0 of the response horizon); 0 means no response, 1 means that the impulse leads to a significant increase in the response variable, -1 means that the impulse leads to a significant decrease in the response variable.
print(c("effect of", Z1$impulse, "on", Z1$response), quote=F)
print(IRF1)
print(c("effect of", Z2$impulse, "on", Z2$response), quote=F)
print(IRF2)
print(c("effect of", Z3$impulse, "on", Z3$response), quote=F)
print(IRF3)
print(c("effect of", Z4$impulse, "on", Z4$response), quote=F)
print(IRF4)
print(c("effect of", Z5$impulse, "on", Z5$response), quote=F)
print(IRF5)
print(IRF_mean) #prints the irf coefficients
print(IRF_upper) #prints the irf upper confidence interval values
print(IRF_lower) #prints the irf lower confidence interval values
print(IRF1_mean)
print(IRF1_upper)
print(IRF1_lower)
print(IRF2_mean)
print(IRF2_upper)
print(IRF2_lower)
print(IRF3_mean)
print(IRF3_upper)
print(IRF3_lower)
print(IRF4_mean)
print(IRF4_upper)
print(IRF4_lower)
print(IRF5_mean)
print(IRF5_upper)
print(IRF5_lower)
