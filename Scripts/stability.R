

var_stability<- function(X,periods=12){

# loading functions


# Creation of the data matrix
ts_matrix<-X

varmodel_list<-list() #Defines a list in which all the models with the rolling windows are later saved.
vardegree_list<-list()
varirf_list<-list()

Portmanteau<-NULL #Defines a vector in which the results of the Ljung-Box Portmanteau test will be collected.
Normality<-NULL
ARCH<-NULL

lm<-length(ts_matrix[,1])


eval_var<-list.eval_var[[I(ncol(X)-1)]]

for (j in 0:I(periods-1)){
  
  X<-ts_matrix[1:I(lm-j),] #This truncates the data matrix to form the rolling windows.
  
  var<-eval_var(X)
  var_model<-var$var.fit
  
  
  Portmanteau[j+1]<-serial.test(var_model, lags.pt=12, type="PT.asymptotic")$serial$p.value #Tests the independence of residuals using an asymptotic Portmanteau test; if this is significant, the model is invalid. In this case another lag for the VAR model may be tried, to be set in the control point against very high lags. 
  #Note that this sets a uniform lag for all rolling windows.
  Normality[j+1]<- normality.test(var_model)$jb.mul$JB$p.value#Test the normality of residuals using a multivariate Jarque-Bera test, a multivariate skewness test and a multivariate kurtosis test; best if not significant, but not crucial for valid models. 
  #Be prepared that it will be significant in models with a high number of variables.
  ARCH[j+1]<-arch.test(var_model, lags.multi = 5)$arch.mul$p.value #Shows the autoregressive conditional heteroskedasticity with a Lagrange multiplier; 
  #best if not significant, but not crucial for valid models.
  
  # The results of the test may also be plotted by plot(arch.test(var_model, lags.multi = 5)) , but this makes the run needing continuous user attention and generates a lot of extra plots. The same applied to stability tests plot(stability(var_macicarb), plot.type=c("s")) .
  #Interpreting the models using impulse-response functions.
  #variables Z, Z1, etc. is the impulse response function for impulse in quotation marks affecting the response variable in quotation marks. The impulse and the response variables should be the names from the model, i.e. those in the data matrix. Z numbers can be extended as desired by copying the command lines below. Alternatively, both impulse and response may be defined by column numbers of the data matrix (X), by replacing =”timeseries1” with =colnames(X)[a], where [a] is the number of the variable in the data matrix.
  varmodel_list[[j+1]]<-coeftest(var_model)
  vardegree_list[[j+1]]<-var$degree
  varirf_list[[j+1]]<-var$impulse.response
  
}

output_vars<-list(Portmanteau,
Normality,
ARCH,
varmodel_list,
vardegree_list, 
varirf_list
)

}