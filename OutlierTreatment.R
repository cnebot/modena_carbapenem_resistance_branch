

library(rstudioapi)
rstudioapi::getActiveDocumentContext
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


OutlierTreatment<- function(y,ld=5.rad=3){
  
#Generic program for outlier detection in Time series based on  
#Structural Change Chow test detection, stationary subseries 

#Parameters
level.detection=ld
r=rad  # radius for mean around outlier position for substitution 
 # which column of data is the time series to be considered 

 
 #Loading data
#input data file to be clean Data.csv 
 Data.OUT<-read.csv("Data.csv", header = T,sep=";",dec=".")

 name<-colnames(Data.OUT)[column]


 l<-length(Data.OUT[,column])
 y<-Data.OUT[,column]
 plot.ts(y)

#packages
#install.packages("strucchange")
library(strucchange)
#install.packages("car")
library(car)
#install.packages("tseries")
library(tseries)

# Looks for breakpoints

bp.y <- breakpoints(y ~ 1)
summary(bp.y)
split<-bp.y$breakpoints
ls<- length(split)

# Split into stationary subseries

split.components<-c(1,split,l+1)

par(mfrow=c(2,2))
plot.ts(y,main="All y")
for(i in 1:I(ls+1)){
  yy<-y[split.components[i]:I(split.components[i+1]-1)]
  nam <- paste("yy", i, sep = "")
  assign(nam, yy)

  plot.ts(yy,main = paste("Split y", i, sep = ""))

  #Augmented Dickey Fuller test
  nam.df <- paste("adf.yy", i, sep = "")
  k<-round((length(yy)-1)^(1/3))
  assign(nam.df, adf.test(yy,k=k))

}

 ## Print Augmented Dickey fuller test for subseries 
for(i in 1:I(ls+1)){
  nam.df <- paste("adf.yy", i, sep = "")
  print(get( nam.df))
}

par(mfrow=c(1,1))


# Outlier detection by subseries and substitute in ynew the mean value 


ynew<-y


for ( i in 1:I(ls+1)){
subseries<-get(paste("yy",i,sep=""))

mod<-lm(subseries~1)
outlierTest(mod)
cooksd1 <- cooks.distance(mod)
plot(cooksd1, pch="·", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = level.detection*mean(cooksd1, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd1)+1, y=cooksd1, labels=ifelse(cooksd1>level.detection*mean(cooksd1, na.rm=T),names(cooksd1),""), col="red")  # add labels

outliers1 <- cooksd1[cooksd1>I(level.detection*mean(cooksd1, na.rm=T))]
outlier.pos1<-as.vector(as.numeric(names(outliers1)))

if(length(outliers1)!=0){
  
  for(k in 1:length(outliers1)){
    j<-outlier.pos1[k]+(split.components[i]-1)
    ball<-y[I(j-r):I(j+r)]
    numball<-2*r-(sum(is.na(ball)))
    ynew[j]<-(sum(ball,na.rm=T)-y[j])/(numball)
  }
 }


}

par(mfrow=c(1,1))
plot.ts(y,ylim=c(0,max(y)))
plot.ts(ynew,ylim=c(0,max(y)))

file.name=paste("new",name,"(",level.detection,").csv",sep="")

changes<-1-(y==ynew) 
result<-cbind(y,ynew,changes )
colnames<-c("y","ynew", "changed")
write.table(result, file.name,col.names=TRUE,dec=".",sep=";",row.names = FALSE)

