#==============================================================#
# Generic program for outlier detection in Time series based on  
# Structural Change Chow test detection, stationary subseries 
# y must be a 
#==============================================================#

# 1. SETUP  
library(rstudioapi)
rstudioapi::getActiveDocumentContext
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(warn=-1)



# 2. Packages
library(ggplot2)       # install.packages("ggplot2")
library(strucchange)   # install.packages("strucchange")
library(car)           # install.packages("car")
library(tseries)       # install.packages("tseries")



OutlierTreatment2 <- function(y,freq=12,start.date,level.detection=5,radius=2, adf.tests=FALSE,only.clean=TRUE){

if(!is.ts(y)) y<-ts(y,frequency=freq, start=start.date)  
  
  
# Parameters
ld<-level.detection         # Cook's distance to determine Outlier 
r<-radius                   # radius for mean around outlier position for substitution 

# Characteristics of the series
l<-length(y)

# Looking for breakpoints for main segments

bp.y <- breakpoints(y ~ 1)
summary(bp.y)
breakpoints(bp.y)
split<-bp.y$breakpoints

if(!is.na(split)){
breakdates(bp.y)
splitd<-breakdates(bp.y)
ci.y <- confint(bp.y)

ls<- length(split)

fm0 <- lm(y ~ 1)
fm1 <- lm(y ~ breakfactor(bp.y, breaks = length(split)))

## confidence interval
ci.y <- breakdates(confint(bp.y))


# Split into subseries

split.components<-c(1,split,l+1)

for(i in 1:I(ls+1)){
  yy<-y[split.components[i]:I(split.components[i+1]-1)]
  nam <- paste("yy", i, sep = "")
  assign(nam, yy)

    #Augmented Dickey Fuller test
  nam.df <- paste("adf.yy", i, sep ="")
  k<-round((length(yy)-1)^(1/3))
  assign(nam.df, adf.test(yy,k=k))
}

if(adf.tests==TRUE){
## Print Augmented Dickey fuller test for subseries 
for(i in 1:I(ls+1)){
  nam.df <- paste("adf.yy", i, sep = "")
  print(get(nam.df))
  }
 }
}else{ls<-0; yy1<-y;split.components<-c(1,l+1)}

# Outlier detection by subseries and substitute in ynew the mean value 
ynew<-y

for (i in 1:I(ls+1)){
subseries<-get(paste("yy",i,sep=""))

mod<-lm(subseries~1)
outlierTest(mod)
cooksd1 <- cooks.distance(mod)


outliers1 <- cooksd1[cooksd1>I(ld*mean(cooksd1, na.rm=T))]
outlier.pos1<-as.vector(as.numeric(names(outliers1)))

if(length(outliers1)!=0){
  
  for(k in 1:length(outliers1)){
    j<-outlier.pos1[k]+(split.components[i]-1)
    ball<-y[max(I(j-r),1):I(j+r)]
    numball<-2*r-(sum(is.na(ball)))
    ynew[j]<-(sum(ball,na.rm=T)-y[j])/(numball)
    }
  }


}

plot<-autoplot(cbind(y,ynew), xlab = "Time", ylab = deparse(substitute(object)))
               
changes<-1-(y==ynew) 

result<-list("y"=y,"ynew"=ynew,"changes"=changes,"plot"=plot)
if(only.clean==TRUE) result<-as.vector(ynew)

return(result)
}
