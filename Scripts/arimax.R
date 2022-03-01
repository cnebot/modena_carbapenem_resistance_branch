
## ENTERS

data<- # name of the dataframe or ts or tsibble object
vars <- # a alist of variables such as c("var1", "var2",....)
data<-ts_p2
vars<-c("carbapenem", "ecoli_esbl_all_clinical")


arimax <- function(data, vars, startyear=2008, startmonth=01, frequency=12){

# create temporary data object
  date<- make_yearmonth(year = startyear, month = startmonth)
  dat <- data %>% as_tsibble() %>% filter_index(~ date) %>%
    dplyr::select(vars)

# pivot to long format for plotting purposes
  datlong<- dat %>%
    pivot_longer(cols=everything(),
    names_to="variable",values_to="value")
  
  datlongd<- datlong %>% as_tibble () %>% 
    group_by(variable) %>% 
    mutate(diff_var = difference(value))
  
# transform to a "ts" object to facilitate further analysis
datts<- ts(dat, start=c(as.numeric(start[[1]]),as.numeric(start[[2]])), frequency=frequency)
  
# Plots
  tsplot<-ggplot(datlongd, aes(date, value, color=variable))+
    geom_line()+
    facet_grid(variable~., scales="free_y")+
    scale_y_continuous("Original variables")+
    theme_pnaat()+
    theme(panel.grid=element_blank())
  
  tsdplot <- ggplot(datlongd, aes(date, 
                                diff_var, color=variable))+
    geom_line()+
    geom_hline(yintercept=0, linetype="dashed",
               color="grey")+
    facet_grid(variable~., scales="free_y")+
    scale_y_continuous("First differenced variables")+
    theme_pnaat()+
    theme(panel.grid=element_blank())
  
  plot<- tsplot+tsdplot+plot_layout(guides="collect")
  
  #----------------------------------------
  # Check for stationarity
  
  
  #----------------------------------------
  # Check for stationarity
  
  adf.matrix<-rep(0,2) 
  
  adfs_and_diffs <-function(X){
    for(i in 1:length(X[1,])){
      endo_var<- X[,i]
      adf1<- adfTest(endo_var, type="c")@test$p.value 
      adf2<- adfTest(endo_var, type="ct")@test$p.value 
      adfs <-cbind(adf1,adf2)
      adf.matrix <-rbind(adf.matrix, adfs) 
    }
    
    adf.matrix<-adf.matrix[-1,]
    nam<-c(paste0(colnames(X)))
    adftab<-cbind(nam,adf.matrix) %>% as_tibble() 
    adftab <-adftab %>% mutate(conclusion=
                                 if_else(adf1<0.05,"stationary in level",
                                         if_else(adf2<0.05, "stationary in trend","nonstationary")))
  }
  
  adfs <- adfs_and_diffs(X)
  adfs
  
  #