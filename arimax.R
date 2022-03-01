
arimax <- function(data, vars){
  X <- data
  vars<- c("date", "gmneg_all_clinical", "fq",
           "cephs", "pseudo_3gc","noncarb","other_monitor_abx")
  
  X <- X %>% as_tsibble() %>% dplyr::select(vars)
  
  Xlong<- X %>% as_tsibble() %>% pivot_longer(cols=(!date), names_to="variable",values_to="value")
  
  Xlongd<- Xlong %>% as_tibble () %>% 
    group_by(variable) %>% 
    mutate(diff_var = difference(value))
  
  X<- X %>% ts()
  #----------------
  # Plots
  tsplot<-ggplot(Xlongd, aes(date, value, color=variable))+
    geom_line()+
    facet_grid(variable~., scales="free_y")+
    scale_y_continuous("Original variables")+
    theme_pnaat()+
    theme(panel.grid=element_blank())
  
  tsdplot <- ggplot(Xlongd, aes(date, 
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