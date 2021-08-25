####################################################################################
## Script includes function defining general ggplot2 theme for Cochrane DTA review##
####################################################################################

# Function
theme_pnaat <- function(){ 
  font <- "Source Sans Pro"   #assign font family up front
  
  theme_grey() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.minor = element_blank(),    #strip major gridlines
      panel.background = element_rect(fill = "#fafafa", colour = NA),
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 10,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 10),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 10,                 #font size
        hjust = 1),               #right align
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 10),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 9),                #font size

      legend.text = element_text(
        family= font,
        size=10)
    )
}
