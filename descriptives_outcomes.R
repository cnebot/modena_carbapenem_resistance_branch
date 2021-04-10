#################################################################################
# MULTI-MODAL CONTROL OF CARBAPENEM RESISTANT ORGANISMS (CRO) IN MODENA, ITALY.##
# SCRIPT (1): DESCRIPTIVE STATISTICS AND PLOTS FOR OUTCOMES                    ##
#################################################################################

### 1. SET-UP ###
## 1.1 Load packages required (use install.packages if not already installed) ##
library(here)
library(tidyverse)
library(lubridate)
library(timeSeries)
library(forcats)
library(ggplot2)
library(viridis)
library(ggfortify)
library(lme4)
library(lmtest)


## 1.2 check that your working directory is 'modena_carbapenem_resistance' ##
# recommend open project in new session to open 'modena_carbapeme_resistance' pjt from #
# location where stored the zip file download from github #
here () #check returned as [1] "C:/Users/User/Documents/R/modena_carbapenem_resistance"#

## 1.3 read in data files ##
outcome <- read_csv("outcomets.csv",
                    col_types = cols(
                      phases = col_factor (),
                      year = col_integer(),
                      month=col_integer(),
                      day=col_integer(),
                      date=col_date(format="%d/%m/%Y")))

parse_date("01/02/2015","%d/%m/%Y")


### 4. FIGURES ###
p <- ggplot(outcome, aes(x=date, y=cspa_ho_clinical_incid)) +
  geom_line(mapping=aes(x=date, y=cspa_ho_clinical_incid), color="darkgrey") +
  geom_point(,colour="darkgrey") +
  xlab("Year")+ ylab("CSPsA cases per 10 000 OBDs")
p
