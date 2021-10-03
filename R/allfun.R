source("ini_value.R") #setting of number of subjects, etc.
source("gpz.R")       #survival, density and random number of Gompertz     
source("f.R")         #Parametric Maximum Likelihood Estimation
source("T1_func.R")   #Parametric survival and density function 
source("dlr.R")       #Nonparametric estimator for survival (Turnbull)
source("dlrb.R")      #Turnbull for bootstrap samples
source("Turnbull3.R") #Turnbull function in ReIns package 
source("getsurv2.R")  #Turnbull function in interval package
source("inverse.R")   #obtain i.i.d. LB sample from Turnbull estimator

source("lamda_nbw.R") #bootstrap value of the test statistic when H0 is weibull
source("lamda_nbu.R") #bootstrap value of the test statistic when H0 is uniform
source("lamda_nbl.R") #bootstrap value of the test statistic when H0 is lognormal
source("lamda_nbg.R") #bootstrap value of the test statistic when H0 is gompertz

source("myfun.R")    #main test procedure without combination not specifying True
source("myfun1.R")    #LB test procedure without combination when H0 is weibull
source("myfun2.R")    #LB test procedure without combination when H0 is gompertz
source("myfun3.R")    #LB test procedure without combination when H0 is lognormal
source("myfun4.R")    #LB test procedure without combination when H0 is uniform
source("myfunc1.R")   #main test procedure with combination specifying True +5
source("myfunc2.R")   #main test procedure with combination specifying True +10
source("myfunc3.R")   #main test procedure with combination specifying True +15
source("myfunc4.R")   #main test procedure with combination specifying True +20
#source("bupa.R")     #re-get estimation of parameters if forgotten in myfun when True is weibull
#source("bupa2.R")    #re-get estimation of parameters if forgotten in myfun when True is gompertz
#source("bupa3.R")    #re-get estimation of parameters if forgotten in myfun when True is lognormal
#source("bupa4.R")    #re-get estimation of parameters if forgotten in myfun when True is uniform
#=============================library===========================
library(survival)
library(maxLik)
library(readr)        #write.csv
library(magrittr)     #%>%
library(pracma)       #Root finding in T1_func
library(ReIns)        #Turnbull3 and Kaplan3
library(tictoc)
library(interval)