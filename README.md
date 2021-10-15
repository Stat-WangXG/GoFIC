# Goodness_of_fit-Inteval_Censored-cure
Supplement to 'Goodness-of-fit test for a parametric mixture cure model under partly interval-censored data'.

R code of data analysis for Ziqi Geng, Jialiang Li, Yi Niu and Xiaogung Wang. Goodness-of-fit test for a parametric mixture cure model under partly interval-censored data. Biometrical Journal. 

Here is a general explanation of what is contained in the folders.

## R
gof_ic.R:     Main R application for testing the null hypothesis.  
gqz.R:        Density, distribution function and random generation for the Gompertz distribution with parameters shape and scale.  
f.R:          Log-likelihood function.  
dlrb_n.R:     Generation of bootstrap samples.  
lamda_nbwB.R: The bootstrap test statistic calculated by assuming Weibull-distributed original data.  
lamda_nbgB.R: The bootstrap test statistic calculated by assuming Gompertz-distributed original data.  
lamda_nblB.R: The bootstrap test statistic calculated by assuming lognormal-distributed original data.  
Turnbull3.R:  A slight transformation of 'Turnbull' function from 'ReIns' package for the purpose of vectorization.

## data
sc.dat: A PIC version of the data (7.6% exact observations) similarly to Gao et al.(2017) by treating observations with interval lengths less than or equal to 0.962 years as exact observations at the geometric mid-point of the interval to obtain 17 exact observations.  
actg.dat: Interval-censored observations with interval lengths less than 2 months are treated as exact observations at the geometric mid-point of the interval to obtain 46 exact observations.
