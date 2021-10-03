##basic constants for partly interval censored data with a cured fraction
n <- 5e2                   #sample size
m <- 20                    #total visit times for an individual
alpha1 <- 0.30             #percentage of exact observations
phi <- 0.80                #percentage of susceptible subjects
#Be careful that phi should be larger than alpha1 because
#cured ones cannot have exact observations

##sample size of various compositions of data
n1 <- round(n*alpha1)       #exact observations
n2 <- round(n*(1-alpha1))   #censored observations
n3 <- round(n*phi)          #susceptible subjects
n4 <- round(n*(1-phi))      #cured subjects

##true value of parameters for various distribution functions to be tested 
min_unif <- 0; max_unif <- 4              #Uniform
shape_weib <- 1.50; scale_weib <- 2.25    #Weibull
meanlog = 0.50; sdlog = 0.60              #Lognormal
labda = 0.16; gama = 0.75                 #Gompertz

##basic parameters for classical bootstrap
B <- 5e2                   #bootstrap times
classb <- rep(1:B,each=n)  #index for class of bootstrap sample

##basic constants for leveraged bootstrap
alpha <- 0.05
C_a <- .46136
e <- 0.05616
rou <- alpha/2

alpha_x <- 0.1
C_ax <- .34730
ex <- .03601
roux <- alpha_x/2

eta <- 0.10
epsilon <- 0.02
gaa <- 1/3
m1 <- n^gaa

##
tt <- 1e3                  #replicating times