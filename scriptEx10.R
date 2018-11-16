rm(list = ls())  # clears global environment

setwd("~/Biocomp-Fall2018-181109-Exercise10")  # sets working directory to exercise file

library(ggplot2)  # reads in the ggplot2 library
library(deSolve)  # reads in deSolve library for ode


# QUESTION 1

data <- read.table("data.txt", header = TRUE, sep = ",")  # reads in data.txt file

nllLinear <- function(p, x, y){  # creates likelihood function for a linear model 
  B0 = p[1]
  B1 = p[2]
  sigma = exp(p[3])
  expected = B0 + B1*x
  
  nll = -sum(dnorm(x = y, mean = expected, sd = sigma, log = TRUE))
  
  return(nll)
}

nllQuadratic <- function (p, x, y){  # creates likelihood function for a quadratic model
  B0 = p[1]
  B1 = p[2]
  B2 = p[3]
  sigma = exp(p[4])
  expected = B0 + B1*x + B2*x^2
  
  nll = -sum(dnorm(x = y, mean = expected, sd = sigma, log = TRUE))
  
  return(nll)
}

initialGuess = c(1, 1, 1)  # (B0, B1, sigma)
linearFit = optim(par = initialGuess, fn = nllLinear, x = data$x, y = data$y)  # optimization for linear model parameters

print(linearFit)  # prints linear fit model

initialGuess2 = c(1, 1, 1, 1)  # (B0, B1, B2, sigma)
quadraticFit = optim(par = initialGuess2, fn = nllQuadratic, x = data$x, y = data$y)  # optimization for quadratic model parameters

print(quadraticFit)  # prints quadratic fit model

testStat = 2 * (linearFit$value - quadraticFit$value)  # test statistic for models
df = length(quadraticFit$par) - length(linearFit$par)  # degrees of freedom for models

1 - pchisq(testStat, df)  # 1 - chi squared p-value to determine the likelihood of accepting the null hypothesis

# value is to rest between 0 and 1 (in this case, value is close to 1)
# value > 0.5 (linear is better), value < 0.5 (quadratic is better)
# quadratic model is bad fit in this case since the value is greater than 0.5


# QUESTION 2

# Lotka-Volterra model describes competition between two species

simFunc <- function(t, y, p){  # creates function that uses the Lotka-Volterra model for competition using two diff equations
  N1 = y[1]
  N2 = y[2]
  
  R1 = p[1]
  alpha11 = p[2]
  alpha12 = p[3]
  
  R2 = p[4]
  alpha22 = p[5]
  alpha21 = p[6]
  
  dN1dt = R1 * (1 - N1*alpha11 - N2*alpha12) * N1  # first diff equation
  dN2dt = R2 * (1 - N2*alpha22 - N2*alpha21) * N2  # second diff equation
  
  return(list(c(dN1dt, dN2dt)))
}

params = c(0.52, 0.05, 0.009, 0.52, 0.05, 0.007)  # alpha12 < alpha11, alpha21 < alpha22
params2 = c(0.43, 0.05, 0.001, 0.43, 0.05, 0.002)  # alpha12 < alpha11, alpha21 < alpha22
params3 = c(0.26, 0.005, 0.01, 0.26, 0.005, 0.02) # alpha12 > alpha11, alpha21 > alpha22
params4 = c(0.33, 0.005, 0.03, 0.33, 0.005, 0.002)  # alpha12 > alpha11, alpha21 < alpha22
times = 1:100  # time values used over the course of model, held constant among all models
y0 = c(1, 1)  # initial values for populations, held constant among all models

modelSim = ode(y = y0, times = times, func = simFunc, parms = params)  # first model simulation, uses ode function to optimize function and parameters
modelOutput = data.frame(time = modelSim[, 1], N1 = modelSim[, 2], N2 = modelSim[, 3])  # puts model output into data frame

ggplot(modelOutput, aes(x = time, y = N1)) +  # creates ggplot that uses time as x axis and N1 as y axis
  geom_line() + 
  geom_line(modelOutput, mapping = aes(x = time, y = N2), col = 'red') + 
  theme_classic()

modelSim2 = ode(y = y0, times = times, func = simFunc, parms = params2) # second model simulation, uses ode function to optimize function and parameters
modelOutput2 = data.frame(time = modelSim2[, 1], N1 = modelSim2[, 2], N2 = modelSim2[, 3])  # puts model output into data frame

ggplot(modelOutput2, aes(x = time, y = N1)) +  # creates ggplot that uses time as x axis and N1 as y axis
  geom_line() + 
  geom_line(modelOutput2, mapping = aes(x = time, y = N2), col = 'green') +
  theme_classic()

modelSim3 = ode(y = y0, times = times, func = simFunc, parms = params3)  # third model simulation, uses ode function to optimize function and parameters
modelOutput3 = data.frame(time = modelSim3[, 1], N1 = modelSim3[, 2], N2 = modelSim3[, 3])  # puts model output into data frame

ggplot(modelOutput3, aes(x = time, y = N1)) +  # creates ggplot that uses time as x axis and N1 as y axis
  geom_line() + 
  geom_line(modelOutput3, mapping = aes(x = time, y = N2), col = 'blue') +
  theme_classic()

modelSim4 = ode(y = y0, times = times, func = simFunc, parms = params4)  # fourth model simulation, uses ode function to optimize function and parameters
modelOutput4 = data.frame(time = modelSim4[, 1], N1 = modelSim4[, 2], N2 = modelSim4[, 3])  # puts model output into data frame

ggplot(modelOutput4, aes(x = time, y = N1)) +  # creates ggplot that uses time as x axis and N1 as y axis
  geom_line() + 
  geom_line(modelOutput3, mapping = aes(x = time, y = N2), col = 'orange') +
  theme_classic()

# The L-V model holds that if alpha11 < alpha12 and alpha21 < alpha22, then the two populations can coexist.
# In this case, modelSim and modelSim2 are both able to coexist as they follow the rules. Their graphs follow each other fairly well.
# However, modelSim3 and modelSim4 show dominant populations with modelSim3 showcasing a population that levels out early (blue), 
# and modelSim4 showcasing a population (black) that is torn apart to extinction.