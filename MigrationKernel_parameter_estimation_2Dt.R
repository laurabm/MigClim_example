##############################################################
####  Scripts for migration kernel parameterisation in the manuscript:
####  Migration-based scenarios for Canadian trees show limited range tracking under climate change  

###https://cran.r-project.org/web/packages/rootSolve/index.html
install.packages("./rootSolve_1.8.2.3.tar.gz", repos=NULL, type="source")  ##Need Rtools to compile
install.packages("readxl", dependencies=T)

library(readxl)
library(rootSolve)

#### Aim of prog:
## To study the 2Dt distribution

#### Clear memory and graphs
rm(list = ls())
graphics.off()


setwd("") ##Set wd here

#### Define the distribution from the median and the mean:
#		Here, we must solve the system that link a and b (the parameters) with the median
#		and mean (Eq. 1 and 2). Unfortunately, the function gamma in Eq. 1 is non-linear
#		and there is no explicit reverse function for gamma. Therefore, we use a Newton-Raphson method to find our parameters


#### Comment:
#### Let us solve the equation μ = m Sqrt[π (2a + 1)]/(2 Γ[a]) Γ[a - 0.5] numerically (using multiroot function from rootSolve package)
## Define the function f(a) = m Sqrt[π (2a + 1)]/(2 Γ[a]) Γ[a - 0.5] - μ
f = function(a, parms)
{
  median = parms[1]
  mu = parms[2]
  return (median * sqrt(pi / (2^(1/a) - 1)) * gamma(a - 0.5)/(2 * gamma(a)) - mu)
}


# The problem with rootSolve is that it might not converge to the solution for 3 reasons:
#	1) R cannot evaluate the function gamma (after 172, R return Inf value)
#	2) multiroot is for unconstrained optimisation, but we have the constraint a > 0.5, and sometimes R goes to negative values, producing NaN
#	3) there is no solution...

## Test with all the species, using the median. Dispersal_metrics must be loaded
MigrationKernel = read.csv("./MigrationKernel.csv")  ### Input file already has estimates for a and b.

dispersal_metrics <- MigrationKernel[,1:6]
colnames(dispersal_metrics) = c("SpeciesCode", "SpeciesName", "disp_median", "disp_mean", "LDD_lwrCL", "LDD_uppCL")

list_index_problem = integer()
result_parameters = data.frame(SpeciesCode = dispersal_metrics$SpeciesCode, SpeciesName = dispersal_metrics$SpeciesName, a = numeric(nrow(dispersal_metrics)), b = numeric(nrow(dispersal_metrics)))

for (i in 1:nrow(dispersal_metrics))
{
  isThereResult = 1.065*as.numeric(dispersal_metrics[i, "disp_median"]) - as.numeric(dispersal_metrics[i, "disp_mean"]) < 0
  if (isThereResult)
  {
    print(paste0("running for ", dispersal_metrics[i, "SpeciesName"], "; median = ", dispersal_metrics[i, "disp_median"], " mean = ", dispersal_metrics[i, "disp_mean"]))
    median = as.numeric(dispersal_metrics[i, "disp_median"])
    mu = as.numeric(dispersal_metrics[i, "disp_mean"])
    a_estimate = multiroot(f = f, start = 0.63, maxiter = 1e3, useFortran = TRUE, parms = c(median, mu))$root
    b_estimate = median^2/(2^(1/a_estimate) - 1)
    print(paste0("a = ", a_estimate, "; b = ", b_estimate))
    result_parameters[i, "a"] = a_estimate
    result_parameters[i, "b"] = b_estimate
  }
  else
  {
    print(paste0("***   no result for ", dispersal_metrics[i, "SpeciesName"], "   ***"))
    list_index_problem = c(list_index_problem, i)
  }
}

## Trying in case of for species that have a problem
for (i in list_index_problem)
{
  print(paste0("running for ", dispersal_metrics[i, "SpeciesName"], "; median = ", dispersal_metrics[i, "disp_median"], " mean = ", dispersal_metrics[i, "disp_mean"]))
  median = as.numeric(dispersal_metrics[i, "disp_median"])
  mu = as.numeric(dispersal_metrics[i, "disp_mean"])
  a_estimate = multiroot(f = f, start = 1, maxiter = 1e3, useFortran = TRUE, parms = c(mode, mu))$root
  b_estimate = median^2/(2^(1/a_estimate) - 1)
  print(paste0("a = ", a_estimate, "; b = ", b_estimate))
}

###If satisfied, bind columns with parameters for use in next modelling step 
dispersal_metrics<-cbind(dispersal_metrics, result_parameters[,3:4])

write.csv(dispersal_metrics, file="./MigrationKernel.csv")   ##Write outputs to location of choice

#######################################################################################

####### Obtain probability of dispersal at least certain distance ####
## Calibrate kernel with estimated a and b for a specific species 
#### Example here using Abies balsamea, ABIBAL


set.seed(19690818) # Woodstock seed

## Parameters (from Greene2004)
a = dispersal_metrics[dispersal_metrics$SpeciesCode=="ABIBAL",7]
b = dispersal_metrics[dispersal_metrics$SpeciesCode=="ABIBAL",8]

LDD_lwrCL<-as.numeric(dispersal_metrics[dispersal_metrics$SpeciesCode=="ABIBAL",5])  ###Model estimate for family used instead of order
LDD_uppCL<-as.numeric(dispersal_metrics[dispersal_metrics$SpeciesCode=="ABIBAL",6])  ###Model estimate for family used instead of order


F2 = function(parms)
{
  a = parms[1]
  b = parms[2]
  f2 = function(x)
    2*a*x/b*(1 + x^2/b)^(-a-1)
  
  return(f2)
}

## Probability density function (pdf) of the 2Dt distribution
pdf_2Dt = function(x, a, b)
  return (2*a*x/b*(1 + x^2/b)^(-a-1))

## Cumulative distribution function (cdf) of the 2Dt distribution
cdf_2Dt = function(x, a, b)
  return (1 - (b/(b + x^2))^a)

## The distribution from Greene2004, with a slight correction. The power should be -S-1 but the graph is correct
distrib_2Dt = function(x, a, b)
  return (a/(pi*b)*(1 + x^2/b)^(-a-1))



### Calculate theoretical mean
theoretical_mean = a*sqrt(pi*b)/(2*gamma(1 + a))*gamma(a - 0.5)

#Proof that theoretical mean equals the observed mean
cat("Theoretical mean:", theoretical_mean)
cat("Mean observed:", as.numeric(result_parameters[result_parameters$SpeciesCode=="ABIBAL",4]))

## Theoretical median
#Proof that theoretical median equals the observed median
theoretical_meadiane = sqrt(b*(2^(1/a) - 1))
print(cdf_2Dt(theoretical_meadiane, a, b))  ##Proof that the probability is 0.5

cat("Theoretical median:", theoretical_meadiane)
cat("Median observed:", as.numeric(dispersal_metrics[dispersal_metrics$SpeciesCode=="ABIBAL",3]))


#### Plot distribution curve 
par(mfrow=c(3,1))
curve(expr = distrib_2Dt(x, a, b), from = 0, to = 1500, xlab = "Distance", ylab = "f(Distance)")
abline(v = theoretical_meadiane, col = "blue", lwd = 2)
abline(v = theoretical_mean, col = "red", lwd = 2)
#abline(h = 0.5, col = "red", lwd = 2)
abline(v = LDD_lwrCL, col = "grey60", lwd = 2)
abline(v = LDD_uppCL, col = "grey60", lwd = 2, lty = 2)

# Plot pdf 2Dt (should be like the histogram that follow in next part)
curve(expr = pdf_2Dt(x, a, b), from = 0, to = 1500, xlab = "Distance", ylab = "Probability f(Distance)")
abline(v = theoretical_meadiane, col = "blue", lwd = 2)
abline(v = theoretical_mean, col = "red", lwd = 2)
#abline(h = 0.5, col = "red", lwd = 2)
abline(v = LDD_lwrCL, col = "grey60", lwd = 2)
abline(v = LDD_uppCL, col = "grey60", lwd = 2, lty = 2)

# Cumulative distribution function which should tend to 0 when x ->0 and to 1 when x -> infinity
curve(expr = cdf_2Dt(x, a, b), from = 0, to = 1500, xlab = "Distance", ylab = "Cumulative probability f(Distance)")
abline(v = theoretical_meadiane, col = "blue", lwd = 2)
abline(v = theoretical_mean, col = "red", lwd = 2)
#abline(h = 0.5, col = "red", lwd = 2)
abline(v = LDD_lwrCL, col = "grey60", lwd = 2)
abline(v = LDD_uppCL, col = "grey60", lwd = 2, lty = 2)


###### Obtain probabilities for use in modelling ######
#### Scenario 1 - Probabilities calculated at distances that match grid resolution, up to set distance
##### Probabilties calculated for typical dispersal are calculated up to the lower CI of the maximum estimated dispersal
##### LDD_lwrCL obtained from the dispeRsal model (Tamme et al. 2014)
#### Modelling conducted at 25m resolution. 

LDD_lwrCL_25<-floor(as.numeric(dispersal_metrics[dispersal_metrics$SpeciesCode=="ABIBAL",5])/25)   #divide by 25 to match grid resolution


dist<-seq(0,LDD_lwrCL_25*25-1,by=25)
prob<-0

dist.prob<-cbind(dist,prob)

int2 = function (parms, minDist) 
  integrate(F2(parms),  lower = minDist, upper = Inf)$value

for(k in 1:length(dist.prob[,1]))
{
  dist.prob[k,2]<-int2(c(a,b), as.numeric(dist.prob[k,1])) 
}

View(dist.prob)

disp.kernel.spp<-dist.prob[,2]   #Extract vector of probabilities to use in modelling

disp.kernel.spp[1]<-1  ###In Migclim.migrate() function, dispKernel=disp.kernel.spp

##### Scenario 2 - Probabilities can also be calculated between set distances
##### MigClim uses minLDD and maxLDD to bound LDD, so we calculated a probability between these distances.
## LDD_lwrCL and LDD_uppCL obtained from the dispeRsal model (Tamme et al. 2014)

LDD_lwrCL_25<-floor(as.numeric(dispersal_metrics[dispersal_metrics$SpeciesCode=="ABIBAL",5])/25)   #divide by 25 to match grid resolution
LDD_uppCL_25<-ceiling(as.numeric(dispersal_metrics[dispersal_metrics$SpeciesCode=="ABIBAL",6])/25) #divide by 25 to match grid resolution
 
LDD.prob <- integrate(F2(c(a,b)),  lower = LDD_lwrCL_25*25, upper = LDD_uppCL_25*25)$value    ###In Migclim.migrate() function, lddFreq=LDD.prob argument



################################
####		END SCRIPT		####
################################
