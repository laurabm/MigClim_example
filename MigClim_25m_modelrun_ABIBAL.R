#install.packages ("MigClim", dependencies = TRUE)

install.packages(c("sp","raster","rasterVis", "maptools","rgeos","rgdal"), dependencies = T)
install.packages("https://cran.r-project.org/src/contrib/Archive/SDMTools/SDMTools_1.1-221.1.tar.gz", repo=NULL, type="source")   ###RTools needs to be installed first
install.packages("https://cran.r-project.org/src/contrib/Archive/MigClim/MigClim_1.6.tar.gz", repo=NULL, type="source")           ###RTools needs to be installed first

library(MigClim)
library(SDMTools)
library(sp)
library(raster)
library(rasterVis)
library(maptools)
library(rgeos)
library(rgdal)


############################################
###     MigClim example script 
### Example provided: Abies balsamea, ABIBAL
### Initial distribution and SDM based on fictitious distribution
############################################

###########START HERE#####################
setwd("") ##Set wd here

###### 1 - SETTING MODEL PARAMETERS ##########
MigKernel<-read.csv("./MigrationKernel.csv")
SeedProd<-read.csv("./SeedProduction.csv")

####Seed production#####

iniSeed.Age = as.numeric(SeedProd[SeedProd$SpeciesCode=="ABIBAL",3])    ###From trait data, ASM
optSeed.Age = as.numeric(SeedProd[SeedProd$SpeciesCode=="ABIBAL",4])  ###From trait data, SDOPT


iniSeed.Age<-round(iniSeed.Age)

x<-seq(0,optSeed.Age-iniSeed.Age, 1)

sigmoid = function(params, x) {
  params[1] / (1 + exp(-params[2] * (x - params[3])))
}

##ABIBAL: Slow growth
params_shape<-as.numeric(MigClim_parameters[MigClim_parameters$SpeciesCode=="ABIBAL",5])
ageInflection<-as.numeric(MigClim_parameters[MigClim_parameters$SpeciesCode=="ABIBAL",6])

params<-c(1,params_shape,ageInflection)   ## a = max prob; b = shape (inflection point); c = nb years passed at inflection (prob=0.5) point. Use plot to check


age.prob<-sigmoid(params,x)

plot(x, sigmoid(params,x), col='blue')


##sigmoid(params, 31)

### Migration kernel ####
 
###### Obtain probabilities for use in modelling ######
F2 = function(parms)
{
  a = parms[1]
  b = parms[2]
  f2 = function(x)
    2*a*x/b*(1 + x^2/b)^(-a-1)
  
  return(f2)
}

int2 = function (parms, minDist) 
  integrate(F2(parms),  lower = minDist, upper = Inf)$value


##### a and b are provided

a = as.numeric(MigKernel[MigKernel$SpeciesCode=="ABIBAL",7])
b = as.numeric(MigKernel[MigKernel$SpeciesCode=="ABIBAL",8])

LDD_lwrCL<-as.numeric(MigKernel[MigKernel$SpeciesCode=="ABIBAL",5])  ###Model estimate for family used instead of order
LDD_uppCL<-as.numeric(MigKernel[MigKernel$SpeciesCode=="ABIBAL",6])  ###Model estimate for family used instead of order

#### Probabilities calculated for typical dispersal distances


LDD_lwrCL_25<-floor(LDD_lwrCL/25)   #divide by 25 to match grid resolution (25m) and round down

dist<-seq(0,LDD_lwrCL_25*25-1,by=25)  ### Calculate probabilities for cell up to one before LDDmin. NOTE: dist = 0, prob = probability of falling between 0 and 25m... up to dist = LDD_lwrCL25 - 1, prob = probability of falling between (LDD_lwrCL_25-1)*25 and LDD_lwrCL_25*25
prob<-0

dist.prob<-cbind(dist,prob)


for(k in 1:length(dist.prob[,1]))
{
  dist.prob[k,2]<-int2(c(a,b), as.numeric(dist.prob[k,1])) 
}

View(dist.prob)

disp.kernel<-dist.prob[,2]   #Extract vector of probabilities to use in modelling

disp.kernel[1]<-1  ###MigClim needs a max dispersal probability of 1 (not 1.00000000...)

##### Probability of LDD between set distances. LDD_lwrCL and LDD_uppCL obtained from the dispeRsal model (Tamme et al. 2014)

LDD_lwrCL_25<-floor(LDD_lwrCL/25)   #divide by 25 to match grid resolution. Round down to near 25m. To be most inclusive as possible of LDD
LDD_uppCL_25<-ceiling(LDD_uppCL/25) #divide by 25 to match grid resolution. Round up to near 25m. To be most inclusive as possible of LDD

LDD.prob <- integrate(F2(c(a,b)),  lower = LDD_lwrCL_25*25, upper = LDD_uppCL_25*25)$value


####### 2 - RUN MIGCLIM ####################
library(MigClim)
simulName<-"ABIBAL_test"   ###Include name of climate suitability scenario. 

  MigClim.ABIBAL<-function(iniDist, hsMap,simulName){
    N<-MigClim.migrate(iniDist=iniDist,
                       hsMap=hsMap, rcThreshold=rc.cutoff.sp,  #rc.cutoff.sp = cutoff for climatic suitability/unsuitability. Varies between ANUCLIM and Maxent
                       envChgSteps=3, dispSteps=30, dispKernel=disp.kernel.ABIBAL,
                       barrier="LakesCan_25", barrierType="strong",
                       iniMatAge=iniSeed.Age, propaguleProd=age.prob,
                       lddFreq=LDD.prob, lddMinDist=LDD_lwrCL_25+1, lddMaxDist=LDD_uppCL_25,  #LDD_lwrCL_25+1 is to get next cell after last one in typical dispersal distance. In this case, next cell is between 275m and 300m from source 
                       simulName=simulName, replicateNb=2, overWrite=TRUE,
                       testMode=FALSE, fullOutput=FALSE, keepTempFiles=FALSE)
  }
  
  print(Sys.time())
  try(MigClim.ABIBAL("iniDist_25","hsMap_25_",simulName))
  print(Sys.time())


###### 3 - OPTIONAL: RECLASSIFY OUTPUTS TO YEAR COLONIZED 1 to 90 ####
##### Create year reclassification matrix ####
env.step1<-as.data.frame(seq(100,129, by=1))
env.step2<-as.data.frame(seq(200,229, by=1))
env.step3<-as.data.frame(seq(300,329, by=1))
colnames(env.step1)<-"x"
colnames(env.step2)<-"x"
colnames(env.step3)<-"x"

year<-as.data.frame(seq(1:90))

env.step11<-as.data.frame(seq(101,130, by=1))
env.step22<-as.data.frame(seq(201,230, by=1))
env.step33<-as.data.frame(seq(301,330, by=1))
colnames(env.step11)<-"x"
colnames(env.step22)<-"x"
colnames(env.step33)<-"x"

envstep<-rbind(env.step1,env.step2,env.step3)

envstep1<-rbind(env.step11,env.step22,env.step33)
rename<-cbind(envstep,envstep1,year)

iniDist<-c(0,1,-998)  ###IniDist value = -998
not.col<-c(29999,30001,NA)
rename<-rbind(iniDist,rename,not.col)

colnames(rename)<-c("start","finish","step")


rename<-as.matrix(rename)

##Reclassify raster. Add replicate number 
spp_proj<-raster("./ABIBAL_test/ABIBAL_test1_raster.asc")

proj_reclass_year<-reclassify(spp_proj, rename)   ###This is the final output

writeRaster(proj_reclass_year, "./ABIBAL_test/ABIBAL_test1_raster_reclass.tif")
