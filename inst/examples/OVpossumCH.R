#OVpossum example.
#Adapted from an example given in 
#https://www.otago.ac.nz/density/pdfs/openCR-examples.pdf

library(openpopscr)
library(secr)
library(rgdal)
library(RcppParallel)
setThreadOptions(defaultNumThreads() - 1)

data <- OVpossumCH
#Days between sessions according to pdf
intervals(data) <- c(123, 91, 154, 118, 85)/365

#Join multiple sessions into one
data <- join(data)

#Ignore Recoveries
data[data < 0] <- 1

#Extract the detectors
detectors <- traps(data)

#Set up primary occasions
primary <- c(rep(1:6, each = 5))

#Number of days between sessions.
intervals <- c(123, 91, 154, 118, 85)
cintervals <- cumsum(intervals)

##Time 

#how should it be included?

#Sessions took place Feb, June, Sep on two consecutive years

#No time included gives phi = 0.85

#time(months) # gives phi>0.9
time1 <- c(2, 6, 9, 14, 18, 21)

#or should it be this gives phi = 0.87
time2 <- rep(c(2, 6, 9, 14, 18, 21), each = 5)

#time(days) # phi = 0.99
time3 <- rep(1:5, 6) + rep(c(0,cintervals), each =5)

#time (months) # phi = 0.804
time4 <- time3 / 30

#time(weeks) #phi = 0.95
time5 <- time3 / 7

#time(years) #phi = 0.07
time6 <- time3 / 365

#Create a season covariate
season <- as.factor(rep(rep(1:3, each = 5), 2))

#Create year covariate
year <- as.factor(rep(1:2, each = 15))

##Code for mask
#Copied from OpenCR package examples
datadir <- system.file("extdata", package = "secr") 
OVforest <- rgdal::readOGR(dsn = datadir, layer = "OVforest")
ovmask <- make.mask(detectors, buffer = 120, type = 'trapbuffer', 
                    poly = OVforest[1:2,], spacing = 15, keep.poly = FALSE) 
ovmask2 <- make.mask(detectors, buffer = 200, type = 'trapbuffer', 
                     poly = OVforest[1:2,], spacing = 15, keep.poly = FALSE)

##Create ScrData object
scrdat <- ScrData$new(data, ovmask, primary = primary, time = time2)

#Add covs
scrdat$add_covariate("season", season, "k")
scrdat$add_covariate("year", year, "k")


#Find Starting values
start  <- get_start_values(scrdat, model = "JsModel")

##Constant across primary occasions
#Formula for the model
form <- list(lambda0 ~ 1, 
                  sigma ~ 1, 
                  phi ~ 1, 
                  beta ~ 1)

#Create Model object
jsmod <- JsModel$new(form, scrdat, start = start)

#Fit Model
jsmod$fit()

#Varying across primary occasions
#Formula for the model
form_prim <- list(lambda0 ~ primary, 
            sigma ~ primary, 
            phi ~ primary, 
            beta ~ primary)

#Create Model object
jsmod_prim <- JsModel$new(form_prim, scrdat, start)

#Fit Model
jsmod_prim$fit()

##Varying across seasons
#Formula for the model
form_seas <- list(lambda0 ~ season, 
                  sigma ~ season, 
                  phi ~ season, 
                  beta ~ season)

#Create Model object
jsmod_seas <- JsModel$new(form_seas, scrdat, start)

#Fit Model
jsmod_seas$fit()

##Varying across years
#Formula for the model
form_year <- list(lambda0 ~ year, 
                  sigma ~ year, 
                  phi ~ year, 
                  beta ~ year)

#Create Model object
jsmod_year <- JsModel$new(form_year, scrdat, start)

#Fit Model
jsmod_year$fit()

##Varying across years and seasons
#Formula for the model
form_seasyear <- list(lambda0 ~ season + year, 
                  sigma ~ season + year, 
                  phi ~ season + year, 
                  beta ~ season + year)

#Create Model object
jsmod_seasyear <- JsModel$new(form_seasyear, scrdat, start)

#Fit Model
jsmod_seasyear$fit()

##Transient Model

mean(scrdat$encrange(each = TRUE))
scrdat$encrange()

#There is evidence of transience

trans_form <- list(beta ~ 1, lambda0 ~ 1, sigma ~ 1, phi ~ 1, sd ~ 1)
trans_start <- get_start_values(scrdat, "JsTransientModel")
trans_obj <- JsTransientModel$new(form, scrdat, trans_start, print = TRUE)
trans_obj$fit()


##Compare with openCR

library(openCR)

spltcap <- split(scrdat$capthist(), as.factor(primary), byoccasion = TRUE)
fit.sp <- openCR.fit(spltcap, type = "JSSAsecrb", mask = scrdat$mesh(), trace = TRUE)

baseargs <- list(capthist = 'OVpossumCH', mask = ovmask, type = 'JSSAsecrb', 
                 sessioncov = c(1,2,3,1,2,3)) 
args <- rep(list(baseargs), 3) 
names(args) <- c('null','seasonal.phi.sigma','seasonal.lambda0.phi.sigma') 
args[[1]]$model <- list(lambda0~1, phi~1, sigma~1) 
args[[2]]$model <- list(lambda0~1, phi~scov, sigma~scov) 
args[[3]]$model <- list(lambda0~scov, phi~scov, sigma~scov) # byseason 
fits <- par.openCR.fit(args, ncores = 3, trace = TRUE)

fits <- readRDS("fits")
saveRDS(fits, "fits")

