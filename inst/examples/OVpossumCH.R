#OVpossum example.
#Adapted from an example given in 
#https://www.otago.ac.nz/density/pdfs/openCR-examples.pdf

library(openpopscr)
library(secr)
library(rgdal)
library(RcppParallel)
setThreadOptions(defaultNumThreads() - 1)

data <- OVpossumCH

#extract the detectors
detectors <- traps(data[[1]])

#Extract the capture history for each session
dat1 <- data$`49`
dat2 <- data$`50`
dat3 <- data$`51`
dat4 <- data$`52`
dat5 <- data$`53`
dat6 <- data$`54`

#Code for mask
#Copied from OpenCR package examples
datadir <- system.file("extdata", package = "secr") 
OVforest <- rgdal::readOGR(dsn = datadir, layer = "OVforest")
ovmask <- make.mask(detectors, buffer = 120, type = 'trapbuffer', poly = OVforest[1:2,], spacing = 15, keep.poly = FALSE) 
ovmask2 <- make.mask(detectors, buffer = 200, type = 'trapbuffer', poly = OVforest[1:2,], spacing = 15, keep.poly = FALSE)

#Create a ScrData object for each session
scrdat1 <- ScrData$new(dat1, ovmask)
scrdat2 <- ScrData$new(dat2, ovmask)
scrdat3 <- ScrData$new(dat3, ovmask)
scrdat4 <- ScrData$new(dat4, ovmask)
scrdat5 <- ScrData$new(dat5, ovmask)
scrdat6 <- ScrData$new(dat6, ovmask)

#Find Starting values
start  <- get_start_values(scrdat1, model = "JsModel")

#Formula for the model
form <- list(lambda0 ~ 1, 
            sigma ~ 1, 
            phi ~ 1, 
            beta ~ 1)

#Create model object for first session
scrmod1 <- JsModel$new(form, scrdat1, start)

start <- list(lambda0=2, sigma = 20, phi = 0.7, beta=0.2, D= 2000)

#Model objects for other sessions.
scrmod2 <- JsModel$new(form, scrdat2, start)
scrmod3 <- JsModel$new(form, scrdat3, start)
scrmod4 <- JsModel$new(form, scrdat4, start)
scrmod5 <- JsModel$new(form, scrdat5, get_start_values(scrdat5, "JsModel"))
scrmod6 <- JsModel$new(form, scrdat6, get_start_values(scrdat6, "JsModel"))

#Fit models where no parameters are shared across sessions.
scrmods <- list(scrmod1, scrmod2, scrmod3, scrmod4, scrmod5, scrmod6)
for(i in scrmods){i$fit()}
sum(sapply(1:6, FUN = function(i)AIC(scrmods[[i]])))

#A model where all parameters are shared across sessions
scrdat <- list(s1 = scrdat1, s2 = scrdat2, s3 = scrdat3, s4 = scrdat4,
               s5 = scrdat5, s6 = scrdat6)
shared_form <- list(beta ~ 1, lambda0 ~ 1, phi ~ 1, sigma ~ 1)
private_form <- list(s1=list(), s2=list(), s3=list(), s4=list(), s5=list(), s6 = list())
obj <- StrataModel$new(scrdat, "JsModel", shared_form, private_form, start)
obj$fit()
AIC(obj)


#Transient Model

form <- list(beta ~ 1, lambda0 ~ 1, sigma ~ 1, phi ~ 1, sd ~ 1)
trans_start1 <- get_start_values(scrdat1, "JsTransientModel")
trans_obj <- JsTransientModel$new(form, scrdat1, trans_start1, print = TRUE)
trans_obj$fit()

trans_obj <- readRDS("ovpossumtrans.rds")
trans_obj$fit()
saveRDS(trans_obj, "OvpossumTrans.rds")
