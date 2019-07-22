library(openpopscr)
library(secr)


# set true parameters 
true_par <- list(D = 1000, lambda0 = 2, sigma = 20)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

scrdat <- simulate_scr(true_par, n_occasions, detectors, mesh, seed = 15483)

# set each parameter to be a constant 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1)

# get some starting values based on data 
start <- get_start_values(scrdat)


mod <- ScrModel$new(form, scrdat, start)

mod$fit()


# set seed so random covariate values will always be the same 
set.seed(58285)
# Simulate, for each detector, an age class: old or new. In reality, this
# covariate would be observed. 
age <- as.factor(sample(c("old", "new"), size = scrdat$n_traps(), replace = T))

# Add the covariate to the ScrData 
scrdat$add_covariate("age", age, "j")

# copy formulae from previous model
form_detage <- form
# update formula for lambda0 to include age 
form_detage[[1]]<- lambda0 ~ age

# create new model object 
mod_detage <- ScrModel$new(form_detage, scrdat, start)
# fit model
mod_detage$fit()

# set truth 
true_par <- list(D = 1000, lambda0 = 2, sigma = 10, sd = 20)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 5 

# simulate ScrData 
scrtransdat <- simulate_scr(true_par, n_occasions, detectors, mesh, move = TRUE, seed = 13854)

# set formulae and start values 
stat_form <- list(lambda0 ~ 1, sigma ~ 1)
start <- get_start_values(scrdat)

# create model object 
stat <- ScrModel$new(stat_form, scrtransdat, start)
# fit model 
stat$fit()

# specify formulas for parameters 
form <- list(lambda0 ~ 1, 
sigma  ~ 1, 
sd ~ 1)

start <- get_start_values(scrtransdat, model = "ScrTransientModel")

trans <- ScrTransientModel$new(form, scrtransdat, start)

# fit model 
trans$fit()
