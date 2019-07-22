library(openpopscr)
library(secr)

true_par <- list(lambda0 = 2, sigma = 20, phi = 0.7)

# set number of individuals to track
N <- 100

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 10 

scrdat <- simulate_cjs_openscr(true_par, N, n_occasions, detectors, mesh, seed = 19539)


# set each parameter to be a constant 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1,
             phi ~ 1)

# get some starting values based on data 
start <- get_start_values(scrdat, model = "CjsModel")


mod <- CjsModel$new(form, scrdat, start)

mod$fit()

# simulate a disturbance covariate (% deforestation for example)
set.seed(385638)
disturb <- runif(scrdat$n_occasions() - 1) 
# add covariate 
scrdat$add_covariate("disturb", disturb, "k")
# copy formulae from previous model
form_surv <- form
# update formula for phi to include disturbance 
form_surv[[3]]<- phi ~ disturb 

# create new model object 
mod_surv <- CjsModel$new(form_surv, scrdat, start)
# fit model
mod_surv$fit()


# set truth 
true_par <- list(lambda0 = 2, sigma = 20, phi = 0.7, sd = 20)
# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")
# set number of occasions to simulate
n_occasions <- 5
# set number of individuals tracked
N <- 100
# simulate ScrData 
scrtransdat <- simulate_cjs_openscr(true_par, 
N, 
n_occasions, 
detectors, 
mesh, 
move = TRUE, 
seed = 95811)

# set formulae and start values 
stat_form <- list(lambda0 ~ 1, sigma ~ 1, phi ~ 1)
start <- get_start_values(scrtransdat, model = "CjsModel")

# create model object 
stat <- CjsModel$new(stat_form, scrtransdat, start)
# fit model 
stat$fit()

form <- list(lambda0 ~ 1, 
sigma  ~ 1,
phi ~ 1,
sd ~ 1)

start <- get_start_values(scrtransdat, model = "CjsTransientModel")

trans <- CjsTransientModel$new(form, scrtransdat, start)

# fit model 
trans$fit()
