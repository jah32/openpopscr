library(openpopscr)

# set true parameters 
true_par <- list(lambda0 = 2, sigma = 20, phi = 0.8, beta = 0.2, D = 1000)

# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")

# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")

# set number of occasions to simulate
n_occasions <- 10

scrdat <- simulate_js_openscr(true_par, n_occasions, detectors, mesh, seed = 52381)

# set each parameter to be a constant 
form <- list(lambda0 ~ 1, 
             sigma  ~ 1,
             phi ~ 1, 
             beta ~ 1)

# get some starting values based on data 
start <- get_start_values(scrdat, model = "JsModel")


mod <- JsModel$new(form, scrdat, start)

mod$fit()

# simulate a wet/dry covariate 
set.seed(385638)
season <- as.factor(sample(c("wet", "dry"), size = scrdat$n_occasions(), replace = TRUE))

scrdat$add_covariate("season", season, "k")
# copy formulae from previous model
form_seas <- form
# update formula for beta to include season
form_seas[[4]]<- beta ~ season
# create new model object 
mod_seas <- JsModel$new(form_seas, scrdat, start)
# fit model
mod_seas$fit()

# set truth 
true_par <- list(D = 1000, lambda0 = 2, sigma = 10, phi = 0.8, beta = 0.2, sd = 20)
# make detectors array 
detectors <- make.grid(nx = 7, ny = 7, spacing = 20, detector = "count")
# make mesh 
mesh <- make.mask(detectors, buffer = 100, nx = 64, ny = 64, type = "trapbuffer")
# set number of occasions to simulate
n_occasions <- 10
# simulate ScrData 
scrtransdat <- simulate_js_openscr(true_par, n_occasions, detectors, mesh, move = TRUE, seed = 34583)

# set formulae and start values 
stat_form <- list(lambda0 ~ 1, sigma ~ 1, phi ~ 1, beta ~ 1)
start <- get_start_values(scrtransdat, model = "JsModel")
# create model object 
stat <- JsModel$new(stat_form, scrtransdat, start)
# fit model 
stat$fit()

# specify formulas for parameters 
form <- list(lambda0 ~ 1, 
sigma  ~ 1,
phi ~ 1,
beta ~ 1,
sd ~ 1)

start <- get_start_values(scrtransdat, model = "JsTransientModel")

trans <- JsTransientModel$new(form, scrtransdat, start)

# fit model 
trans$fit()
