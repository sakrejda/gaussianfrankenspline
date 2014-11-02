library(ggplot2)
library(lubridate)
library(MASS)
library(rstan)

source("package_dir/GFS/R/basics.R")
source("package_dir/GFS/R/spline_functions.R")
source("package_dir/GFS/R/smooth_functions.R")

n_points <- 2000
theta_mu <- 13.55
mu <- rep(theta_mu, n_points)
n_steps <- n_points-1
time <- cumsum(rexp(n=n_points, rate=0.2))
time <- ymd('2010-01-01') + days(round(time))
day_of_year <- yday(time)
angle <- day_to_circle(day_of_year)

n_knots <- 30 
knot_weights_mu <- rep(3,n_knots)
last_knot_day <- 366 - 366/n_knots
day_of_knot <- 366/n_knots * (0:(n_knots-1)) #seq(from=0, to=last_knot_day, length.out=n_knots)
angle_of_knot <- day_to_circle(day_of_knot)

## We construct some arbitrary knot weights to plot the spline first:
knot_weights <- c(1:(n_knots/2-1),(n_knots/2+1):1)

spline_position <- circular_spline(
	x=angle, knot_points=angle_of_knot,
	knot_weights=knot_weights, knot_scale=pi/n_knots)

## Shows where the spline goes based on some knots:
qplot(day_of_year, spline_position, geom='line', ylim=c(0,1.1*max(spline_position)))

theta_eta_sq <- 10
theta_rho_sq <- 1
theta_sigma_sq <- 1

knot_weights <- GP(1, angle_of_knot, knot_weights_mu, theta_eta_sq, theta_rho_sq, theta_sigma_sq, TRUE)
qplot(angle_of_knot, knot_weights)



spline_position <- circular_spline(
	x=angle, knot_points=angle_of_knot,
	knot_weights=knot_weights, knot_scale=pi/n_knots, normalize=FALSE)

## Shows where the spline goes based on some knots:
qplot(time, spline_position, geom='line', ylim=c(0,1.1*max(spline_position)))
qplot(day_of_year, spline_position, geom='line') + 
	geom_point(data=data.frame(x=day_of_knot, y=knot_weights), aes(x=x, y=y)) 



data <- list(
	n_points = n_points,
	positions = spline_position + rnorm(length(spline_position),0,.05),
	day_of_year = day_of_year,
	n_knots = n_knots 
)

with(data=data, expr=stan_rdump(list=ls(), file='circular-rw.data'))


initial_fun <- function(id) {
	list(
		theta_eta_sq = theta_eta_sq * runif(1,0.9,1.1),
		theta_rho_sq = theta_rho_sq * runif(1,0.9,1.1),
		theta_sigma_sq = theta_sigma_sq * runif(1,0.9,1.1),
		theta_mu = 13,
		knot_weights = knot_weights,
		positions_sigma = 0.01
	)
}

inits <- initial_fun(1)

with(data=inits, expr=stan_rdump(list=ls(), file='circular-rw.inits'))

truth <- list(
	theta_mu = theta_mu,
	theta_eta_sq = theta_eta_sq,
	theta_rho_sq = theta_rho_sq,
	theta_sigma_sq = theta_sigma_sq,
	knot_weights = knot_weights
)

with(data=truth, expr=stan_rdump(list=ls(), file='circular-rw.truth'))



