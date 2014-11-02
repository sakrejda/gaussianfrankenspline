
data {
	int<lower=2> n_points;
	vector[n_points] positions;
	vector[n_points] day_of_year;

	int n_knots;   // Must be big enough to let you estimates GP parameters.
}

transformed data {
	vector[n_knots] knot_points;
	real<lower=0.0> knot_scale;

	knot_points <- make_circle_knots(n_knots);
	knot_scale <- pi()/n_knots;
}

parameters {
	real<lower=0> theta_eta_sq;
	real<lower=0> theta_rho_sq;
	real<lower=0> theta_sigma_sq;
	real theta_mu;

	vector[n_knots] knot_weights; 
}

model {
	matrix[n_knots, n_knots] theta_Sigma;
	vector[n_knots] mu;
	vector[n_points] positions_mu;
	
	// Priors on hyperparameters of GP:
	theta_mu ~ normal(0,10); 
	theta_eta_sq ~ cauchy(0,5);
	theta_rho_sq ~ cauchy(0,5);
	theta_sigma_sq ~ cauchy(0,5);
	theta_Sigma <- gp_circ_generalized_sq_exp(theta_eta_sq, theta_rho_sq, theta_sigma_sq, knot_points);

	knot_weights ~ multi_normal(rep_vector(theta_mu,n_knots), theta_Sigma);

	positions ~ normal(yday_circular_spline(day_of_year, knot_points, knot_weights, knot_scale), 1);
}

