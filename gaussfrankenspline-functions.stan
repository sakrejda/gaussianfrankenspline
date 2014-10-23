functions {
	real sq_distance(real x, real y) {
		return pow(x-y,2);
	}

	real sq_circ_distance(real x, real y) {
		real d;
		vector[2] xy;
		if (x > y) {
			xy[1] <- x-y;
			xy[2] <- (2.0*pi()-x)+y;
			d <- min(xy);
		} else {
			xy[1] <- y-x;
			xy[2] <- (2.0*pi()-y)+x;
			d <- min(xy);
		}
		return pow(d,2);
	}

	matrix gp_generalized_sq_exp(
		real eta_sq, real rho_sq, real sigma_sq, vector x
	) {
		matrix[dims(x)[1], dims(x)[1]] Sigma;
		// off-diagonal for Sigma
		for (i in 1:(dims(x)[1])) {
			for (j in i:(dims(x)[1])) {
				Sigma[i,j] <- eta_sq * exp(-rho_sq * sq_distance(x[i],x[j]));
				Sigma[j,i] <- Sigma[i,j];
			}
		}

		// diagonal elements
		for (k in 1:(dims(x)[1])) {
			Sigma[k,k] <- eta_sq + sigma_sq;
		}

		return Sigma;
	}

	matrix gp_circ_generalized_sq_exp(
		real eta_sq, real rho_sq, real sigma_sq, vector x
	) {
		matrix[dims(x)[1], dims(x)[1]] Sigma;
		// off-diagonal for Sigma
		for (i in 1:(dims(x)[1])) {
			for (j in i:(dims(x)[1])) {
				Sigma[i,j] <- eta_sq * exp(-rho_sq * sq_circ_distance(x[i],x[j]));
				Sigma[j,i] <- Sigma[i,j];
			}
		}

		// diagonal elements
		for (k in 1:(dims(x)[1])) {
			Sigma[k,k] <- eta_sq + sigma_sq;
		}

		return Sigma;
	}

	
	vector make_knots(int n_knots, real lower_knot, real upper_knot) {
		vector[n_knots] knots;
		real interknot_length;
		interknot_length <- (upper_knot - lower_knot) / (n_knots-1);

		if (n_knots > 1) {
			for ( i in 1:(n_knots)) {
				knots[i] <- interknot_length*(i-1);
			}
		} else {
			knots[1] <- (upper_knot - lower_knot) / 2.0;
		}
		if (n_knots < 1) {
			print("'n_knots' < 1 is non-sensical. Assumed 1 knot.");
		}
		return knots;
	}

	vector make_circle_knots(int n_knots) {
		vector[n_knots] knots;
		if (n_knots > 1) {
			for ( i in 1:(n_knots)) {
				knots[i] <- 2.0*pi()/n_knots*(i-1);
			}
		} else {
			knots[1] <- 0;
		}
		if (n_knots < 1) {
			print("'n_knots' < 1 is non-sensical. Assumed 1 knot.");
		}
		return knots;
	}

	real spline(
		real x,
		vector knot_points, vector knot_weights,
		real knot_scale
	) {
		real y;
		int J;
		y <- 0.0;
		J <- dims(knot_points)[1];
		for ( j in 1:J ) {
			y <- y + knot_weights[j] * normal_log(x,knot_points[j],knot_scale);
		}
		return y;
	}

	real circular_spline(
		real theta,
		vector knot_points, vector knot_weights, 
		real knot_scale
	) {
		real y;
		int J;
		y <- 0.0;
		J <- dims(knot_points)[1];
		for ( j in 1:J ) {
			y <- y + knot_weights[j] * 1/(knot_scale*pow(2.0*pi(),.5)) * (
            exp(-(pow(sq_circ_distance(theta,knot_points[j]),.5)+2.0*pi())^2/(2.0*pow(knot_scale,2))) +
            exp(-(pow(sq_circ_distance(theta,knot_points[j]),.5)         )^2/(2.0*pow(knot_scale,2))) +
            exp(-(pow(sq_circ_distance(theta,knot_points[j]),.5)+2.0*pi())^2/(2.0*pow(knot_scale,2)))
      );
		}
		return y;
	}

	real day_to_angle(real day) {
		real angle;
		angle <- (day/366)*2*pi();
		return angle;
	}

	real yday_circular_spline(
		real yday,
		vector knot_points, vector knot_weights, 
		real knot_scale
	) {
		real theta;
		real y;
		theta <- day_to_angle(yday);
		y <- circular_spline(theta, knot_points, knot_weights, knot_scale);
		return y;
	}
}

