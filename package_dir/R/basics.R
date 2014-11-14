`%out%` <- function (x, y) !(x %in% y)
pow <- function(x,y) x^y

day_to_circle <- function(yday) yday/366*2*pi


sq_distance <- function(x, y) {
  return(pow(x-y,2));
}
  
sq_circ_distance <- function(x, y) {
    d <- ifelse(x > y,
			{a <- x-y; b <- (2.0*pi-x)+y; mapply(min,a,b)},
			{a <- y-x; b <- (2.0*pi-y)+x; mapply(min,a,b)}
    )
    return(pow(d,2))
}
  
make_knots <- function(n_knots, lower_knot, upper_knot) {
  knots <- vector(mode='numeric', length=n_knots)
  interknot_length <- (upper_knot - lower_knot) / (n_knots-1)                                      

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
	
	return(knots)
}

make_circle_knots <- function(n_knots) {
	knots <- vector(mode='numeric', length=n_knots)
  if (n_knots > 1) {
    for ( i in 1:(n_knots)) {
      knots[i] <- 2.0*pi/n_knots*(i-1);
    }
  } else {
    knots[1] <- 0;
  }
  if (n_knots < 1) {
    print("'n_knots' < 1 is non-sensical. Assumed 1 knot.");
  }
	return(knots)
}

