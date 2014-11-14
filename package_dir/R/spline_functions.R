
## @knitr generic-radial-spline
# Generic radial spline:
radial_spline <- function(x, f, knot_points, knot_weights, knot_scale, normalize=FALSE) {
	if (normalize) {
		normalize_height <- f(x=knot_points[1], center=knot_points[1], scale=knot_scale[1])
	} else {
		normalize_height <- 1
	}
  if (length(x) == 1) {
    knot_contributions <- f(x=x, center=knot_points, scale=knot_scale)*knot_weights/normalize_height
    return(sum(knot_contributions))
  } else {
    return(sapply(
			X=x,
			FUN=radial_spline, 
			f=f, knot_points=knot_points, knot_weights=knot_weights, knot_scale=knot_scale,normalize=normalize))
  }
}

## @knitr gaussian-radial-spline
# Generic radial spline specialized to Gaussian.
gaussian_spline <- function(x, knot_points, knot_weights, knot_scale) {
  f <- function(x, center, scale) dnorm(x=x, mean=center, sd=scale)
  return(radial_spline(x=x, f=f, knot_points=knot_points, knot_weights=knot_weights, knot_scale=knot_scale))
}


## @knitr gaussian-spline-integral
# Integral of Gaussian radial spline.
gaussian_spline_integral <- function(knot_points, knot_weights, knot_scale, a, b) {
  contributions <- mapply(
    FUN=function(point, weight, scale, a, b) {
      weight * (pnorm(q=b, mean=point, sd=scale) - pnorm(q=a, mean=point, sd=scale))
    },
    point = knot_points, weight = knot_weights,
    MoreArgs = list(scale=knot_scale, a=a, b=b)
  )
  return(sum(contributions))
}

## @knitr circular-normal function for spline, only k=-1,0,1
# Circular normal.
circular_normal_density <- function(x, center, scale) {
  o <- 1/(scale*sqrt(2*pi)) * (
 		exp(-(sq_circ_distance(x,center)^.5-2*pi)^2/(2*scale^2)) +
    exp(-(sq_circ_distance(x,center)^.5     )^2/(2*scale^2)) +
    exp(-(sq_circ_distance(x,center)^.5+2*pi)^2/(2*scale^2))
  )
}

circular_spline <- function(x, knot_points, knot_weights, knot_scale, normalize=FALSE) {
	f <- circular_normal_density
  o <- radial_spline(x=x, f=f, knot_points=knot_points, knot_weights=knot_weights, knot_scale=knot_scale, normalize)
	return(o)
}
	







