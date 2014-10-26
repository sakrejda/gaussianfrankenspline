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
  


