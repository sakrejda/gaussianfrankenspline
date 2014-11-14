
gp_generalized_sq_exp <- function(eta_sq,  rho_sq,  sigma_sq,  x) {
    Sigma <- matrix(data=0, nrow=length(x), ncol=length(x)) 
    for (i in 1:length(x)) {
      for (j in i:length(x)) {
        Sigma[i,j] <- eta_sq * exp(-rho_sq * sq_distance(x[i],x[j]));
        Sigma[j,i] <- Sigma[i,j];
      }
    }
    
    for (k in 1:length(x)) {
      Sigma[k,k] <- eta_sq + sigma_sq;
    }
    
    return(Sigma)
  }
  
gp_circ_generalized_sq_exp <- function(eta_sq, rho_sq, sigma_sq, x) {
  Sigma <- matrix(data=0, nrow=length(x), ncol=length(x)) 
  for (i in 1:length(x)) {
    for (j in i:length(x)) {
      Sigma[i,j] <- eta_sq * exp(-rho_sq * sq_circ_distance(x[i],x[j]));
    	Sigma[j,i] <- Sigma[i,j];
  	}
  }
    
  for (k in 1:length(x)) {
  	Sigma[k,k] <- eta_sq + sigma_sq;
  }
    
	return(Sigma);
}

GP <- function(n=1, x, mu=rep(0,length(x)), eta, rho, sigma, circular=FALSE) {
  if (circular) {
    Sigma <- gp_circ_generalized_sq_exp(eta, rho, sigma, x)
  } else {
    Sigma <- gp_generalized_sq_exp(eta, rho, sigma, x)
  }
  y <- mvrnorm(n=n, mu=mu, Sigma=Sigma)
  return(y)
}




#angle <- (yday(time)/366-0.5)*2*pi
#
#theta_eta_sq <- 10
#theta_rho_sq <- 1
#theta_sigma_sq <- 1
#theta_Sigma <- gp_circ_generalized_sq_exp(theta_eta_sq, theta_rho_sq, theta_sigma_sq, angle)
#positions <- mvrnorm(n=1, mu=mu, Sigma=theta_Sigma)
#qplot(angle, positions)

