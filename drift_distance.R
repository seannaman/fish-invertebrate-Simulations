## Function for the drift settlement model from the supplement.The Model predicts the number 
## of animals remaining in suspension at X length downstrem of a given release point. 
## This is assumed to follow the negative exponential distribution:
## Nx = N0*e^RX

# X is the length of a given upstream riffle, 
# N0 is the initial number in the drift, 
# R is the settlment rate - defined as 1/drift distance (D).
# D is a function of water velocity (Anderson et al. 2013. Ecological Modelling)

## This function also incorporates uncertainty and variance in drift distances by randomly drawing from an 
## exponential distribution around a given D value.

Nx_estimate <- function(V, X, N0){ 
  D <- vector()
  Nx <- vector()
  mean_Nx <- vector()
  sd_Nx <- vector()
  for(i in 1:10000){
    D <- rexp(N0, (1/(0.1346*V + 0.7442))) # Drawing a random set of drift distances from exp. distribution
    Nx[i] <- length(D[D > X]) ## The length of cases within drift distances > specified length
  }
  mean_Nx <- mean(Nx) # Mean value from 
  se_Nx <- sd(Nx)/sqrt(length(Nx))
  return(cbind(mean_Nx,se_Nx)) ##End
} 