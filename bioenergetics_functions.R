### Function to estimate consumption based on observed growth. This code is essentially obselete
### with the release of the R-based Bioenergetics 4.0 by Deslauries et al. (2017) Fisheries.
## Nevertheless, I'll provide the code I wrote for reproducability purposes. The approach is to model growth 
## with different values of pCmax, then use pCmax values that most closely match modeled growth to observed
## growth [based on final mass]. The input data required are final and initial fish mass (as columns) and 
## daily temperatures [a dataframe called temperature].The function returns pCmax, cumulative consumption, and predicted final mass for each fish.


## Set bioenergetics parameters - Here for coho salmon
## Species specific parameters can be found in the original bioenergetics documentation or the Bioenergetics 4.0 manual
## f(T) = Ka*Kb
CQ <- 5
CK1 <- 0.36
CTO <- 15
CTM <- 18
CTL <- 24
CK4 <- 0.01
## Cmax function parameters
CA <- 0.303
CB <- -0.275
# W <- fish mass (from observed growth dataframe)

## Respiration
RA <- 0.0046
RB <- -0.217
RQ <- 2.1
RTO <- 18
RTM <- 26

## Egestion and excretion  
FA <- 0.212
FB <- -0.522
FG <- 0.631
UA <- 0.0214
UB <- 0.380
UG <- -0.299

## SDA
SDA <- 0.172

## Energy density
alpha1 <- 4111
beta1 <- 155
alpha2 <- 7602
beta2 <- 0.5266
cutoff <- 10 ## switches from alpha1 and beta1 to alpha2 beta2 at a mass of 10g


# Set up bioenergetics function ####
## Inputs are initial and final mass for each fish. Function calculates growth and consumption at a daily time step
## Temperature is mean daily temp from ibutton data loggers. In data frame "temp"
bio_sim <- function(Mass, final.mass){
  p_Cmax <- seq(0,1,0.05) # Sequence of pCmax values to run through during simulations
  consumption <- vector() # daily consumption rate
  consumption.g <- vector() # total consumption in grams
  Cmax <- vector() # Maximum ration
  R <- vector()
  egestion <- vector()
  cum_consumption <- vector()
  growth_fit <- vector() # Growth over each interval
  end.weight <- vector() ## fitted final mass
  W <- vector()
  f_T <- vector()
  R <- vector()
  U <- vector()
  S <- vector()
  temp <- vector()
  f_TR <- vector()
  final.mass <- final.mass
  Mass <- Mass
  for(i in 1:length(p_Cmax)){ 
    for(j in 1:length(temperature$Day)){
      temp[j] <- temperature$Temp[temperature$Day==j] 
      W[j] <- ifelse(j==1, Mass, W[j-1]+growth_fit[j-1]) 
      ### Set up bioenergetics equations 
      ## Consumption at given p_Cmax
      # Temperature dependence function
      G1 <-(1/(CTO-CQ)) * log((0.98*(1-CK1))/(CK1 * 0.02)) 	#1st
      L1 <- exp(G1*(temp[j] - CQ))	#2nd 
      Ka <- (CK1 *L1) / (1 + CK1 * (L1-1))	#3rd 
      G2 <-(1/(CTL-CTM)) * log((0.98 * (1-CK4)) / (CK4 * 0.02))	#4th 
      L2 <-exp(G2*(CTL - temp[j])) 
      Kb <-(CK4*L2) / (1 + CK4* (L2-1)) 
      f_T[j] <- Ka*Kb  
      # Consumption equation
      Cmax[j] <- CA * W[j]^(CB)       
      consumption[j] <- Cmax[j]*f_T[j]*p_Cmax[i] 
      consumption.g[j] <- consumption[j]*W[j] 
      ## Respiration
      # Temperature dependence function
      Y = log(RQ)*(RTM-RTO+2) 
      Z = log(RQ)*(RTM-RTO) 
      X = ((Z^2 *(1+40/Y)^0.5)^2)/400
      V = (RTM-temp[j])/(RTM-RTO)
      f_TR[j] = V^X*exp(X*(1-V))
      # Respiration equation with ACT=2
      R[j] = RA*W[j]^RB*f_TR[j]*2
      ## Egestion
      egestion[j] <- FA*temp[j]^FB*exp(FG*p_Cmax[i])*consumption[j]
      ## Excretion
      U[j] <- UA * temp[j]^(UB) * exp(UG * p_Cmax[i]) * (consumption[j]-egestion[j])
      ## Specific dynamic action as proportion of total consumption
      S[j] <- SDA*(consumption[j]-egestion[j])
      ### Then estimate daily growth for given ration (p_Cmax)
      growth_fit[j] <- consumption[j]- (R[j]+S[j]+egestion[j]+U[j])
    }
    ### final mass and total consumption after simulation period of 28 days
    end.weight[i] <- W[[28]]
    cum_consumption[i] <- sum(consumption.g)
    
  }
  ### bind desired output into dataframe 
  output <- data.frame(cbind(p_Cmax,end.weight,cum_consumption))
  ### Find which fitted final mass is closest to the observed final mass then retain p_Cmax and consumption
  p_best <- output[which(abs(output$end.weight-final.mass)==min(abs(output$end.weight-final.mass))),]
  return(list(p_best)) 
} 


##########################
##########################
## Bioenergetics function to estimate potential fish growth from consuming additional prey ####
## Consumption here is known from fitted values and is an input. Growth (g per day) is an output
## Additional consumption can be simulated by iterating this function over user defined values 
## of consumption.gday

consumption_sim <- function(consumption.gday, Mass){
  p_Cmax <- vector()
  consumption <- vector() # daily consumption rate
  consumption.gday <- consumption.gday # total consumption in grams
  Cmax <- vector() # Maximum ration
  R <- vector()
  egestion <- vector()
  cum_consumption <- vector()
  growth_fit <- vector() # Growth over each interval
  end.weight <- vector() ## fitted final mass
  W <- vector()
  f_T <- vector()
  R <- vector()
  U <- vector()
  S <- vector()
  temp <- vector()
  f_TR <- vector()
  Mass <- Mass
  for(j in 1:length(temperature$Day)){
    temp[j] <- temperature$Temp[temperature$Day==j] 
    W[j] <- ifelse(j==1, Mass, W[j-1]+growth_fit[j-1]) 
    ### Set up bioenergetics equations 
    ## Convert fitted consumption to g per g of fish
    consumption[j] <- consumption.gday/W[j]
    # Temperature dependence function
    G1 <-(1/(CTO-CQ)) * log((0.98*(1-CK1))/(CK1 * 0.02)) 	#1st
    L1 <- exp(G1*(temp[j] - CQ))	#2nd 
    Ka <- (CK1 *L1) / (1 + CK1 * (L1-1))	#3rd 
    G2 <-(1/(CTL-CTM)) * log((0.98 * (1-CK4)) / (CK4 * 0.02))	#4th 
    L2 <-exp(G2*(CTL - temp[j])) 
    Kb <-(CK4*L2) / (1 + CK4* (L2-1)) 
    f_T[j] <- Ka*Kb  
    # Consumption equation. This time solve for p_Cmax at given consumption value 
    Cmax[j] <- CA * W[j]^(CB)       
    p_Cmax[j] <- Cmax[j]*f_T[j]*consumption[j] 
    ifelse(p_Cmax[j] >1,1,p_Cmax[j]) # limit pCmax to 1, i.e., fish can't consume beyond max daily ration
    ## Respiration
    # Temperature dependence function
    Y = log(RQ)*(RTM-RTO+2) 
    Z = log(RQ)*(RTM-RTO) 
    X = ((Z^2 *(1+40/Y)^0.5)^2)/400
    V = (RTM-temp[j])/(RTM-RTO)
    f_TR[j] = V^X*exp(X*(1-V))
    # Respiration equation with ACT=2
    R[j] = RA*W[j]^RB*f_TR[j]*2
    ## Egestion
    egestion[j] <- FA*temp[j]^FB*exp(FG*p_Cmax[j])*consumption[j]
    ## Excretion
    U[j] <- UA * temp[j]^(UB) * exp(UG * p_Cmax[j]) * (consumption[j]-egestion[j])
    ## Specific dynamic action as proportion of total consumption
    S[j] <- SDA*(consumption[j]-egestion[j])
    ### Then estimate daily growth for given ration (p_Cmax)
    growth_fit[j] <- consumption[j]- (R[j]+S[j]+egestion[j]+U[j])
  }
  ### final mass and total consumption after simulation period of 28 days
  end.weight <- W[[28]]
  mass.change <- end.weight - Mass
  ### bind desired output into dataframe 
  output <- data.frame(cbind(end.weight,mass.change))
  return(list(output)) 
} 




