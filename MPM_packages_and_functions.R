                 # -------------------------------------- #
                 #                FUNCTIONS               #
                 #    Roe deer matrix population model    #
                 # -------------------------------------- #

# ------------------------------------------------------------------------ #
#                     Functions used in the simulation                     #
#                     for roe deer population dynamics                     #
# ------------------------------------------------------------------------ #


# ----- Install packages needed -----

# - Distribution packages:
# Install gtools for the Dirichlet Distribution
if(!require(gtools)){
  install.packages("gtools")
  library(gtools)
}

# Install truncnorm for the Trunc Normal Distribution
if(!require(truncnorm)){
  install.packages("truncnorm")
  library(truncnorm)
}

# - Tidy plot and data 
# Install tidyverse for ggplot2, tidyr and dplyr
if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}
# Combine separate ggplots into the same graphic
if(!require(patchwork)){
  install.packages("patchwork")
  library(patchwork)
}

# Time your simulation
if(!require(tictoc)){
  install.packages("tictoc")
  library("tictoc")
}



# ----- *Functions* -----

# - Stochastic initial density ------------------------------------------------

Stoch.init <- function(density, name, area, variation = 0.1, stoch = T){
  # ---
  # This function make the initial structure and density to vary (stochastic)
  # ---
  # density is the initial density of the population
  # name is the name of each stage of the population
  # area is the study area in square kilometres
  # variation is the proportion variation around the initial density
  # stoch permit the turn ON/OFF the stochasticity
  # ---
  if (isTRUE(stoch)){
    # Number of individuals in every stage at t0 with 10% variation
    N.pop <- round(runif(1, density*area - (density*area*variation), density*area + (density*area*variation)))
    # Population abundance vector
    N.init <- VectorPop(round(rdirichlet(1, c(30, 20, 15, 30, 30, 15, 30))
                              * N.pop), name)
  } else { # Start at the initial state from fixed initial density in the parameter file
    source("MPM_load_parameters.R")
    N.init <- N.roe.init
  }
  return(N.init)
}



# - Population abundance vector -----------------------------------------------

VectorPop <- function(vec.abund, vec.name = NULL){
  # ---
  # This function create the abundance vector for the population projection model
  # ---
  # vec.abund is the abundance of each stages of the population
  # vec.name is the name of all the stages in the population
  # ---
  N.pop.init <- matrix(0, length(vec.abund), ncol=1)
  for (i in 1:(length(vec.abund))){
    N.pop.init[i] <- vec.abund[i]
  }
  colnames(N.pop.init) <- "Abundance"
  rownames(N.pop.init) <- vec.name
  return(N.pop.init)
}



# - Reproduction function -----------------------------------------------------

Repro <- function(repro.type, sex, k, h=1, a=25, Nb.female, Nb.male){
  if (repro.type == "HFF"){
    # ---
    # Harmonic Fertility Function
    # This function calculate the contribution of each sex to the reproduction
    # used in Caswell, 2001 - Gerber & White, 2014
    # ---
    # Sex allows to differentiate between male or female fertility function
    # k is the average litter size
    # h is the number of female per male for reproduction
    # Nb.female and Nb.male is the number of individuals in each sex
    # ---
    Ff <- (k*Nb.male)/(Nb.male + Nb.female/h)
    Fm <- (k*Nb.female)/(Nb.male + Nb.female/h)
    if (Nb.female == 0 || Nb.male == 0 || is.na(Nb.female) || is.na(Nb.male)){
      Ff <- 0
      Fm <- 0
    }
    if (sex == "F"){
      return(Ff)
    } else if (sex == "M") {
      return(Fm)
    } else {
      stop("Error in sex reproduction parameter")
    }
    
  } else if (repro.type == "FemaleFF"){
    # ---
    # Female Fertility Function
    # This function calculate female fertility
    # considering the sex ratio of male/female in the model
    # allowing for polygamous reproduction
    # ---
    # k is the average litter size
    # h is the number of female per male for reproduction, h=1 is monogamous
    # Nb.female and Nb.male is the number of individuals in each sex
    # ---
    Fertility <- k * min(1, Nb.male/(Nb.female/h))
    
    if (Nb.female == 0 || Nb.male == 0 || is.na(Nb.female) || is.na(Nb.male)){
      Fertility <- 0
    }
    if (sex == "F"){
      return(Fertility)
    } else if (sex == "M") {
      return(0)
    } else {
      stop("Error in sex reproduction parameter")
    }
    
  } else if (repro.type == "ExpFF"){
    # ---
    # Exponential Female Fertility Function
    # This function calculate female fertility
    # considering the sex ratio of female/male in the model
    # allowing for polygamous reproduction
    # ---
    # k is the average litter size
    # a is the strength parameter of the exponential 
    # Nb.female and Nb.male is the number of individuals in each sex
    # ---
    Fertility <- k * (1-exp(-a*(Nb.male/Nb.female)))
    
    if (Nb.female == 0 || Nb.male == 0 || is.na(Nb.female) || is.na(Nb.male)){
      Fertility <- 0
    }
    if (sex == "F"){
      return(Fertility)
    } else if (sex == "M") {
      return(0)
    } else {
      stop("Error in sex reproduction parameter")
    }
    
  } else {
   stop("Reproduction function not available") 
  }
}



# - Survival functions ---------------------------------------------------------

# Survival rate
Srate <- function(mortality.rates){
  # ---
  # Function calculating the survival rate from different sources of mortality
  # ---
  # mortality.rate is a vector of different mortality rates
  # ---
  Sx <- 1
  for (i in 1:length(mortality.rates)){
    if (mortality.rates[i] > 1){
      mortality.rates[i] <- 1
    } else if (mortality.rates[i] < 0){
      mortality.rates[i] <- 0
    }
    Sx <- Sx - mortality.rates[i]
  }
  if (Sx <= 0){
    Sx <- 0.01 # Not everyone die at once
  }
  return(Sx)
}



# Lynx predation mortality
LynxPred <- function(kill.rate = NULL, months.pred = NULL, lynx.hunter = NULL, roe.tot = NULL, roe.select = NULL, lynx.pred){
  # ---
  # Function using kill rate of lynx per 30 days
  # ---
  # kill.rate is the number of roe deer killed by lynx per month
  # months.pred is the number of months in the time step simulation
  # (assume that lynx prey on roe deer the same way all over the year)
  # lynx.hunter is the stage of lynx that are preying in roe deer
  # roe.tot is the total number of roe deer
  # roe.select is the predation proportion for each stage of roe deer
  # lynx.pred is the type of predation from lynx on roe deer
  # ("Predator" for the dynamic function, "selective" or "none")
  # ---
  if (lynx.pred == "predator") {
    pred.tot <- 0
    if (is.null(lynx.hunter)){
      stop("Lynx population not found")
    }
    if (is.null(kill.rate) || is.null(months.pred) || is.null(roe.tot)){
      stop("kill.rate, months.pred and roe.tot arguments need to be filled")
    }
    for (i in 1:length(lynx.hunter)){
      pred.tot <- pred.tot + kill.rate[i]*months.pred*lynx.hunter[i]
    }
    pred.rate <- (pred.tot/roe.tot)
    
    if (pred.tot <= 0 || roe.tot <= 0 || is.na(roe.tot) || is.na(pred.tot)){
      pred.rate <- 0
    }
    if (pred.rate >= 1) {
      pred.rate <- 0.999
    }
    pred.rate[is.na(pred.rate)] <- 0
    return(pred.rate) # Proportion of roe deer killed by lynx in the population
    
  } else if (lynx.pred == "selective") {
    if (is.null(roe.select)){
      stop("Predation rate not found")
    }
    pred.rate <- roe.select
    return(pred.rate)
  } else if (lynx.pred == "none") {
    pred.rate <- 0
    return(pred.rate)
  } else {
    stop("Lynx predation not available, predation type unknown")
  }
}


# Hunting mortality
HuntingHarvest <- function(sex = NULL, pop = NULL, pop.goal = NULL, pop.threshold = 0, prop.hunt = NULL, select.hunt = NULL, fem.prop = NULL, hunting.type){
  # ---
  # The function aims to add an extra mortality to imitate hunting quota.
  # When the population is approaching the quota, hunters will start to hunt 10% of the population.
  # When the population is exceeding the quota, 30% of the population will be hunted.
  # It returns the additive proportion of individuals that will be hunted.
  # ---
  # sex is the sex of the stage
  # pop is the number of individuals in the population at time t
  # pop.goal is the maximum number of individuals in the population they wants to keep
  # pop.threshold is the minimum population abundance where hunter start stopping/slowing hunting
  # prop.hunt is the fixed proportion of the population that gets hunted
  # select.hunt is the selective hunting proportion for each stage
  # fem.prop = Nb.fem.hunt / Nb.male of the same age
  # hunting.type to choose the type of hunting within:
    # "biased" hunting is hunting the same number of individuals but only males,
    # "goal" harvest is to keep the population at a certain level of abundance,
    # "relaxed" hunting is hunting only males when the population gets to a very low density,
    # "stopped" hunting is to stop hunting completely when the population gets to a very low density,
    # "fixed" proportion is hunting with a of the whole population,
    # "selective" hunting is using proportion for each stage, 
    # "none" is for no hunting
  # ---
  if (hunting.type == "biased"){
    if (is.null(select.hunt) || is.null(fem.prop)){
      stop("select.hunt and fem.prop arguments need to be filled")
    }
    if (sex == "M"){
      if (is.na(fem.prop) | is.infinite(fem.prop)){
        fem.prop <- 0
      }
      hunting.rate <- select.hunt + fem.prop
      return(hunting.rate)
    } else if (sex == "F"){
      hunting.rate <- 0
      return(hunting.rate)
    } else {
      stop("Hunting not available, sex is missing")
    }
    
  } else if (hunting.type == "goal"){
    if (is.null(pop.goal)){
      stop("Population goal argument is missing")
    }
    if (pop >= pop.goal - 0.25*pop.goal){
      hunting.rate <- 0.1
      if (pop >= pop.goal + 0.25*pop.goal){
        hunting.rate <- 0.3
      }
    } else {
      hunting.rate <- 0
    }
    return(hunting.rate)

  } else if (hunting.type == "relaxed"){
    if (is.null(select.hunt) || is.null(pop.threshold)){
      stop("select.hunt and pop.threshold arguments need to be filled")
    }
    if (pop <= pop.threshold){
      if (sex == "F"){
        hunting.rate <- 0
      } else {
        hunting.rate <- select.hunt
      }
    } else {
      hunting.rate <- select.hunt
    }
    return(hunting.rate)
    
  } else if (hunting.type == "stopped"){
    if (is.null(select.hunt) || is.null(pop.threshold)){
      stop("select.hunt and pop.threshold arguments need to be filled")
    }
    if (pop <= pop.threshold){
      hunting.rate <- 0
    } else {
      hunting.rate <- select.hunt
    }
    return(hunting.rate)
      
  } else if (hunting.type == "fixed"){
    if (is.null(prop.hunt)){
      stop("Hunting proportion is missing")
    }
    hunting.rate <- prop.hunt
    return(hunting.rate)
    
  } else if (hunting.type == "selective"){
    if (is.null(select.hunt)){
      stop("select.hunt argument is missing")
    }
    hunting.rate <- select.hunt
    return(hunting.rate)
    
  } else if (hunting.type == "none"){
    hunting.rate <- 0
    return(hunting.rate)
    
  } else {
    stop("Hunting not available, hunting type unknown")
  }
}



# Red fox predation mortality
RedFox <- function(redfox = NULL, redfoxHL = NULL, predation){
  # ---
  # This function is used to change the magnitude of the red fox predation on fawns
  # From low red foxes predation to high red foxes predation
  # ---
  # redfox is the mean predation rate
  # redfoxHL is the vector with Low and High red fox predation rate
  # predation is the argument to set the type of predation you want to use
  # ---
  if (predation == "fixed"){
    if (is.null(redfox)){
      stop("Red fox predation rate is missing")
    }
    rate <- redfox
    return(rate)
    
  } else if (predation == "low"){
    if (is.null(redfox)){
      stop("High and low red fox predation rates argument is missing")
    }
    rate <- min(redfoxHL)
    return(rate)
    
  } else if (predation == "high"){
    if (is.null(redfox)){
      stop("High and low red fox predation rates argument is missing")
    }
    rate <- max(redfoxHL)
    return(rate)
    
  } else if (predation == "none"){
    rate <- 0
    return(rate)
    
  } else {
    stop("Red fox predation not available, predation type unknown")
  }
}


# Snowfall mortality (can be turned off)
SnowEvent <- function(snow.rate, snow.freq, snow.event = F){
  # ---
  # This function simulate heavy snowfall by adding an additional mortality
  # ---
  # snow.rate is the mortality induced by this heavy snowfall
  # snow.freq is the frequency of those event happening
  # snow.event is to turn ON/OFF the snow mortality
  # ---
  if (isTRUE(snow.event)){ # If snow event is turned ON
    snow <- rbinom(1, 1, snow.freq)*snow.rate
  } else { # Otherwise, there is no effect of snow
    snow <- 0
  }
  return(snow)
}



# Stop hunting function
StopHunt <- function(snow.mortality, hunting.type, stop.hunting = F){
  # ---
  # Stop hunting after a snow catastrophic event the year before:
  # ---
  # snow.mortality is the magnitude of the mortality induced by snow event
  # hunting.type is the type of hunting running during the simulation
  # stop.hunting allows to switch ON/OFF the function
  # ---
  if (isTRUE(stop.hunting) & snow.mortality != 0){ # if a catastrophic snow event happened
    hunting <- "none" # don't hunt this year !
  } else {
    hunting <- hunting.type # keep hunting them...
  }
  return(hunting)
}



# - Density dependence function -----------------------------------------------
# Roe deer density dependence
DDroe <- function(N.density){
  # ---
  # This DD function will affect the survival rate of juveniles
  # by adding an additional mortality on fawns
  # ---
  # N.density is the number of roe deer per square kilometre
  # ---
  m.juv <- 0.015*N.density - 0.15
  return(max(0, m.juv, na.rm = T))
}



# - Matrix projection model ---------------------------------------------------

# (Stage structured - Lefkovitch)
MPM <- function(name.stages, rho, Fx, NDx, Sx){
  # ---
  # This function will create the matrix model out of the vital rate for the population projection model 
  # ---
  # name.stages is the name of all the stages in the population
  # rho is the proportion of female at birth
  # Fx is the fertility rate of each stage
  # NDx is the dispersal rate of each stage
  # Sx is the survival rate of each stage
  # ---
  matrix.pop <- matrix(0, length(name.stages), length(name.stages))
  colnames(matrix.pop) <- name.stages
  rownames(matrix.pop) <- name.stages
  # Fecundity and recruitment ---
  matrix.pop[1,2] <- rho*Fx[1] # Primiparous female producing juvenile female
  matrix.pop[1,3] <- rho*Fx[2] # First year adult female producing juvenile female
  matrix.pop[1,4] <- rho*Fx[3] # Adult female producing juvenile female
  if (length(Fx) == 5) {
    matrix.pop[1,6] <- rho*Fx[4] # First year adult male producing juvenile female
    matrix.pop[1,7] <- rho*Fx[5] # Adult male producing juvenile female
  }
  matrix.pop[5,2] <- (1-rho)*Fx[1] # Primiparous female producing juvenile male
  matrix.pop[5,3] <- (1-rho)*Fx[2] # First year adult female producing juvenile male
  matrix.pop[5,4] <- (1-rho)*Fx[3] # Adult female producing juvenile male
  if (length(Fx) == 5){
    matrix.pop[5,6] <- (1-rho)*Fx[4] # First year adult male producing juvenile male
    matrix.pop[5,7] <- (1-rho)*Fx[5] # Adult male producing juvenile male
  }
  # Natal Dispersal ---
  matrix.pop[2,2] <- NDx[1] # Yearling female dispersal
  matrix.pop[6,6] <- NDx[2] # Yearling male dispersal
  # Survival rate ---
  matrix.pop[2,1] <- Sx[1] # Juvenile female survival
  matrix.pop[3,2] <- Sx[2] # Yearling female survival
  matrix.pop[4,3] <- Sx[3] # Primiparous female survival
  matrix.pop[4,4] <- Sx[4] # Adult female survival
  matrix.pop[6,5] <- Sx[5] # Juvenile male survival
  matrix.pop[7,6] <- Sx[6] # Yearling male survival
  matrix.pop[7,7] <- Sx[7] # Adult male survival
  
  return(matrix.pop)
}



# - Data preparation function -------------------------------------------------


# From array to dataframe
ArrayToDataframe <- function(simulation){
  # ---
  # 
  # ---
  pop <- as.data.frame(ftable(simulation))
  colnames(pop) <- c("stage","time","simulation","value")
  pop <- pop[, c("simulation","time","stage","value")]
  pop$time <- as.numeric(pop$time)-1
  return(pop)
}



# Create the population sum
PopSum <- function(simulation){
  # ---
  # 
  # ---
  sum.pop <- as.data.frame(ftable(colSums(simulation)))
  # if a population reach quasi extinction in the simulation,
  # sum population to 0:
  sum.pop[is.na(sum.pop)] <- 0 
  colnames(sum.pop) <- c("time","simulation","value")
  sum.pop <- sum.pop[, c("simulation","time","value")]
  sum.pop$time <- as.numeric(sum.pop$time)-1
  return(sum.pop)
}


# Create population quantiles
PopQuantile <- function(simulation){
  # ---
  # 
  # ---
  sum.pop <- PopSum(simulation)
  pop.quant <- data.frame()
  for (i in 0:max(sum.pop$time)){
    pop.quant[i+1, 1] <- i # time
    # Nb of simulations running each time:
    pop.quant[i+1, 2] <- table(!is.na(sum.pop[which(sum.pop$time == i), "value"]))["TRUE"]
    # Quantiles:
    pop.quant[i+1, 3] <- quantile(sum.pop[which(sum.pop$time == i), "value"],
                                  c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[1] #0%
    pop.quant[i+1, 4] <- quantile(sum.pop[which(sum.pop$time == i), "value"],
                                  c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[2] #5%
    pop.quant[i+1, 5] <- quantile(sum.pop[which(sum.pop$time == i), "value"],
                                  c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[3] #25%
    pop.quant[i+1, 6] <- quantile(sum.pop[which(sum.pop$time == i), "value"],
                                  c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[4] #50%
    pop.quant[i+1, 7] <- quantile(sum.pop[which(sum.pop$time == i), "value"],
                                  c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[5] #75%
    pop.quant[i+1, 8] <- quantile(sum.pop[which(sum.pop$time == i), "value"],
                                  c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[6] #95%
    pop.quant[i+1, 9] <- quantile(sum.pop[which(sum.pop$time == i), "value"],
                                  c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[7] #100%
    pop.quant[i+1, 10] <- mean(sum.pop[which(sum.pop$time == i), "value"], na.rm = T) # Mean
  }
  colnames(pop.quant) <- c("time", "simulations", "q0", "q5", "q25", "q50", "q75", "q95", "q100", "mean")
  return(pop.quant)
}


# Growth rate
GrowthRate <- function(simulation){
  # ---
  # 
  # ---
  GR <- PopSum(simulation)
  for (i in (2:nrow(GR)-1)){
    # Lambda
    GR[i+1, 4] <- GR[i+1, 3]/GR[i, 3]
    # Growth rate percent
    GR[i+1, 5] <- ((GR[i+1, 3]-GR[i, 3])/GR[i, 3])*100
    # if the population went extinct (no individuals) 
    if (is.na(GR[i, 3]) | GR[i, 3] == 0){
      GR[i, 4] <- NA # NA or 0
      GR[i, 5] <- NA # NA or -100
    }
    # No growth rate at t=0
    if (GR$time[i] == 0){
      GR[i, 4] <- NA
      GR[i, 5] <- NA
    }
  }
  colnames(GR) <- c("simulation","time","abundance","lambda","growth.rate")
  return(GR)
}


# Growth rate quantiles
GrowthRateQuantile <- function(simulation){
  # ---
  # 
  # ---
  GR <- GrowthRate(simulation)
  GR.quant <- data.frame()
  for (i in 1:max(GR$time)){
    GR.quant[i, 1] <- i # time
    # GR.quant[i, 2] <- table(!is.na(GR[which(GR$time == i), "abundance"]))["TRUE"] # Nb of simulations running each time
    if (is.na(table(GR[which(GR$time == i), "abundance"] != 0)["TRUE"])){
      GR.quant[i, 2] <- 0
    } else {
      GR.quant[i, 2] <- table(GR[which(GR$time == i), "abundance"] != 0)["TRUE"] # Nb of simulations running each time
    }
    
    # Lambda
    GR.quant[i, 3] <- quantile(GR[which(GR$time == i), "lambda"],
                                   c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[1] #0%
    GR.quant[i, 4] <- quantile(GR[which(GR$time == i), "lambda"],
                                   c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[2] #5%
    GR.quant[i, 5] <- quantile(GR[which(GR$time == i), "lambda"],
                                   c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[3] #25%
    GR.quant[i, 6] <- quantile(GR[which(GR$time == i), "lambda"],
                                   c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[4] #50%
    GR.quant[i, 7] <- quantile(GR[which(GR$time == i), "lambda"],
                                   c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[5] #75%
    GR.quant[i, 8] <- quantile(GR[which(GR$time == i), "lambda"],
                                   c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[6] #95%
    GR.quant[i, 9] <- quantile(GR[which(GR$time == i), "lambda"],
                                   c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[7] #100%
    GR.quant[i, 10] <- mean(GR[which(GR$time == i), "lambda"], na.rm = T) # Mean
    
    # Growth rate percentage increase
    GR.quant[i, 11] <- quantile(GR[which(GR$time == i), "growth.rate"],
                                   c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[1] #0%
    GR.quant[i, 12] <- quantile(GR[which(GR$time == i), "growth.rate"],
                                    c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[2] #5%
    GR.quant[i, 13] <- quantile(GR[which(GR$time == i), "growth.rate"],
                                    c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[3] #25%
    GR.quant[i, 14] <- quantile(GR[which(GR$time == i), "growth.rate"],
                                    c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[4] #50%
    GR.quant[i, 15] <- quantile(GR[which(GR$time == i), "growth.rate"],
                                    c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[5] #75%
    GR.quant[i, 16] <- quantile(GR[which(GR$time == i), "growth.rate"],
                                    c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[6] #95%
    GR.quant[i, 17] <- quantile(GR[which(GR$time == i), "growth.rate"],
                                    c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = T)[7] #100%
    GR.quant[i, 18] <- mean(GR[which(GR$time == i), "growth.rate"], na.rm = T) # Mean
  }
  colnames(GR.quant) <- c("time", "simulations",
                          "l.q0", "l.q5", "l.q25", "l.q50",
                          "l.q75", "l.q95", "l.q100", "l.mean",
                          "p.q0", "p.q5", "p.q25", "p.q50",
                          "p.q75", "p.q95", "p.q100", "p.mean")
  return(GR.quant)
}


# Proportion of populations at quasi-extinction and irruption
StateProp <- function(simulation, Ne, Ni){
  # ---
  # 
  # ---
  sum.pop <- PopSum(simulation)
  state <- data.frame()
  for (i in 0:max(sum.pop$time)){
    state[i+1, 1] <- i
    state[i+1, 2] <- prop.table(table(is.na(sum.pop[which(sum.pop$time == i), 3]) | sum.pop[which(sum.pop$time == i), 3] < Ne))["TRUE"]
    state[i+1, 3] <- prop.table((table(sum.pop[which(sum.pop$time == i), 3] > Ni)))["TRUE"]
    state[is.na(state)] <- 0
    state[i+1, 4] <- 1 - state[i+1, 3] - state[i+1, 2]
  }
  colnames(state) <- c("time", "extinct", "irruption", "stable")
  state <- state[, c("time", "extinct", "stable", "irruption")]
  state.pivot <- pivot_longer(state, 2:4, names_to = "state", values_to = "proportion", cols_vary = "slowest")
  state.pivot$state <- factor(state.pivot$state, levels = c("irruption", "stable", "extinct"))
  return(state.pivot)
}



# - Plot function -------------------------------------------------------------

DensityPlot <- function(data, study.area, species, x.years = NULL, y.lim = NULL, carrying.capacity = NULL){
  # Prepare data
  quant.data <- PopQuantile(data)
  if (species == "Lynx"){quant.data[, 3:10] <- quant.data[, 3:10]*100}
  # Plot it
  ggplot(quant.data) +
    geom_ribbon(aes(x=time, ymin=q95/study.area, ymax=q5/study.area), fill="darkolivegreen", alpha=0.5) +
    geom_ribbon(aes(x=time, ymin=q95/study.area, ymax=q5/study.area), fill="grey", alpha=0.5) +
    geom_ribbon(aes(x=time, ymin=q75/study.area, ymax=q25/study.area), fill="darkolivegreen", alpha=0.5) +
    geom_line(aes(x=time, y=q0/study.area), color="darkgrey", linetype="dashed", size=0.5) +
    geom_line(aes(x=time, y=mean/study.area), color="black", size=1) +
    geom_line(aes(x=time, y=q100/study.area), color="darkgrey", linetype="dashed", size=0.5) +
    {if (!is.null(carrying.capacity))geom_hline(yintercept = carrying.capacity)} +
    theme(legend.key.height= unit(2.5, 'cm'),
          text = element_text(size = 15),
          axis.text = element_text(size = 15)) +
    xlim(0, min(x.years, max(quant.data$time))) +
    ylim(0, max(y.lim, max(quant.data[1:min(x.years, max(quant.data$time)), 3:10]/study.area))) +
    labs(title = paste(species, "population density"), x = "Time (Years)",
         y = case_when(species == "Lynx" ~ "Population density (N/100km2)",
                       species == "Roe deer" ~ "Population density (N/km2)"))
}



GrowthRatePlot <- function(data, study.area, species, x.years = NULL){
  # Prepare data
  GR.quant.data <- GrowthRateQuantile(data)
  # Plot it
  ggplot(GR.quant.data) +
    geom_ribbon(aes(x=time, ymin=l.q95, ymax=l.q5), fill="darkolivegreen", alpha=0.5) +
    geom_ribbon(aes(x=time, ymin=l.q95, ymax=l.q5), fill="grey", alpha=0.5) +
    geom_ribbon(aes(x=time, ymin=l.q75, ymax=l.q25), fill="darkolivegreen", alpha=0.5) +
    geom_line(aes(x=time, y=l.mean), color="black", size=1) +
    geom_line(aes(x=time, y=l.q0), color="darkgrey", size=0.5, linetype="dashed") +
    geom_line(aes(x=time, y=l.q100), color="darkgrey", size=0.5, linetype="dashed") +
    geom_hline(yintercept=1, linetype="dotdash", color="black", size=0.75) +
    theme(legend.key.height= unit(2.5, 'cm'),
          text = element_text(size = 15),
          axis.text = element_text(size = 15)) +
    xlim(0, min(x.years, max(GR.quant.data$time))) +
    ylim(0, 2) +
    labs(title = paste(species, "population growth rate"), x = "Time (Years)", y = "Population growth rate (λ)")
}



StatePropPlot <- function(data, study.area, species, Ne, Ni, x.years = NULL){
  # Prepare data
  data.state <- StateProp(data, Ne=Ne, Ni=Ni)
  # Plot it
  ggplot(data.state, aes(x=time, y=proportion*100, fill=state)) +
    geom_col(position = "stack", width = 0.75) +
    scale_fill_manual("State", values = c("irruption" = "deepskyblue4", "stable" = "snow3", "extinct" = "darkred")) +
    theme(legend.key.height= unit(2.5, 'cm'),
          text = element_text(size = 15),
          axis.text = element_text(size = 15)) +
    xlim(0, min(x.years, max(data.state$time))) +
    labs(title = "State proportion of populations", x = "Time (Years)", y = "Proportion (%)")
}


