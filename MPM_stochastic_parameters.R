# ------------------------------- #
#          - Chapter I -          #
#             SETTINGS            #
#     Packages and parameters     #
# ------------------------------- #



# ----- Load default parameters values -----
source("MPM_load_parameters.R")

# --------------------------------.
### --- ROE DEER POPULATION ------
# - Stochastic parameters values -.
# --------------------------------.

# ----- Reproduction and recruitment parameters -----
# Average litter size: Nb of offspring per 2 year-old female
k1.r <- max(1, rnorm(1, k1.r, 0.38)) # k1.r <- runif(1, k1.r-0.45, k1.r+0.45)
# Average litter size: Nb of offspring per adult female
k2.r <- max(1, rnorm(1, k2.r, 0.73)) # k2.r <- runif(1, k2.r-0.83, k2.r+0.83)
k.r <- c(k1.r, k2.r) # average litter size vector
# Proportion of female that reproduce
pr1.r <- runif(1, pr1.r-0.1385, pr1.r+0.1385) # First reproduction
pr2.r <- runif(1, pr2.r-0.0708, pr2.r+0.0708) # Adult reproduction
pr.r <- c(pr1.r, pr2.r) # Breeding proportion vector


# ----- Survival parameters -----
# Mortality - Used in survival rate function calculation

if (exists("Simulation.roe")){
  if (stoch.distrib == "norm"){
    # Normal distribution, estimated from 95% CI (Melis et al., 2013)
    # print("Normal distrib")
    # --- Traffic, natural and unknown mortality:
    m.jroef <- rtruncnorm(1, a=0, b=1, mean = m.jroef, sd=0.056)
    m.yroef <- rtruncnorm(1, a=0, b=1, mean = m.yroef, sd=0.093)
    m.proef <- rtruncnorm(1, a=0, b=1, mean = m.proef, sd=0.044)
    m.aroef <- rtruncnorm(1, a=0, b=1, mean = m.aroef, sd=0.044)
    
    m.jroem <- rtruncnorm(1, a=0, b=1, mean = m.jroem, sd=0.066)
    m.yroem <- rtruncnorm(1, a=0, b=1, mean = m.yroem, sd=0.089)
    m.aroem <- rtruncnorm(1, a=0, b=1, mean = m.aroem, sd=0.049)
    
    # --- Red fox predation mortality:
    m.redfox <- rtruncnorm(1, a=0, b=1, mean = m.redfox, sd=0.047)
    
    # --- Lynx predation mortality:
    m.lynxpred <- rtruncnorm(1, a=0, b=1, mean = m.lynxpred, sd=0.027)
    
    pred.jroef <- rtruncnorm(1, a=0, b=1, mean = pred.jroef, sd=0.049)
    pred.yroef <- rtruncnorm(1, a=0, b=1, mean = pred.yroef, sd=0.083)
    pred.proef <- rtruncnorm(1, a=0, b=1, mean = pred.proef, sd=0.035)
    pred.aroef <- rtruncnorm(1, a=0, b=1, mean = pred.aroef, sd=0.035)
    
    pred.jroem <- rtruncnorm(1, a=0, b=1, mean = pred.jroem, sd=0.051)
    pred.yroem <- rtruncnorm(1, a=0, b=1, mean = pred.yroem, sd=0.113)
    pred.aroem <- rtruncnorm(1, a=0, b=1, mean = pred.aroem, sd=0.048)
    
    # --- Hunting mortality:
    m.hunting.r <- rtruncnorm(1, a=0, b=1, mean = m.hunting.r, sd=0.030)
    
    hunt.jroef <- rtruncnorm(1, a=0, b=1, mean = hunt.jroef, sd=0.051)
    hunt.yroef <- rtruncnorm(1, a=0, b=1, mean = hunt.yroef, sd=0.121)
    hunt.proef <- rtruncnorm(1, a=0, b=1, mean = hunt.proef, sd=0.030)
    hunt.aroef <- rtruncnorm(1, a=0, b=1, mean = hunt.aroef, sd=0.030)
    
    hunt.jroem <- rtruncnorm(1, a=0, b=1, mean = hunt.jroem, sd=0.047)
    hunt.yroem <- rtruncnorm(1, a=0, b=1, mean = hunt.yroem, sd=0.118)
    hunt.aroem <- rtruncnorm(1, a=0, b=1, mean = hunt.aroem, sd=0.057)
    
  } else if (stoch.distrib == "unif"){
    # Uniform distribution
    # print("Uniform distrib")
    # --- Traffic, natural and unknown mortality:
    m.jroef <- runif(1, 0.082, 0.267)
    m.yroef <- runif(1, 0.088, 0.393)
    m.proef <- runif(1, 0.170, 0.316)
    m.aroef <- runif(1, 0.170, 0.316)
    
    m.jroem <- runif(1, 0.159, 0.377)
    m.yroem <- runif(1, 0.000, 0.267)
    m.aroem <- runif(1, 0.081, 0.243)
    
    # --- Red fox predation mortality:
    m.redfox <- runif(1, 0.183, 0.343)
    
    # --- Lynx predation mortality:
    m.lynxpred <- runif(1, 0.098, 0.189)
    
    pred.jroef <- runif(1, 0.060, 0.221)
    pred.yroef <- runif(1, 0.000, 0.270)
    pred.proef <- runif(1, 0.082, 0.198)
    pred.aroef <- runif(1, 0.082, 0.198)
  
    pred.jroem <- runif(1, 0.052, 0.220)
    pred.yroem <- runif(1, 0.027, 0.400)
    pred.aroem <- runif(1, 0.073, 0.232)
    
    # --- Hunting mortality:
    m.hunting.r <- runif(1, 0.090, 0.187)
  
    hunt.jroef <- runif(1, 0.033, 0.188)
    hunt.yroef <- runif(1, 0.128, 0.518)
    hunt.proef <- runif(1, 0.028, 0.124)
    hunt.aroef <- runif(1, 0.028, 0.124)
  
    hunt.jroem <- runif(1, 0.039, 0.191)
    hunt.yroem <- runif(1, 0.028, 0.418)
    hunt.aroem <- runif(1, 0.137, 0.325)
    
  } else {
    stop("Stochasticity distribution not available")
  }
}





# --------------------------------.
### --- LYNX POPULATION ----------
# - Stochastic parameters values -.
# --------------------------------.

# ----- Reproduction and recruitment parameters -----
# Average litter size: Nb of offspring per 2 year-old female
k1.l <- max(1, rnorm(1, k1.l, 0.59)) # k1.l <- runif(1, k1.l-0.96, k1.l+1.78)
# k1.l <- rnorm(1, k1.l, 2.31) # SEM calculated in Excel sheet
# Average litter size: Nb of offspring per adult female
k2.l <- max(1, rnorm(1, k2.l, 0.19)) # k2.l <- runif(1, k2.l-0.27, k2.l+0.31)
# assign("k2.l", rnorm(1, k2.l, 0.49) # SEM calculated in Excel sheet
k.l <- c(k1.l, k2.l) # average litter size vector
# Proportion of female that reproduce
pr1.l <- runif(1, 0.09, 0.81) # First reproduction
pr2.l <- runif(1, 0.53, 0.81) # Adult reproduction
pr.l <- c(pr1.l, pr2.l) # Breeding proportion vector


# ----- Survival parameters -----
# Mortality - Used in survival rate calculation
# --- Traffic kills, natural and unknown mortality: (+ poaching)
m.jlynf <- rtruncnorm(1, a=0, b=1, mean = m.jlynf, sd=0.076)
m.ylynf <- rtruncnorm(1, a=0, b=1, mean = m.ylynf, sd=0.113)
m.plynf <- rtruncnorm(1, a=0, b=1, mean = m.plynf, sd=0.039)
m.alynf <- rtruncnorm(1, a=0, b=1, mean = m.alynf, sd=0.039)

m.jlynm <- rtruncnorm(1, a=0, b=1, mean = m.jlynm, sd=0.059)
m.ylynm <- rtruncnorm(1, a=0, b=1, mean = m.ylynm, sd=0.104)
m.alynm <- rtruncnorm(1, a=0, b=1, mean = m.alynm, sd=0.057)

# m.jlynf <- runif(1, max(0, m.jlynf-0.126), min(m.jlynf+0.126, 1))
# m.ylynf <- runif(1, max(0, m.ylynf-0.187), min(m.ylynf+0.187, 1))
# m.plynf <- runif(1, max(0, m.plynf-0.05), min(m.plynf+0.065, 1))
# m.alynf <- runif(1, max(0, m.alynf-0.05), min(m.alynf+0.065, 1))
# 
# m.jlynm <- runif(1, max(0, m.jlynm-0.097), min(m.jlynm+0.097, 1))
# m.jlynm <- runif(1, max(0, m.ylynm-0.170), min(m.ylynm+0.170, 1))
# m.alynm <- runif(1, max(0, m.alynm-0.05), min(m.alynm+0.257, 1))


# ----- Predator behaviour ??
# pred.lynFS <- 2.71 # Number of roe deer killed per month per lynx solitary female
# pred.lynFK <- 6.23 # Number of roe deer killed per month per lynx female with kittens
# pred.lynM <- 4.85 # Number of roe deer killed per month per lynx male
# 
# 
# pred.lyn <- c(pred.lynFS, pred.lynFK, pred.lynM)
# 
# half.sat <- 0.37 #0.33 # Half saturation density of prey (Nilsen et al., 2009)
# kill.a <- 10.14 # Asymptotic kill rate (roe deer / 100 days)

