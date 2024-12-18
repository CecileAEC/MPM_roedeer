# ------------------------------- #
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
    # m.lynxpred <- rtruncnorm(1, a=0, b=1, mean = m.lynxpred, sd=0.027)
    
    pred.jroef <- rtruncnorm(1, a=0, b=1, mean = pred.jroef, sd=0.049)
    pred.yroef <- rtruncnorm(1, a=0, b=1, mean = pred.yroef, sd=0.083)
    pred.proef <- rtruncnorm(1, a=0, b=1, mean = pred.proef, sd=0.035)
    pred.aroef <- rtruncnorm(1, a=0, b=1, mean = pred.aroef, sd=0.035)
    
    pred.jroem <- rtruncnorm(1, a=0, b=1, mean = pred.jroem, sd=0.051)
    pred.yroem <- rtruncnorm(1, a=0, b=1, mean = pred.yroem, sd=0.113)
    pred.aroem <- rtruncnorm(1, a=0, b=1, mean = pred.aroem, sd=0.048)
    
    # --- Hunting mortality:
    # m.hunting.r <- rtruncnorm(1, a=0, b=1, mean = m.hunting.r, sd=0.030)
    
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
    # m.lynxpred <- runif(1, 0.098, 0.189)
    
    pred.jroef <- runif(1, 0.060, 0.221)
    pred.yroef <- runif(1, 0.000, 0.270)
    pred.proef <- runif(1, 0.082, 0.198)
    pred.aroef <- runif(1, 0.082, 0.198)
  
    pred.jroem <- runif(1, 0.052, 0.220)
    pred.yroem <- runif(1, 0.027, 0.400)
    pred.aroem <- runif(1, 0.073, 0.232)
    
    # --- Hunting mortality:
    # m.hunting.r <- runif(1, 0.090, 0.187)
  
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

