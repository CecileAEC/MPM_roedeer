                    # ------------------------------- #
                    #             SETTINGS            #
                    #    Parameters default values    #
                    # ------------------------------- #

# --------------------------------------------------------------------------- #
#                              Here make changes                              #
#        !!! Don't forget to save your changes before simulation !!!          #
# --------------------------------------------------------------------------- #


# ----- Parameters - Set default values -----

# ----- Simulation (Set up in the simulation file)
# Set number of simulations:
# nb.sim <- 500
# Set time simulation (years):
# nYears <- 25

#  ----- Stochasticity in simulation True/False
initial.stoch <- T # Initial density stochasticity
init.var <- 0.1 # Percentage variation around the initial density
envir.stoch <- T # Parameter stochasticity
# Draw stochasticity from uniform or normal distribution for mortality rates
stoch.distrib <- "norm" # "unif" or "norm"

# ----- Reproduction type
# Select the type of reproduction function to use
R.type <- "FemaleFF" # "FemaleFF" or "ExpFF" or "HFF"


# ----- Scenarios
Redfox.predation <- c("fixed", "none") # "low", "high", 
Lynx.predation <- c("selective", "none")
# Hunting type:
# "biased", "goal", "relaxed", "stopped", "fixed", "selective", "none"
Hunting.type.roe <- c("selective", "relaxed", "stopped", "none")
Hunting.type.lynx <- c("goal", "selective", "none")
Snow.ON.OFF <- F # ON:True or OFF:False
Stop.hunt.ON.OFF <- F # ON:True or OFF:False

# ----- Study site
study.area <- 10000 # Study area in square kilometre

# -------------------------------.
### --- ROE DEER POPULATION -----
# -------------------------------.

# ----- Set roe deer stage population name
nameRoeDeer <- c("Jroe.F", # Juvenile roe deer female
                 "Yroe.F", # Yearling roe deer female
                 "Proe.F", # Two year-old roe deer female
                 "Aroe.F", # Adult roe deer female
                 "Jroe.M", # Juvenile roe deer male
                 "Yroe.M", # Yearling roe deer male
                 "Aroe.M") # Adult roe deer male

# ----- Set initial density -----
roe.density <- 10 # Nb of roe deer per square kilometre

# - Initial population structure (set abundance from density)
# Number of roe deer in every stage at t0 from population density
# - Population abundance vector started based on sex and age distribution
# of a population of 213 roe deer in Kalø, Denmark
N.roe.init <- VectorPop(round(c(0.237, 0.105, 0.045, 0.147, 0.237, 0.082, 0.147)
                              * (roe.density*study.area)), nameRoeDeer)

# Set a quasi extinction threshold:
Ne.roe <- 0.03*study.area # From density i.e. 300 individuals



# ----- Reproduction and recruitment parameters -----
rho.r <- 0.5 # Proportion of female roe deer at birth
k1.r <- 2.14 # Average litter size: Nb of offspring per 2 year-old female
k2.r <- 2.18 # Average litter size: Nb of offspring per adult female
k.r <- c(k1.r, k2.r) # average litter size vector
pr1.r <- 0.778 # Proportion of female that first reproduce
pr2.r <- 0.846 # Proportion of adult female that reproduce
pr.r <- c(pr1.r, pr2.r) # Breeding proportion vector
h.r <- 10 # Number of female that a male can reproduce with during the season
alpha.r <- 25 # Strength parameter of the exponential penalty function



# ----- Survival parameters -----
# Mortality - Used in survival rate function calculation
# Basic mortality:
m.jroef <- 0.175
m.yroef <- 0.241
m.proef <- 0.243
m.aroef <- 0.243

m.jroem <- 0.268
m.yroem <- 0.116
m.aroem <- 0.162

# Additional mortality factor:
m.redfox <- 0.26 # Red fox predation mortality rate (Melis et al., 2013)
redfox.HL <- c(0.183, 0.343) # High and low redfox predation

m.lynxpred <- 0.143 # Fixed proportion of the population (Melis et al., 2013) # Removed from function!
pred.jroef <- 0.141 # Predation proportion of each stage (Melis et al., 2013)
pred.yroef <- 0.131 # Predation proportion of each stage (Melis et al., 2013)
pred.proef <- 0.140 # Predation proportion of each stage (Melis et al., 2013)
pred.aroef <- 0.140 # Predation proportion of each stage (Melis et al., 2013)
pred.jroem <- 0.136 # Predation proportion of each stage (Melis et al., 2013)
pred.yroem <- 0.213 # Predation proportion of each stage (Melis et al., 2013)
pred.aroem <- 0.152 # Predation proportion of each stage (Melis et al., 2013)


m.hunting.r <- 0.138 # Fixed proportion of the population harvested (Melis et al., 2013)
hunt.jroef <- 0.111 # Hunting proportion of each stage (selection) (Melis et al., 2013)
hunt.yroef <- 0.323 # Hunting proportion of each stage (selection) (Melis et al., 2013)
hunt.proef <- 0.076 # Hunting proportion of each stage (selection) (Melis et al., 2013)
hunt.aroef <- 0.076 # Hunting proportion of each stage (selection) (Melis et al., 2013)
hunt.jroem <- 0.115 # Hunting proportion of each stage (selection) (Melis et al., 2013)
hunt.yroem <- 0.223 # Hunting proportion of each stage (selection) (Melis et al., 2013)
hunt.aroem <- 0.231 # Hunting proportion of each stage (selection) (Melis et al., 2013)

# snowevent.r <- 0 # For when it is fixed
snow.frequency <- 0.08 # 0.04, 0.08 or 0.12 # How often the heavy snow happens
snow.magnitude <- 0.30 # How much it impact roe deer survival


# ----- Population density and dispersal rate -----
NDf.r <- 0 # Natal dispersal rate female (rate of income/outcome)
NDm.r <- 0 # Natal dispersal rate male (rate of income/outcome)
