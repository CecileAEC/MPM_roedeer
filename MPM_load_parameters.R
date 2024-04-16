                    # ------------------------------- #
                    #          - Chapter I -          #
                    #             SETTINGS            #
                    #    Parameters default values    #
                    # ------------------------------- #

# --------------------------------------------------------------------------- #
#                              Here make changes                              #
#                  !!! Don't forget to save your changes !!!                  #
# --------------------------------------------------------------------------- #


# ----- Parameters - Set default values -----

# ----- Simulation (Also set up in the simulation file)
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
# Range from 5000 to 10000 km2
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
# roe.density <- seq(0, 40, by=0.5) # or a gradient in the simulation file

# - Initial state (set default abundance "density")
# Number of roe deer in every stage at t0 from population density

# - Based on sex and age distribution of a population of 213 roe deer
# in Kalø, Denmark
# Jroef_init <- round(0.2159624*N_Roe_pop) #46
# Yroef_init <- round(0.1126761*N_Roe_pop) #24
# Proef_init <- round(0.09859155*N_Roe_pop) #21
# Aroef_init <- round(0.1455399*N_Roe_pop) #31
# 
# Jroem_init <- round(0.2112676*N_Roe_pop) #45
# Yroem_init <- round(0.07981221*N_Roe_pop) #17
# Aroem_init <- round(0.1361502*N_Roe_pop) #29
# N.roe.init <- VectorPop(round(c(0.216, 0.113, 0.099, 0.145, 0.211, 0.080, 0.136)
#                               * (roe.density*study.area)), nameRoeDeer)

# - Population abundance vector started at stable stage (???)
# Stable stage:
# sum(c(0.237, 0.105, 0.045, 0.147, 0.237, 0.082, 0.147)) #???
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

# Storfosna data 1992
# m.jroef <- 1-0.70
# m.yroef <- 1-0.90
# m.proef <- 1-0.91
# m.aroef <- 1-0.91
# 
# m.aroef <- 1-0.73
# m.yroem <- 1-0.89
# m.aroem <- 1-0.93

# Storfosna data 1993
# m.jroef <- 1-0.61
# m.yroef <- 1-0.64
# m.proef <- 1-0.94
# m.aroef <- 1-0.94
# 
# m.aroef <- 1-0.62
# m.yroem <- 1-0.72
# m.aroem <- 1-0.71

# Additional mortality factor:
# roadkill.r <- 0.15 # Included in the basic mortality? (from Grøtan et al., 2005)
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
  # stagepred <- c(pred.jroef, pred.yroef, pred.proef, pred.aroef, pred.jroem, pred.yroem, pred.aroem)
  
# roe.goal <- 20 # Goal population size (in density)
roe.threshold <- 50 # Percent of the population left before stopping hunting
m.hunting.r <- 0.138 # Fixed proportion of the population harvested (Melis et al., 2013)
  hunt.jroef <- 0.111 # Hunting proportion of each stage (selection) (Melis et al., 2013)
  hunt.yroef <- 0.323 # Hunting proportion of each stage (selection) (Melis et al., 2013)
  hunt.proef <- 0.076 # Hunting proportion of each stage (selection) (Melis et al., 2013)
  hunt.aroef <- 0.076 # Hunting proportion of each stage (selection) (Melis et al., 2013)
  hunt.jroem <- 0.115 # Hunting proportion of each stage (selection) (Melis et al., 2013)
  hunt.yroem <- 0.223 # Hunting proportion of each stage (selection) (Melis et al., 2013)
  hunt.aroem <- 0.231 # Hunting proportion of each stage (selection) (Melis et al., 2013)
  # select.hunt.r <- c(hunt.jroef, hunt.yroef, hunt.proef, hunt.aroef, hunt.jroem, hunt.yroem, hunt.aroem)

# snowevent.r <- 0 # For when it is fixed
snow.frequency <- 0.08 # 0.04, 0.08 or 0.12 # How often the heavy snow happens
snow.magnitude <- 0.30 # How much it impact roe deer survival


# ----- Population density and dispersal rate -----
# Find the parameter (kernel distribution?) for dispersal
# Is it Density dependent?
NDf.r <- 0 # Drate(0) # Natal dispersal rate female (rate of income/outcome)
NDm.r <- 0 # Drate(0) # Natal dispersal rate male (rate of income/outcome)





# ---------------------------.
### --- LYNX POPULATION -----
# ---------------------------.

# ----- Set lynx stages population name
nameLynx <- c("Jlyn.F", # Juvenile lynx female
              "Ylyn.F", # Yearling lynx female
              "Plyn.F", # Two year-old lynx female
              "Alyn.F", # Adult lynx female
              "Jlyn.M", # Juvenile lynx male
              "Ylyn.M", # Yearling lynx male
              "Alyn.M") # Adult lynx male


# ----- Set initial density -----
lyn.density <- 0.3 # Nb of lynx per 100 square kilometre
# lyn.density <- seq(0.1, 6, by=0.1) # or a gradient in the simulation file

# - Initial state (set default abundance "density")
# Number of lynx in every stage at t0 from population density

# - Based on Number of individuals, Andrén et al., 2002 - Table 1 (Hedmark data)
# Jlynf_init <- round(0.182*N_Lyn_pop) # 12
# Ylynf_init <- round(0.167*N_Lyn_pop) # 11
# Plynf_init <- round(0.060*N_Lyn_pop) # 4
# Alynf_init <- round(0.136*N_Lyn_pop) # 9 (#13 tot)
# 
# Jlynm_init <- round(0.182*N_Lyn_pop) # 12
# Ylynm_init <- round(0.152*N_Lyn_pop) # 10
# Alynm_init <- round(0.121*N_Lyn_pop) # 8
# N.lyn.init <- VectorPop(round(c(0.182, 0.167, 0.060, 0.136, 0.182, 0.152, 0.121)
#                               * ((lyn.density/100)*study.area)), nameLynx)

# - Population abundance vector started at stable stage (???)
# Stable stage:
# sum(c(0.182, 0.167, 0.060, 0.136, 0.182, 0.152, 0.121)) #???
N.lyn.init <- VectorPop(round(c(0.213, 0.112, 0.042, 0.174, 0.213, 0.059, 0.187)
                              * ((lyn.density/100)*study.area)), nameLynx)

# Set a quasi extinction threshold:
Ne.lyn <- (0.05/100)*study.area # From density



# ----- Reproduction and recruitment parameters -----
rho.l <- 0.5 # Proportion of female lynx at birth
k1.l <- 2.09 # Average litter size: Nb of offspring per 2 year-old female
k2.l <- 2.10 # Average litter size: Nb of offspring per adult female
k.l <- c(k1.l, k2.l) # average litter size vector
pr1.l <- 0.40 # Proportion of female that first reproduce
pr2.l <- 0.69 # Proportion of adult female that reproduce
pr.l <- c(pr1.l, pr2.l) # Breeding proportion vector
h.l <- 5  # Number of female that a male can reproduce with during the season
alpha.l <- 15 # Strength parameter of the exponential penalty function



# ----- Survival parameters -----
# Mortality - Used in survival rate calculation
# Basic mortality:
m.jlynf <- 0.407
m.ylynf <- 0.571
m.plynf <- 0.056
m.alynf <- 0.056

m.jlynm <- 0.635
m.ylynm <- 0.271
m.alynm <- 0.03 # corrected 

# Additional mortality factor:
# roadkill.l <- 0.05 # Increase with density?
# poaching.l <- 0.12 # Get some additional poaching?
# poaching and traffic road kill ?? in addition to the one included in natural mortality??

lyn.goal <- 0.5 # Goal population size (in density, indiv/100km2)
m.hunting.l <- 0.15 # Fixed proportion of the population harvested
  hunt.jlynf <- 0.03 # Proportion of each stage harvested (Andrén et al., 2006 - Table 1) - Corrected
  hunt.ylynf <- 0.03 # Proportion of each stage harvested (Andrén et al., 2006 - Table 1) - Corrected
  hunt.plynf <- 0.083 # Proportion of each stage harvested (Andrén et al., 2006 - Table 1)
  hunt.alynf <- 0.083 # Proportion of each stage harvested (Andrén et al., 2006 - Table 1)
  hunt.jlynm <- 0.067 # Proportion of each stage harvested (Andrén et al., 2006 - Table 1)
  hunt.ylynm <- 0.146 # Proportion of each stage harvested (Andrén et al., 2006 - Table 1)
  hunt.alynm <- 0.117 # Proportion of each stage harvested (Andrén et al., 2006 - Table 1)
  # select.hunt.l <- c(hunt.jlynf, hunt.ylynf, hunt.plynf, hunt.alynf, hunt.jlynm, hunt.ylynm, hunt.alynm)


# ----- Population density and dispersal rate -----
# Find the parameter (kernel distribution?) for dispersal
# Is it Density dependent?
NDf.l <- 0 # Drate(0) # Natal dispersal rate female (rate of income/outcome)
NDm.l <- 0 # Drate(0) # Natal dispersal rate male (rate of income/outcome)


# ----- Predator prey parameter -----
pred.lynFS <- 2.71 # Number of roe deer killed per month per lynx solitary female
pred.lynFK <- 6.23 # Number of roe deer killed per month per lynx female with kittens
pred.lynM <- 4.85 # Number of roe deer killed per month per lynx male


pred.lyn <- c(pred.lynFS, pred.lynFK, pred.lynM)

half.sat <- 0.37 #0.33 # Half saturation density of prey (Nilsen et al., 2009)
kill.a <- 10.14 # Asymptotic kill rate (roe deer / 100 days)

