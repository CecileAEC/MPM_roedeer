                    # ------------------------------- #
                    #          - Chapter I -          #
                    #           SIMULATIONS           #
                    #        Lynx and Roe deer        #
                    #           populations           #
                    # ------------------------------- #



# -------------------------------.
### --- ROE DEER POPULATION -----
# -------------------------------.

rm(list = ls())
# Source the functions and parameters files
source("MPM_packages_and_functions.R")
source("MPM_load_parameters.R")

# --- Simulations, Roe deer population -----

# --- Two loop simulation:
# The first one for the different simulations
# (i.e. to analyse different results with stochasticity)
# The second one for the time years simulation
# (i.e. how long the simulation last)


# Set number of simulations:
nb.sim <- 500
# Set time simulation (years):
nYears <- 25
# # Stochasticity in simulation True/False
# initial.stoch <- T # Initial density stochasticity
# init.var <- 0.1 # Percentage variation around the initial density
# envir.stoch <- T # Parameter stochasticity
# # Draw stochasticity from uniform or normal distribution for mortality rates
# stoch.distrib <- "norm" # "unif" or "norm"
# # Select the type of reproduction function
# R.type <- "FemaleFF" # "FemaleFF" or "ExpFF" or "HFF"
# Initialise list of stochastic matrices
mat.roe.list <- list()
# Initialise scenario
scena <- expand.grid("redfox" = Redfox.predation,
                     "lynx" = Lynx.predation,
                     "hunt" = Hunting.type.roe)
## - Select only some scenario...
scena <- rbind(scena[5, ], scena[9, ])
# scena <- scena[5, ]

# Initialise the simulation array:
Simulation.roe <- array(rep(NA, length(N.roe.init)*nYears*nb.sim),
                        c(length(N.roe.init), nYears+1, nb.sim),
                        list(nameRoeDeer, seq(0, nYears),
                             paste('sim', seq(1, nb.sim))))

# --- Start simulation
# tic("simulation")
for (i in 1:nrow(scena)){ # For each combination of scenario

  # - Initialise scenarios:
  Redfox.pred <- scena[i, "redfox"]
  Lynx.pred <- scena[i, "lynx"]
  Hunting.roe <- scena[i, "hunt"]
  
  # - Initialise the simulation array for each scenario:
  Simulation.roe <- array(rep(NA, length(N.roe.init)*nYears*nb.sim),
                          c(length(N.roe.init), nYears+1, nb.sim),
                          list(nameRoeDeer, seq(0, nYears),
                               paste('sim', seq(1, nb.sim))))
  
  for (s in 1:(nb.sim)){ # For each simulation
    
    # - Stochastic population structure and abundance
    N.roe.init <- Stoch.init(roe.density, nameRoeDeer, study.area, init.var, initial.stoch)

    # - Initialise population vector
    N.roe <- N.roe.init
    Simulation.roe[, 1, s] <- N.roe
    
    # - Initialise survival rate at t0
    # No catastrophic snow event before simulation started
    m.snow <- 0

    # Run the simulations over the years...
    for (t in 2:(nYears+1)){ # For each time step
      
      # ----- Environmental stochasticity
      if (isTRUE(envir.stoch)){
        # Reload the parameter file with stochasticity
        source("MPM_stochastic_parameters.R")
      }
      
      # ----- Density dependence function:
      m.dens <- DDroe(sum(N.roe)/study.area)
      
      # ----- Vital rates at each time step
      # Changes related to population abundance
      
      # --- Survival rates ---
      # - Snow catastrophic event:
      # Stop hunting after a snow event the year before
      Hunting.roe <- StopHunt(m.snow, scena[i, "hunt"], Stop.hunt.ON.OFF)
      # Snow event mortality rate:
      m.snow <- SnowEvent(snow.rate = snow.magnitude, snow.freq = snow.frequency, Snow.ON.OFF)
      # - Survival function from mortality rates:
      S.Jroef <- Srate(c(m.jroef, m.dens,
                         RedFox(m.redfox, redfox.HL, Redfox.pred),
                         LynxPred(roe.select = pred.jroef, lynx.pred = Lynx.pred),
                         HuntingHarvest("F", pop = sum(N.roe), pop.goal = roe.goal*study.area, prop.hunt = m.hunting.r, select.hunt = hunt.jroef, hunting.type = Hunting.roe),
                         m.snow))
      S.Yroef <- Srate(c(m.yroef,
                         LynxPred(roe.select = pred.yroef, lynx.pred = Lynx.pred),
                         HuntingHarvest("F", pop = sum(N.roe), pop.goal = roe.goal*study.area, prop.hunt = m.hunting.r, select.hunt = hunt.yroef, hunting.type = Hunting.roe),
                         m.snow))
      S.Proef <- Srate(c(m.proef,
                         LynxPred(roe.select = pred.proef, lynx.pred = Lynx.pred),
                         HuntingHarvest("F", pop = sum(N.roe), pop.goal = roe.goal*study.area, prop.hunt = m.hunting.r, select.hunt = hunt.proef, hunting.type = Hunting.roe),
                         m.snow))
      S.Aroef <- Srate(c(m.aroef,
                         LynxPred(roe.select = pred.aroef, lynx.pred = Lynx.pred),
                         HuntingHarvest("F", pop = sum(N.roe), pop.goal = roe.goal*study.area, prop.hunt = m.hunting.r, select.hunt = hunt.aroef, hunting.type = Hunting.roe),
                         m.snow))
  
      S.Jroem <- Srate(c(m.aroef, m.dens,
                         RedFox(m.redfox, redfox.HL, Redfox.pred),
                         LynxPred(roe.select = pred.jroem, lynx.pred = Lynx.pred),
                         HuntingHarvest("M", pop = sum(N.roe), pop.goal = roe.goal*study.area, prop.hunt = m.hunting.r, select.hunt = hunt.jroem, fem.prop = (N.roe[1] * hunt.jroef)/N.roe[5], hunting.type = Hunting.roe),
                         m.snow))
      S.Yroem <- Srate(c(m.yroem,
                         LynxPred(roe.select = pred.yroem, lynx.pred = Lynx.pred),
                         HuntingHarvest("M", pop = sum(N.roe), pop.goal = roe.goal*study.area, prop.hunt = m.hunting.r, select.hunt = hunt.yroem, fem.prop = (N.roe[2] * hunt.yroef)/N.roe[6], hunting.type = Hunting.roe),
                         m.snow))
      S.Aroem <- Srate(c(m.aroem,
                         LynxPred(roe.select = pred.aroem, lynx.pred = Lynx.pred),
                         HuntingHarvest("M", pop = sum(N.roe), pop.goal = roe.goal*study.area, prop.hunt = m.hunting.r, select.hunt = hunt.aroem, fem.prop = (N.roe[3] * hunt.proef + N.roe[4] * hunt.aroef)/N.roe[7], hunting.type = Hunting.roe),
                         m.snow))
      # Vector of all the survival rates
      S.roe <- c(S.Jroef, S.Yroef, S.Proef, S.Aroef, S.Jroem, S.Yroem, S.Aroem)
      
      
      
      # --- Reproduction rates ---
      # # Reproductive individuals in each stage
      # N.Yf <- round(N.roe[2]*S.Yroef)
      # N.Pf <- round(N.roe[3]*S.Proef)
      # N.Af <- round(N.roe[4]*S.Aroef)
      # N.Ym <- round(N.roe[6]*S.Yroem)
      # N.Am <- round(N.roe[7]*S.Aroem)
      # Reprofemale <- N.Yf + N.Pf + N.Af # all reproductive females
      # Repromale <- N.Ym + N.Am # all reproductive males
      Reprofemale <- N.roe[2] + N.roe[3] + N.roe[4] # all reproductive females
      Repromale <- N.roe[6] + N.roe[7] # all reproductive males
      
      # - Reproduction function:
      F.Yroef <- S.Yroef*pr.r[1]*Repro(R.type, "F", k.r[1], h.r, alpha.r,
                                       Nb.female = Reprofemale,
                                       Nb.male = Repromale)

      F.Proef <- S.Proef*pr.r[2]*Repro(R.type, "F", k.r[2], h.r, alpha.r,
                                       Nb.female = Reprofemale,
                                       Nb.male = Repromale)

      F.Aroef <- S.Aroef*pr.r[2]*Repro(R.type, "F", k.r[2], h.r, alpha.r,
                                       Nb.female = Reprofemale,
                                       Nb.male = Repromale)

      prm.r <- ifelse(Reprofemale == 0, 0, (pr.r[1]*N.roe[2] + pr.r[2]*N.roe[3] + pr.r[2]*N.roe[4])/Reprofemale)
      km.r <- (k.r[1]*N.roe[2] + k.r[2]*N.roe[3] + k.r[2]*N.roe[4])/Reprofemale
      F.Yroem <- prm.r*Repro(R.type, "M", km.r, h.r, alpha.r,
                             Nb.female = Reprofemale,
                             Nb.male = Repromale)
      F.Aroem <- prm.r*Repro(R.type, "M", km.r, h.r, alpha.r,
                             Nb.female = Reprofemale,
                             Nb.male = Repromale)
      # Vector of all the fertility rate
      F.roe <- c(F.Yroef, F.Proef, F.Aroef, F.Yroem, F.Aroem)

      # Sex ratio at birth
      j.roe <- round(sum(F.roe * c(N.roe[2], N.roe[3], N.roe[4], N.roe[6], N.roe[7])))
      rho.r <- rbinom(1, max(j.roe, Ne.roe, na.rm = T), 0.5)/max(j.roe, Ne.roe, na.rm = T)
      
      

      # --- Dispersal rate ---
      # NDf.r <- 0 # Natal dispersal rate female
      # NDm.r <- 0 # Natal dispersal rate male
      ND.roe <- c(NDf.r, NDm.r)
      
      
      
      # ----- Update the matrix:
      mat.roe <- MPM(nameRoeDeer, rho.r, F.roe, ND.roe, S.roe)
      # --- Save stochastic matrices:
      mat.roe.list <- c(mat.roe.list, list(mat.roe))
      
      # ----- Run the MPM:
      Simulation.roe[, t, s] <- pmax(round(mat.roe %*% Simulation.roe[, t-1, s], digit = 0), 0)
      
      # Save population abundance
      # or stop running if population reaches quasi extinction:
      if (sum(Simulation.roe[, t, s], na.rm = T) > Ne.roe){
        N.roe[, 1] <- Simulation.roe[, t, s]
      } else {
        N.roe[, 1] <- rep(0, 7)
        break
      }
    }
    # - Name the scenario simulation
    scenario_name <- paste0("roe", nYears, "_", nb.sim, "-",
                            "fox_", Redfox.pred, "-",
                            "lynx_", Lynx.pred, "-",
                            "hunt_", Hunting.roe)
    # - Save matrices for each scenario
    assign(paste0("mat.list.", scenario_name), mat.roe.list)
    # --- Save simulation for each scenario
    saveRDS(Simulation.roe,
            file = paste0(if(isTRUE(Stop.hunt.ON.OFF)){"SH_"},
                          if(isTRUE(Snow.ON.OFF)){paste0("SNOW_", snow.frequency*100, "_", snow.magnitude*100, "-")},
                          # "Repro_", R.type, "-stoch_", stoch.distrib, "-",
                          scenario_name, ".rds"))
  }
}
# toc()


# --- Popbio analysis ---------------------------------------------------------

# --- From popbio package ------------------.
# https://cran.r-project.org/web/packages/popbio/popbio.pdf

# ! ################################## ! #
# !!! Stochastic analysis are needed !!! #
# ! ################################## ! #

# stable.stage(mat.roe) # Different at each time step with stochasticity...
# reproductive.value(mat.roe) # Relative contribution to the next generation 
# lambda(mat.roe) # Growth rate from the (last) matrix...
# 
# elasticity(mat.roe) # Same for elasticity and sensitivity
# sensitivity(mat.roe) #Default = F
# sens <- sensitivity(mat.roe, T)
# sens <- elasticity(mat.roe)

# stoch.sens(mat.list) # !!!!!
# mat.roe.list <- `mat.list.roe25_500-fox_fixed-lynx_selective-hunt_none`
length(mat.roe.list) # 12500 = 25 years * 500 runs
sens <- stoch.sens(mat.roe.list, 500)

## IMAGE plot with smaller boxes
image2(sens$elasticities, mar=c(1,3.5,5,1), box.offset=.1)
title("Elasticity matrix for scenario with red fox, lynx and hunting", line = 3)
      # no red fox predation, no lynx predation and no hunting", cex = 1, line=3)

image2(sens$sensitivities, mar=c(1,3.5,5,1), box.offset=.1)
title("Sensitivity matrix for scenario with red fox, lynx and hunting", line = 3)
# no red fox predation, no lynx predation and no hunting", cex = 1, line=3)

## MATPLOT
matplot2(sens$elasticities, log='y', type='b', yaxt='n', ltitle="Fate",
         ylab=expression(paste("Sensitivity of ",lambda)),
         main="Sensitivity matrix using matplot2")
pwrs <- -4:1
axis(2, 10^pwrs, parse(text=paste("10^", pwrs, sep = "")), las=1)

generation.time(mat.roe)


# ---------------------------.
### --- LYNX POPULATION -----
# ---------------------------.

rm(list = ls())
# Source the functions and parameters files
source("MPM_packages_and_functions.R")
source("MPM_load_parameters.R")

# --- Simulations, Lynx population -----

# --- Two loop simulation:
# The first one for the multiple simulations
# (i.e. start with slightly different values of initial states: Stochasticity)
# The second one for the time years simulation
# (i.e. how long the simulation last)


# Set number of simulations:
nb.sim <- 500
# Set time simulation (years):
nYears <- 25
# # Stochasticity in simulation True/False
# initial.stoch <- T # Initial density stochasticity
# init.var <- 0.1 # Percentage variation around the initial density
# envir.stoch <- T # Parameter stochasticity
# # Select the type of reproduction function
# R.type <- "FemaleFF" # "FemaleFF" or "ExpFF" or "HFF"
# Initialise list of stochastic matrices
mat.lyn.list <- list()
# Fix hunting scenario:
Hunting.type.lynx <- "goal"

# Initialise the simulation array:
Simulation.lyn <- array(rep(NA, length(N.lyn.init)*nYears*nb.sim),
                        c(length(N.lyn.init), nYears+1, nb.sim),
                        list(nameLynx, seq(0, nYears), paste('sim', seq(1, nb.sim))))


# --- Start simulation
# tic("simulation")
for (h in 1:length(Hunting.type.lynx)){
  # ----- Hunting type scenario
  Hunting.lynx <- Hunting.type.lynx[h]
  
# for (l in 1:3){
  l <- 0 # Lag in hunting decisions
  
  # Initialise the simulation array for each scenario:
  Simulation.lyn <- array(rep(NA, length(N.lyn.init)*nYears*nb.sim),
                        c(length(N.lyn.init), nYears+1, nb.sim),
                        list(nameLynx, seq(0, nYears), paste('sim', seq(1, nb.sim))))
  
  for (s in 1:(nb.sim)){ # For each simulation
    
    # Stochastic population structure and abundance
    N.lyn.init <- Stoch.init(lyn.density, nameLynx, study.area, init.var, initial.stoch)
    
    # Initialise population vector
    N.lyn <- N.lyn.init
    Simulation.lyn[, 1, s] <- N.lyn
    
    # run the simulation over the years...
    for (t in 2:(nYears+1)){ # For each time step
      
      # ----- Environmental stochasticity
      if (isTRUE(envir.stoch)){
        # Reload the parameter file with stochasticity
        source("MPM_stochastic_parameters.R")
      }
      
      # ----- Density dependence function:
      pr.l <- DDlynx((sum(N.lyn)/study.area)*100, pr.l)
      
      # ----- Vital rates at each time step
      # Changes related to population abundance
      
      # --- Survival rates ---
      # - Survival function from mortality rates:
      S.Jlynf <- Srate(c(m.jlynf,
                         HuntingHarvest("F", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.jlynf, hunting.type = Hunting.lynx)))
      S.Ylynf <- Srate(c(m.ylynf, #poaching.l,
                         HuntingHarvest("F", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.ylynf, hunting.type = Hunting.lynx)))
      S.Plynf <- Srate(c(m.plynf, #poaching.l,
                         HuntingHarvest("F", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.plynf, hunting.type = Hunting.lynx)))
      S.Alynf <- Srate(c(m.alynf, #poaching.l,
                         HuntingHarvest("F", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.alynf, hunting.type = Hunting.lynx)))
      
      S.Jlynm <- Srate(c(m.jlynm,
                         HuntingHarvest("M", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.jlynm, hunting.type = Hunting.lynx)))
      S.Ylynm <- Srate(c(m.ylynm, #poaching.l,
                         HuntingHarvest("M", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.ylynm, hunting.type = Hunting.lynx)))
      S.Alynm <- Srate(c(m.alynm,
                         HuntingHarvest("M", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.alynm, hunting.type = Hunting.lynx)))
      # Vector of all the survival rates
      S.lyn <- c(S.Jlynf, S.Ylynf, S.Plynf, S.Alynf, S.Jlynm, S.Ylynm, S.Alynm)
      
      
      
      # --- Reproduction rates ---
      # # Reproductive individuals in each stage
      # N.Yf <- round(N.lyn[2]*S.Ylynf)
      # N.Pf <- round(N.lyn[3]*S.Plynf)
      # N.Af <- round(N.lyn[4]*S.Alynf)
      # N.Ym <- round(N.lyn[6]*S.Ylynm)
      # N.Am <- round(N.lyn[7]*S.Alynm)
      # Reprofemale <- N.Yf + N.Pf + N.Af # all reproductive females
      # Repromale <- N.Ym + N.Am # all reproductive males
      Reprofemale <- N.lyn[2] + N.lyn[3] + N.lyn[4] # all reproductive females
      Repromale <- N.lyn[6] + N.lyn[7] # all reproductive males
      
      # - Reproduction function:
      F.Ylynf <- S.Ylynf*pr.l[1]*Repro(R.type, "F", k.l[1], h.l, alpha.l,
                                       Nb.female = Reprofemale,
                                       Nb.male = Repromale)

      F.Plynf <- S.Plynf*pr.l[2]*Repro(R.type, "F", k.l[2], h.l, alpha.l,
                                       Nb.female = Reprofemale,
                                       Nb.male = Repromale)

      F.Alynf <- S.Alynf*pr.l[2]*Repro(R.type, "F", k.l[2], h.l, alpha.l,
                                       Nb.female = Reprofemale,
                                       Nb.male = Repromale)

      prm.l <- ifelse(Reprofemale == 0, 0, (pr.l[1]*N.lyn[2] + pr.l[2]*N.lyn[3] + pr.l[2]*N.lyn[4])/Reprofemale)
      km.l <- (k.l[1]*N.lyn[2] + k.l[2]*N.lyn[3] + k.l[2]*N.lyn[4])/Reprofemale
      F.Ylynm <- prm.l*Repro(R.type, "M", km.l, h.l, alpha.l,
                             Nb.female = Reprofemale,
                             Nb.male = Repromale)
      F.Alynm <- prm.l*Repro(R.type, "M", km.l, h.l, alpha.l,
                             Nb.female = Reprofemale,
                             Nb.male = Repromale)
      # Vector of all the fertility rate
      F.lyn <- c(F.Ylynf, F.Plynf, F.Alynf, F.Ylynm, F.Alynm)
      
      # Sex ratio at birth
      j.lyn <- round(sum(F.lyn * c(N.lyn[2], N.lyn[3], N.lyn[4], N.lyn[6], N.lyn[7])))
      rho.l <- rbinom(1, max(j.lyn, 100, na.rm = T), 0.5)/max(j.lyn, 100, na.rm = T)
      
      

      # --- Dispersal rate ---
      # NDf.l <- 0 # Natal dispersal rate female
      # NDm.l <- 0 # Natal dispersal rate male
      ND.lyn <- c(NDf.l, NDm.l)
      
      
      
      # ----- Update the matrix:
      mat.lyn <- MPM(nameLynx, rho.l, F.lyn, ND.lyn, S.lyn)
      # --- Save stochastic matrices:
      mat.lyn.list <- c(mat.lyn.list, list(mat.lyn))
      
      # Run the MPM:
      Simulation.lyn[, t, s] <- pmax(round(mat.lyn %*% Simulation.lyn[, t-1, s], digit = 0), 0)
      
      # Save population abundance
      # or stop running if population reaches quasi extinction:
      if (sum(Simulation.lyn[, t, s], na.rm = T) > Ne.lyn){
        N.lyn[, 1] <- Simulation.lyn[, t, s]
      } else {
        N.lyn[, 1] <- rep(0, 7)
        break
      }
    }
    # Name the scenario simulation
    scenario_name <- paste0("lynx", nYears, "_", nb.sim, "-",
                            "hunt_", Hunting.lynx,
                            if (Hunting.lynx == "goal"){paste0(lyn.goal)},
                            if (l != 0){paste0("_lag_", l)}
                            )
    # Save matrices for each scenario
    assign(paste0("mat.list.", scenario_name), mat.lyn.list)
    # Save simulation for each scenario
    saveRDS(Simulation.lyn,
            file = paste0(scenario_name, ".rds"))
            # file = paste0("TEST-", scenario_name, ".rds"))
  }
}
# toc()
# Simulation.lyn <- readRDS("lynx50_500-hunt_none.rds")


# --- Popbio analysis ---------------------------------------------------------

stable.stage(mat.lyn) # Different for each scenario and initial density...
lambda(mat.lyn)

elasticity(mat.lyn) # Same for elasticity and sensitivity
sensitivity(mat.lyn) #Default = F
sens <- sensitivity(mat.lyn, T)
sens <- elasticity(mat.lyn)
# --- From popbio package ------------------.
## IMAGE plot with smaller boxes
image2(sens, mar = c(1, 3.5, 5, 1), box.offset = .1)
title("Sensitivity matrix using image2", line = 2.5)
## MATPLOT
matplot2(sens, log='y', type='b', yaxt='n', ltitle="Fate",
         ylab=expression(paste("Sensitivity of ",lambda)),
         main="Sensitivity matrix using matplot2")
pwrs <- -4:1
axis(2, 10^pwrs, parse(text=paste("10^", pwrs, sep = "")), las=1)


# -----------------------.
# --- PREDATOR PREY -----
# -----------------------.

rm(list = ls())
# Source the functions and parameters files
source("MPM_packages_and_functions.R")
source("MPM_load_parameters.R")

# Set number of simulations:
nb.sim <- 500
# Set time simulation (years):
nYears <- 50
# # Stochasticity in simulation True/False
# initial.stoch <- T # Initial density stochasticity
# init.var <- 0.1 # Percentage variation around the initial density
# envir.stoch <- T # Parameter stochasticity
# # Draw stochasticity from uniform or normal distribution for mortality rates
# stoch.distrib <- "norm" # "unif" or "norm"
# # Select the type of reproduction function
# R.type <- "FemaleFF" # "FemaleFF" or "ExpFF" or "HFF"

# Initialise scenario
Redfox.pred <- "fixed"
Lynx.pred <- "predator"
# Hunting.roe <- "selective"
# Hunting.lynx <- "none"
# Several scenario
hunt.scena <- expand.grid("hunt.roe" = Hunting.type.roe,
                     "hunt.lyn" = Hunting.type.lynx)
## - Select only some scenario...
# hunt.scena <- rbind(hunt.scena[2:3, ], hunt.scena[6:7, ], hunt.scena[10:11, ])
# hunt.scena <- hunt.scena[6:7, ]

# Initialise the simulation array:
# for roe deer and lynx populations
Simulation.roe <- array(rep(NA, length(N.roe.init)*nYears*nb.sim),
                        c(length(N.roe.init), nYears+1, nb.sim),
                        list(nameRoeDeer, seq(0, nYears), paste('sim', seq(1, nb.sim))))
Simulation.lyn <- array(rep(NA, length(N.lyn.init)*nYears*nb.sim),
                        c(length(N.lyn.init), nYears+1, nb.sim),
                        list(nameLynx, seq(0, nYears), paste('sim', seq(1, nb.sim))))

tic()
for (h in 1:nrow(hunt.scena)){ # For each combination of scenario

  # - Initialise scenarios:
  Hunting.roe <- hunt.scena[h, "hunt.roe"]
  Hunting.lynx <- hunt.scena[h, "hunt.lyn"]
  
  # - Initialise the simulation array for each scenario:
  Simulation.roe <- array(rep(NA, length(N.roe.init)*nYears*nb.sim),
                          c(length(N.roe.init), nYears+1, nb.sim),
                          list(nameRoeDeer, seq(0, nYears),
                               paste('sim', seq(1, nb.sim))))
  Simulation.lyn <- array(rep(NA, length(N.lyn.init)*nYears*nb.sim),
                          c(length(N.lyn.init), nYears+1, nb.sim),
                          list(nameLynx, seq(0, nYears),
                               paste('sim', seq(1, nb.sim))))
  # --- Start simulation
  for (s in 1:(nb.sim)){ # For each simulation
    
    # Stochastic population structure and abundance
    N.roe.init <- Stoch.init(roe.density, nameRoeDeer, study.area, init.var, initial.stoch)
    N.lyn.init <- Stoch.init(lyn.density, nameLynx, study.area, init.var, initial.stoch)
    
    # Initialise population vector (Stochastic)
    N.roe <- N.roe.init
    Simulation.roe[, 1, s] <- N.roe
    # ---
    N.lyn <- N.lyn.init
    Simulation.lyn[, 1, s] <- N.lyn
    
    # run the simulation over the years...
    for (t in 2:(nYears+1)){ # For each time step
      
      # ----- Environmental stochasticity
      if (isTRUE(envir.stoch)){
        # Reload the parameter file with stochasticity
        source("MPM_stochastic_parameters.R")
      }
      
      # ----- Scenarios (when change during simulations)
      # Redfox.pred <- Redfox.predation[rp]
      # Lynx.pred <- Lynx.predation.lr[lp]
      # Hunting.roe <- Hunting.type.roe[hr]
      # Hunting.lynx <- Hunting.type.lynx[hl]
      l <- 0 #lag in goal hunting
  
      
      # ----- Density dependence function:
      m.dens <- DDroe(sum(N.roe)/study.area)
      # ---
      pr.l <- DDlynx((sum(N.lyn)/study.area)*100, pr.l, sum(N.roe)/study.area)
  
      
      
      # ----- Vital rates at each time step
      # Changes related to population abundance
      
      # --- Survival rates ---
      # - Lynx kill rate:
      FS <- TypeIIKillRate(pred.lynFS, half.sat, sum(N.roe)/study.area, 30)
      FK <- TypeIIKillRate(pred.lynFK, half.sat, sum(N.roe)/study.area, 30)
      M <- TypeIIKillRate(pred.lynM, half.sat, sum(N.roe)/study.area, 30)
      pred.lyn <- c(FS, FK, M)
      lynx.kill <- c(N.lyn[2] + (1-pr.l[1])*N.lyn[3] + (1-pr.l[2])*N.lyn[4], pr.l[1]*N.lyn[3] + pr.l[2]*N.lyn[4], N.lyn[6] + N.lyn[7])
      
      # - Snow catastrophic event: (roe deer)
      # Snow event mortality rate: (Roe deer)
      m.snow <- SnowEvent(snow.rate = snow.magnitude, snow.freq = snow.frequency, Snow.ON.OFF)
      # --- Survival rates: Roe deer ---
      # - Survival function from mortality rates:
      S.Jroef <- Srate(c(m.jroef, m.dens,
                         RedFox(m.redfox, redfox.HL, Redfox.pred),
                         LynxPred(kill.rate = pred.lyn, months.pred = 3,
                                  lynx.hunter = lynx.kill, roe.tot = sum(N.roe),
                                  lynx.pred = Lynx.pred),
                         HuntingHarvest("F", pop = sum(N.roe), pop.goal = roe.goal*study.area, pop.threshold = sum(N.roe.init)/(100/roe.threshold), prop.hunt = m.hunting.r, select.hunt = hunt.jroef, hunting.type = Hunting.roe),
                         m.snow))
      S.Yroef <- Srate(c(m.yroef,
                         LynxPred(kill.rate = pred.lyn, months.pred = 12,
                                  lynx.hunter = lynx.kill, roe.tot = sum(N.roe),
                                  lynx.pred = Lynx.pred),
                         HuntingHarvest("F", pop = sum(N.roe), pop.goal = roe.goal*study.area, pop.threshold = sum(N.roe.init)/(100/roe.threshold), prop.hunt = m.hunting.r, select.hunt = hunt.yroef, hunting.type = Hunting.roe),
                         m.snow))
      S.Proef <- Srate(c(m.proef,
                         LynxPred(kill.rate = pred.lyn, months.pred = 12,
                                  lynx.hunter = lynx.kill, roe.tot = sum(N.roe),
                                  lynx.pred = Lynx.pred),
                         HuntingHarvest("F", pop = sum(N.roe), pop.goal = roe.goal*study.area, pop.threshold = sum(N.roe.init)/(100/roe.threshold), prop.hunt = m.hunting.r, select.hunt = hunt.proef, hunting.type = Hunting.roe),
                         m.snow))
      S.Aroef <- Srate(c(m.aroef,
                         LynxPred(kill.rate = pred.lyn, months.pred = 12,
                                  lynx.hunter = lynx.kill, roe.tot = sum(N.roe),
                                  lynx.pred = Lynx.pred),
                         HuntingHarvest("F", pop = sum(N.roe), pop.goal = roe.goal*study.area, pop.threshold = sum(N.roe.init)/(100/roe.threshold), prop.hunt = m.hunting.r, select.hunt = hunt.aroef, hunting.type = Hunting.roe),
                         m.snow))
  
      S.Jroem <- Srate(c(m.aroef, m.dens,
                         RedFox(m.redfox, redfox.HL, Redfox.pred),
                         LynxPred(kill.rate = pred.lyn, months.pred = 3,
                                  lynx.hunter = lynx.kill, roe.tot = sum(N.roe),
                                  lynx.pred = Lynx.pred),
                         HuntingHarvest("M", pop = sum(N.roe), pop.goal = roe.goal*study.area, pop.threshold = sum(N.roe.init)/(100/roe.threshold), prop.hunt = m.hunting.r, select.hunt = hunt.jroem, hunting.type = Hunting.roe),
                         m.snow))
      S.Yroem <- Srate(c(m.yroem,
                         LynxPred(kill.rate = pred.lyn, months.pred = 12,
                                  lynx.hunter = lynx.kill, roe.tot = sum(N.roe),
                                  lynx.pred = Lynx.pred),
                         HuntingHarvest("M", pop = sum(N.roe), pop.goal = roe.goal*study.area, pop.threshold = sum(N.roe.init)/(100/roe.threshold), prop.hunt = m.hunting.r, select.hunt = hunt.yroem, hunting.type = Hunting.roe),
                         m.snow))
      S.Aroem <- Srate(c(m.aroem,
                         LynxPred(kill.rate = pred.lyn, months.pred = 12,
                                  lynx.hunter = lynx.kill, roe.tot = sum(N.roe),
                                  lynx.pred = Lynx.pred),
                         HuntingHarvest("M", pop = sum(N.roe), pop.goal = roe.goal*study.area, pop.threshold = sum(N.roe.init)/(100/roe.threshold), prop.hunt = m.hunting.r, select.hunt = hunt.aroem, hunting.type = Hunting.roe),
                         m.snow))
      # Vector of all the survival rates
      S.roe <- c(S.Jroef, S.Yroef, S.Proef, S.Aroef, S.Jroem, S.Yroem, S.Aroem)
      # ---
      # --- Survival rates: Lynx ---
      # ---
      S.Jlynf <- Srate(c(m.jlynf,
                         HuntingHarvest("F", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.jlynf, hunting.type = Hunting.lynx)))
      S.Ylynf <- Srate(c(m.ylynf, #poaching.l,
                         HuntingHarvest("F", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.ylynf, hunting.type = Hunting.lynx)))
      S.Plynf <- Srate(c(m.plynf, #poaching.l,
                         HuntingHarvest("F", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.plynf, hunting.type = Hunting.lynx)))
      S.Alynf <- Srate(c(m.alynf, #poaching.l,
                         HuntingHarvest("F", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.alynf, hunting.type = Hunting.lynx)))
      
      S.Jlynm <- Srate(c(m.jlynm,
                         HuntingHarvest("M", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.jlynm, hunting.type = Hunting.lynx)))
      S.Ylynm <- Srate(c(m.ylynm, #poaching.l,
                         HuntingHarvest("M", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.ylynm, hunting.type = Hunting.lynx)))
      S.Alynm <- Srate(c(m.alynm,
                         HuntingHarvest("M", pop = sum(Simulation.lyn[, max(1, (t-1)-l), s]), pop.goal = (lyn.goal/100)*study.area, prop.hunt = m.hunting.l, select.hunt = hunt.alynm, hunting.type = Hunting.lynx)))
      # Vector of all the survival rates
      S.lyn <- c(S.Jlynf, S.Ylynf, S.Plynf, S.Alynf, S.Jlynm, S.Ylynm, S.Alynm)
      
      
      
      # --- Reproduction rates ---
      # - Roe deer
      # # Reproductive individuals in each stage
      Reprofemale.r <- N.roe[2] + N.roe[3] + N.roe[4] # all reproductive females
      Repromale.r <- N.roe[6] + N.roe[7] # all reproductive males
      
      # - Reproduction function:
      F.Yroef <- S.Yroef*pr.r[1]*Repro(R.type, "F", k.r[1], h.r, alpha.r,
                                       Nb.female = Reprofemale.r,
                                       Nb.male = Repromale.r)
  
      F.Proef <- S.Proef*pr.r[2]*Repro(R.type, "F", k.r[2], h.r, alpha.r,
                                       Nb.female = Reprofemale.r,
                                       Nb.male = Repromale.r)
  
      F.Aroef <- S.Aroef*pr.r[2]*Repro(R.type, "F", k.r[2], h.r, alpha.r,
                                       Nb.female = Reprofemale.r,
                                       Nb.male = Repromale.r)
  
      prm.r <- ifelse(Reprofemale.r == 0, 0, (pr.r[1]*N.roe[2] + pr.r[2]*N.roe[3] + pr.r[2]*N.roe[4])/Reprofemale.r)
      km.r <- (k.r[1]*N.roe[2] + k.r[2]*N.roe[3] + k.r[2]*N.roe[4])/Reprofemale.r
      F.Yroem <- prm.r*Repro(R.type, "M", km.r, h.r, alpha.r,
                             Nb.female = Reprofemale.r,
                             Nb.male = Repromale.r)
      F.Aroem <- prm.r*Repro(R.type, "M", km.r, h.r, alpha.r,
                             Nb.female = Reprofemale.r,
                             Nb.male = Repromale.r)
      # Vector of all the fertility rate
      F.roe <- c(F.Yroef, F.Proef, F.Aroef, F.Yroem, F.Aroem)
  
      # Sex ratio at birth
      j.roe <- round(sum(F.roe * c(N.roe[2], N.roe[3], N.roe[4], N.roe[6], N.roe[7])))
      rho.r <- rbinom(1, max(j.roe, Ne.roe, na.rm = T), 0.5)/max(j.roe, Ne.roe, na.rm = T)
      
      # ---
      # - Lynx
      # # Reproductive individuals in each stage
      Reprofemale.l <- N.lyn[2] + N.lyn[3] + N.lyn[4] # all reproductive females
      Repromale.l <- N.lyn[6] + N.lyn[7] # all reproductive males
      
      # - Reproduction function:
      F.Ylynf <- S.Ylynf*pr.l[1]*Repro(R.type, "F", k.l[1], h.l, alpha.l,
                                       Nb.female = Reprofemale.l,
                                       Nb.male = Repromale.l)
  
      F.Plynf <- S.Plynf*pr.l[2]*Repro(R.type, "F", k.l[2], h.l, alpha.l,
                                       Nb.female = Reprofemale.l,
                                       Nb.male = Repromale.l)
  
      F.Alynf <- S.Alynf*pr.l[2]*Repro(R.type, "F", k.l[2], h.l, alpha.l,
                                       Nb.female = Reprofemale.l,
                                       Nb.male = Repromale.l)
  
      prm.l <- ifelse(Reprofemale.l == 0, 0, (pr.l[1]*N.lyn[2] + pr.l[2]*N.lyn[3] + pr.l[2]*N.lyn[4])/Reprofemale.l)
      km.l <- (k.l[1]*N.lyn[2] + k.l[2]*N.lyn[3] + k.l[2]*N.lyn[4])/Reprofemale.l
      F.Ylynm <- prm.l*Repro(R.type, "M", km.l, h.l, alpha.l,
                             Nb.female = Reprofemale.l,
                             Nb.male = Repromale.l)
      F.Alynm <- prm.l*Repro(R.type, "M", km.l, h.l, alpha.l,
                             Nb.female = Reprofemale.l,
                             Nb.male = Repromale.l)
      # Vector of all the fertility rate
      F.lyn <- c(F.Ylynf, F.Plynf, F.Alynf, F.Ylynm, F.Alynm)
      
      # Sex ratio at birth
      j.lyn <- round(sum(F.lyn * c(N.lyn[2], N.lyn[3], N.lyn[4], N.lyn[6], N.lyn[7])))
      rho.l <- rbinom(1, max(j.lyn, 100, na.rm = T), 0.5)/max(j.lyn, 100, na.rm = T)
      
      
  
      # --- Dispersal rate ---
      ND.roe <- c(NDf.r, NDm.r) # Roe deer
      # ---
      ND.lyn <- c(NDf.l, NDm.l) # Lynx
      
      # Update the matrix:
      mat.roe <- MPM(nameRoeDeer, rho.r, F.roe, ND.roe, S.roe)
      mat.lyn <- MPM(nameLynx, rho.l, F.lyn, ND.lyn, S.lyn)
      
      # Run the MPM:
      Simulation.roe[, t, s] <- pmax(round(mat.roe %*% Simulation.roe[, t-1, s], digit = 0), 0)
      Simulation.lyn[, t, s] <- pmax(round(mat.lyn %*% Simulation.lyn[, t-1, s], digit = 0), 0)
      
      # Save population abundance:
      # or stop running if population reaches quasi extinction:
      if (sum(Simulation.roe[, t, s], na.rm = T) > Ne.roe){
        N.roe[, 1] <- Simulation.roe[, t, s]
      } else {
        N.roe[, 1] <- rep(0, 7)
        Simulation.roe[, t, s] <- rep(NA, 7)
      }
      # ---
      if (sum(Simulation.lyn[, t, s], na.rm = T) > Ne.lyn){
        N.lyn[, 1] <- Simulation.lyn[, t, s]
      } else {
        N.lyn[, 1] <- rep(0, 7)
        Simulation.lyn[, t, s] <- rep(NA, 7)
      }
    }
    # Name the scenario simulation
    scenario_name <- paste0(nYears, "_", nb.sim, "-",
                            "roe_", Hunting.roe,
                            if (Hunting.roe == "stopped" | Hunting.roe == "relaxed" ){paste0(roe.threshold)},
                            "-lynx_", Hunting.lynx,
                            if (Hunting.lynx == "goal"){paste0(lyn.goal)})
                            # if (l != 0){paste0("_lag_", l)})
    
    # Save simulation for each scenario
    saveRDS(Simulation.roe,
            file = paste0("pred_prey-roe", scenario_name, ".rds"))
    # ---
    saveRDS(Simulation.lyn,
            file = paste0("pred_prey-lyn", scenario_name, ".rds"))
  }
}
toc()
Simulation.roe
Simulation.lyn






# The end -----


