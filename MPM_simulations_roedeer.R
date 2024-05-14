                    # ------------------------------- #
                    #           SIMULATIONS           #
                    #       Roe deer population       #
                    # ------------------------------- #



# -------------------------------.
### --- ROE DEER POPULATION -----
# Single simulation study -------
# -------------------------------.

rm(list = ls())
# Source the functions and parameters files
source("MPM_packages_and_functions.R")
source("MPM_load_parameters.R")

# --- Simulations, roe deer population -----

# --- Three loop simulation:
# The first one for the different scenario
# (i.e. set up which mortality to account for in the scenario)
# The second one for the different simulations
# (i.e. to analyse different results with stochasticity)
# The third one for the time years simulation
# (i.e. how long the simulation last)


# Set number of simulations:
nb.sim <- 500
# Set time simulation (years):
nYears <- 25

# Initialise scenario
scena <- expand.grid("redfox" = Redfox.predation,
                     "lynx" = Lynx.predation,
                     "hunt" = Hunting.type.roe)

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

    # --- Save simulation for each scenario
    saveRDS(Simulation.roe,
            file = paste0(if(isTRUE(Stop.hunt.ON.OFF)){"SH_"},
                          if(isTRUE(Snow.ON.OFF)){paste0("SNOW_", snow.frequency*100, "_", snow.magnitude*100, "-")},
                          scenario_name, ".rds"))
  }
}
# toc()






# -------------------------------.
### --- ROE DEER POPULATION -----
# Theoretical gradient study ----
# -------------------------------.

rm(list = ls())
# Source the functions and parameters files
source("MPM_packages_and_functions.R")
source("MPM_load_parameters.R")

# --- Simulations, roe deer population -----

# --- Four loop simulation:
# The first one for the different scenario
# (i.e. set up which mortality to account for in the scenario)
# The second one for each combination of the gradient
# (i.e. yearling.grad and adult.grad)
# The third one for the different simulations
# (i.e. to analyse different results with stochasticity)
# The fourth one for the time years simulation
# (i.e. how long the simulation last)


# Set number of simulations:
nb.sim <- 500
# Set time simulation (years):
nYears <- 25

# Initialise gradient
yearling.grad <- seq(0, 0.03*sum(N.roe.init), by=150) # 0-0.5, by 1000 = 51
adult.grad <- seq(0, 100, by=5) # 0-100, by 2 = 51
grad <- expand.grid("yd" = yearling.grad, "ar" = adult.grad)

# Initialise scenario
scena <- expand.grid("redfox" = Redfox.predation,
                     "lynx" = Lynx.predation,
                     "hunt" = Hunting.type.roe)

# Initialise the simulation array:
Simulation.roe <- array(rep(NA, length(N.roe.init)*nYears*nb.sim),
                        c(length(N.roe.init), nYears+1, nb.sim),
                        list(nameRoeDeer, seq(0, nYears),
                             paste('sim', seq(1, nb.sim))))

library("tictoc")
# --- Start simulation
tic("scenario")
for (sc in 1:nrow(scena)){ # For each combination of scenario

  # - Initialise scenarios:
  Redfox.pred <- scena[sc, "redfox"]
  Lynx.pred <- scena[sc, "lynx"]
  Hunting.roe <- scena[sc, "hunt"]
  
  # - Initialise the simulation array for each scenario:
  Simulation.roe <- array(rep(NA, length(N.roe.init)*nYears*nb.sim),
                          c(length(N.roe.init), nYears+1, nb.sim),
                          list(nameRoeDeer, seq(0, nYears),
                               paste('sim', seq(1, nb.sim))))
  # --- Start simulation
  Roe.rate <- data.frame()
  tic("simulation")
  for (g in 1:nrow(grad)){ # For each step in the gradient
    
    # Each step in the gradient:
    yd <- grad[g, "yd"]
    ar <- grad[g, "ar"]
    
    # Initialise the simulation array for each gradient combination:
    Simulation.roe <- array(rep(NA, length(N.roe.init)*nYears*nb.sim),
                            c(length(N.roe.init), nYears+1, nb.sim),
                            list(nameRoeDeer, seq(0, nYears),
                                 paste('sim', seq(1, nb.sim))))
    
    for (s in 1:(nb.sim)){ # For each simulation
      
      # Stochastic population structure and abundance
      N.roe.init <- Stoch.init(roe.density, nameRoeDeer, study.area, init.var, initial.stoch)
      
      # Initialise population vector
      N.roe <- N.roe.init
      Simulation.roe[,1,s] <- N.roe
      
      # - Initialise survival rate at t0
      # No catastrophic snow event before simulation started
      m.snow <- 0
      
      # run the simulation over the years...
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
        Hunting.roe <- StopHunt(m.snow, Hunting.roe, Stop.hunt.ON.OFF)
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
                           LynxPred(roe.select = pred.aroef, lynx.pred = Lynx.pred)-LynxPred(roe.select = pred.aroef, lynx.pred = Lynx.pred)*(ar/100),
                           HuntingHarvest("F", pop = sum(N.roe), pop.goal = roe.goal*study.area, prop.hunt = m.hunting.r, select.hunt = hunt.aroef, hunting.type = Hunting.roe),
                           m.snow))
    
        S.Jroem <- Srate(c(m.aroef, m.dens,
                           RedFox(m.redfox, redfox.HL, Redfox.pred),
                           LynxPred(roe.select = pred.jroem, lynx.pred = Lynx.pred),
                           HuntingHarvest("M", pop = sum(N.roe), pop.goal = roe.goal*study.area, prop.hunt = m.hunting.r, select.hunt = hunt.jroem, hunting.type = Hunting.roe),
                           m.snow))
        S.Yroem <- Srate(c(m.yroem,
                           LynxPred(roe.select = pred.yroem, lynx.pred = Lynx.pred),
                           HuntingHarvest("M", pop = sum(N.roe), pop.goal = roe.goal*study.area, prop.hunt = m.hunting.r, select.hunt = hunt.yroem, hunting.type = Hunting.roe),
                           m.snow))
        S.Aroem <- Srate(c(m.aroem,
                           LynxPred(roe.select = pred.aroem, lynx.pred = Lynx.pred)-LynxPred(roe.select = pred.aroem, lynx.pred = Lynx.pred)*(ar/100),
                           HuntingHarvest("M", pop = sum(N.roe), pop.goal = roe.goal*study.area, prop.hunt = m.hunting.r, select.hunt = hunt.aroem, hunting.type = Hunting.roe),
                           m.snow))
        # Vector of all the survival rates
        S.roe <- c(S.Jroef, S.Yroef, S.Proef, S.Aroef, S.Jroem, S.Yroem, S.Aroem)
        
        
        
        # --- Reproduction rates ---
        # # Reproductive individuals in each stage
        Reprofemale <- N.roe[2] + N.roe[3] + N.roe[4] # all reproductive females
        Repromale <- N.roe[6] + N.roe[7] # all reproductive males
        
        # - Reproduction function:
        F.Yroef <- S.Yroef*pr.r[1]*Repro(R.type, "F", k.r[1], h.r, alpha,
                                         Nb.female = Reprofemale,
                                         Nb.male = Repromale)
  
        F.Proef <- S.Proef*pr.r[2]*Repro(R.type, "F", k.r[2], h.r, alpha,
                                         Nb.female = Reprofemale,
                                         Nb.male = Repromale)
  
        F.Aroef <- S.Aroef*pr.r[2]*Repro(R.type, "F", k.r[2], h.r, alpha,
                                         Nb.female = Reprofemale,
                                         Nb.male = Repromale)
  
        prm.r <- ifelse(Reprofemale == 0, 0, (pr.r[1]*N.roe[2] + pr.r[2]*N.roe[3] + pr.r[2]*N.roe[4])/Reprofemale)
        km.r <- (k.r[1]*N.roe[2] + k.r[2]*N.roe[3] + k.r[2]*N.roe[4])/Reprofemale
        F.Yroem <- prm.r*Repro(R.type, "M", km.r, h.r,
                               Nb.female = Reprofemale,
                               Nb.male = Repromale)
        F.Aroem <- prm.r*Repro(R.type, "M", km.r, h.r,
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
        
        # ----- Run the MPM:
        Simulation.roe[, t, s] <- pmax(round(mat.roe %*% Simulation.roe[, t-1, s], digit = 0), 0)
        # Additional dispersal...
        Simulation.roe[2, t, s] <- Simulation.roe[2, t, s] + yd/2 #Yearling female
        Simulation.roe[6, t, s] <- Simulation.roe[6, t, s] + yd/2 #Yearling male
        
        # Save population abundance
        # or stop running if population reaches quasi extinction:
        if (sum(Simulation.roe[, t, s], na.rm = T) > Ne.roe){
          N.roe[, 1] <- Simulation.roe[, t, s]
        } else {
          N.roe[, 1] <- rep(NA, 7)
          break
        }
        
      }
      
      Roe.state <- data.frame()
      for (i in 1:(nb.sim)){
        for (t in seq(1, 7, by = 2)){
          # State of the population at tend (25), t20, t15 and t10
          time <- c(nYears, 20, 15, 10)
          Roe.state[i, t] <- sum(Simulation.roe[, time[ceiling(t/2)], i])
    
          if (is.na(Roe.state[i, t]) | Roe.state[i, t] <= sum(N.roe.init)/20){
            Roe.state[i, t+1] <- "extinction"
          } else if (Roe.state[i, t] > 4*sum(N.roe.init)){ #or directly by saying K = 40 ind/km2
            Roe.state[i, t+1] <- "irruption"
          } else {
            Roe.state[i, t+1] <- "stable"
          }
        }
      }
      colnames(Roe.state) <- c("sumpop25", "state25",
                             "sumpop20", "state20",
                             "sumpop15", "state15",
                             "sumpop10", "state10")
    }
    
    print(paste("yd:", yd/1000, "ar:", ar, "| Done:", round((g/nrow(grad))*100), "%"))
    Roe.rate[g, 1] <- yd
    Roe.rate[g, 2] <- ar
    for (t in seq(1, 7, by = 2)){
      # --- Get extinction proportion at each time 25, 20, 15 and 10:
      if (any(Roe.state[, t+1] == "extinction")){
        Roe.rate[g, t*2+1] <- table(Roe.state[which(Roe.state[, t+1] == "extinction"), t+1])/nrow(Roe.state)
      } else {
        Roe.rate[g, t*2+1] <- 0
      }
      # --- Get stable proportion:
      if (any(Roe.state[, t+1] == "stable")){
        Roe.rate[g, t*2+2] <- table(Roe.state[which(Roe.state[, t+1] == "stable"), t+1])/nrow(Roe.state)
      } else {
        Roe.rate[g, t*2+2] <- 0
      }
      # --- Get irruption proportion:
      if (any(Roe.state[, t+1] == "irruption")){
        Roe.rate[g, t*2+3] <- table(Roe.state[which(Roe.state[, t+1] == "irruption"), t+1])/nrow(Roe.state)
      } else {
        Roe.rate[g, t*2+3] <- 0
      }
      # --- Sum up proportions:
      if (Roe.rate[g, t*2+2] == 1) { # If it is stable:
        Roe.rate[g, t*2+4] <- 0
      } else if (Roe.rate[g, t*2+3] != 0) { # If it is irruption:
        Roe.rate[g, t*2+4] <- Roe.rate[g, t*2+3]
      } else if (Roe.rate[g, t*2+1] != 0) { # If it is extinction:
        Roe.rate[g, t*2+4] <- - Roe.rate[g, t*2+1]
      }
    }
    colnames(Roe.rate) <- c("yd", "ar", "ext25", "stable25", "irrup25", "prop25",
                            "ext20", "stable20", "irrup20", "prop20",
                            "ext15", "stable15", "irrup15", "prop15",
                            "ext10", "stable10", "irrup10", "prop10")
    
    # --- Save roe deer state rate:
    saveRDS(Roe.rate,
            file = paste0("Roedeer-gradientrate-", nYears, "_", nb.sim, "-",
                          "fox_", Redfox.pred, "-",
                          "lynx_", Lynx.pred, "-",
                          "hunt_", Hunting.roe, "-",
                          "var", init.var*100,
                          "-", nrow(grad), ".rds"))
  
  }
  toc("simulation")
}
toc("scenario")

# --- Plots -----
# Load data:
Roerate <- Roe.rate

Roerate <- readRDS("Roedeer-gradientrate-25_500-fox_none-lynx_selective-hunt_none-var10-441.rds")
Roerate <- readRDS("Roedeer-gradientrate-25_500-fox_fixed-lynx_selective-hunt_none-var10-441.rds")
Roerate <- readRDS("Roedeer-gradientrate-25_500-fox_fixed-lynx_selective-hunt_selective-var10-441.rds")
# Select part of the gradient only
# Roerate <- Roe.rate[which(Roe.rate$yd >= 20000 & Roe.rate$ar >= 80),]


# Plot
ggplot(Roerate, aes(x=ar, y=(yd/sum(N.roe.init))*100, fill=prop25)) +
  geom_tile() +
  # scale_fill_gradient2(high="lightblue", mid="white", low="darkred",
  #                      midpoint=0, limits=range(-1, 1)) +
  # scale_fill_gradientn(colours = c("darkred","red4", "darkorange2", "orange1",
  #                                  "white",
  #                                  "powderblue", "cadetblue3", "cadetblue", "paleturquoise4"),
  scale_fill_gradientn(colours = c("darkred", "orange",
                                   "white"),
                                   #"lightblue", "deepskyblue4"),
                       breaks=c(-1, -0.75, -0.5, -0.25, 0), #0.25, 0.5, 0.75, 1),
                       labels=c("quasi-extinct", "75%", "50%", "25%","stable"), #"25%", "50%", "75%","irruption"),
                       limits=c(-1, 0)) +
  labs(title = "Roe deer population state proportion
- Red fox and lynx predation",
# - Lynx predation only", - Red fox and lynx predation", - Hunting, red fox and lynx predation",
       x = "Adult refuge
(% of refuge from lynx predation)",
       y = "Yearling immigration
(% of initial population density)",
       fill = "State
proportion") +
  theme(legend.key.height= unit(3, 'cm'),
        text = element_text(size = 25),
        axis.text = element_text(size = 25))


ggplot(Roerate, aes(x=ar, y=(yd/sum(N.roe.init))*100, fill=prop25)) +
  geom_tile() +
  # scale_fill_gradient2(high="lightblue", mid="white", low="darkred",
  #                      midpoint=0, limits=range(-1, 1)) +
  # scale_fill_gradientn(colours = c("darkred","red4", "darkorange2", "orange1",
  #                                  "white",
  #                                  "powderblue", "cadetblue3", "cadetblue", "paleturquoise4"),
  scale_fill_gradientn(colours = c("darkred", "firebrick3",
                                   "white",
                                   "lightblue", "deepskyblue4"),
                       breaks=c(-1, -0.75, -0.5, -0.25, 0),
                       labels=c("quasi-extinct", "75%", "50%", "25%","stable"),
                       limits=c(-1, 0)) +
  labs(title = "Roe deer population state proportion
- Hunting, red fox and lynx predation only",
       x = "Adult refuge
(% of refuge from lynx predation)",
       y = "Yearling immigration
(% of initial population density)",
       fill = "State
proportion") +
  ylim(0, 3) +
  theme(legend.key.height= unit(3, 'cm'),
        text = element_text(size = 25),
        axis.text = element_text(size = 25))

# plot100 <- ggplot(Roerate, aes(x=ar, y=(yd/sum(N.roe.init))*100, fill=prop25)) +
#   geom_tile() +
#   scale_fill_gradientn(colours = c("darkred", "orange",
#                                    "white",
#                                    "lightblue", "deepskyblue4"),
#                        breaks=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
#                        labels=c("quasi-extinct", "75%", "50%", "25%","stable", "25%", "50%", "75%","irruption"),
#                        limits=c(-1, 1)) +
#   labs(title = "Stability rate of roe deer population
# with gradient of yearling immigration dispersal and adult refugee",
#        x = "Adult refugees (% of additional survival)",
#        y = "Yearling immigration
#        (% of additional yearling individuals)",
#        fill = "Stability rate") +
#   theme(legend.key.height= unit(2, 'cm'),
#         text = element_text(size = 15),
#         axis.text = element_text(size = 15))
# 
# plot500 <- ggplot(Roe.rate, aes(x=ar, y=(yd/sum(N.roe.init))*100, fill=prop25)) +
#   geom_tile() +
#   scale_fill_gradientn(colours = c("darkred", "orange",
#                                    "white",
#                                    "lightblue", "deepskyblue4"),
#                        breaks=c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
#                        labels=c("quasi-extinct", "75%", "50%", "25%","stable", "25%", "50%", "75%","irruption"),
#                        limits=c(-1, 1)) +
#   labs(title = "Stability rate of roe deer population
# with gradient of yearling immigration dispersal and adult refugee",
#        x = "Adult refugees (% of additional survival)",
#        y = "Yearling immigration
#        (% of additional yearling individuals)",
#        fill = "Stability rate") +
#   theme(legend.key.height= unit(2, 'cm'),
#         text = element_text(size = 15),
#         axis.text = element_text(size = 15))
# 
# grid.arrange(plot100, plot500, ncol=2)







