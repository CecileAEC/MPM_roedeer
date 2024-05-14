                    # ------------------------------- #
                    #           SIMULATIONS           #
                    #       Roe deer population       #
                    # ------------------------------- #



# -------------------------------.
### --- ROE DEER POPULATION -----
# -------------------------------.

rm(list = ls())
# Source the functions and parameters files
source("MPM_packages_and_functions.R")
source("MPM_load_parameters.R")

# --- Simulations, Roe deer population -----

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
    # - Save matrices for each scenario
    assign(paste0("mat.list.", scenario_name), mat.roe.list)
    # --- Save simulation for each scenario
    saveRDS(Simulation.roe,
            file = paste0(if(isTRUE(Stop.hunt.ON.OFF)){"SH_"},
                          if(isTRUE(Snow.ON.OFF)){paste0("SNOW_", snow.frequency*100, "_", snow.magnitude*100, "-")},
                          scenario_name, ".rds"))
  }
}
# toc()











