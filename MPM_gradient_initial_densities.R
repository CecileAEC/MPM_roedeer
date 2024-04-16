# ------------------------------- #
#          - Chapter I -          #
#      GRADIENT  SIMULATIONS      #
#        Lynx and Roe deer        #
#           populations           #
# ------------------------------- #


rm(list=ls())
# Source the functions and parameters files
source("MPM_packages_and_functions.R")
source("MPM_load_parameters.R")

# --- Set simulation parameters -----

# Set number of simulations:
nb.sim <- 500
# Set time simulation (years):
nYears <- 25
# Stochasticity in simulation True/False
# initial.stoch <- T # Initial density stochasticity
# init.var <- 0.1 # Percentage variation around the initial density
# envir.stoch <- T # Parameter stochasticity
# # Draw stochasticity from uniform or normal distribution for mortality rates
# stoch.distrib <- "norm" # "unif" or "norm"
# # Select the type of reproduction function
# R.type <- "FemaleFF" # "FemaleFF" or "ExpFF" or "HFF"

# --- Initialise gradient and scenario -----

## Initialise gradient:
# roe.dens.grad <- seq(0.3, 30, by=0.99)
# lyn.dens.grad <- seq(0.3, 3, by=0.09)
## Quicker gradient check:
roe.dens.grad <- seq(0.3, 30, by=3.3)
lyn.dens.grad <- seq(0.3, 3, by=0.3)
grad.init <- expand.grid("roe" = roe.dens.grad, "lynx" = lyn.dens.grad)

## Initialise scenario
Redfox.pred <- "fixed"
Lynx.pred <- "predator"
Hunting.roe <- "selective"
Hunting.lynx <- "goal"
## Multiple scenario 
pred.prey.hunting <- expand.grid("roe" = c("selective", "relaxed", "stopped"), # + "none",
                                 "lynx" = Hunting.type.lynx)
pred.prey.hunting <- rbind(pred.prey.hunting[5:6, ], pred.prey.hunting[8:9, ])

## Initialise output dataframe
sim.density.roe <- data.frame()
sim.density.lyn <- data.frame()

# --- Start gradient simulation -----

tic()
for (ht in 1:nrow(pred.prey.hunting)){
  # Source the functions and parameters files
  source("MPM_packages_and_functions.R")
  source("MPM_load_parameters.R")
  
  # Initialise output dataframe
  sim.density.roe <- data.frame()
  sim.density.lyn <- data.frame()

  # scenario hunting lynx
  Hunting.roe <- pred.prey.hunting[ht, 1]
  Hunting.lynx <- pred.prey.hunting[ht, 2]
  for (g in 1:nrow(grad.init)){ # For each step in the initial density gradient
    # Initialise the simulation array:
    # for roe deer and lynx populations
    Simulation.roe <- array(rep(NA, length(N.roe.init)*nYears*nb.sim),
                            c(length(N.roe.init), nYears+1, nb.sim),
                            list(nameRoeDeer, seq(0, nYears), paste('sim', seq(1, nb.sim))))
    Simulation.lyn <- array(rep(NA, length(N.lyn.init)*nYears*nb.sim),
                            c(length(N.lyn.init), nYears+1, nb.sim),
                            list(nameLynx, seq(0, nYears), paste('sim', seq(1, nb.sim))))
  
    # Where are you in the loop?
    print(paste("roe density:", grad.init[g,"roe"], "| lynx density:", grad.init[g,"lynx"], "| Done:", round((g/nrow(grad.init))*100), "%"))
    
    # --- Start simulation
    for (s in 1:(nb.sim)){ # For each simulation
    
      # Stochastic population structure and abundance
      N.roe.init <- Stoch.init(grad.init[g, "roe"], nameRoeDeer, study.area, 0.1, initial.stoch)
      N.lyn.init <- Stoch.init(grad.init[g, "lynx"], nameLynx, study.area, 0.1, initial.stoch)
    
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
        }
        # ---
        if (sum(Simulation.lyn[, t, s], na.rm = T) > Ne.lyn){
          N.lyn[, 1] <- Simulation.lyn[, t, s]
        } else {
          N.lyn[, 1] <- rep(0, 7)
        }
      }
    }
    # Name the scenario simulation
    scenario_name <- paste0(nYears, "_", nb.sim, "-",
                            "roe_", Hunting.roe,
                            if (Hunting.roe == "stopped" | Hunting.roe == "relaxed"){paste0(roe.threshold)},
                            "-",
                            "lynx_", Hunting.lynx,
                            if (Hunting.lynx == "goal"){paste0(lyn.goal)})
                            # if (l != 0){paste0("_lag_", l)})
    
    # Save simulation for each scenario
    # saveRDS(Simulation.roe,
    #         file = paste0("GRAD-pred_prey-roe", scenario_name,
    #                       "-grad_", grad.init[g,"roe"], "_", grad.init[g,"lynx"], ".rds"))
    # ---
    # saveRDS(Simulation.lyn,
    #         file = paste0("GRAD-pred_prey-lyn", scenario_name,
    #                       "-grad_", grad.init[g,"roe"], "_", grad.init[g,"lynx"], ".rds"))
    
    # Make new outputs!!!
    Roe.sumpop <- PopSum(Simulation.roe)
    Roe.quant <- PopQuantile(Roe.sumpop)
    
    sim.density.roe[g, 1:2] <- grad.init[g, ]
    sim.density.roe[g, 3:4] <- Roe.quant[which(Roe.quant$time == 25), 1:2]
    sim.density.roe[g, 5:12] <- Roe.quant[which(Roe.quant$time == 25), 3:10]/study.area
    sim.density.roe[is.na(sim.density.roe)] <- 0
    # Average growth rate over the 25 years simulation
    for (i in (1:nrow(sim.density.roe))){
      sim.density.roe[i, 13] <- sim.density.roe[i, 12]/sim.density.roe[i, 1]
    }
    colnames(sim.density.roe)[13] <- "lambda"
    # --> All lambda? How much variance?
    # ---
    Lyn.sumpop <- PopSum(Simulation.lyn)
    Lyn.quant <- PopQuantile(Lyn.sumpop)
    
    sim.density.lyn[g, 1:2] <- grad.init[g, ]
    sim.density.lyn[g, 3:4] <- Lyn.quant[which(Lyn.quant$time == 25), 1:2]
    sim.density.lyn[g, 5:12] <- (Lyn.quant[which(Lyn.quant$time == 25), 3:10]*100)/study.area
    sim.density.lyn[is.na(sim.density.lyn)] <- 0
    # Average growth rate over the 25 years simulation
    for (i in (1:nrow(sim.density.lyn))){
      sim.density.lyn[i, 13] <- sim.density.lyn[i, 12]/sim.density.lyn[i, 2]
    }
    colnames(sim.density.lyn)[13] <- "lambda"
  
    # Save output
    saveRDS(sim.density.roe,
            file = paste0("Roe-pred_prey", scenario_name,
                          "-grad_", length(roe.dens.grad), "x", length(lyn.dens.grad),
                          ".rds"))
    
    saveRDS(sim.density.lyn,
            file = paste0("Lyn-pred_prey", scenario_name,
                          "-grad_", length(roe.dens.grad), "x", length(lyn.dens.grad),
                          ".rds"))
  }
}
toc()

# Roe <- readRDS("Roe-pred_prey25_500-roe_none-lynx_none-grad_31x31.rds")
# # Lynx <- readRDS("Lyn-pred_prey25_500-roe_none-lynx_none-grad_31x31.rds")
# Roe <- readRDS("Roe-Pred_prey25_500-roe_selective-lynx_none-grad_31x31.rds")
# # Lynx <- readRDS("Lyn-Pred_prey25_500-roe_selective-lynx_none-grad_31x31.rds")
# Roe <- readRDS("Roe-Pred_prey25_500-roe_selective-lynx_selective-grad_31x31.rds")
# Lynx <- readRDS("Lyn-Pred_prey25_500-roe_selective-lynx_selective-grad_31x31.rds")

Roe.data <- sim.density.roe
Lynx.data <- sim.density.lyn

# No hunting on both
Roe.nohunt <- readRDS("Roe-pred_prey25_500-roe_none-lynx_none-grad_10x10.rds")
Lynx.nohunt <- readRDS("Lyn-pred_prey25_500-roe_none-lynx_none-grad_10x10.rds")

# Hunting roe deer only
Roe.propnone <- readRDS("Roe-pred_prey25_500-roe_selective-lynx_none-grad_10x10.rds")
Lynx.propnone <- readRDS("Lyn-pred_prey25_500-roe_selective-lynx_none-grad_10x10.rds")
Roe.SH20none <- readRDS("Roe-pred_prey25_500-roe_stopped20-lynx_none-grad_10x10.rds") # 20%
Lynx.SH20none <- readRDS("Lyn-pred_prey25_500-roe_stopped20-lynx_none-grad_10x10.rds") # 20%
Roe.RH20none <- readRDS("Roe-pred_prey25_500-roe_relaxed20-lynx_none-grad_10x10.rds") # 20%
Lynx.RH20none <- readRDS("Lyn-pred_prey25_500-roe_relaxed20-lynx_none-grad_10x10.rds") # 20%
Roe.SH50none <- readRDS("Roe-pred_prey25_500-roe_stopped50-lynx_none-grad_10x10.rds") # 50%
Lynx.SH50none <- readRDS("Lyn-pred_prey25_500-roe_stopped50-lynx_none-grad_10x10.rds") # 50%
Roe.RH50none <- readRDS("Roe-pred_prey25_500-roe_relaxed50-lynx_none-grad_10x10.rds") # 50%
Lynx.RH50none <- readRDS("Lyn-pred_prey25_500-roe_relaxed50-lynx_none-grad_10x10.rds") # 50%

# Hunting on both
Roe.data <- readRDS("Roe-pred_prey25_500-roe_selective-lynx_selective-grad_10x10.rds")
Lynx.data <- readRDS("Lyn-pred_prey25_500-roe_selective-lynx_selective-grad_10x10.rds")
Roe.data <- readRDS("Roe-pred_prey25_500-roe_stopped20-lynx_selective-grad_10x10.rds") # 20%
Lynx.data <- readRDS("Lyn-pred_prey25_500-roe_stopped20-lynx_selective-grad_10x10.rds") # 20%
Roe.data <- readRDS("Roe-pred_prey25_500-roe_relaxed20-lynx_selective-grad_10x10.rds") # 20%
Lynx.data <- readRDS("Lyn-pred_prey25_500-roe_relaxed20-lynx_selective-grad_10x10.rds") # 20%
Roe.data <- readRDS("Roe-pred_prey25_500-roe_stopped50-lynx_selective-grad_10x10.rds") # 50%
Lynx.data <- readRDS("Lyn-pred_prey25_500-roe_stopped50-lynx_selective-grad_10x10.rds") # 50%
Roe.data <- readRDS("Roe-pred_prey25_500-roe_relaxed50-lynx_selective-grad_10x10.rds") # 50%
Lynx.data <- readRDS("Lyn-pred_prey25_500-roe_relaxed50-lynx_selective-grad_10x10.rds") # 50%

Roe.data <- readRDS("Roe-pred_prey25_500-roe_selective-lynx_goal1-grad_10x10.rds")
Lynx.data <- readRDS("Lyn-pred_prey25_500-roe_selective-lynx_goal1-grad_10x10.rds")
Roe.data <- readRDS("Roe-pred_prey25_500-roe_stopped20-lynx_goal1-grad_10x10.rds") # 20% - 1
Lynx.data <- readRDS("Lyn-pred_prey25_500-roe_stopped20-lynx_goal1-grad_10x10.rds") # 20% - 1
Roe.data <- readRDS("Roe-pred_prey25_500-roe_relaxed20-lynx_goal1-grad_10x10.rds") # 20% - 1
Lynx.data <- readRDS("Lyn-pred_prey25_500-roe_relaxed20-lynx_goal1-grad_10x10.rds") # 20% - 1

Roe.data <- readRDS("Roe-pred_prey25_500-roe_selective-lynx_goal0.5-grad_10x10.rds")
Lynx.data <- readRDS("Lyn-pred_prey25_500-roe_selective-lynx_goal0.5-grad_10x10.rds")
Roe.data <- readRDS("Roe-pred_prey25_500-roe_stopped50-lynx_goal0.5-grad_10x10.rds") # 50% - 0.5
Lynx.data <- readRDS("Lyn-pred_prey25_500-roe_stopped50-lynx_goal0.5-grad_10x10.rds") # 50% - 0.5
Roe.data <- readRDS("Roe-pred_prey25_500-roe_relaxed50-lynx_goal0.5-grad_10x10.rds") # 50% - 0.5
Lynx.data <- readRDS("Lyn-pred_prey25_500-roe_relaxed50-lynx_goal0.5-grad_10x10.rds") # 50% - 0.5

Roe.data <- list()
Lynx.data <- list()

# --- Plots -----

RoePop <- GradientPlotPopDensity(Roe.data, "Roe deer")
RoeGR <- GradientPlotGR(Roe.data, "Roe deer")

LynxPop <- GradientPlotPopDensity(Lynx.data, "Lynx")
LynxGR <- GradientPlotGR(Lynx.data, "Lynx")

RoePop + LynxPop + RoeGR + LynxGR +
  plot_annotation(title = "Predator prey simulation with hunting both: Relaxed hunting roe deer + Lynx hunting goal 0.5", theme = theme(plot.title = element_text(hjust = 0.5, size = 20)))


# Roe deer
roepop <- ggplot(Roe, aes(x=roe, y=lynx, fill=mean)) +
  geom_tile() +
  # scale_fill_gradient2(high="lightblue", mid="white", low="darkred",
  #                      midpoint=0, limits=range(-1, 1)) +
  scale_fill_gradientn(colours = c("gainsboro", "darkolivegreen", "black"),
                       limits = c(0, 30)) + #max(Roe$mean)
                         # c("darkred","red4", "darkorange2", "orange1",
                         #           "white",
                         #           "powderblue", "cadetblue3", "cadetblue", "paleturquoise4"),
                       
  # xlim(10, 30) +
  # ylim() +
  # stat_contour(color = "white", linewidth = .7, bins = 5) +
  labs(title = "Roe deer population density", #after 25 years simulation
       x = "Roe deer initial density (N/km2)",
       y = "Lynx initial density (N/100km2)",
       fill = "Population\ndensity") +
  theme(legend.key.height= unit(1.5, 'cm'),
        text = element_text(size = 15),
        axis.text = element_text(size = 15))

roegrow <- ggplot(Roe, aes(x=roe, y=lynx, fill=lambda)) +
  geom_tile() +
  # scale_fill_gradient2(high="deepskyblue4", mid="white", low="darkred",
  #                      midpoint=1, limits=c(0, 6)) + #max(Roe$lambda)
  scale_fill_gradientn(colours = c("darkred", #"red4", "darkorange2", "orange1",
                                   "white",
                                   "powderblue", "cadetblue3", "cadetblue", "paleturquoise4", "darkslategray"),
                       limits = c(0, 6)) + #max(Roe$lambda)
  # xlim(10, 30) +
  # ylim() +
  labs(title = "Roe deer growth rate", #for the 25 years simulation
       x = "Roe deer initial density (N/km2)",
       y = "Lynx initial density (N/100km2)",
       fill = "Lambda") +
  theme(legend.key.height= unit(1.5, 'cm'),
        text = element_text(size = 15),
        axis.text = element_text(size = 15))

# Lynx
lynxpop <- ggplot(Lynx, aes(x=roe, y=lynx, fill=mean)) +
  geom_tile() +
  # scale_fill_gradient2(high="lightblue", mid="white", low="darkred",
  #                      midpoint=0, limits=range(-1, 1)) +
  scale_fill_gradientn(colours = c("gainsboro", "darkolivegreen", "black"),
                       limits = c(0, 3)) + #max(Lynx$mean)
                       # c("darkred","red4", "darkorange2", "orange1",
                       #             "white",
                       #             "powderblue", "cadetblue3", "cadetblue", "paleturquoise4"),
                       
  # xlim(10, 30) +
  # ylim() +
  labs(title = "Lynx population density", #after 25 years simulation
       x = "Roe deer initial density (N/km2)",
       y = "Lynx initial density (N/100km2)",
       fill = "Population\ndensity") +
  theme(legend.key.height= unit(1.5, 'cm'),
        text = element_text(size = 15),
        axis.text = element_text(size = 15))

lynxgrow <- ggplot(Lynx, aes(x=roe, y=lynx, fill=lambda)) +
  geom_tile() +
  # scale_fill_gradient2(high="deepskyblue4", mid="white", low="darkred",
  #                      midpoint=1, limits=c(0, 6)) + #max(Lynx$lambda)
  scale_fill_gradientn(colours = c("darkred", #"red4", "darkorange2", "orange1",
                                   "white",
                                   "powderblue", "cadetblue3", "cadetblue", "paleturquoise4", "darkslategray"),
                       limits = c(0, 6)) + #max(Lynx$lambda)
  # xlim(10, 30) +
  # ylim() +
  labs(title = "Lynx growth rate", #for the 25 years simulation
       x = "Roe deer initial density (N/km2)",
       y = "Lynx initial density (N/100km2)",
       fill = "Lambda") +
  theme(legend.key.height= unit(1.5, 'cm'),
        text = element_text(size = 15),
        axis.text = element_text(size = 15))



roepop + lynxpop + roegrow + lynxgrow +
  plot_annotation(title = "Predator prey simulation with goal hunting 1", theme = theme(plot.title = element_text(hjust = 0.5)))

                                                        