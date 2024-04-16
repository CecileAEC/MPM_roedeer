                    # ------------------------------- #
                    #          - Chapter I -          #
                    #         Simulation plots        #
                    #        Lynx and Roe deer        #
                    #           populations           #
                    # ------------------------------- #


# -
# ----- Data preparation functions -----
# -
source("MPM_packages_and_functions.R")

# -
# ----- Load default parameters -----
# -
source("MPM_load_parameters.R") # Get the default parameters!
# from this one you just need:
# study.area <- 10000 # Study area in square kilometre
# roe.density <- 10 # Nb of roe deer per square kilometre
# nameRoeDeer <- c("Jroe.F", # Juvenile roe deer female
#                  "Yroe.F", # Yearling roe deer female
#                  "Proe.F", # Two year-old roe deer female
#                  "Aroe.F", # Adult roe deer female
#                  "Jroe.M", # Juvenile roe deer male
#                  "Yroe.M", # Yearling roe deer male
#                  "Aroe.M") # Adult roe deer male
# 
# N.roe.init <- VectorPop(round(c(0.237, 0.105, 0.045, 0.147, 0.237, 0.082, 0.147)
#                               * (roe.density*study.area)), nameRoeDeer)

# -
# ----- Roe deer - DATA -----
# -
# Explore different scenarios
Roe.heaven.normFe <- readRDS("Repro_FemaleFF-stoch_norm-roe25_500-fox_none-lynx_none-hunt_none.rds")
Roe.fox.normFe <- readRDS("Repro_FemaleFF-stoch_norm-roe25_500-fox_fixed-lynx_none-hunt_none.rds")
Roe.lyn.normFe <- readRDS("Repro_FemaleFF-stoch_norm-roe25_500-fox_none-lynx_selective-hunt_none.rds")
Roe.pred.normFe <- readRDS("Repro_FemaleFF-stoch_norm-roe25_500-fox_fixed-lynx_selective-hunt_none.rds")
Roe.hell.normFe <- readRDS("Repro_FemaleFF-stoch_norm-roe25_500-fox_fixed-lynx_selective-hunt_selective.rds")
NormFe.list <- list(Roe.heaven.normFe, Roe.fox.normFe, Roe.lyn.normFe, Roe.pred.normFe, Roe.hell.normFe)

# Explore Snow event
Sim.roe.snow.heaven <- readRDS("SNOW_8_30-Repro_FemaleFF-stoch_norm-roe25_500-fox_none-lynx_none-hunt_none.rds")
Sim.roe.snow.fox <- readRDS("SNOW_8_30-Repro_FemaleFF-stoch_norm-roe25_500-fox_fixed-lynx_none-hunt_none.rds")
Sim.roe.snow.lynx <- readRDS("SNOW_8_30-Repro_FemaleFF-stoch_norm-roe25_500-fox_none-lynx_selective-hunt_none.rds")
Sim.roe.snow.pred <- readRDS("SNOW_8_30-Repro_FemaleFF-stoch_norm-roe25_500-fox_fixed-lynx_selective-hunt_none.rds")
Sim.roe.snow.hell <- readRDS("SNOW_8_30-Repro_FemaleFF-stoch_norm-roe25_500-fox_fixed-lynx_selective-hunt_selective.rds")
Roe.snow08.list <- list(Sim.roe.snow.heaven, Sim.roe.snow.fox, Sim.roe.snow.lynx, Sim.roe.snow.pred, Sim.roe.snow.hell)

# Explore snow event with stop hunting the next year
Simulation.roe <- readRDS("SH_SNOW_8_30_Repro_FemaleFF-stoch_norm-roe25_500-fox_fixed-lynx_selective-hunt_selective.rds")
# Harvesting male only
Simulation.roe <- readRDS("Repro_FemaleFF-stoch_norm-roe25_500-fox_fixed-lynx_selective-hunt_biased-var10.rds")
# Simulation.roe <- readRDS("Repro_FemaleFF-stoch_norm-roe25_500-fox_fixed-lynx_selective-hunt_biased.M13.8percent.rds")



# -
# ----- Lynx - DATA -----
# -
Sim.lyn.goal01 <- readRDS("lynx25_500-hunt_goal0.1.rds")
Sim.lyn.goal05 <- readRDS("lynx25_500-hunt_goal0.5.rds")
Sim.lyn.goal1 <- readRDS("lynx25_500-hunt_goal1.rds")
Sim.lyn.select <- readRDS("lynx25_500-hunt_selective.rds")
Sim.lyn.none <- readRDS("lynx25_500-hunt_none.rds")
Lynx.list <- list(Sim.lyn.none, Sim.lyn.select, Sim.lyn.goal01, Sim.lyn.goal05, Sim.lyn.goal1)

# Simulation.lyn0 <- readRDS("lynx25_500-hunt_goal1-var10.rds")
# Simulation.lyn1 <- readRDS("lynx25_500-hunt_goal1_lag_1-var10.rds")
# Simulation.lyn2 <- readRDS("lynx25_500-hunt_goal1_lag_2-var10.rds")
# Simulation.lyn3 <- readRDS("lynx25_500-hunt_goal1_lag_3-var10.rds")

# -
# ----- Predator-prey - DATA -----
# -
Predprey.roe.nonenone <- readRDS("pred_prey-roe50_500-roe_none-lynx_none.rds")

Predprey.roe.huntnone <- readRDS("pred_prey-roe50_500-roe_selective-lynx_none.rds")
Predprey.roe.huntSH20none <- readRDS("pred_prey-roe50_500-roe_stopped20-lynx_none.rds")
Predprey.roe.huntSH50none <- readRDS("pred_prey-roe50_500-roe_stopped50-lynx_none.rds")
Predprey.roe.huntRH20none <- readRDS("pred_prey-roe50_500-roe_relaxed20-lynx_none.rds")
Predprey.roe.huntRH50none <- readRDS("pred_prey-roe50_500-roe_relaxed50-lynx_none.rds")

Predprey.roe.hunthunt <- readRDS("pred_prey-roe50_500-roe_selective-lynx_selective.rds")
Predprey.roe.huntSH20hunt <- readRDS("pred_prey-roe50_500-roe_stopped20-lynx_selective.rds")
Predprey.roe.huntRH20hunt <- readRDS("pred_prey-roe50_500-roe_relaxed20-lynx_selective.rds") # ???
Predprey.roe.huntSH50hunt <- readRDS("pred_prey-roe50_500-roe_stopped50-lynx_selective.rds")
Predprey.roe.huntRH50hunt <- readRDS("pred_prey-roe50_500-roe_relaxed50-lynx_selective.rds") # ???

Predprey.roe.huntgoal1 <- readRDS("pred_prey-roe50_500-roe_selective-lynx_goal1.rds")
Predprey.roe.huntgoal05 <- readRDS("pred_prey-roe50_500-roe_selective-lynx_goal0.5.rds")

Predprey.roe.huntSH20goal1 <- readRDS("pred_prey-roe50_500-roe_stopped20-lynx_goal1.rds")
Predprey.roe.huntSH50goal1 <- readRDS("pred_prey-roe50_500-roe_stopped50-lynx_goal1.rds")
Predprey.roe.huntSH20goal05 <- readRDS("pred_prey-roe50_500-roe_stopped20-lynx_goal0.5.rds")
Predprey.roe.huntSH50goal05 <- readRDS("pred_prey-roe50_500-roe_stopped50-lynx_goal0.5.rds")

Predprey.roe.huntRH20goal1 <- readRDS("pred_prey-roe50_500-roe_relaxed20-lynx_goal1.rds")
Predprey.roe.huntRH50goal1 <- readRDS("pred_prey-roe50_500-roe_relaxed50-lynx_goal1.rds")
Predprey.roe.huntRH20goal05 <- readRDS("pred_prey-roe50_500-roe_relaxed20-lynx_goal0.5.rds")
Predprey.roe.huntRH50goal05 <- readRDS("pred_prey-roe50_500-roe_relaxed50-lynx_goal0.5.rds")

roe.predprey <- list(Predprey.roe.nonenone, #1
                     Predprey.roe.huntnone, #2
                     Predprey.roe.huntSH20none, #3
                     Predprey.roe.huntSH50none, #4
                     Predprey.roe.huntRH20none, #5
                     Predprey.roe.huntRH50none, #6
                     Predprey.roe.hunthunt, #7
                     Predprey.roe.huntSH20hunt, #8
                     Predprey.roe.huntRH20hunt, #9
                     Predprey.roe.huntSH50hunt, #10
                     Predprey.roe.huntRH50hunt, #11
                     Predprey.roe.huntgoal1, #12
                     Predprey.roe.huntgoal05, #13
                     Predprey.roe.huntSH20goal1, #14
                     Predprey.roe.huntSH50goal1, #15
                     Predprey.roe.huntSH20goal05, #16
                     Predprey.roe.huntSH50goal05, #17
                     Predprey.roe.huntRH20goal1, #18
                     Predprey.roe.huntRH50goal1, #19
                     Predprey.roe.huntRH20goal05, #20
                     Predprey.roe.huntRH50goal05) #21

Predprey.lyn.nonenone <- readRDS("pred_prey-lyn50_500-roe_none-lynx_none.rds")

Predprey.lyn.huntnone <- readRDS("pred_prey-lyn50_500-roe_selective-lynx_none.rds")
Predprey.lyn.huntSH20none <- readRDS("pred_prey-lyn50_500-roe_stopped20-lynx_none.rds")
Predprey.lyn.huntSH50none <- readRDS("pred_prey-lyn50_500-roe_stopped50-lynx_none.rds")
Predprey.lyn.huntRH20none <- readRDS("pred_prey-lyn50_500-roe_relaxed20-lynx_none.rds")
Predprey.lyn.huntRH50none <- readRDS("pred_prey-lyn50_500-roe_relaxed50-lynx_none.rds")

Predprey.lyn.hunthunt <- readRDS("pred_prey-lyn50_500-roe_selective-lynx_selective.rds")
Predprey.lyn.huntSH20hunt <- readRDS("pred_prey-lyn50_500-roe_stopped20-lynx_selective.rds")
Predprey.lyn.huntRH20hunt <- readRDS("pred_prey-lyn50_500-roe_relaxed20-lynx_selective.rds") # ???
Predprey.lyn.huntSH50hunt <- readRDS("pred_prey-lyn50_500-roe_stopped50-lynx_selective.rds")
Predprey.lyn.huntRH50hunt <- readRDS("pred_prey-lyn50_500-roe_relaxed50-lynx_selective.rds") # ???

Predprey.lyn.huntgoal1 <- readRDS("pred_prey-lyn50_500-roe_selective-lynx_goal1.rds")
Predprey.lyn.huntgoal05 <- readRDS("pred_prey-lyn50_500-roe_selective-lynx_goal0.5.rds")

Predprey.lyn.huntSH20goal1 <- readRDS("pred_prey-lyn50_500-roe_stopped20-lynx_goal1.rds")
Predprey.lyn.huntSH50goal1 <- readRDS("pred_prey-lyn50_500-roe_stopped50-lynx_goal1.rds")
Predprey.lyn.huntSH20goal05 <- readRDS("pred_prey-lyn50_500-roe_stopped20-lynx_goal0.5.rds")
Predprey.lyn.huntSH50goal05 <- readRDS("pred_prey-lyn50_500-roe_stopped50-lynx_goal0.5.rds")

Predprey.lyn.huntRH20goal1 <- readRDS("pred_prey-lyn50_500-roe_relaxed20-lynx_goal1.rds")
Predprey.lyn.huntRH50goal1 <- readRDS("pred_prey-lyn50_500-roe_relaxed50-lynx_goal1.rds")
Predprey.lyn.huntRH20goal05 <- readRDS("pred_prey-lyn50_500-roe_relaxed20-lynx_goal0.5.rds")
Predprey.lyn.huntRH50goal05 <- readRDS("pred_prey-lyn50_500-roe_relaxed50-lynx_goal0.5.rds")

lyn.predprey <- list(Predprey.lyn.nonenone, #1
                     Predprey.lyn.huntnone, #2
                     Predprey.lyn.huntSH20none, #3
                     Predprey.lyn.huntSH50none, #4
                     Predprey.lyn.huntRH20none, #5
                     Predprey.lyn.huntRH50none, #6
                     Predprey.lyn.hunthunt, #7
                     Predprey.lyn.huntSH20hunt, #8
                     Predprey.lyn.huntRH20hunt, #9
                     Predprey.lyn.huntSH50hunt, #10
                     Predprey.lyn.huntRH50hunt, #11
                     Predprey.lyn.huntgoal1, #12
                     Predprey.lyn.huntgoal05, #13
                     Predprey.lyn.huntSH20goal1, #14
                     Predprey.lyn.huntSH50goal1, #15
                     Predprey.lyn.huntSH20goal05, #16
                     Predprey.lyn.huntSH50goal05, #17
                     Predprey.lyn.huntRH20goal1, #18
                     Predprey.lyn.huntRH50goal1, #19
                     Predprey.lyn.huntRH20goal05, #20
                     Predprey.lyn.huntRH50goal05) #21


# ---
# --- Roe deer - PLOTS -----
# ---

# --- Data preparation ---
# Roepop <- ArrayToDataframe(Simulation.roe)
# Simulation.roe <- NormFe.list[[4]]
Roe.sumpop <- PopSum(Simulation.roe)
Roe.quant <- PopQuantile(Simulation.roe)
GR.roe <- GrowthRate(Simulation.roe)
GR.roe.quant <- GrowthRateQuantile(Simulation.roe)
Roe.state <- StateProp(Simulation.roe, Ne=sum(N.roe.init)/20, Ni=40*study.area)
# Roe.state.pivot <- pivot_longer(Roe.state, 2:4, names_to = "state", values_to = "proportion", cols_vary = "slowest")

# --- Plots ---
# # Roe deer population density over time - QUANTILES
# Roe.quant <- PopQuantile(Roe.sumpop)
ggplot(Roe.quant) +
  geom_ribbon(aes(x=time, ymin=q95/study.area, ymax=q5/study.area), fill="darkolivegreen", alpha=0.5) +
  geom_ribbon(aes(x=time, ymin=q95/study.area, ymax=q5/study.area), fill="grey", alpha=0.5) +
  geom_ribbon(aes(x=time, ymin=q75/study.area, ymax=q25/study.area), fill="darkolivegreen", alpha=0.5) +
  geom_line(aes(x=time, y=q0/study.area), color="darkgrey", linetype="dashed", size=0.5) +
  geom_line(aes(x=time, y=q50/study.area), color="black", size=1) +
  geom_line(aes(x=time, y=q100/study.area), color="darkgrey", linetype="dashed", size=0.5) +
  # geom_hline(yintercept = 40) + # Carrying capacity ~40
  theme(legend.key.height= unit(5, 'cm'),
        text = element_text(size = 50),
        axis.text = element_text(size = 50)) +
  ylim(0, 65) + #80 or 17.5
  labs(title = "Roe deer population density", x = "Time (Years)", y = "Population density (N/km2)")



# # Lambda growth rate - QUANTILES
# GR.roe.quant <- GrowthRateQuantile(GR.roe)
ggplot(GR.roe.quant) +
  geom_ribbon(aes(x=time, ymin=l.q95, ymax=l.q5), fill="darkolivegreen", alpha=0.5) +
  geom_ribbon(aes(x=time, ymin=l.q95, ymax=l.q5), fill="grey", alpha=0.5) +
  geom_ribbon(aes(x=time, ymin=l.q75, ymax=l.q25), fill="darkolivegreen", alpha=0.5) +
  geom_line(aes(x=time, y=l.q50), color="black", size=1) +
  geom_line(aes(x=time, y=l.q0), color="darkgrey", size=0.5, linetype="dashed") +
  geom_line(aes(x=time, y=l.q100), color="darkgrey", size=0.5, linetype="dashed") +
  geom_hline(yintercept=1, linetype="dotdash", color="black", size=0.75) +
  theme(legend.key.height= unit(5, 'cm'),
        text = element_text(size = 50),
        axis.text = element_text(size = 50)) +
  ylim(0, 2) +
  labs(title = "Roe deer population growth rate", x = "Time (Years)", y = "Population growth rate (λ)")


# # Roe deer state - PROPORTIONS
ggplot(Roe.state, aes(x=time, y=proportion*100, fill=state)) +
  geom_col(position = position_stack(reverse = T), width = 0.75) +
  scale_fill_manual("State", values = c("irruption" = "deepskyblue4", "stable" = "snow3", "extinct" = "darkred")) +
  theme(legend.key.height= unit(5, 'cm'),
        text = element_text(size = 50),
        axis.text = element_text(size = 50)) +
  labs(title = "State proportion of populations", x = "Time (Years)", y = "Proportion (%)")



# --- Plot all at once!
# grid.arrange(Pop.1, Pop.2, Pop.3, Pop.4,
#              GR.1, GR.2, GR.3, GR.4,
#              Ext.1, Ext.2, Ext.3, Ext.4,
#              ncol=4, top=textGrob("Roe deer populations - 500 runs"))
# Use patchwork instead!


# ---
# --- Lynx - PLOTS -----
# ---

# --- Data preparation ---
# Simulation.lyn <- Lynx.list[[1]]
Lyn.sumpop <- PopSum(Simulation.lyn)
Lyn.quant <- PopQuantile(Simulation.lyn)
GR.lyn <- GrowthRate(Simulation.lyn)
GR.lyn.quant <- GrowthRateQuantile(Simulation.lyn)
Lyn.state <- StateProp(Simulation.lyn, Ne=sum(N.lyn.init)/50, Ni=(5.5/100)*study.area) # Ni???


# Lynx population density
ggplot(Lyn.quant) +
  geom_ribbon(aes(x=time, ymin=q95*100/study.area, ymax=q5*100/study.area), fill="darkolivegreen", alpha=0.5) +
  geom_ribbon(aes(x=time, ymin=q95*100/study.area, ymax=q5*100/study.area), fill="grey", alpha=0.5) +
  geom_ribbon(aes(x=time, ymin=q75*100/study.area, ymax=q25*100/study.area), fill="darkolivegreen", alpha=0.5) +
  geom_line(aes(x=time, y=q0*100/study.area), color="darkgrey", linetype="dashed", size=0.5) +
  geom_line(aes(x=time, y=q50*100/study.area), color="black", size=1) +
  geom_line(aes(x=time, y=q100*100/study.area), color="darkgrey", linetype="dashed", size=0.5) +
  # geom_hline(yintercept = 36.5) + # Carrying capacity ~47
  theme(legend.key.height= unit(5, 'cm'),
        text = element_text(size = 50),
        axis.text = element_text(size = 50)) +
  # ylim(0, 4) + # 5 or 2
  labs(title = "Lynx population density", x = "Time (Years)", y = "Population density (N/100km2)")


# Lambda growth rate
ggplot(GR.lyn.quant) +
  geom_ribbon(aes(x=time, ymin=l.q95, ymax=l.q5), fill="darkolivegreen", alpha=0.5) +
  geom_ribbon(aes(x=time, ymin=l.q95, ymax=l.q5), fill="grey", alpha=0.5) +
  geom_ribbon(aes(x=time, ymin=l.q75, ymax=l.q25), fill="darkolivegreen", alpha=0.5) +
  geom_line(aes(x=time, y=l.q50), color="black", size=1) +
  geom_line(aes(x=time, y=l.q0), color="darkgrey", size=0.5, linetype="dashed") +
  geom_line(aes(x=time, y=l.q100), color="darkgrey", size=0.5, linetype="dashed") +
  geom_hline(yintercept=1, linetype="dotdash", color="black", size=0.75) +
  theme(legend.key.height= unit(5, 'cm'),
        text = element_text(size = 50),
        axis.text = element_text(size = 50)) +
  ylim(0, 2) +
  labs(title = "Lynx population growth rate", x = "Time (Years)", y = "Population growth rate (λ)")


# # Lynx state - PROPORTIONS
ggplot(Lyn.state.pivot, aes(x=time, y=proportion*100, fill=state)) +
  # geom_bar(position = "fill", stat = "identity") +
  geom_col(position = position_stack(reverse = T), width = 0.75) +
  scale_fill_manual("State", values = c("irruption" = "deepskyblue4", "stable" = "snow3", "extinct" = "darkred")) +
  theme(legend.key.height= unit(5, 'cm'),
        text = element_text(size = 50),
        axis.text = element_text(size = 50)) +
  labs(title = "State proportion of populations", x = "Time (Years)", y = "Proportion (%)")



# ---
# --- Predator prey - PLOTS -----
# ---

Simulation.roe <- roe.predprey[[20]]
Simulation.lyn <- lyn.predprey[[20]]

Roe1 <- DensityPlot(Simulation.roe, study.area, "Roe deer", x.years = 25, carrying.capacity = (sum(N.roe.init)/(100/20))/study.area)
Roe2 <- GrowthRatePlot(Simulation.roe, study.area, "Roe deer", x.years = 25)
Roe3 <- StatePropPlot(Simulation.roe, study.area, "Roe deer", Ne=sum(N.roe.init)/20, Ni=40*study.area, x.years = 25)

Lyn1 <- DensityPlot(Simulation.lyn, study.area, "Lynx", x.years = 25, y.lim = 2.75) #y.lim = 2.75
Lyn2 <- GrowthRatePlot(Simulation.lyn, study.area, "Lynx", x.years = 25)
Lyn3 <- StatePropPlot(Simulation.lyn, study.area, "Lynx", Ne=sum(N.lyn.init)/20, Ni=40*study.area, x.years = 25)

(Roe1 + Roe2 + Roe3) / (Lyn1 + Lyn2 + Lyn3)


# Prepare data
ratio.data <- PredPreyQuantile(Simulation.lyn, Simulation.roe)
# Plot it
ggplot(ratio.data) +
  geom_ribbon(aes(x=time, ymin=ratio.q95, ymax=ratio.q5), fill="darkolivegreen", alpha=0.5) +
  geom_ribbon(aes(x=time, ymin=ratio.q95, ymax=ratio.q5), fill="grey", alpha=0.5) +
  geom_ribbon(aes(x=time, ymin=ratio.q75, ymax=ratio.q25), fill="darkolivegreen", alpha=0.5) +
  geom_line(aes(x=time, y=ratio.q0), color="darkgrey", linetype="dashed", size=0.5) +
  geom_line(aes(x=time, y=ratio.mean), color="black", size=1) +
  geom_line(aes(x=time, y=ratio.q100), color="darkgrey", linetype="dashed", size=0.5) +
  theme(legend.key.height= unit(2.5, 'cm'),
        text = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  xlim(0, max(ratio.data$time)) +
  ylim(0, max(ratio.data$ratio.q100)) +
  labs(title = ("Predator prey ratio"),
       x = "Time (Years)",
       y = "Ratio log(Roe deer/Lynx)")


data <- PredPreyRatioGR(Simulation.lyn, Simulation.roe)
# plot(x = data$ratio, y = data$lambda.prey)
ggplot(data) + #aes(x = ratio, y = lambda.prey)
  geom_ribbon(aes(x = ratio, ymin = min(lambda.prey, na.rm = T), ymax = max(lambda.prey, na.rm = T)), fill="darkolivegreen", alpha=0.5) +
  # xlim(0, max(ratio.data$time)) +
  ylim(-2, 2) +
  labs(title = ("Predator prey growth rate ~ ratio"),
       x = "Ratio log(Roe deer/Lynx)",
       y = "Growth rate (lambda)")
  # geom_ribbon(aes(x = ratio, ymin = min(lambda.pred), ymax = max(lambda.pred)), fill="grey", alpha=0.5) +

  
#### Experiment? (Plots in progress) ####
# Create both population data frame
# Sim.pops <- merge(Roe.sumpop, Lyn.sumpop, by = c("simulation", "time"))
# colnames(Sim.pops) <- c("simulation", "time", "roe", "lyn")
# Two populations, one simulation
ggplot(AllTimeSim.pops, aes(x=time), group=simulation) +
  geom_line(aes(y=roe/study.area, color= "dark green")) +
  geom_line(aes(y=(lyn*100)/study.area, color= "orange")) +
  labs(title = "Predator prey population densities", x = "Time (Years)", y = "Population density (N/100 km2)")

#### Obsolete plots ####

# # Lynx population all stages over time for each simulation
Lynpop <- ArrayToDataframe(Lynx.list[[1]])
ggplot(Lynpop) +
  geom_line(aes(x=time, y=value, group=interaction(stage,simulation), color=stage), size=1) +
  scale_color_brewer("Stages", palette = "Paired") +
  scale_x_continuous(breaks = seq(0, nYears, by = 2)) +
  scale_y_continuous(limits = c(0, max(Lynpop$value))) +
  labs(title = "Lynx population", x = "Time (Years)", y = "Population size (N)")
# 
# # Lynx population density over time for each simulation
# ggplot(Lyn.sumpop) +
#   geom_line(aes(x=time, y=(value/study.area)*100, alpha=simulation), size=1) +
#   scale_x_continuous(breaks = seq(0, nYears, by = 2)) +
#   theme(legend.position = "none") +
#   labs(title = "Lynx population density simulation", x = "Time (Years)", y = "Population density (N/100 km2)")
# 
# # Lynx population density over time QUANTILES
# ggplot(Lyn.quant) +
#   geom_ribbon(aes(x=time, ymin=(q100/study.area)*100, ymax=(q0/study.area)*100), fill="grey", alpha=0.5) +
#   geom_line(aes(x=time, y=(q0/study.area)*100), color="dark grey", size=0.5) +
#   geom_line(aes(x=time, y=(q25/study.area)*100), color="orange", linetype="dashed", size=1) +
#   geom_line(aes(x=time, y=(q50/study.area)*100), color="dark green", size=1) +
#   geom_line(aes(x=time, y=(q75/study.area)*100), color="orange", linetype="dashed", size=1) +
#   geom_line(aes(x=time, y=(q100/study.area)*100), color="dark grey", size=0.5) +
# 
#   labs(title = "Lynx population density - quantiles", x = "Time (Years)", y = "Population density (N/100 km2)")
# 
# # The annual population growth rate for each simulation
# # Lambda growth rate
# ggplot(GR.lyn) +
#   geom_line(aes(x=time, y=lambda, alpha=simulation), size=1) +
#   geom_hline(yintercept=1, linetype="dashed", color="orange", size=1 ) +
#   scale_x_continuous(breaks = seq(0, nYears, by = 2)) +
#   # scale_y_continuous(limits = c(0, ceiling(max(GR.lyn[, 4], na.rm = T)))) +
#   theme(legend.position = "none") +
#   labs(title = "Lynx population growth rate", x = "Time (Years)", y = "Population growth rate (λ)")
# 
# # Lambda growth rate - QUANTILES
# ggplot(GR.lyn.quant) +
#   geom_ribbon(aes(x=time, ymin=l.q100, ymax=l.q0), fill="grey", alpha=0.5) +
#   geom_line(aes(x=time, y=l.q0), color="dark grey", size=0.5) +
#   geom_line(aes(x=time, y=l.q25), color="orange", linetype="dashed", size=1) +
#   geom_line(aes(x=time, y=l.q50), color="dark green", size=1) +
#   geom_line(aes(x=time, y=l.q75), color="orange", linetype="dashed", size=1) +
#   geom_line(aes(x=time, y=l.q100), color="dark grey", size=0.5) +
#   geom_hline(yintercept=1, linetype="twodash", color="dark red", size=1) +
#   labs(title = "Roe deer population growth rate - quantiles", x = "Time (Years)", y = "Population growth rate (λ)")
# 
# # Growth rate in %
# ggplot(GR.lyn) +
#   geom_line(aes(x=time, y=growth.rate, alpha=simulation), size=1) +
#   geom_hline(yintercept=0, linetype="dashed", color="orange", size=1 ) +
#   scale_x_continuous(breaks = seq(0, nYears, by = 2)) +
#   scale_y_continuous(breaks = seq(round(min(GR.lyn[, 5], na.rm = T)), round(max(GR.lyn[, 5], na.rm = T)), by = 4)) +
#   theme(legend.position = "none") +
#   labs(title = "Lynx population growth rate", x = "Time (Years)", y = "Population growth rate (%)")
# 
# # Growth rate in % QUANTILES
# ggplot(GR.lyn.quant) +
#   geom_ribbon(aes(x=time, ymin=p.q100, ymax=p.q0), fill="grey", alpha=0.5) +
#   geom_line(aes(x=time, y=p.q0), color="dark grey", size=0.5) +
#   geom_line(aes(x=time, y=p.q25), color="orange", linetype="dashed", size=1) +
#   geom_line(aes(x=time, y=p.q50), color="dark green", size=1) +
#   geom_line(aes(x=time, y=p.q75), color="orange", linetype="dashed", size=1) +
#   geom_line(aes(x=time, y=p.q100), color="dark grey", size=0.5) +
#   geom_hline(yintercept=0, linetype="twodash", color="dark red", size=1) +
#   labs(title = "Roe deer population growth rate - quantiles", x = "Time (Years)", y = "Population growth rate (%)")
