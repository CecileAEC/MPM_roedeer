                    # ------------------------------- #
                    #         Population plots        #
                    #       Roe deer populations      #
                    # ------------------------------- #


# -
# ----- Data preparation and plot functions -----
# -
source("MPM_packages_and_functions.R")

# -
# ----- Load default parameters -----
# -
# You just need:
study.area <- 10000 # Study area in square kilometre
roe.density <- 10 # Nb of roe deer per square kilometre
nameRoeDeer <- c("Jroe.F", # Juvenile roe deer female
                 "Yroe.F", # Yearling roe deer female
                 "Proe.F", # Two year-old roe deer female
                 "Aroe.F", # Adult roe deer female
                 "Jroe.M", # Juvenile roe deer male
                 "Yroe.M", # Yearling roe deer male
                 "Aroe.M") # Adult roe deer male

N.roe.init <- VectorPop(round(c(0.237, 0.105, 0.045, 0.147, 0.237, 0.082, 0.147)
                              * (roe.density*study.area)), nameRoeDeer)



# -
# ----- Roe deer population - DATA -----
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



# ---
# --- Roe deer population - PLOTS -----
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
Simulation.roe <- roe.predprey[[20]]

Roe1 <- DensityPlot(Simulation.roe, study.area, "Roe deer", x.years = 25, carrying.capacity = (sum(N.roe.init)/(100/20))/study.area)
Roe2 <- GrowthRatePlot(Simulation.roe, study.area, "Roe deer", x.years = 25)
Roe3 <- StatePropPlot(Simulation.roe, study.area, "Roe deer", Ne=sum(N.roe.init)/20, Ni=40*study.area, x.years = 25)





# -
# ----- Roe deer gradient - DATA -----
# -
Roerate <- Roe.rate

Roerate <- readRDS("Roedeer-gradientrate-25_500-fox_none-lynx_selective-hunt_none-var10-441.rds")
Roerate <- readRDS("Roedeer-gradientrate-25_500-fox_fixed-lynx_selective-hunt_none-var10-441.rds")
Roerate <- readRDS("Roedeer-gradientrate-25_500-fox_fixed-lynx_selective-hunt_selective-var10-441.rds")

# -
# ----- Roe deer gradient - PLOTS -----
# -
ggplot(Roerate, aes(x=ar, y=(yd/sum(N.roe.init))*100, fill=prop25)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred", "orange",
                                   "white"),
                       breaks=c(-1, -0.75, -0.5, -0.25, 0),
                       labels=c("quasi-extinct", "75%", "50%", "25%","stable"),
                       limits=c(-1, 0)) +
  labs(title = "Roe deer population state proportion
- Red fox and lynx predation",
# Choose title depending on sceanario:
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

(Roe1 + Roe2 + Roe3)


