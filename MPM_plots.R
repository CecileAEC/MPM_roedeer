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
# From:
# source("MPM_load_parameters.R")
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
Roe.favourable <- readRDS("roesim-fox_none-lynx_none-hunt_none.rds")
Roe.fox <- readRDS("roesim-fox_fixed-lynx_none-hunt_none.rds")
Roe.lynx <- readRDS("roesim-fox_none-lynx_selective-hunt_none.rds")
Roe.predation <- readRDS("roesim-fox_fixed-lynx_selective-hunt_none.rds")
Roe.allmortality <- readRDS("roesim-fox_fixed-lynx_selective-hunt_selective.rds")
Roesim.list <- list(Roe.favourable, Roe.fox, Roe.lynx, Roe.predation, Roe.allmortality)

# Explore Snow event
Roe.snow.favourable <- readRDS("roesim-fox_none-lynx_none-hunt_none-SNOW.rds")
Roe.snow.fox <- readRDS("roesim-fox_fixed-lynx_none-hunt_none-SNOW.rds")
Roe.snow.lynx <- readRDS("roesim-fox_none-lynx_selective-hunt_none-SNOW.rds")
Roe.snow.predation <- readRDS("roesim-fox_fixed-lynx_selective-hunt_none-SNOW.rds")
Roe.snow.allmortality <- readRDS("roesim-fox_fixed-lynx_selective-hunt_selective-SNOW.rds")
Roesim.snow.list <- list(Roe.snow.favourable, Roe.snow.fox, Roe.snow.lynx, Roe.snow.predation, Roe.snow.allmortality)


# ---
# --- Roe deer population - PLOTS -----
# ---

# --- Data preparation ---
Simulation.roe <- Roesim.list[[1]] # Select which scenario you want to plot from the lists above
Roe.sumpop <- PopSum(Simulation.roe)
Roe.quant <- PopQuantile(Simulation.roe)
GR.roe <- GrowthRate(Simulation.roe)
GR.roe.quant <- GrowthRateQuantile(Simulation.roe)
Roe.state <- StateProp(Simulation.roe, Ne=sum(N.roe.init)/20, Ni=40*study.area)

# --- Single plots ---
# # Roe deer population density over time - QUANTILES
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
  ylim(0, 65) +
  labs(title = "Roe deer population density", x = "Time (Years)", y = "Population density (N/km2)")



# # Lambda growth rate - QUANTILES
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
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual("State", values = c("irruption" = "deepskyblue4", "stable" = "snow3", "extinct" = "darkred")) +
  theme(legend.key.height= unit(5, 'cm'),
        text = element_text(size = 50),
        axis.text = element_text(size = 50)) +
  labs(title = "State proportion of populations", x = "Time (Years)", y = "Proportion (%)")



# --- Three plots from one scenario at once:
Simulation.roe <- Roesim.list[[1]] # Select which scenario you want to plot from the lists of data

Roe1 <- DensityPlot(Simulation.roe, study.area, "Roe deer", x.years = 25, y.lim = 65)
Roe2 <- GrowthRatePlot(Simulation.roe, study.area, "Roe deer", x.years = 25)
Roe3 <- StatePropPlot(Simulation.roe, study.area, "Roe deer", Ne=sum(N.roe.init)/20, Ni=40*study.area, x.years = 25)

Roe1 + Roe2 + Roe3



# -
# ----- Roe deer gradient - DATA -----
# -

Roerate1 <- readRDS("roesim-gradientrate-fox_none-lynx_selective-hunt_none.rds")
Roerate2 <- readRDS("roesim-gradientrate-fox_fixed-lynx_selective-hunt_none.rds")
Roerate3 <- readRDS("roesim-gradientrate-fox_fixed-lynx_selective-hunt_selective.rds")

# -
# ----- Roe deer gradient - PLOTS -----
# -
Roe1 <- ggplot(Roerate1, aes(x=ar, y=(yd/sum(N.roe.init))*100, fill=prop25)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred", "orange",
                                   "white"),
                       breaks=c(-1, -0.75, -0.5, -0.25, 0),
                       labels=c("100%", "75%", "50%", "25%","0%"),
                       limits=c(-1, 0)) +
  labs(title = "Roe deer population state proportion
- Lynx predation only",
       x = "Adult refuge
(% of refuge from lynx predation)",
       y = "Yearling immigration
(% of initial population density)",
       fill = "Extinction
proportion") +
  theme(legend.key.height= unit(3, 'cm'),
        text = element_text(size = 25),
        axis.text = element_text(size = 25))

Roe2 <- ggplot(Roerate2, aes(x=ar, y=(yd/sum(N.roe.init))*100, fill=prop25)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred", "orange",
                                   "white"),
                       breaks=c(-1, -0.75, -0.5, -0.25, 0),
                       labels=c("100%", "75%", "50%", "25%","0%"),
                       limits=c(-1, 0)) +
  labs(title = "Roe deer population state proportion
- Red fox and lynx predation",
       x = "Adult refuge
(% of refuge from lynx predation)",
       y = "Yearling immigration
(% of initial population density)",
       fill = "Extinction
proportion") +
  theme(legend.key.height= unit(3, 'cm'),
        text = element_text(size = 25),
        axis.text = element_text(size = 25))

Roe3 <- ggplot(Roerate3, aes(x=ar, y=(yd/sum(N.roe.init))*100, fill=prop25)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred", "orange",
                                   "white"),
                       breaks=c(-1, -0.75, -0.5, -0.25, 0),
                       labels=c("100%", "75%", "50%", "25%","0%"),
                       limits=c(-1, 0)) +
  labs(title = "Roe deer population state proportion
 - Hunting, red fox and lynx predation",
       x = "Adult refuge
(% of refuge from lynx predation)",
       y = "Yearling immigration
(% of initial population density)",
       fill = "Extinction
proportion") +
  theme(legend.key.height= unit(3, 'cm'),
        text = element_text(size = 25),
        axis.text = element_text(size = 25))

Roe1 + Roe2 + Roe3


