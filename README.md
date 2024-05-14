# R code and files to run and plot the roe deer matrix population model

This repository is attached to the manuscript entitled "Understanding Roe Deer Population Dynamics in Boreal Environments: Exploring the Impact of Predation, Harvest Pressures and environmental variability on Population Stability and Trends."

This repository contains all R scripts and data files to run and plot the results of the manuscript.

- The R file called "MPM_packages_and_functions.R" contain all the necessary packages and the functions used in the simulations and plots.
- The R file called "MPM_load_parameters.R" contain all the parameters value to run the model.
- The R file called "MPM_stochastic_parameters.R" is used to add stochastcity around the parameter value of the model.
The stochasticity can be turned off in the "MPM_load_parameters.R" file to get a deterministic run with initial parameter values.
- The R file called "MPM_simulations_roedeer.R" contain the model to run the matrix population model as well as the gradient study model.
- The R file called "MPM_plots.R" load the data file .RDS for the different scenarios.

The RDS file are called after their scenario names:
- Scenario without hunting and predation:
"roesim-fox_none-lynx_none-hunt_none.rds"
- Scenario with red fox predation:
"roesim-fox_fixed-lynx_none-hunt_none.rds"
- Scenario with lynx predation:
"roesim-fox_none-lynx_selective-hunt_none.rds"
- Scenario with both red fox and lynx predation:
"roesim-fox_fixed-lynx_selective-hunt_none.rds"
- Scenario with predation and hunting:
"roesim-fox_fixed-lynx_selective-hunt_selective.rds"

The same scenario including snow event are called:
- "roesim-fox_none-lynx_none-hunt_none-SNOW.rds"
- "roesim-fox_fixed-lynx_none-hunt_none-SNOW.rds"
- "roesim-fox_none-lynx_selective-hunt_none-SNOW.rds"
- "roesim-fox_fixed-lynx_selective-hunt_none-SNOW.rds"
- "roesim-fox_fixed-lynx_selective-hunt_selective-SNOW.rds"

The gradient study data are also called after their scenario name:
- "roesim-gradientrate-fox_none-lynx_selective-hunt_none.rds"
- "roesim-gradientrate-fox_fixed-lynx_selective-hunt_none.rds"
- "roesim-gradientrate-fox_fixed-lynx_selective-hunt_selective.rds"

