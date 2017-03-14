# HPV-regional-vaccine-heterogeneity-model
Mathematical model of HPV-16 transmission in the context of heterogeneous vaccination uptake (example of Switzerland)

In this depository you will find the following files:

1: 
N_01_params.R : Script which codes the parameters used for the model

It requires the following data: 
cant_vacc_cov_pmob_ordered.csv <- cantonal/regional vaccination coverages
weightsYA_pmob_ordered.csv     <- weight for each canton, according to its population size (18-24y.o men and women, www.BFS.ch      2013)
pmob.csv                        <- population mobility matrix, received from the Swiss Federal Office for Spatial Development (ARE), represents the sum of commuters (average per day over one working week, 2010) between all pairs of cantons, by public transport and personal motorised vehicles)

2. 
N_02_models.R
Contains all function required to run the core script

3.
Core_script.R

This script is sub-divided into three parts:

Part 1: One population model
System with only one population. It gives you the prevalence under different vaccination coverages and the R0 (1.a), sensibility analysis of the per-partnership transmissibility parameter with different sexual behaviour assumptions (Swiss based or UK based) (1.b), and shows the relative reduction in prevalence over different vaccination coverages at different time points (years after vaccination onset) and compared it the published data from a systematic review (Drolet et al, 2015 Lancet).

Part 2: Two sub-populations model
System with two sub-population. Used to investigate the effect of heterogeneity in a simple model with two sub-population in which the vaccination coverage ranges from 0 to 100% AND investigate the effect of sexual mixing between the sub-population under different vaccination coverage situations. 

Part 3: 26 sub-populations/cantons model
System representing Switzerland and it's 26 cantons. 
3.a 26 cantons with different vaccination coverages (based on the Federal Office of Public Health vaccination survey of 16y. old women, 2011-2013, www.bag.ch). Scenario 1 (fully assortative, no sexual mixing occurs between cantons) and scenario 2 (partial random mixing, 80% of sexual contacts occur within the own canton, 20% randomly distributed among all other cantons according to population size).
3.b Same as 3.a but for scenario 3 (80% of sexual contacts occur within the own canton (overall weighted average), but the 20% are distributed according to the population mobility matrix. 
3.c Gives the number of cantons achieving certain number of relative risk reduction over different sexual mixing levels (0= completely random to 1= completely assortative) and at different time points (number of years after vaccination onset).
3.d Migration plot or chord diagram of the cantonal interconnectivity, according to each scenario.

