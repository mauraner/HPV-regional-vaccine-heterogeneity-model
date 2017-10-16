require(deSolve)

# Args from array in bash
args=(commandArgs(TRUE))
i=as.numeric(unlist(args))

# set number of cantons
k <- 26

#define the number of simulations wanted for the sensitivity analysis
sens_rep <- 1

#run the parameters script and models over again
source("N_01_params_rev.R") # N_01_params runs the parameters of the model
source("N_02_models.R") # N_02_models runs the two main SIRV models: STI_1 for scenario 1&2 and STI_2 for scenario 3
source("N_01b_params_rev.R") # This creates a sens_rep amount of sets of parameters

# pset <- matrix(ncol=14, nrow=sens_rep)
# for( i in 1:sens_rep){
#   pset[i,] <- c(mu=mu, omega=omega[i], gamma=gamma[i], 
#                 epsilon=epsilon[i], intercantmix=intercantmix, m=m[i], p=mp, C=C, N=N, n=n, beta=beta[i])
# }
# colnames(pset) <- c("mu", "omega", "gamma", "epsilon", "intercantmix", "m", "p", "C1", "C2", "N1", "N2", "n1", "n2", "beta" )

#simulation time
time <- seq(0, 400, 1)

# set mixing between population (not relevant for single canton)
intercantmix <- 1
# intercantmix <- 0.8
# intercantmix <- weighted.mean(diag(sigma), w)


# pb <- txtProgressBar(1, sens_rep, style=3)
#run the model for 26 cantons (heterogeneous situation)
init <- as.vector(matrix(unlist(inits(k)), ncol = k, byrow = TRUE))
sto_sim_het <- array(NA, c(length(time), ((r*8)*k)+1, sens_rep ) )

  C <- cbind(rep(sens_clow, k), rep(sens_chigh, k))
  n <- cbind(rep(sens_nlow,k), rep(sens_nhigh,k))
  N <- cbind(rep(sens_nlow,k)*w, rep(sens_nhigh,k)*w)

  params <- c(mu=mu, omega=omega, gamma=gamma, 
              epsilon=epsilon, intercantmix=intercantmix, m=m, p=p, C=C, N=N, n=n, beta=beta)
  parms <- params
  simulation_26k<- as.data.frame(ode(init, time, STI_1, parms=params)) #for the assortative mixing
  # simulation_1k<- as.data.frame(ode(init, time, STI_2, parms=params)) #for the mobility informed mixing
  simulation_26k <- as.matrix(simulation_26k)
  sto_sim_het[,,] <- simulation_26k
  # setTxtProgressBar(pb,i)


# save(sto_sim_het, file = "sens100sim26k_eps1.RData")



k <- 1
w <- 1
p <- mp
init <- as.vector(matrix(unlist(inits(k)), ncol = k, byrow = TRUE))
sto_sim_hom <- array(NA, c(length(time), ((r*8)*k)+1, sens_rep ) )

  C <- cbind(rep(sens_clow, k), rep(sens_chigh, k))
  n <- cbind(rep(sens_nlow,k), rep(sens_nhigh,k))
  N <- cbind(rep(sens_nlow,k)*w, rep(sens_nhigh,k)*w)
  # rho_ki <- rho_ki_f(epsilon, intercantmix, k,C)
  # rho_ki <- rho_sra(epsilon[i], intercantmix, C, N, k)
  params <- c(mu=mu, omega=omega, gamma=gamma, 
              epsilon=epsilon, intercantmix=intercantmix, m=m, p=p, C=C, N=N, n=n, beta=beta)
  parms <- params
  simulation_1k<- as.data.frame(ode(init, time, STI_1, parms=params))
  simulation_1k <- as.matrix(simulation_1k)
  sto_sim_hom[,,] <- simulation_1k


# save(sto_sim_hom, file = "sens1sim1k_eps1.RData")


pset1 <- matrix(ncol=14, nrow=sens_rep)

  C <- cbind(rep(sens_clow, k), rep(sens_chigh, k))
  n <- cbind(rep(sens_nlow,k), rep(sens_nhigh,k))
  N <- cbind(rep(sens_nlow,k)*w, rep(sens_nhigh,k)*w)
  
  pset1<- c(mu=mu, omega=omega, gamma=gamma, 
                 epsilon=epsilon, intercantmix=intercantmix, m=m, p=mp, C=C, N=N, n=n, beta=beta)
  

# t1 <- c(mu=mu, omega=omega, gamma=gamma, epsilon=epsilon, intercantmix=intercantmix, m=m, p=mp, C=C, N=N, n=n, beta=beta)
# names(t1)
# colnames(pset1) <- names(t1)
# save(pset1, file = "pset1k1sim.RData")

res = list(r1=sto_sim_het, r2=sto_sim_hom, r3=pset1)

save(res, file="p.RData")

# save(res, file=paste("p",i,".RData", sep=""))
