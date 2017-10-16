# Sensitivity analyse, 

########################### parameters that do not change
# 
# #### Load vaccination coverage file (cantons are always ordered in the same way)
# kv <-  read.csv("cant_vacc_cov_pmob_ordered.csv")
# kv <- as.numeric(kv[,2])
# 
# #vaccination coverage for each canton
# p <- kv/100
# 
# ####. Load the population size weights (w_k = N_k/ mean(N_k) and N_k is the population size of 18-24 y.o in 2013, men and women) )
# kw <-read.csv("weightsYA_pmob_ordered.csv")
# kw <- as.numeric(kw[,2])
# 
# 
# w <- if(k==26)  kw else (rep(1,k))
# 
# 
# # weighted mean vaccination coverage
# mp <- weighted.mean(p, kw)
# 
# r <- 2
# mu <- 0.143

#sens_rep <- 1000

#summary(sens_clow)
########################### parameters that can vary

#vaccine efficacy
ve <-runif(sens_rep, min=0.908531, max=0.9556387)
# ve mean is : 1-0.0637

#sexual activity
# sens_clow <-rnorm(sens_rep, mean= 0.3728684, sd=0.02237743 )
sens_clow <-rnorm(sens_rep, mean= 0.171, sd=0.01388649 )



# sens_chigh <-rnorm(sens_rep, mean= 3.5960606, sd=0.20415769 )
sens_chigh <-rnorm(sens_rep, mean= 2.409, sd=0.12000095 )

# 
# sens_clow <-rep(0.171 ,sens_rep, )
# sens_chigh <-rep( 2.409,sens_rep, )

# sens_nlow <-rnorm(sens_rep, mean= 0.853, sd=0.0354 )
sens_nlow <-rnorm(sens_rep, mean= 0.853, sd=0.01021346 )

sens_nhigh <- 1-sens_nlow

# sens_nlow <-rep(0.853,sens_rep )
# sens_nhigh <- 1-sens_nlow


# C_ch <- cbind(rep(sens_clow[i], k), rep(sens_chigh[i], k))
# n <- cbind(rep(sens_nlow[i],k), rep(sens_nhigh[i],k))
# N <- cbind(rep(sens_nlow[i],k)*w, rep(sens_nhigh[i],k)*w)


#other
m <- runif(sens_rep, min=0, max=1)
# m <- rep(1,sens_rep)
epsilon <- runif(sens_rep, min=0, max=1)
# epsilon <- rep(0.5,sens_rep)

#immunity duration (1/omega= duration for becoming susceptible again)
omega <- runif(sens_rep, min=0.011, max=0.032)
# omega <- rep( 0.024,sens_rep)

#infection clearance rate (1/gamma = duration of the infectiousness)
gamma <- runif(sens_rep, min=0.146, max=1.159)
# gamma <- rep( 0.542,sens_rep)

#Transmission probability "matrix"
beta <- runif(sens_rep, min=0.60, max=0.99)
# beta <- rep( 0.8,sens_rep)

##initial params: from normal model after 200y (without vacc)
# init_2 <- read.csv("O:/Test_R/R/2016_Septembre/init_2.csv")
# 
# init_2 <- init_2[2:17]
# 
# # itest <- as.numeric(init_2)
# 0.01981521/N[1]
# 0.01361767/N[2]



