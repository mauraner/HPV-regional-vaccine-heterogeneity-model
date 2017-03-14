# 10.02.2017
# 
# 01_Params


#### Load vaccination coverage file (cantons are always ordered in the same way)
kv <-  read.csv("cant_vacc_cov_pmob_ordered.csv")
kv <- as.numeric(kv[,2])

####. Load the population size weights (w_k = N_k/ mean(N_k) and N_k is the population size of 18-24 y.o in 2013, men and women) )
kw <-read.csv("weightsYA_pmob_ordered.csv")
kw <- as.numeric(kw[,2])


w <- if(k==26)  kw else (rep(1,k))
# w <- rep(1,k)

#vaccination coverage for each canton
p <- kv/100
# weighted mean vaccination coverage
mp <- weighted.mean(p, kw)


S <- 680168    # S is the total population size in Switzerland of the targeted age group
# meanksize <- S/26  
# w1 <- w*meanksize  # This gives back again the N_k

#### Number cantons
# k <- 26
# k <- 2 # if k is not 26, set w to  1 (k times)
# k <- 1

#### Number group risks
r <- 2

#### contacts low in cantons 
Cl_ch <- rep(0.171, k) # According to Swiss data from SIR survey
Cl_uk <- rep(0.373, k) # According to British and Scottland data from NATSAL-3 survey

#contacts high in cantons 
Ch_ch <- rep(2.409, k) #CH
Ch_uk <- rep(3.596, k) 


#Pop low in cantons
nl_1 <- rep(0.853,k)

#Pop high in cantons
nh_1 <- rep(0.147,k)

#proportion in each risk group
n <- cbind(nl_1,nh_1)

#Weighted Pop low in cantons
nl <- nl_1 *w
#Weighted Pop high in cantons
nh <- nh_1 * w

C_ch <- cbind(Cl_ch, Ch_ch) # depends what parameters wanted
C_uk <- cbind(Cl_uk, Ch_uk)
N <- cbind(nl,nh) 


#people changing group risk
m <- 1

#mixing assortativity between risk groups:
epsilon <- 0.5

#mixing for mobility informed matrix scenarios


#loads the population mobility matrix (sum of motorized and public transports mean persons per day)
pmob <- read.csv("pmob.csv")
pmob_1 <- pmob
pmob_1 <- as.data.frame(pmob_1)
pmob_1 <- apply(pmob_1[2:27],2,as.numeric)
pmob_1[is.na(pmob_1)] <- 0
diag(pmob_1) <- 0

# s is the scaling factor which allows to scale the population movement matrix troughout all cantons in
# order to set the proportion of contacts wihtin the same cantons in an wished range
s <- 0.0345 # mean (weighted by population sizes) of 80% of all contacts being made within the own cantons

fsigma <- function(pmob_1, s, w){
  pmob_2 <- sweep(s*pmob_1, 2, (S/k)*w, "/") #this divides each col by the vector w which corresponds to the population weight
  pmob_3 <- c()
  for (i in 1:k) {
    pmob_3[i] <- 1 - (sum(pmob_2[,i]))
  }
  diag(pmob_2) <- pmob_3 #put pmob3 as diagnonal in the pmob2 matrix
  
  return(pmob_2)
}

sigma <- fsigma(pmob_1, s, w)

sigma <- t(sigma)

diag(sigma) #this represents the percentage of contacts made within each ind. cantons
# weighted.mean(diag(sigma), w)


#########################
######################### HPV-specific PARAMETERS
#########################


#rate of people entering and exiting the group / time period observed
## 1/mu = time period
mu <- 0.143

#immunity duration (1/omega= duration for becoming susceptible again)
omega <- 0.024

#infection clearance rate (1/gamma = duration of the infectiousness)
gamma <- 0.542


#Transmission probability "matrix"
beta <- 0.8


#initial params
I0_1 <- rep(0.01, k)
I0_2 <- rep(0.1, k)

I0 <- as.matrix(cbind(I0_1,I0_2))





