# 10.02.2017
# 
# 02_Models

require("deSolve")

#### Vector of Initial values, depending number of cantons in the model (k) (I0, see N_01_params)
#women low
ISlf <- c(); IIlf <- c() ; IRlf <- c(); IVlf <- c()
#women high
IShf <- c(); IIhf <- c(); IRhf <- c(); IVhf <- c()
#men low
ISlm <- c(); IIlm <- c(); IRlm <- c(); IVlm <- c()
#men high
IShm <- c(); IIhm <- c(); IRhm <- c(); IVhm <- c()

inits <- function(k){
  for(i in 1:k){
    ISlf[i] <- N[i,1]-I0[i,1]*N[i,1]
    IIlf[i] <- I0[i,1]*N[i,1]
    IRlf[i] <- 0
    IVlf[i] <- 0  
    
    IShf[i] <- N[i,2]-I0[i,2]*N[i,2]
    IIhf[i] <- I0[i,2]*N[i,2]
    IRhf[i] <- 0
    IVhf[i] <- 0
    
    ISlm[i] <- N[i,1]-I0[i,1]*N[i,1]
    IIlm[i] <- I0[i,1]*N[i,1]
    IRlm[i] <- 0
    IVlm[i] <- 0  
    
    IShm[i] <- N[i,2]-I0[i,2]*N[i,2]
    IIhm[i] <- I0[i,2]*N[i,2]
    IRhm[i] <- 0
    IVhm[i] <- 0
  }
  return(list(c(ISlf,IIlf,IRlf,IVlf, IShf,IIhf,IRhf,IVhf, ISlm, IIlm,IRlm,IVlm, IShm,IIhm,IRhm,IVhm)))
}

#########################
######################### FOR SCENARIO 1&2 , function for rho_kk'rr'
#########################



#matrix with number of available contacts of the target group on the total contacts
#the equation is separated in two parts, first part for epsilon1, mixing between the high and low risk groups and second part
# for epsilon1, micing between the cantons. The two parts will be put together at the end 
#####  #########################
#- 1B
#sort of Kronecker delta for mutliple cantons and two risk groups
##### This only works if R=2!!!!
# - 2B*
#sort of Kronecker delta for mutliple cantons and two risk groups
#kronecker delta for epsilon1 

rho_ki_f <- function(epsilon, epsilon1, k, C){
  PA1_1 <- matrix(0,nrow=k,ncol=r)
  PA1_2 <- c()
  PAB1_f <- matrix(0,nrow=k*r,ncol=k*r)
  PA2_1 <- c()
  PA2_1a <- matrix(NA,nrow=k,ncol=r)
  PA2_2 <- c()
  PAB2_f <- matrix(0,nrow=k*r,ncol=k*r)
  #- final matrix 1 and 2 toghether
  rho_ki <- matrix(0,nrow=k*r,ncol=k*r)
  P1B_k_delta <- matrix(rep(c( rep(c(epsilon,0),k), rep(c(0,epsilon),k) ), k), ncol=k*r, nrow=k*r)
  P2B_k_delta<- kronecker(diag(epsilon1, k), matrix(rep(1,4), ncol = r))
  for (i in 1:k){
    for (j in 1:r){
      PA1_1[i,j]<- (1- epsilon)*(C[i,j]*N[i,j]/sum(C[i,1]*N[i,1],C[i,2]*N[i,2]))
      PA1_2 <- as.vector(t(PA1_1))
      PAB1_f <-  P1B_k_delta + t(matrix(rep(PA1_2, k*r), ncol=k*r, nrow=k*r)) #add part b (for risk groups)
      
      PA2_1[i]<- (1- epsilon1)*sum(C[i,1]*N[i,1],C[i,2]*N[i,2])/sum(C*N) 
      PA2_1a<- as.matrix(rep(t(as.vector(PA2_1)), each=r), nrow=k, ncol=r)
      PA2_2 <- as.vector(t(PA2_1a))  
      PAB2_f <- P2B_k_delta + t(matrix(rep(PA2_2, k*r), ncol=k*r, nrow=k*r))   
      #multiply the two parts of the equations together
      rho_ki <- PAB1_f*PAB2_f 
    }
  }
  return(as.matrix(rho_ki))
}

#########################
######################### FOR SCENARIO 3 , function for rho_rr'
#########################
rho_f <- function(epsilon){
  rho_m <- matrix(0, nrow=r, ncol=r)
  rholl <- epsilon*1 + (1- epsilon) * (C[1,1]*n[1,1])/((C[1,1]*n[1,1])+(C[1,2]*n[1,2]))
  rholh <- epsilon*0 + (1- epsilon) * (C[1,2]*n[1,2])/((C[1,1]*n[1,1])+(C[1,2]*n[1,2]))
  rhohl <- epsilon*0 + (1- epsilon) * (C[1,1]*n[1,1])/((C[1,1]*n[1,1])+(C[1,2]*n[1,2]))
  rhohh <- epsilon*1 + (1- epsilon) * (C[1,2]*n[1,2])/((C[1,1]*n[1,1])+(C[1,2]*n[1,2]))
  
  rho_m <- t(matrix(c(rholl, rholh, rhohl, rhohh), nrow=2, ncol=2))
  return(as.matrix(rho_m))
} 




############################################################################################################
############ Meta-pop. model 1: for with assortative or random mixing between cantons (scenario 1 & 2)
############################################################################################################


STI_1<- function(t, X, params, verbose=FALSE)
{
  with(as.list(params), {
    dX <- numeric()
    pt <- c()
    epsilon1 <- intercantmix
    rho_ki<- rho_ki_f(epsilon, epsilon1, k,C)
    for (i in 1:k){ 
      Blf<- c()
      Bhf<- c()
      Blm<- c()
      Bhm<- c()
      for (z in 1:k){
        #infectious acquirance for women (from men)
        # low risk group
        Blf[z]<-(sum(beta*rho_ki[(i*2)-1,(z*2)-1]* X[(z*16)-06]/N[z,1]
                     ,beta*rho_ki[(i*2)-1, (z*2)]* X[(z*16)-02]/N[z,2]))
        Blf_a <- sum(Blf)
        # high risk group
        Bhf[z]<-(sum(beta*rho_ki[(i*2),(z*2)] * X[(z*16)-02] /N[z,2]
                     ,beta*rho_ki[(i*2), (z*2)-1] * X[(z*16)-06]/N[z,1]))
        Bhf_a <- sum(Bhf)
        #infectious acquirance for men (from women)
        # low risk group
        Blm[z]<-(sum(beta*rho_ki[(i*2)-1,(z*2)-1]* X[(z*16)-14]/N[z,1]
                     ,beta*rho_ki[(i*2)-1, (z*2)]* X[(z*16)-10]/N[z,2]))
        Blm_a <- sum(Blm)
        # high risk group
        Bhm[z]<-(sum(beta*rho_ki[(i*2),(z*2)] * X[(z*16)-10] /N[z,2]
                     ,beta*rho_ki[(i*2), (z*2)-1] * X[(z*16)-14]/N[z,1]))
        Bhm_a <- sum(Bhm)     
      }       
      #SIRS two level equations
      pt[i]<- ifelse(t>200, p[i], 0)
      #################### FOR WOMEN 
      ####LOW RISK GROUP
      #S
      dX[(i*16)-15] <- (((1-pt[i])*mu*N[i,1])  
                        - C[i,1]*X[(i*16)-15] * Blf_a
                        - (mu*X[(i*16)-15]) + (omega*X[(i*16)-13])
                        - m*X[(i*16)-15] + m*n[i,1]*(X[(i*16)-15]+X[(i*16)-11]))
      #I
      dX[(i*16)-14] <-  (C[i,1]*X[(i*16)-15]* Blf_a
                         - mu*X[(i*16)-14]- gamma*X[(i*16)-14] - m*X[(i*16)-14] + m*n[i,1]*(X[(i*16)-14]+X[(i*16)-10]))
      #R
      dX[(i*16)-13] <-  (gamma*X[(i*16)-14] - mu*X[(i*16)-13] - (omega*X[(i*16)-13])- m*X[(i*16)-13] + m*n[i,1]*(X[(i*16)-13]+X[(i*16)-09]))
      #V
      dX[(i*16)-12] <-  (pt[i]*mu*N[i,1] - mu*X[(i*16)-12] - m*X[(i*16)-12] + m*n[i,1]*(X[(i*16)-12]+X[(i*16)-08]))
      ####High RISK GROUP
      #S
      dX[(i*16)-11] <- (((1-pt[i])*mu*N[i,2])  
                        - C[i,2]*X[(i*16)-11] * Bhf_a
                        - (mu*X[(i*16)-11]) + (omega*X[(i*16)-1])
                        - m*X[(i*16)-11] + m*n[i,2]*(X[(i*16)-15]+X[(i*16)-11]))
      #I
      dX[(i*16)-10] <-  (C[i,2]*X[(i*16)-11] * Bhf_a
                         - mu*X[(i*16)-10]- gamma*X[(i*16)-10] - m*X[(i*16)-10] + m*n[i,2]*(X[(i*16)-14]+X[(i*16)-10]))
      #R
      dX[(i*16)-09] <-  (gamma*X[(i*16)-10] - mu*X[(i*16)-09] - (omega*X[(i*16)-09]) - m*X[(i*16)-09] + m*n[i,2]*(X[(i*16)-13]+X[(i*16)-09]))
      #V
      dX[(i*16)-08] <-  (pt[i]*mu*N[i,2] - mu*X[(i*16)-08] - m*X[(i*16)-08] + m*n[i,2]*(X[(i*16)-12]+X[(i*16)-08]))
      
      #################### FOR MEN 
      ####LOW RISK GROUP
      #S
      dX[(i*16)-07] <- (mu*N[i,1]  
                        - C[i,1]*X[(i*16)-07] * Blm_a
                        - (mu*X[(i*16)-07]) + (omega*X[(i*16)-05])
                        - m*X[(i*16)-07] + m*n[i,1]*(X[(i*16)-07]+X[(i*16)-03]))
      #I
      dX[(i*16)-06] <-  (C[i,1]*X[(i*16)-07]* Blm_a
                         - mu*X[(i*16)-06]- gamma*X[(i*16)-06] - m*X[(i*16)-06] + m*n[i,1]*(X[(i*16)-06]+X[(i*16)-02]))
      #R
      dX[(i*16)-05] <-  (gamma*X[(i*16)-06] - mu*X[(i*16)-05] - (omega*X[(i*16)-05])- m*X[(i*16)-05] + m*n[i,1]*(X[(i*16)-05]+X[(i*16)-01]))
      #V
      dX[(i*16)-04] <-  (0*N[i,1] - mu*X[(i*16)-04] - m*X[(i*16)-04] + m*n[i,1]*(X[(i*16)-04]+X[(i*16)-00]))
      ####High RISK GROUP
      #S
      dX[(i*16)-03] <- (mu*N[i,2]  
                        - C[i,2]*X[(i*16)-03] * Bhm_a
                        - (mu*X[(i*16)-03]) + (omega*X[(i*16)-01])
                        - m*X[(i*16)-03] + m*n[i,2]*(X[(i*16)-07]+X[(i*16)-03]))
      #I
      dX[(i*16)-02] <-  (C[i,2]*X[(i*16)-03] * Bhm_a
                         - mu*X[(i*16)-02]- gamma*X[(i*16)-02] - m*X[(i*16)-02] + m*n[i,2]*(X[(i*16)-06]+X[(i*16)-02]))
      #R
      dX[(i*16)-01] <-  (gamma*X[(i*16)-02] - mu*X[(i*16)-01] - (omega*X[(i*16)-01]) - m*X[(i*16)-01] + m*n[i,2]*(X[(i*16)-05]+X[(i*16)-01]))
      #V
      dX[(i*16)-00] <-  (0*N[i,2] - mu*X[(i*16)-00] - m*X[(i*16)-00] + m*n[i,2]*(X[(i*16)-04]+X[(i*16)-00]))
    }
    list(dX)                  
  })
}


############################################################################################################
############ Meta-pop. model 2: with mobility informed matrix mixing between cantons 3
############################################################################################################


STI_2<- function(t, X, params, verbose=FALSE)
{ 
  with(as.list(params), {
    dX <- numeric()
    pt <- c()
    rho <- rho_f(epsilon)
    # s <- intercantmix
    sigma <- t(fsigma(pmob_1, s, w))   
    for (i in 1:k){ 
      lambdal_ll1f <- c();lambdal_ll2f <- c();lambdal_lh1f <- c();lambdal_lh2f <- c()
      lambdah_hh1f <- c();lambdah_hh2f <- c();lambdah_hl1f <- c();lambdah_hl2f <- c()
      lambdalf <- c();lambdahf <- c()
      
      lambdal_ll1m <- c();lambdal_ll2m <- c();lambdal_lh1m <- c();lambdal_lh2m <- c()
      lambdah_hh1m <- c();lambdah_hh2m <- c();lambdah_hl1m <- c();lambdah_hl2m <- c()
      lambdalm <- c();lambdahm <- c()
      for (z in 1:k){  
        #infectious acquirance for women (from men)
        lambdal_ll1f[z] <- (sigma[i,z] * beta * rho[1,1] * X[(z*16)-06]/N[z,1] )
        lambdal_ll2f    <- sum(lambdal_ll1f) # for low with low
        
        lambdal_lh1f[z] <- (sigma[i,z] * beta * rho[1,2] * X[(z*16)-02]/N[z,2] )
        lambdal_lh2f <-  sum(lambdal_lh1f) #for low with high   
        
        lambdalf <-  C[i,1]* (lambdal_ll2f + lambdal_lh2f) 
        
        lambdah_hh1f[z] <- (sigma[i,z] * beta * rho[2,2] * X[(z*16)-02]/N[z,2] )
        lambdah_hh2f    <- sum(lambdah_hh1f)
        
        lambdah_hl1f[z] <- (sigma[i,z] * beta * rho[2,1] * X[(z*16)-06]/N[z,1] ) 
        lambdah_hl2f    <- sum(lambdah_hl1f)
        
        lambdahf <-  C[i,2]* (lambdah_hh2f+lambdah_hl2f)
        
        #infectious acquirance for men (from women)
        lambdal_ll1m[z] <- (sigma[i,z] * beta * rho[1,1] * X[(z*16)-14]/N[z,1] )
        lambdal_ll2m    <- sum(lambdal_ll1m) # for low with low
        
        lambdal_lh1m[z] <- (sigma[i,z] * beta * rho[1,2] * X[(z*16)-10]/N[z,2] )
        lambdal_lh2m <-  sum(lambdal_lh1m) #for low with high   
        
        lambdalm <-  C[i,1]* (lambdal_ll2m + lambdal_lh2m) 
        
        lambdah_hh1m[z] <- (sigma[i,z] * beta * rho[2,2] * X[(z*16)-10]/N[z,2] )
        lambdah_hh2m    <- sum(lambdah_hh1m)
        
        lambdah_hl1m[z] <- (sigma[i,z] * beta* rho[2,1] * X[(z*16)-14]/N[z,1] ) 
        lambdah_hl2m    <- sum(lambdah_hl1m)
        
        lambdahm <-  C[i,2]* (lambdah_hh2m+lambdah_hl2m)
      }       
      #SIRS two level equations
      pt[i]<- ifelse(t>200, p[i], 0)
      #  pt[i] <- p[i]
      #################### FOR WOMEN 
      ####LOW RISK GROUP
      #S
      dX[(i*16)-15] <- (((1-pt[i])*mu*N[i,1])  
                        - X[(i*16)-15] * lambdalf
                        - (mu*X[(i*16)-15]) + (omega*X[(i*16)-13])
                        - m*X[(i*16)-15] + m*n[i,1]*(X[(i*16)-15]+X[(i*16)-11]))
      #I
      dX[(i*16)-14] <-  (X[(i*16)-15]* lambdalf
                         - mu*X[(i*16)-14]- gamma*X[(i*16)-14] - m*X[(i*16)-14] + m*n[i,1]*(X[(i*16)-14]+X[(i*16)-10]))
      #R
      dX[(i*16)-13] <-  (gamma*X[(i*16)-14] - mu*X[(i*16)-13] - (omega*X[(i*16)-13])- m*X[(i*16)-13] + m*n[i,1]*(X[(i*16)-13]+X[(i*16)-09]))
      #V
      dX[(i*16)-12] <-  (pt[i]*mu*N[i,1] - mu*X[(i*16)-12] - m*X[(i*16)-12] + m*n[i,1]*(X[(i*16)-12]+X[(i*16)-08]))
      ####High RISK GROUP
      #S
      dX[(i*16)-11] <- (((1-pt[i])*mu*N[i,2])  
                        - X[(i*16)-11] * lambdahf
                        - (mu*X[(i*16)-11]) + (omega*X[(i*16)-1])
                        - m*X[(i*16)-11] + m*n[i,2]*(X[(i*16)-15]+X[(i*16)-11]))
      #I
      dX[(i*16)-10] <-  (X[(i*16)-11] * lambdahf
                         - mu*X[(i*16)-10]- gamma*X[(i*16)-10] - m*X[(i*16)-10] + m*n[i,2]*(X[(i*16)-14]+X[(i*16)-10]))
      #R
      dX[(i*16)-09] <-  (gamma*X[(i*16)-10] - mu*X[(i*16)-09] - (omega*X[(i*16)-09]) - m*X[(i*16)-09] + m*n[i,2]*(X[(i*16)-13]+X[(i*16)-09]))
      #V
      dX[(i*16)-08] <-  (pt[i]*mu*N[i,2] - mu*X[(i*16)-08] - m*X[(i*16)-08] + m*n[i,2]*(X[(i*16)-12]+X[(i*16)-08]))
      
      #################### FOR MEN 
      ####LOW RISK GROUP
      #S
      dX[(i*16)-07] <- (mu*N[i,1]  
                        - X[(i*16)-07] * lambdalm
                        - (mu*X[(i*16)-07]) + (omega*X[(i*16)-05])
                        - m*X[(i*16)-07] + m*n[i,1]*(X[(i*16)-07]+X[(i*16)-03]))
      #I
      dX[(i*16)-06] <-  (X[(i*16)-07]* lambdalm
                         - mu*X[(i*16)-06]- gamma*X[(i*16)-06] - m*X[(i*16)-06] + m*n[i,1]*(X[(i*16)-06]+X[(i*16)-02]))
      #R
      dX[(i*16)-05] <-  (gamma*X[(i*16)-06] - mu*X[(i*16)-05] - (omega*X[(i*16)-05])- m*X[(i*16)-05] + m*n[i,1]*(X[(i*16)-05]+X[(i*16)-01]))
      #V
      dX[(i*16)-04] <-  (0*mu*N[i,1]- mu*X[(i*16)-04] - m*X[(i*16)-04] + m*n[i,1]*(X[(i*16)-04]+X[(i*16)-00]))
      ####High RISK GROUP
      #S
      dX[(i*16)-03] <- (mu*N[i,2]  
                        - X[(i*16)-03] * lambdahm
                        - (mu*X[(i*16)-03]) + (omega*X[(i*16)-01])
                        - m*X[(i*16)-03] + m*n[i,2]*(X[(i*16)-07]+X[(i*16)-03]))
      #I
      dX[(i*16)-02] <-  (X[(i*16)-03] * lambdahm
                         - mu*X[(i*16)-02]- gamma*X[(i*16)-02] - m*X[(i*16)-02] + m*n[i,2]*(X[(i*16)-06]+X[(i*16)-02]))
      #R
      dX[(i*16)-01] <-  (gamma*X[(i*16)-02] - mu*X[(i*16)-01] - (omega*X[(i*16)-01]) - m*X[(i*16)-01] + m*n[i,2]*(X[(i*16)-05]+X[(i*16)-01]))
      #V
      dX[(i*16)-00] <-  ( 0*mu*N[i,2] - mu*X[(i*16)-00] - m*X[(i*16)-00] + m*n[i,2]*(X[(i*16)-04]+X[(i*16)-00]))
    }
    list(dX)                  
  })
}


############################################################################################################
############ Meta-pop. model 3: Different assumption for vaccination (fig.A8)
############################################################################################################

# Here a proportion p of the susceptible are getting vaccinated at each time step.

STI_3 <- function(t, X, params, verbose=FALSE)
{
  
  with(as.list(params), {
    dX <- numeric()
    pt <- c()
    rho_ki<- rho_ki_f(epsilon, epsilon1, k,C)
    for (i in 1:k){ 
      Blf<- c()
      Bhf<- c()
      Blm<- c()
      Bhm<- c()
      for (z in 1:k){
        #infectious acquirance for women (from men)
        # low risk group
        Blf[z]<-(sum(beta*rho_ki[(i*2)-1,(z*2)-1]* X[(z*16)-06]/N[z,1]
                     ,beta*rho_ki[(i*2)-1, (z*2)]* X[(z*16)-02]/N[z,2]))
        Blf_a <- sum(Blf)
        # high risk group
        Bhf[z]<-(sum(beta*rho_ki[(i*2),(z*2)] * X[(z*16)-02] /N[z,2]
                     ,beta*rho_ki[(i*2), (z*2)-1] * X[(z*16)-06]/N[z,1]))
        Bhf_a <- sum(Bhf)
        
        #infectious acquirance for men (from women)
        # low risk group
        Blm[z]<-(sum(beta*rho_ki[(i*2)-1,(z*2)-1]* X[(z*16)-14]/N[z,1]
                     ,beta*rho_ki[(i*2)-1, (z*2)]* X[(z*16)-10]/N[z,2]))
        Blm_a <- sum(Blm)
        # high risk group
        Bhm[z]<-(sum(beta*rho_ki[(i*2),(z*2)] * X[(z*16)-10] /N[z,2]
                     ,beta*rho_ki[(i*2), (z*2)-1] * X[(z*16)-14]/N[z,1]))
        Bhm_a <- sum(Bhm)     
      }       
      #SIRS two level equations
      pt[i]<- ifelse(t>200, mp[i], 0)
      #  pt[i] <- p[i]
      #################### FOR WOMEN 
      ####LOW RISK GROUP
      #S
      dX[(i*16)-15] <- (mu*N[i,1] -  pt[i]*X[(i*16)-15]
                        - C[i,1]*X[(i*16)-15] * Blf_a
                        - (mu*X[(i*16)-15]) + (omega*X[(i*16)-13])
                        - m*X[(i*16)-15] + m*n[i,1]*(X[(i*16)-15]+X[(i*16)-11]))
      #I
      dX[(i*16)-14] <-  (C[i,1]*X[(i*16)-15]* Blf_a
                         - mu*X[(i*16)-14]- gamma*X[(i*16)-14] - m*X[(i*16)-14] + m*n[i,1]*(X[(i*16)-14]+X[(i*16)-10]))
      #R
      dX[(i*16)-13] <-  (gamma*X[(i*16)-14] - mu*X[(i*16)-13] - (omega*X[(i*16)-13])- m*X[(i*16)-13] + m*n[i,1]*(X[(i*16)-13]+X[(i*16)-09]))
      #V
      dX[(i*16)-12] <-  (pt[i]*X[(i*16)-15] - mu*X[(i*16)-12] - m*X[(i*16)-12] + m*n[i,1]*(X[(i*16)-12]+X[(i*16)-08]))
      ####High RISK GROUP
      #S
      dX[(i*16)-11] <- (mu*N[i,2]  - pt[i]*X[(i*16)-11]
                        - C[i,2]*X[(i*16)-11] * Bhf_a
                        - (mu*X[(i*16)-11]) + (omega*X[(i*16)-1])
                        - m*X[(i*16)-11] + m*n[i,2]*(X[(i*16)-15]+X[(i*16)-11]))
      #I
      dX[(i*16)-10] <-  (C[i,2]*X[(i*16)-11] * Bhf_a
                         - mu*X[(i*16)-10]- gamma*X[(i*16)-10] - m*X[(i*16)-10] + m*n[i,2]*(X[(i*16)-14]+X[(i*16)-10]))
      #R
      dX[(i*16)-09] <-  (gamma*X[(i*16)-10] - mu*X[(i*16)-09] - (omega*X[(i*16)-09]) - m*X[(i*16)-09] + m*n[i,2]*(X[(i*16)-13]+X[(i*16)-09]))
      #V
      dX[(i*16)-08] <-  (pt[i]*X[(i*16)-11] - mu*X[(i*16)-08] - m*X[(i*16)-08] + m*n[i,2]*(X[(i*16)-12]+X[(i*16)-08]))
      
      #################### FOR MEN 
      ####LOW RISK GROUP
      #S
      dX[(i*16)-07] <- (mu*N[i,1]  
                        - C[i,1]*X[(i*16)-07] * Blm_a
                        - (mu*X[(i*16)-07]) + (omega*X[(i*16)-05])
                        - m*X[(i*16)-07] + m*n[i,1]*(X[(i*16)-07]+X[(i*16)-03]))
      #I
      dX[(i*16)-06] <-  (C[i,1]*X[(i*16)-07]* Blm_a
                         - mu*X[(i*16)-06]- gamma*X[(i*16)-06] - m*X[(i*16)-06] + m*n[i,1]*(X[(i*16)-06]+X[(i*16)-02]))
      #R
      dX[(i*16)-05] <-  (gamma*X[(i*16)-06] - mu*X[(i*16)-05] - (omega*X[(i*16)-05])- m*X[(i*16)-05] + m*n[i,1]*(X[(i*16)-05]+X[(i*16)-01]))
      #V
      dX[(i*16)-04] <-  (0 - mu*X[(i*16)-04] - m*X[(i*16)-04] + m*n[i,1]*(X[(i*16)-04]+X[(i*16)-00]))
      ####High RISK GROUP
      #S
      dX[(i*16)-03] <- (mu*N[i,2]  
                        - C[i,2]*X[(i*16)-03] * Bhm_a
                        - (mu*X[(i*16)-03]) + (omega*X[(i*16)-01])
                        - m*X[(i*16)-03] + m*n[i,2]*(X[(i*16)-07]+X[(i*16)-03]))
      #I
      dX[(i*16)-02] <-  (C[i,2]*X[(i*16)-03] * Bhm_a
                         - mu*X[(i*16)-02]- gamma*X[(i*16)-02] - m*X[(i*16)-02] + m*n[i,2]*(X[(i*16)-06]+X[(i*16)-02]))
      #R
      dX[(i*16)-01] <-  (gamma*X[(i*16)-02] - mu*X[(i*16)-01] - (omega*X[(i*16)-01]) - m*X[(i*16)-01] + m*n[i,2]*(X[(i*16)-05]+X[(i*16)-01]))
      #V
      dX[(i*16)-00] <-  (0*N[i,2] - mu*X[(i*16)-00] - m*X[(i*16)-00] + m*n[i,2]*(X[(i*16)-04]+X[(i*16)-00]))
    }
    list(dX)                  
  })
}






#############
############ Output 1: function that outputs prev in women in single cantons and in overall pop.
#############



#function that makes a vector out of two vectors/4
outp <- function(sim, k ){
  kprev <- c()
  kprev_b <- c()
  kprev1 <- matrix(0,ncol=k, nrow=sim$t)
  kprev1_b <- matrix(0,ncol=k, nrow=sim$t)
  ### only women
  for(i in 1:k){
    kprev[i]   <- ( (sim[(i*16)-14+1]+ sim[(i*16)-10+1] ) )/w[i] #combines low and high infectious compartement ONLY WOMEN
    kprev_b[i] <-  (sim[(i*16)-14+1]+ sim[(i*16)-10+1] ) 
  }
  
  kprev1<- do.call(cbind, kprev)   #proportion infectious in each canton (N=1*k)
  kprev1_b <- do.call(cbind, kprev_b) 
  kprev1_a <- rowSums(kprev1_b)/k  # proportion in total population 
  #add time in the matrix and mean prevalence het as first column
  kprev2<- cbind(sim$time,  kprev1, kprev1_a)
  #kprev2<- cbind(simulation1$time,  kprev1,  kprev1_a)
  kprev2 <- as.data.frame(kprev2)
  return(kprev2)
}



#############
############ Output 2: function that outputs prev in women, overall, in a two sub-populations model 
#                       (varying vaccination coverages between 0% and 100% in both populations)
#############



prev.2k.diffvaccf2 <-  function(diff.eps, diff.vacc, timeseq){
  outp <- array(NA,  dim=c( length(diff.eps), length(diff.vacc[,1]), 1) )
  pb <- txtProgressBar(0, length(diff.vacc[,1]), style=3)
  di <-  seq(1,length(diff.vacc[,1]),1)
 
   for(y in 1:length(diff.eps)){
    for(z in 1:length(diff.vacc[,1])){
      
      p <- diff.vacc[z,]
      p <- as.vector(p)
      params <- c(p=p, intercantmix=diff.eps[y])
      time <- seq(0,timeseq, 1)
      simulation_A <- as.matrix(ode(init, time, STI_1, parms=params))
      kprev_b <- c()
      kprev1_global <- c()
      ### only women
      for(x in 1:k){
        kprev_b[x] <-  simulation_A[timeseq,(x*16)-14+1]+ simulation_A[timeseq,(x*16)-10+1] #combines low and high infectious compartement ONLY WOMEN
      }
      # vector with prevalence in the global pop
      kprev1_global <- sum(kprev_b)/k
      outp[y,z,] <- kprev1_global
      setTxtProgressBar(pb, z)
    }}
  return( outp)
}



#############
############ Output 3: functions that outputs the number of canton achieving different thresholds of RR reduction
#                       using STI_1, over different assortativity indexes
#############

#first function which gives prevalence overall and in each ind. canton
# for different epsilon1 (different assortivity indexes)
# takes time
prev.free.k_fn <- function(mix){
  outp <- array(dim=c(length(mix), length(time)-200, k ))
  pb <- txtProgressBar(0, length(mix), style=3)
  for(i in 1:length(mix)){  
    params <- c(intercantmix=mix[i])
    simulation_A <- as.data.frame(ode(init, time, STI_1, params))
    simulationA1 <- simulation_A[simulation_A$time>=200,]
    kprev <- c()
    kprev_b <- c()
    kprev1 <- matrix(0,ncol=k, nrow=simulationA1$t)
    kprev1_b <- matrix(0,ncol=k, nrow=simulationA1$t)
    kprev1_global <- c()
    ### only women
    for(j in 1:k){
      kprev[j]   <- ( (simulationA1[(j*16)-14+1]+ simulationA1[(j*16)-10+1] ) )/w[j] #combines low and high infectious compartement ONLY WOMEN
      kprev_b[j] <-  (simulationA1[ (j*16)-14+1 ]+ simulationA1[(j*16)-10+1] )
    }
    # matrix with prevalences in ind. cantons
    kprev1<- do.call(cbind, kprev)
    # vector with prevalence in Switzerland (weighted by pop.size)
    kprev_b1<- do.call(cbind, kprev_b)
    kprev1_global <- rowSums( kprev_b1 )/k
    outp[i,, ] <- kprev1
    setTxtProgressBar(pb, i)
  }
  return( outp)
}

## This function is based on the matrix output from function: prev.free.k_fn
# This function gives: number of cantons achieving certain threshold of RR reduction AND
# mean prevalence overall (national) over different epsilon1 (assortativity indexes for sexual mixing between cantons)
# It needs the matrix output from the previous function and the number of years after vaccination onset
# and the different thresholds we are interested in. Rapid output function
prev.free.k_fn2 <- function(arr, y.after.v , diff.tresh ) {
  outp2 <- array(dim=c(length(diff.eps), length(diff.tresh),2 )) 
  pb <- txtProgressBar(0, length(diff.eps), style=3)
  for(i in 1: length(diff.eps)){
    for(j in 1: length(diff.tresh)){ 
      outp2[i,j,] <- c( length(which( arr[i,y.after.v, ]/ max(arr[i, , ]) < (1- diff.tresh[j]) ) ) , #nb cantons achieving rr red according to threshold
                        ( sum(arr[i,y.after.v,] * w ) ) / k ) # here is the mean prev in Switz. 
      
      setTxtProgressBar(pb, i)
    } }
  return(outp2)
}







