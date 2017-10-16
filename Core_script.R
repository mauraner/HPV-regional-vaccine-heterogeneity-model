# 10.02.2017
# 
# Modeling the consequences of regional heterogeneity in human papillomavirus (HPV)
# vaccination uptake on transmission in Switzerland
# 
# Maurane Riesen, Victor Garcia, Nicola Low, Christian L. Althaus
# 
# R script


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 1.

# System with 1 population. to make the following calculations:

# - 1.a: pre-vaccine prevalences with Swiss (CH) or British (UK) sexual activity data
# - 1.b: sensitivity analysis on different per partnership virus transmissibilities (with CH or UK parameter), fig. A.7
# - 1.c: Reduction in prevalence under different vaccination coverages after different times (years after vaccination onset):
#   comparision with published data from post-vaccination prevalence studies (Drolet et al, 2015, Lancet)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

######
### - 1.a: pre-vaccine prevalences with Swiss (CH) or British (UK) sexual activity data
######

# set number of cantons
k <- 1

#run the parameters script and models over again
source("N_01_params_rev.R") # N_01_params runs the parameters of the model
source("N_02_models.R") # N_02_models runs the two main SIRV models: STI_1 for scenario 1&2 and STI_2 for scenario 3

# run the initial vector with the right size (right k)
init <- as.vector(matrix(unlist(inits(k)), ncol = k, byrow = TRUE))

# p <- 0 # set vaccination to 0
p <- mp

# p <- p/2
#simulation time
time <- seq(0, 400, 1)

# choose the source for number of sexual contacts (UK or CH), here CH
C <- C_ch
# C <- C_uk

# set mixing between population (not relevant here, can be any value between 0 and 1)
intercantmix <- 1



params <- c(mu=mu, omega=omega, gamma=gamma, 
            epsilon=epsilon, intercantmix=intercantmix, m=m, p=p, C=C, N=N)

simulation_1k<- as.data.frame(ode(init, time, STI_1, parms=params))

mat <- outp(simulation_1k, k)

#prevalence after 200 years:
mat[216,]


#calculating R0

beta_1 <- matrix(rep(beta, (k*r)*(k*r)), ncol=k*r, nrow=k*r)

#for 2 risk groups
C_1 <- matrix(rep(cbind(t(C),t(C)), k), ncol=k*r, nrow=k*r)
#population division matrix
pd <- matrix(c(1, nh/nl, nl/nh, 1), ncol=2, nrow=2)
rho_ki<- rho_ki_f(epsilon, intercantmix, k,C)
#New infections
f <- (beta_1*C_1*rho_ki*pd)
#Time in compartment
V <-  matrix(c(mu+gamma+m*nh, -m*nh, -m*nl, mu+gamma+m*nl ), ncol=2, nrow=2)
G <- f %*% solve(V)
R0 <- max(eigen(G)$values)
R0
#vaccination threshold
1-1/R0
#single sex vaccination threshold
1-1/R0^2

######
### - 1.b: sensitivity analysis on different per partnership virus transmissibilities (with CH or UK parameter), fig. A.7
######

# this function gives the prevalence in a system with 1 canton/population, for women at the end of the simulation time
# over the different beta. 
p <- 0

prev <- function(beta, C){
  prev_overb <- c()
  for (i in 1:length(beta)) {
    params <- c( beta= beta[i], C=C)
    y <- as.data.frame(lsoda(init, time, STI_1, parms=params))
    ########## total prev of HPV16 (LR group and HR group in women)
    prev_overb[i] <-  y[,3][[length(time)]]+y[,7][[length(time)]]
  }
  return(prev_overb)
}

#Make a vector to see how it behaves with different beta
beta <- seq(0,1,0.05)

C <- C_ch
prevCH <- prev(beta, C)
C <- C_uk
prevUK<- prev(beta, C)


# pdf("sensitivity_betaCHUK_newgamma.pdf", width=6, height=6)
plot(prevCH~ beta, type="b", ylim=c(0,0.1), frame=FALSE, 
     ylab="HPV-16 prevalence", xlab= expression(paste("per partnership transmission probability, ", beta)),lwd=2)
lines(prevUK~ beta, type="b", col="blue", lwd=2)


text(beta, prevCH, labels=ifelse(round(prevCH,3)==0, NA, round(prevCH,3) ), pos=2, cex=0.6)
text(beta, prevUK, labels=ifelse(round(prevUK,3)==0, NA, round(prevUK,3) ), pos=2, cex=0.6, col="blue")

legend(0, 0.1, legend=c("British model", "Swiss model"), cex=1,
       col =c('blue', 'black'),lwd = c(2,2),
       lty=c(1,1))
# dev.off()

######
### - 1.c: Reduction in prevalence under different vaccination coverages after different times (years after vaccination onset):
###        comparision with published data from post-vaccination prevalence studies (Drolet et al, 2015, Lancet)
######


#Here we will use the STI_3 model. In the model STI_3, a proportion p of the susceptible are getting vaccinated at each time step.
# In this way, it require less time until the proportion p is reflected in the proportion that are in the vaccinated compartement
library(plotrix)

#HPV 16-18 on 13-19 and 20-24 (13-24y.o) years old women, data from the Drolet et. al study
# c is the vaccination coverage, r is the relative risk reduction and l and u are the lower and upper confidence interval resp.
c <- c(0.34,0.58,0.62,0.77,0.88,0.89, 0.16,0.16,0.18,0.31,0.60,0.83)
rr <- c(0.50,0.47,0.39,0.38,0.04,0.32, 1.09,0.81,1.07,0.7,0.67,0.25)
l <- c(0.34,0.35,0.19,0.25,0.01,0.12, 0.75,0.48,0.74,0.42,0.61,0.17)
u <- c(0.74,0.63,0.79,0.58,0.15,0.89, 1.59,1.38,1.56,1.16,0.74,0.36)


# this function calculates the relative risk reduction over different vaccination coverage (p)


solve.STI <- function( mp, time)
{
  params <- c( mp=mp)
  init <- init
  return(as.data.frame(ode(init, time, STI_3, parms=params)))
}

prev.pop<- function( time)
{
  p <- seq(0,1,0.01)
  kprev_hom_f<- c()
  rro <- c()
  kvacc <- c()
  m1 <- matrix(0,nrow=length(p), ncol=3)
  for (i in 1:length(p)){ 
    temp <- solve.STI(p[i], time)
    kprev_hom_f[i] <- temp[length(time),((1*16)-14)+1] + temp[length(time),((1*16)-10)+1]
    kvacc[i] <- (temp[length(time),((1*16)-12)+1] + temp[length(time),((1*16)-08)+1 ] ) 
    #kvacc[i] <- p * (1-exp(-mu* length(time)-200 ))
    rro[i]<- kprev_hom_f[i]/ ( temp[190,((1*16)-14)+1] + temp[190,((1*16)-10)+1] )
  }
  #   m1 <- cbind(p, kprev_hom_f, kvacc ) #gives prev
  m1 <- cbind(p, rro, kvacc ) #gives relative risk reduction
  # m1 <- ifelse(m1<tr,0,m1)
  return(m1)
}



# give required parameters (again)
k <- 1
source("N_01_params_rev.R")
C <- C_ch

mp <- NULL
p <- NULL
epsilon1 <- 1

# number of years post-vaccination, given that vaccination starts at year 200
time2y <- seq(0,202,1)
time4y <- seq(0,204,1)
time10y <- seq(0,210,1)
time100y <- seq(0,300,1)


#Results for different numbers of years after vaccination onset
m2y <- prev.pop(time2y)
m4y <- prev.pop(time4y)
m10y <- prev.pop(time10y)
m100y <- prev.pop(time100y)


# pdf("rrOnVaccination_lancet_CH_NEW6_.pdf", width=6, height=6)
par(mar = c(5,5,3,1))

plotCI(c,rr,ui=u,li=l,xlim=c(0,1),ylim=c(0,1.8),xlab="vaccination coverage",ylab="HPV-16 relative risk"
       ,frame=FALSE, pch=19,  cex=1, cex.lab=1.5)
#        ,main="HPV prevalence relative risk decrease, \non vaccination coverage")
pointCI(cb, rb, ui=ub, li=lb, xlim=c(0,1), ylim=c(0,1))
legend(0.6,1.55,legend=c("2 years", "4 years","10 years", "100 years"), col=c("purple","black","grey","blue"),lty=1,
       bty="n", cex=1)
legend(0.3,1.7,legend="time after vaccination onset", col="black",bty="n", cex=1.25)

lines(m2y[,2] ~ m2y[,3] ,  type= "l", col="purple", lwd=1, ylim=c(0,1), xlim=c(0,1))

lines(m4y[,2] ~ m4y[,3] , type = 'l',  col="black", lwd=1, ylim=c(0,1), xlim=c(0,1))
# lines(m5y[,2] ~ m5y[,3] , type = 'l',  col="black", lwd=0.5, ylim=c(0,1), xlim=c(0,1))

lines(m10y[,2] ~ m10y[,3] ,  type= "l", col="grey", lwd=1, ylim=c(0,1), xlim=c(0,1))

lines(m100y[,2] ~ m100y[,3] ,  type= "l", col="blue", lwd=1, ylim=c(0,1), xlim=c(0,1))

# dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 2.

# System with 2 population. to make the following calculations (Fig.3):

# - 2.a: 2 cantons with different vaccination coverages, fully assortative.
# - 2.b: 2 cantons with different vaccination coverages, partial random mixing
# - 2.c: prevalence differences between 2.a and 2.b

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


######
### - 2.a & b & c: cantons with different vaccination coverages, fully assortative.
######

# or use fig3_ext to run on UBELIX

# set number of cantons
k <- 2

#run the parameters script and models over again
source("N_01_params_rev.R") # N_01_params runs the parameters of the model
source("N_02_models.R") # N_02_models runs the two main SIRV models: STI_1 for scenario 1&2 and STI_2 for scenario 3

# run the initial vector with the right size (right k)
init <- as.vector(matrix(unlist(inits(k)), ncol = k, byrow = TRUE))

#simulation time (here 50y after vaccination onset)
time <- seq(0, 350, 1)

# choose the source for number of sexual contacts (UK or CH), here CH
C <- C_ch
# C <- C_uk
p <- NULL

#Choose the different vaccination coverages you want for each cantons
p1 <- seq(0.0,1,0.01)
p2 <- seq(0.0,1,0.01)


p_s <- expand.grid(p1,p2)
p_s <- as.matrix(p_s)


#####Use function  prev.2k.diffvaccf which need a vector for diff.eps and a matrix for diff.vacc(=p_s)
# It results an array with following dimensions:   
#   outp <- array( dim=c(  length(diff.eps), length(diff.vacc[,1]), length(time)-200, k + 1 ) )
# Time is time after vaccination onset, and at the end (last dimension of output) is prev k1, prevk2 and prev overall
#   
timeseq <- c(351)

diff.vacc <- p_s

## choose vector for mixing between pop
diff.eps <- c(1,0.6)

#this takes time
outp <- array(NA,  dim=c( length(diff.eps), length(diff.vacc[,1]), 1) )
pb <- txtProgressBar(0, length(diff.vacc[,1]), style=3)

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


outp1 <- outp


## if ubelix was used load:
# load("outp2.Rdata")

# outp1 <- outp

#Choose the different vaccination coverages you want for each cantons
p1 <- seq(0.0,1,0.01)
p2 <- seq(0.0,1,0.01)

p_s <- expand.grid(p1,p2)
p_s <- as.matrix(p_s)


two_k_eps1[,1:2] <- round(two_k_eps1[,1:2],3)
two_k_eps06[,1:2] <- round(two_k_eps06[,1:2],3)

### Plot a & b & c

superpos <- two_k_eps1 

superpos[,3] <- two_k_eps1[,3] - two_k_eps06[,3]


two_k_eps1[,3] <- ifelse( two_k_eps1[,3] < 1e-6,  1e-6, two_k_eps1[,3] )
two_k_eps06[,3] <- ifelse( two_k_eps06[,3] < 1e-6,  1e-6, two_k_eps06[,3] )

summary(round( superpos[,3]*100,3) )

two_k_eps1[,3] <-  two_k_eps1[,3]*100 
two_k_eps06[,3] <- two_k_eps06[,3]*100 
superpos[,3] <-superpos[,3]*100 

summary(two_k_eps1)
summary(two_k_eps06)
summary(superpos)

two_k_eps1[which.min(superpos[,3]),]
two_k_eps06[which.min(superpos[,3]),]

summary((superpos[,3]>0))
h <- hist(superpos[,3])


heatcols <- heat.colors(7)
cv <- colorRampPalette(c("wheat", "orangered4"))( 7)
cv2 <- colorRampPalette(c("white", "orangered4"))( 12) 
cv4 <- colorRampPalette(c("white", "orangered4"))( 13) 
cv3 <- colorRampPalette(c("white", "orangered4"))( 14) 

p1 <- ggplot(two_k_eps1, aes(x = Var1, y = Var2))   +
  geom_tile(aes(fill=cut(otp1, breaks=c(1e-7, 0.5,1,1.5,2,2.5,3,3.5), include.lowest=TRUE ))) + 
  xlab("vaccination coverage in population 2")+
  ylab("vaccination coverage in population 1")+
  coord_equal()+
  scale_fill_manual(values= cv ,  
                    name="HPV-16 \nprevalence (%)\n" ,labels=list(expression(paste("0.0", ""-"",0.5)), expression(paste(0.5,""-"","1.0")), 
                                                                  expression(paste("1.0",""-"",1.5)) ,expression(paste(1.5,""-"","2.0")),
                                                                  expression(paste("2.0",""-"",2.5)), expression(paste(2.5,""-"","3.0")), 
                                                                  expression(paste("3.0",""-"",3.5)) ))



cv <- c("white", cv)

p1 <- ggplot(two_k_eps1, aes(x = Var1, y = Var2))   +
  geom_tile(aes(fill=cut(otp1, breaks=c(-1e-30, 1e-7, 0.5,1,1.5,2,2.5,3,3.5), include.lowest=TRUE ))) + 
  xlab("vaccination coverage in population 2")+
  ylab("vaccination coverage in population 1")+
  coord_equal()+
  scale_fill_manual(values= cv ,  
                    name="HPV-16 \nprevalence (%)\n" ,labels=list(expression(paste("0.0")),
                                                                  expression(paste("0.0", ""-"",0.5)), expression(paste(0.5,""-"","1.0")), 
                                                                  expression(paste("1.0",""-"",1.5)) ,expression(paste(1.5,""-"","2.0")),
                                                                  expression(paste("2.0",""-"",2.5)), expression(paste(2.5,""-"","3.0")), 
                                                                  expression(paste("3.0",""-"",3.5)) ))


##short term
p2 <- ggplot(two_k_eps06, aes(x = Var1, y = Var2))   +
  geom_tile(aes(fill=cut(otp1, breaks=c(0,.5,1,1.5,2,2.5,3,3.5), include.lowest=TRUE ))) + 
  xlab("vaccination coverage in population 2")+
  ylab("vaccination coverage in population 1")+
  coord_equal()+
  scale_fill_manual(values= cv ,  
                    name="HPV-16 \nprevalence (%)\n",
                    labels=list(expression(paste("0.0", ""-"",0.5)), expression(paste(0.5,""-"","1.0")), 
                                expression(paste("1.0",""-"",1.5)) ,expression(paste(1.5,""-"","2.0")),
                                expression(paste("2.0",""-"",2.5)), expression(paste(2.5,""-"","3.0")), 
                                expression(paste("3.0",""-"",3.5)) ))


##long term
p2 <- ggplot(two_k_eps06, aes(x = Var1, y = Var2))   +
  geom_tile(aes(fill=cut(otp1, breaks=c(-1e-30, 1e-7, 0.5,1,1.5,2,2.5,3,3.5), include.lowest=TRUE ))) + 
  xlab("vaccination coverage in population 2")+
  ylab("vaccination coverage in population 1")+
  coord_equal()+
  scale_fill_manual(values= cv ,  
                    name="HPV-16 \nprevalence (%)\n",
                    labels=list(expression(paste("0.0")),
                                expression(paste("0.0", ""-"",0.5)), expression(paste(0.5,""-"","1.0")), 
                                expression(paste("1.0",""-"",1.5)) ,expression(paste(1.5,""-"","2.0")),
                                expression(paste("2.0",""-"",2.5)), expression(paste(2.5,""-"","3.0")), 
                                expression(paste("3.0",""-"",3.5)) ))


cv3 <- colorRampPalette(c("white", "orangered4"))( 13) 
cv3 <- c(rgb( 230, 230, 230, alpha=255, max=255), cv3)


#replace all negative prevalence difference from vacciantion coverage higher than 0.5% with 0
temp <- superpos[which(superpos[,1]>0.45 & superpos[,2]>0.45  & superpos[,1]!= superpos[,2] ), ]
summary(temp)
superpos[which(superpos[,1]>0.45 & superpos[,2]>0.45 & superpos[,1]!= superpos[,2] ), 3] <- 2e-7


temp <- superpos[which(superpos[,1]== superpos[,2] ), ]
summary(temp)


superpos[which(superpos[,1]==1 & superpos[,2]==1 ), 3] <- 0

cv3 <- colorRampPalette(c("white", "orangered4"))( 9) 
cv3 <- c("lightgrey", cv3)
cv3 <- c(rgb( 230, 230, 230, alpha=255, max=255), cv3)

p3d <- ggplot(superpos, aes(x = Var1, y = Var2))   +
  geom_tile(aes(fill=cut(otp1, breaks=c(-1e-4,-1e-7, 1e-7, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2),
                         include.lowest=TRUE ))) +
  xlab("vaccination coverage in population 2")+
  ylab("vaccination coverage in population 1")+
  coord_equal()+
  scale_fill_manual(values= cv3 ,
                    name= expression(paste(Delta, " prevalence" )) ,
                    labels=list( expression(paste("-1", ""%.%"", 10^-4, ""-"",0 )) , expression(paste("0" )) , 
                                 expression(paste(0, -"", 0.1)) , expression(0.1-0.2), expression(0.2-0.3),
                                 expression(0.3-0.4), expression(0.4-0.5),
                                 expression(0.5-0.6), expression(0.6-0.7), expression(0.7-0.8), expression(0.8-0.9),
                                 expression(paste(0.9,""-"","1.0")), expression(paste("1.0",""-"",1.1)), expression(1.1-1.2) ) )


p1 
p2 
p3d  

require("cowplot")
lp <- list(p1 + custom6, p2 + custom6, p3d +  custom6 )

pdf("allthree2cantonsplot_after15yv_20171005.pdf", width=9, height=22)

plot_grid(plotlist=lp, align="v", nrow=3, ncol=1, labels=c("a)", "b)", "c)"), label_size= 25)

dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 2.

# System with 26 populations (cantons). to make the following calculations (Fig.3):

# - 3.a: 26 cantons with different vaccination coverages, fully assortative and partial random mixing
# - 3.b: 26 cantons with different vaccination coverages, partial-mobility based mixing
# - 3.c: Number of cantons achieving different elimination thresholds
# - 3.d: Migration plots

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


######
### - 3.a & b : 26 cantons with different vaccination coverages, fully assortative and partial random mixing
######

# set number of cantons
k <- 26

#run the parameters script and models over again
source("N_01_params_rev.R") # N_01_params runs the parameters of the model
source("N_02_models.R") # N_02_models runs the two main SIRV models: STI_1 for scenario 1&2 and STI_2 for scenario 3

# run the initial vector with the right size (right k)
init <- as.vector(matrix(unlist(inits(k)), ncol = k, byrow = TRUE))

#simulation time
time <- seq(0, 400, 1)

# choose the source for number of sexual contacts (UK or CH), here CH
C <- C_ch
# C <- C_uk

# set mixing between population 1 for scenario 1, 0.8 for scenario 2
intercantmix <- 0.8
# intercantmix <- 1


params <- c(mu=mu, omega=omega, gamma=gamma, 
            epsilon=epsilon, intercantmix=intercantmix, m=m, p=p, C=C, N=N)


#for scenario 1 and 2: use STI_1, for scenario 3: use STI_2
# simulation_26k<- as.data.frame(ode(init, time, STI_1, parms=params))
# simulation_26k<- as.data.frame(ode(init, time, STI_2, parms=params))

simulation_26k<- as.data.frame(ode(init, time, STI_1, parms=params))

# output for 26 cantons AND heterogeneous mean
mat26t199 <- outp(simulation_26k[300:316,], k)

# output for 1 canton (homogeneous mean)
mat1hom <- outp(simulation_1k[300:316,], 1)

k <- 26
#color vector
c1<- rep("grey", k+1)
c1 <- c(c1, "black")
#lwd vector
c2 <- rep(2, k+1)
c2 <- c(c2, 3)

# pdf("HET_26cantons_scen1_15y.pdf", width=8, height=6)
# pdf("HET_26cantons_scen2_15y.pdf", width=8, height=6)
pdf("HET_26cantons_scen3_15y.pdf", width=8, height=6)
# pdf("HET_26cantons_scen1_15y_SAMEWEIGHT_CHn_SAMEVACCCOV.pdf", width=8, height=6)

par(mar = c(5,5,3,1))

plot(NA, xlim=c(0,15), ylim=c(0, 0.035), # 0.035 for CH, 0.08 for UK
     xlab= "years after vaccination onset",
     ylab="HPV-16 prevalence",
     cex.lab=2, frame.plot=FALSE)

for(i in 2:28){
  lines(seq(-1,15), mat26t199[,i],  type="l", col=c1[i], lwd=c2[i], cex=2)
}

lines(seq(-1,15), mat1hom[,2], col="red", lwd=3 )

# add line for homogeneous vaccination with mean coverage. run from l.27-55 with p <- mp 
legend(-0.5, 0.008, legend=c("prevalence in each individual cantons", #0.008 for CH. 0.015 for UK
                             'mean prevalence with homogeneous vaccination', 
                             "mean prevalence with heterogeneous vaccination"), cex=1.3,
       # col =c( "grey", 'red', 'black'),lwd = c(1,2,2),
       col =c( "grey", 'red', 'black'),lwd = c(2,3,3),
       lty=c(1,1,1),  bty="n")  
mtext( expression(bold("a)")), side=3, cex=1.5, line=+1, adj=0, at=par("usr")[1]-2.5 +0.05*diff(par("usr")[1:2] ) )
# mtext( expression(bold("b)")), side=3, cex=1.5, line=+1, adj=0, at=par("usr")[1]-2.5 +0.05*diff(par("usr")[1:2] ) )
# mtext( expression(bold("c)")), side=3, cex=1.5, line=+1, adj=0, at=par("usr")[1]-2.5 +0.05*diff(par("usr")[1:2] ) )
# text(15, mat1hom[17,2], "j")


dev.off()

#select all cantons prevalences after 15y of vaccination onset,
p15yS1 <- mat26t199[17,2:28]

# p15yS1 <- dw[17,2:28]


min(p15yS1*100)
max(p15yS1*100)

# select the average HPV16 prevalence after 15y of vaccination onset, Scenario 1
p15yS1_meanhet <- p15yS1[, 27]
p15yS1_meanhom <- mat1hom[17,2]

p15yS1_meanhet*100
p15yS1_meanhom*100


# RBI
( (mat26t199[1,28]*100-p15yS1_meanhom*100) - (mat26t199[1,28]*100- p15yS1_meanhet*100) ) / (mat26t199[1,28]*100- p15yS1_meanhet*100)

write.csv(mat26t199, file="scen1_ch_diffw_normp.csv")
# 
# dw <-read.csv("scen3_ch_diffw_normp.csv")
# dw <- dw[,2:29]
# mat26t199 <- dw

# p15yS1_meanhom <- 0.01552484

######
### - 3.c : Number of cantons achieving different elimination thresholds
######

intercantmix <- NULL
# set the vector for the different espilon1 (mixing between cantons) we want to look at (x-axis)
diff.eps <- seq(0,1,0.05)

# calculates the prevalences for each canton for all x years (according to time) for the different epsilon1 (takes time)
p_on_eps1 <- prev.free.k_fn(diff.eps)

p_on_eps2 <- prev.free.k_fn(diff.eps)

saveRDS(p_on_eps1, "p_on_eps171002")
# xt <- readRDS("p_on_eps1")
# xt <- p_on_eps1
# p_on_eps1 <- xt
dim(p_on_eps1)

cbind(300:352, apply(p_on_eps1[18,1:52,],1,sum)/26)


sum(p_on_eps1[21,18,])/26
sum(p_on_eps1[17,17,])/26

p_on_eps1[1,,1]

sum(p_on_eps1[21,18,])/26

p_on_eps1[21,17,2]

## plot after 15 y 

# different thresholds wanted:
diff.tresh <- c( 0.5 , 0.75, 0.8, 0.85, 0.9)

## required for plots
shadesOfGrey <- colorRampPalette(c("grey80", "grey0"))
colv <- shadesOfGrey(length(diff.tresh))

#legend
leg <- c()
for ( i in 1:length(diff.tresh)){
  leg[i] <- paste("threshold of", diff.tresh[i]*100,"% RR reduction")
}

## data frame for given time and tresholds
prev.f.k_15y <- prev.free.k_fn2(p_on_eps1, 16, diff.tresh)


pdf("Sens_epsk_difftres15y_20171002.pdf", width=8, height=10)

par(mar=c(5.1, 4.5, 4.1, 2.1))
par(oma=c(10, 1, 0, 3))

plot(NA, xlim=c(min(diff.eps), max(diff.eps)), ylim=c(0,26),
     xlab= expression(assortativity~index~between~cantons~epsilon[k]),
     ylab="number of cantons", 
     #      main= "after 15 years",
     yaxt="n", cex.lab=2)

axis(side=2, at = seq(0,26,1) )



for(i in 1: length(diff.tresh)){
  lines(diff.eps, prev.f.k_15y[,i,1],  col=colv[i], pch=16, type= "b", cex=2)
}


par(new=T)
plot(diff.eps, prev.f.k_15y[,1,2]*100,  axes=F, xlab=NA, ylab=NA, cex=1.2, ylim=c(0,0.7), lty=2,
     type="l", col="red", lwd=4)
axis(side=4)
mtext(side=4, line=3, "HPV-16 prevalence in Switzerland (%)", cex=2)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend=c(leg, "HPV-16 prevalence in Switzerland") , 
       pch=c( rep(16, length(diff.tresh)), NA ), lty =c(rep(1, length(diff.tresh)), 2), col=c(colv, "red"), 
       cex=1.4 , xpd = TRUE , inset=c(1,0.02),lwd= c( (rep(1,length(diff.tresh), 4)) ))
mtext( expression(bold("a)")), side=3, cex=1.5, line=-2, adj=0 )


dev.off()


## plot after 50 y 

# different thresholds wanted:
diff.tresh <- c(0.90, 0.95, 0.990 ,0.995, 0.999)

#change oy to look at the results after x number of years
oy <- 50

oy <- 100

oy <- c(31, 41, 51, 61, 71, 81, 91, 100)

for ( z in 1:length(oy)){
  ####### Makes data frame
  prev.f.k_50y <- prev.free.k_fn2(p_on_eps1, oy[z], diff.tresh)
  
  ## required for plots
  colv <- shadesOfGrey(length(diff.tresh))
  
  #legend
  leg <- c()
  for ( i in 1:length(diff.tresh)){
    leg[i] <- paste("threshold of", diff.tresh[i]*100,"% RR reduction")
  }
  
  
  pdf(paste(oy[z]-1,"y_Sens_epsk_difftres.pdf"), width=8, height=10)
  
  par(mar=c(5.1, 4.5, 4.1, 2.1))
  par(oma=c(10, 1, 0, 3))
  
  plot(NA, xlim=c(min(diff.eps), max(diff.eps)), ylim=c(0,26),
       xlab= expression(assortativity~index~between~cantons~epsilon[k]),
       ylab="number of cantons", 
       #      main= "after 50 years",
       yaxt="n" , cex.lab=2)
  
  axis(side=2, at = seq(0,26,1) )
  
  for(i in 1: length(diff.tresh)){
    lines(diff.eps, prev.f.k_50y[,i,1],  col=colv[i], pch=16, type= "b" , cex=2 )
  }
  
  
  par(new=T)
  plot(diff.eps, prev.f.k_50y[,1,2]*100,  axes=F, xlab=NA, ylab=NA, cex=1.2, ylim=c(0,0.7), lty=2,
       type="l", col="red", lwd=4)
  axis(side=4)
  mtext(side=4, line=3, "HPV-16 prevalence in Switzerland (%)", cex=2)
  
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend=c(leg, "HPV-16 prevalence in Switzerland") , 
         pch=c( rep(16, length(diff.tresh)), NA ), lty =c(rep(1.5, length(diff.tresh)), 2), col=c(colv, "red"), 
         cex=1.4 , xpd = TRUE , inset=c(1,0.02), lwd= c( (rep(1,length(diff.tresh), 4)) ))
  mtext( expression(bold("b)")), side=3, cex=1.5, line=-2, adj=0 )
  
  dev.off()
  
}

######
### - 3.d : Migration plots
######

k <- 26

#run the parameters script and models over again
source("N_01_params.R") # N_01_params runs the parameters of the model

library("circlize")
require("migest")
require("plyr")
require("reshape2")

#sigma is calculated from the pmob matrix in N_01_params.R
sigma

m_s <- as.matrix(sigma[,1:26])



colnames(m_s) <- pmob[,1]
rownames(m_s) <- pmob[,1]




# m_s <- t(m_s)

colSums(m_s)
rowSums(m_s)

diag(m_s) <- c(0)
m_s[lower.tri(m_s, diag = FALSE)] <- 0

m_s2 <- m_s* (w*(S/k))

# m_s2 <- t(m_s2)

chordDiagram(x = m_s)


chordDiagram(x = m_s2)

# 
# m1 <- melt(m)
# head(m1)
# round(colSums(m_s),2)
# round(rowSums(m_s),2)
# m_s[26,]
# round(m_s[,1],3)
# 
# round(m_s[,26],3)
# 
# sum(round(m_s[,26],2))
# 
# pdf("Migrationcircosigma2.pdf", width=6, height=6)
# chordDiagram(x = m_s)
# dev.off()

## Run from here for the scenario 1 2 and 3
d1 <- melt(m_s2)
head(d1)
colnames(d1) <- c("Ka", "Ka2", "value")


# ju <- d1[d1$Ka=="JU",]
# zh <- d1[d1$Ka=="ZH",]
# zh2 <- d1[d1$Ka2=="ZH",]
# sg <- d1[d1$Ka=="SG",]
# sg2 <- d1[d1$Ka2=="SG",]
# 
# 
# sum(ju[,3])
# 
# sum(zh[,3])
# sum(zh2[,3])
# 
# sum(sg[,3])
# sum(sg2[,3])

# chordDiagram(x = d1)

#### attributing colors according to language regions

d1$sprache<- d1$sprache[d1$Ka=="ZH" &
                          d1$Ka=="BE" &
                          d1$Ka=="LU" &
                          d1$Ka=="UR" &
                          d1$Ka=="SZ" &
                          d1$Ka=="OW" &
                          d1$Ka=="NW" &
                          d1$Ka=="GL" &
                          d1$Ka=="ZG" &
                          d1$Ka=="SO" &
                          d1$Ka=="BS" &
                          d1$Ka=="BL" &
                          d1$Ka=="SH" &
                          d1$Ka=="AR" &
                          d1$Ka=="AI" &
                          d1$Ka=="SG" &
                          d1$Ka=="GR" &
                          d1$Ka=="AG" &
                          d1$Ka=="TG" ]<- "DE"

d1$sprache[d1$Ka=="JU"]<- "FR"
d1$sprache[d1$Ka=="NE"]<- "FR"
d1$sprache[d1$Ka=="FR"]<- "FR"
d1$sprache[d1$Ka=="VD"]<- "FR"
d1$sprache[d1$Ka=="GE"]<- "FR" 
d1$sprache[d1$Ka=="VS"]<- "FR"
d1$sprache[d1$Ka=="TI"] <- "IT"

d1$sprache<- as.factor(d1$sprache)

d1$sprache2<- d1$sprache2[d1$Ka2=="ZH" &
                            d1$Ka2=="BE" &
                            d1$Ka2=="LU" &
                            d1$Ka2=="UR" &
                            d1$Ka2=="SZ" &
                            d1$Ka2=="OW" &
                            d1$Ka2=="NW" &
                            d1$Ka2=="GL" &
                            d1$Ka2=="ZG" &
                            d1$Ka2=="SO" &
                            d1$Ka2=="BS" &
                            d1$Ka2=="BL" &
                            d1$Ka2=="SH" &
                            d1$Ka2=="AR" &
                            d1$Ka2=="AI" &
                            d1$Ka2=="SG" &
                            d1$Ka2=="GR" &
                            d1$Ka2=="AG" &
                            d1$Ka2=="TG" ]<- "DE"

d1$sprache2[d1$Ka2=="JU"]<- "FR"
d1$sprache2[d1$Ka2=="NE"]<- "FR"
d1$sprache2[d1$Ka2=="FR"]<- "FR"
d1$sprache2[d1$Ka2=="VD"]<- "FR"
d1$sprache2[d1$Ka2=="GE"]<- "FR" 
d1$sprache2[d1$Ka2=="VS"]<- "FR"
d1$sprache2[d1$Ka2=="TI"] <- "IT"

d1$sprache2 <- as.factor(d1$sprache2)


# 3 categorise between romands, de, or interlanguages

d1$cat <- ifelse( (d1$sprache == d1$sprache2 & d1$sprache == "DE" ), c("DE"), 
                  ifelse( (d1$sprache == d1$sprache2 & d1$sprache == "FR" ), c("FR"), c("inter.") )  )

d2 <- d1

d2$sprache2 <- d2$sprache <- NULL


head(d2)

col.de <- "tan"
col.fr <-  "steelblue2"
col.it <- "darkred"

g.col.l <- ifelse( d1$cat=="FR", "steelblue2", 
                   ifelse( d1$cat=="inter.", "thistle4",
                           "tan" ) )


grid.col =c(ZH = col.de, BE = col.de, LU = col.de , UR = col.de, SZ= col.de, OW= col.de,
            NW= col.de, GL= col.de, ZG= col.de, SO= col.de, BS= col.de, BL= col.de, SH= col.de,
            AR= col.de, AI= col.de, SG= col.de, GR= col.de, AG= col.de, TG= col.de,
            JU =col.fr, NE = col.fr, FR= col.fr, VD= col.fr, GE= col.fr, VS= col.fr,
            TI =col.it)


scen3 <- d2
# pdf("Migrationcircosigma4_col1_b.pdf", width=6, height=6)
# pdf("Migrationcircos_scenario2.pdf", width=6, height=6)
pdf("Migrationcircos_scenario1.pdf", width=6, height=6)


chordDiagram(d2, grid.col= grid.col, col=g.col.l, annotationTrack="grid" , preAllocateTracks =list(track.height = 0.3) )  


circos.trackPlotRegion(track.index = 1, panel.fun =function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", cex=0.75,
              niceFacing = TRUE, adj =c(0, 0.5))}, bg.border = NA)


dev.off()

# 
# chordDiagram(d2, grid.col = grid.col,
#  annotationTrack = "grid", annotationTrackHeight = 0.15)
# for(si in get.all.sector.index()) {
#   xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
#   ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
#   circos.text(mean(xlim),mean(ylim), si, sector.index = si, track.index = 1,
#     facing = "bending.inside", col = "white")
# }

scen3 <- d2

# for scenario 1 and 2:
source("N_02_models.R") # N_01_params runs the parameters of the model
C <- C_ch

# first for scenario 1&2: rhokkrr

# for scenario 1:
epsilon1 <- 1

# for scenario 2:
epsilon1<- 0.8

rho_ki<- rho_ki_f(epsilon, epsilon1, k,C)

#make the sum of H and L sexual behaviour group in each canton
mig <- matrix(nrow=26, ncol=26)
miga <- matrix(nrow=2*26, ncol=26)


for(i in 1:26){
  miga[,i] <- rho_ki[,(i*2)-1] + rho_ki[,(i*2)]
  
}

for(i in 1:26){
  mig[i,] <- miga[(i*2)-1,] + miga[(i*2),]
  
}

colSums(mig)
rowSums(mig)

mig <- mig/2

m_s <- as.matrix(mig)

colnames(m_s) <- pmob[,1]
rownames(m_s) <- pmob[,1]


# for scenario 1 (or do nothing, it is the same)
m_s[lower.tri(m_s, diag = FALSE)] <- 0

#for scenario 2
m_s[lower.tri(m_s, diag = TRUE)] <- 0

m_s2 <- m_s* (w*(S/k))

# chordDiagram(x = m_s)
# chordDiagram(x = m_s2)

### run from line 667 - 777 downwards

scen1 <- d2
scen2 <- d2  


pdf("Migrationcircos_allthree.pdf", width=6, height=16)
par(mfrow=c(3,1))



pdf("Migrationcircos_scenario1.pdf", width=6, height=6)

chordDiagram(scen1, grid.col= grid.col, col=g.col.l, annotationTrack="grid" , preAllocateTracks =list(track.height = 0.3) )  

circos.trackPlotRegion(track.index = 1, panel.fun =function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", cex=0.75,
              niceFacing = TRUE, adj =c(0, 0.5))}, bg.border = NA)
mtext("assortative sexual mixing", font=2, line=-4)
dev.off()


pdf("Migrationcircos_scenario2.pdf", width=6, height=6)
mtext("Assortative sexual mixing", outer=TRUE, line=-17)
chordDiagram(scen2, grid.col= grid.col, col=g.col.l, annotationTrack="grid" , preAllocateTracks =list(track.height = 0.3) )  
circos.trackPlotRegion(track.index = 1, panel.fun =function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", cex=0.75,
              niceFacing = TRUE, adj =c(0, 0.5))}, bg.border = NA)
mtext("proportional sexual mixing", font=2, line=-4)
dev.off()

pdf("Migrationcircos_scenario3.pdf", width=6, height=6)
chordDiagram(scen3, grid.col= grid.col, col=g.col.l, annotationTrack="grid" , preAllocateTracks =list(track.height = 0.3) )  
circos.trackPlotRegion(track.index = 1, panel.fun =function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", cex=0.75,
              niceFacing = TRUE, adj =c(0, 0.5))}, bg.border = NA)
mtext("mobility-informed sexual mixing", font=2, line=-4)
dev.off()


mtext("Overall Title Row 2", outer=TRUE, line=-17)
mtext("Overall Title Row 3", outer=TRUE, line=-34) 

dev.off()
