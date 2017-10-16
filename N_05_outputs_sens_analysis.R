# 18.08.2017
# Script which shows the output from the sensitivity simulations

sens_rep <- 1000

r <- 2

time <- seq(0, 400, 1)

# save(res, file=paste("p",i,".RData", sep=""))
k <- 26

# dim(sto_sim_het)
sto_sim_het <- array(NA, c(length(time), ((r*8)*k)+1, sens_rep ) )
sto_sim_hom <- array(NA, c(length(time), ((r*8)*1)+1, sens_rep ) )
pmat <- array(NA, c( sens_rep, 15 ) )



for( i in 1:sens_rep){
  # load(paste("C:/Users/mriesen.CAMPUS/Documents/data_revisions/outp_scen1a/p", i, ".RData", sep="")) #eps1 scenario 1
  # load(paste("C:/Users/mriesen.CAMPUS/Documents/data_revisions/outp_scen2a/p", i, ".RData", sep="")) #eps0.8 scenario 2
  # load(paste("C:/Users/mriesen.CAMPUS/Documents/data_revisions/outp_scen2/p", i, ".RData", sep="")) #eps0.8 scenario 2
  # load(paste("C:/Users/mriesen.CAMPUS/Documents/data_revisions/outp_scen3a/p", i, ".RData", sep="")) #eps0.8 scenario 2
  load(paste("C:/Users/mriesen.CAMPUS/Documents/data_revisions/outp_scen1_uk/p", i, ".RData", sep="")) 
  sto_sim_het[,,i] <- res$r1
  sto_sim_hom[,,i] <- res$r2
  pmat[i,] <- res$r3
}

#matrix with national prevalence over time and repetitions
natprevf <- function(sim, k ){
  mat <- array(NA, c(length(time),  sens_rep ) )
  for( i in 1:sens_rep){
    mat[,i] <- apply(sim[,,i], 1,  function(x) sum(  x[(1:k*16)-14+1] + x[(1:k*16)-10+1]) / k    )
  } 
  return(mat)
}

mat_het <- natprevf(sto_sim_het, 26)
mat_hom <- natprevf(sto_sim_hom,1)


sto_sim_het[1,,4]

res$r3

#check if the pre-vaccination state is the same in hom and het
round(mat_het[299,],4) == round(mat_hom[299,], 4)
d1 <- mat_hom[299,] - mat_het[299,] 
summary(d1)

d1 <- mat_hom[399,] - mat_het[399,] 
summary(d1)

mat_het[399,which.min(d1)]
mat_hom[399,which.min(d1)]

dim(mat_het)


plot(NA, xlim=c(0,400), ylim=c(0, 0.3),
     xlab= "years after vaccination onset",
     ylab="HPV-16 prevalence",
     cex.lab=2, frame.plot=FALSE)


#with 26k
for(i in 1:length(mat_het[1,]) ) {
  lines( mat_het[,i] ,  type="l", col="black", lwd=1, cex=2)
}


for(i in 1:length(mat_hom[1,])) {
  lines( mat_hom[,i] ,  type="l", col="red", lwd=1, cex=2)
}

#delete simulations where the prevalence is < 0.001 (<0.01%) before vaccination
mat_het <- mat_het[,apply(mat_het,2,function(x) all(x[299]>=0.001))]
mat_hom <- mat_hom[,apply(mat_hom,2,function(x) all(x[299]>=0.001))]

#number of simulations higher than 0.001 prev
nbsim <-length(mat_het[1,])

##difference between het and hom national prevalence:
diffhh <- array(NA, c(length(time),  length(mat_het[1,]) ) )

for(i in 1:length(mat_het[1,])){
  diffhh[,i] <- mat_het[,i] - mat_hom[,i] 
}



# mat_het[,654] - mat_hom[,654] 
# 
# mat_het[1,]

plot(NA, xlim=c(200,400), ylim=c(0, 0.01),
     xlab= "years after vaccination onset",
     ylab="prevalence difference",
     cex.lab=2, frame.plot=FALSE)

#with 26k
for(i in 1:length(mat_het[1,]) ) {
  lines( diffhh[,i] ,  type="l", col="grey", lwd=1, cex=2)
}


#plot with median difference and 95% a and 50%

med <- c()
CI95 <- matrix(nrow=400, ncol=2)
CI50 <- matrix(nrow=400, ncol=2)

for( i in 1:400){
  med[i] <- median(diffhh[i,])
  CI95[i,]<- round(quantile(diffhh[i,],probs=c(.025,.975)), 4)
  CI50[i,] <- round(quantile(diffhh[i,],probs=c(.25,.75)), 4)
}


# pdf("prev_diff_1000rep_400y.pdf", width=8, height=6)
# pdf("prev_diff_1000rep_400yeps08.pdf", width=8, height=6)
pdf("prev_diff_1000rep_400yscen3.pdf", width=8, height=6)
par(mar = c(5,5,3,1))
plot(NA, xlim=c(300, 400), ylim=c(0, 0.01),
     xlab= "years after vaccination onset",
     ylab="prevalence difference",
     cex.lab=2, frame.plot=FALSE)
lines( med ,  type="l", col="black", lwd=2, cex=2)
polygon(c(200:400, rev(200:400)), c(CI95[200:400,2], rev(CI95[200:400,1]) ),
        border = NA, col=rgb(0.01, 0, 0.1,0.1))
polygon(c(200:400, rev(200:400)), c(CI50[200:400,2], rev(CI50[200:400,1]) ),
        border = NA, col=rgb(0.01, 0, 0.8,0.1))


dev.off()




ymax <- 100
ymaxr <- 400 


# pdf("prev_diff_1000rep.pdf", width=8, height=6)
# pdf("prev_diff_1000rep_eps08.pdf", width=8, height=6)
pdf("prev_diff_1000rep_scen3.pdf", width=8, height=6)
par(mar = c(5,5,3,1))
plot(NA, xlim=c(0, ymax), 
      ylim=c(0.00001, 0.005),
     xlab= "years after vaccination onset",
     ylab="prevalence difference",
     cex.lab=2, frame.plot=FALSE, xaxt="n")
     # , log = "y")

axis(1, at = seq(0, ymax, by = 5), las=2)

#with 26k
lines( med[301:ymaxr] ,  type="l", col="black", lwd=2, cex=2)
# lines( CI95[,1] ,  type="l", col="darkred", lwd=1, cex=2)
# lines( CI95[,2] ,  type="l", col="darkred", lwd=1, cex=2)
# lines( CI50[,1] ,  type="l", col="blue", lwd=1, cex=2)
# lines( CI50[,2] ,  type="l", col="blue", lwd=1, cex=2)

polygon(c(1:ymax, rev(1:ymax)), c(CI95[301:ymaxr,2], rev(CI95[301:ymaxr,1] )),
        border = NA, col=rgb(0.01, 0, 0.1,0.1)) 
polygon(c(1:ymax, rev(1:ymax)), c(CI50[301:ymaxr,2], rev(CI50[301:ymaxr,1] )),
        border = NA, col=rgb(0.01, 0, 0.8,0.1)) 


legend(0, 0.005, legend=c("median", "95% of the simulations", 
                          "50% of the simulations"), cex=1,
       # col =c( "grey", 'red', 'black'),lwd = c(1,2,2),
       col =c("black", rgb(0.01, 0, 0.1,0.1), rgb(0.01, 0, 0.8,0.1)),lwd = c(2,3,3),
       lty=c(1,1,1),  bty="n")  

# lines(x=c(which.max(med[301:ymaxr]),which.max(med[301:ymaxr])) , y= c(0,0.003), lty=2)


dev.off()

ymax <- 15
ymaxr <- 215

# pdf("prev_diff_1000_15y_b.pdf", width=8, height=6)
pdf("prev_diff_1000_15y_b_eps08.pdf", width=8, height=6)
par(mar = c(5,5,3,1))
plot(NA, xlim=c(0, ymax), ylim=c(0, 0.0015),
     xlab= "years after vaccination onset",
     ylab="prevalence difference",
     cex.lab=2, frame.plot=FALSE, xaxt="n")

axis(1, at = seq(0, 15, by = 5), las=2)

#with 26k
lines( med[201:ymaxr] ,  type="l", col="black", lwd=2, cex=2)
# lines( CI95[,1] ,  type="l", col="darkred", lwd=1, cex=2)
# lines( CI95[,2] ,  type="l", col="darkred", lwd=1, cex=2)
# lines( CI50[,1] ,  type="l", col="blue", lwd=1, cex=2)
# lines( CI50[,2] ,  type="l", col="blue", lwd=1, cex=2)

polygon(c(1:ymax, rev(1:ymax)), c(CI95[201:ymaxr,2], rev(CI95[201:ymaxr,1] )),
        border = NA, col=rgb(0.01, 0, 0.1,0.1)) 
polygon(c(1:ymax, rev(1:ymax)), c(CI50[201:ymaxr,2], rev(CI50[201:ymaxr,1] )),
        border = NA, col=rgb(0.01, 0, 0.8,0.1)) 


legend(0, 0.005, legend=c("median", "95% of the simulations", 
                          "50% of the simulations"), cex=1,
       # col =c( "grey", 'red', 'black'),lwd = c(1,2,2),
       col =c("black", rgb(0.01, 0, 0.1,0.1), rgb(0.01, 0, 0.8,0.1)),lwd = c(2,3,3),
       lty=c(1,1,1),  bty="n")  

# lines(x=c(which.max(med[201:ymaxr]),which.max(med[201:ymaxr])) , y= c(0,0.003), lty=2)
dev.off()
# for(i in 1:length(mat_het[1,]) ) {
#   lines( diffhh[201:ymaxr,i] ,  type="l", col="grey", lwd=1, cex=2)
# }

ymax <- 100
ymaxr <- 400 

med100y <- med[300:ymaxr]
CI95100y <- CI95[300:ymaxr,]
CI50100y <- CI50[300:ymaxr,]



# pdf("prev_diff_1000_15y_b.pdf", width=8, height=6)
# pdf("prev_diff_1000_50y_scen2_b.pdf", width=8, height=6)
pdf("prev_diff_1000_50y_scen2_20171002_withinfo.pdf", width=8, height=6)
pdf("prev_diff_1000_50y_scen2_20171002.pdf", width=8, height=6)
par(mar = c(5,5,3,1))
plot(NA, xlim=c(0, 50), ylim=c(0, 0.005),
     xlab= "years after vaccination onset",
     ylab="prevalence difference",
     cex.lab=2, frame.plot=FALSE, xaxt="n")

axis(1, at = seq(0, 50, by = 5), las=2)

length(0:100)
dim(mat_hom)


#with 26k
lines( 0:ymax, med100y,  type="l", col="black", lwd=2, cex=2)

polygon(c(0:ymax, rev(0:ymax)), c(CI95100y[,2], rev(CI95100y[,1] )),
        border = NA, col=rgb(0.01, 0, 0.1,0.1)) 
polygon(c(0:ymax, rev(0:ymax)), c(CI50100y[,2], rev(CI50100y[,1] )),
        border = NA, col=rgb(0.01, 0, 0.8,0.1)) 

legend(0, 0.005, legend=c("median", "95% of the simulations", 
                          "50% of the simulations"), cex=1,
       # col =c( "grey", 'red', 'black'),lwd = c(1,2,2),
       col =c("black", rgb(0.01, 0, 0.1,0.1), rgb(0.01, 0, 0.8,0.1)),lwd = c(2,3,3),
       lty=c(1,1,1),  bty="n")  

lines(x=c(which.max(med100y[0:ymaxr]),which.max(med100y[0:ymaxr])) , y= c(0,0.003), lty=2)

text(which.max(med100y[0:ymaxr])-5,  0.0035, paste("max median: ",round(max(med100y),4),
                                                 " (IQR: ", CI50100y[which.max(med100y[0:ymaxr]),1],
                                                 "-",CI50100y[which.max(med100y[0:ymaxr]),2],
                                                 "), " ),  cex=0.8 )
text(which.max(med100y[0:ymaxr])-10,  0.0033, paste("reached after ", which.max(med100y[0:ymaxr]), " years", sep=""), cex=0.8 )

dev.off()

round(med[395],5) == round(med[400],5)


round(med[350:400],6)

summary(diffhh[400,])

mean(diffhh[400,])
median(diffhh[400,])
max(diffhh[400,])

mat_het[400,]

# at year 400
# d1 <- diffhh[400,]
d1 <- diffhh[316,]
#delete values where prevalence is < 0.1% (<0.001)
h <- hist(d1)
plot(d1)

round( median(d1), 4)
round(quantile(d1,probs=c(.025,.975)), 4)
round(quantile(d1,probs=c(.25,.75)), 4)

############## estimation around the pre-vaccination prevalences
d1a <-  mat_hom[299,]


quantile(d1a,probs=c(.025,.975))

qts <- quantile(d1a,probs=c(.05,.95))
hist(d1a, breaks=100)
abline(v=qts[1],col="red")
abline(v=qts[2],col="red")


round( median(d1a)*100, 3)
round(quantile(d1a,probs=c(.025,.975))*100, 3)
round(quantile(d1a,probs=c(.25,.75))*100, 3)



#############



########## b: Calculating R0
dim(mat_hom)
head(mat_hom)
mat_hom[1,]

pmat <- pmat[apply(mat_hom,2,function(x) all(x[299]>=0.001)),]


k <- 1
r <- 2



# colnames(pmat) <- names(res$r3) 
colnames(pmat) <- c("mu", "omega", "gamma","epsilon","intercantmix", "m", "p", "C1", "C2", "N1","N2",
                    "n1","n2","beta", "ve")        

pset <- pmat
vt <- c()
R0_s <- c()
pset[3,"beta"]

for( i in 1:length(pset[,1])){
  beta_1 <- matrix(rep(pset[i,"beta"], (k*r)*(k*r)), ncol=k*r, nrow=k*r)
  
  C_1 <- matrix(c(pset[i,"C1"],  pset[i,"C2"], pset[i,"C1"], pset[i,"C2"]) , ncol=k*r, nrow=k*r)
  
  #population division matrix
  pd <- matrix(c(1, pset[i,"n2"]/pset[i,"n1"], pset[i,"n1"]/pset[i,"n2"], 1), ncol=2, nrow=2)
  
  C_temp <- cbind(pset[i,"C1"], pset[i,"C2"])
  N <- cbind(pset[i,"N1"], pset[i,"N2"])
  
  rho_ki <- rho_ki_f(pset[i,"epsilon"], pset[i,"intercantmix"], k, C_temp)
  
  #New infections
  f <- (beta_1*C_1*rho_ki*pd)
  #Time in compartment
  V <-  matrix(c(pset[i,"mu"]+pset[i,"gamma"]+pset[i,"m"]*pset[i,"n2"], -pset[i,"m"]*pset[i,"n2"], 
                 -pset[i,"m"]*pset[i,"n1"], pset[i,"mu"]+pset[i,"gamma"]+pset[i,"m"]*pset[i,"n1"] ), ncol=2, nrow=2)
  
  G <- f %*% solve(V)
  
  R0 <- max(eigen(G)$values)
  vt[i] <- 1- 1/R0
  R0_s[i]<- R0
}


round( median(R0_s), 2)
round(quantile(R0_s,probs=c(.025,.975)),2)
round(quantile(R0_s,probs=c(.25,.75)), 2)

round( median( (vt*100)/pset[,15]), 2)
round(quantile( (vt*100)/pset[,15],probs=c(.025,.975)),2)
round(quantile( (vt*100)/pset[,15],probs=c(.25,.75)), 2)

round( median(1- 1/R0_s)*100, 2)
round(quantile((1- 1/R0_s)*100 ,probs=c(.025,.975)),2)
round(quantile( (1- 1/R0_s)*100,probs=c(.25,.75)), 2)

round( median( ( (1- 1/R0_s^2)*100)/pset[,15]), 2)
round( quantile( ((1- 1/R0_s^2)*100)/pset[,15] ,probs=c(.025,.975) ),2)
round(quantile( ( (1- 1/R0_s^2)*100)/pset[,15],probs=c(.25,.75) ), 2)


temp <- 1/R0_s

summary(temp)

############## ############## which params brings < 0 prev diff

which(diffhh[325,]< 0)

pmat[which(diffhh[325,]< 0),]

diffhh[300:350,which(diffhh[325,]< 0)]

summary(pmat[which(diffhh[325,]< 0),])

pmat[which.min(diffhh[325,]),]
pmat[which.max(diffhh[325,]),]


summary(pmat)

############## relative risk reduction

rrred <- array(NA, c(length(201:400),  length(mat_het[1,]) ) )

for(i in 1:length(mat_het[1,])){
  rrred[,i] <-  ( (mat_hom[201,i]-mat_hom[201:400,i]) -   (mat_het[201,i] - mat_het[201:400,i]) ) / (mat_het[201,i] - mat_het[201:400,i])
}


plot(NA, xlim=c(0,200), ylim=c(0, 0.2),
     xlab= "years after vaccination onset",
     ylab="prevalence difference",
     cex.lab=2, frame.plot=FALSE)


#with 26k
for(i in 1:length(mat_het[1,]) ) {
  lines( rrred[,i] ,  type="l", col="grey", lwd=1, cex=2)
}


#plot with median difference and 95% a and 50%

med <- c()
CI95 <- matrix(nrow=200, ncol=2)
CI50 <- matrix(nrow=200, ncol=2)

for( i in 1:200){
  med[i] <- median(rrred[i,])
  CI95[i,]<- round(quantile(rrred[i,],probs=c(.025,.975), na.rm=TRUE), 4)
  CI50[i,] <- round(quantile(rrred[i,],probs=c(.25,.75), na.rm=TRUE), 4)
}


# pdf("rrred_perc_1000rep_200y.pdf", width=8, height=6)
pdf("rrred_perc_1000rep_50y_eps08.pdf", width=8, height=6)
par(mar = c(5,5,3,1))
plot(NA, xlim=c(0, 51), ylim=c(0, 10),
     xlab= "years after vaccination onset",
     ylab="relative risk reduction (%)",
     cex.lab=2, frame.plot=FALSE)
lines( med*100 ,  type="l", col="black", lwd=2, cex=2)
polygon(c(1:200, rev(1:200)), c(CI95[0:200,2]*100, rev(CI95[0:200,1]*100) ),
        border = NA, col=rgb(0.01, 0, 0.1,0.1))
polygon(c(1:200, rev(1:200)), c(CI50[0:200,2]*100, rev(CI50[0:200,1]*100) ),
        border = NA, col=rgb(0.01, 0, 0.8,0.1))

# lines(x=c(which.max(med[0:200]),which.max(med[0:200])) , y= c(0,6), lty=2)


dev.off()

length(CI95[0:200,2])


ymax <- 150
ymaxr <- 350
