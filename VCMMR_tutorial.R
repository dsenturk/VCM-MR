## VCMMR_tutorial.R
#############################################################################
## Description: A step-by-step implementation of VCM-MR and the associated  
## procedures described in "Modeling Time-Varying Effects of Multilevel Risk
## Factors of Hospitalizations in Patients on Dialysis".  
#############################################################################
## Functions implemented: 
## VCMMR_simulation.R, VCMMR_estimation.R, VCMMR_bootstrap.R, VCMMR_hypothesistest.R
#############################################################################
## Tutorial Outline:
## 1. Simulate hospitalization outcome data (VCMMR_simulation.R)
## 2. Perform VCMMR estimation (VCMMR_estimation.R)
## 3. Inference on multilevel risk factors via bootstrap (VCMMR_bootstrap.R)
## 4. Hypothesis test of facility-specific fixed effects (VCMMR_hypothesistest.R)
## 5. Prediction of patient- and facility-level risk trajectory
## 6. Visualization of VCMMR results
#############################################################################

# Install missing packages
list.of.packages <- c("MASS", "statmod", "mvtnorm","bisoreg")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 

# Load packages  
library(MASS)
library(statmod)
library(mvtnorm)
library(bisoreg) 


#############################################################################
# 1. Simulate hospitalization outcome data
#############################################################################

# NOTE: Generating one dataset with 100 facilities will take approximately eight minutes.

# Simulate one dataset from the simulation design described in Web Appendix E with 100 facilities
data <- VCMMR_simulation(numF = 100)  # VCMMR_simulation.R

#############################################################################
# 2. Perform VCMMR estimation
#############################################################################

# NOTE: Performing VCMMR estimation with 100 facilities will take approximately four minutes.

VCMMREst <- VCMMR_estimation(data = data, sigma2 = 1, hThetaBeta = .36, hEta = .8, hGamma = c(.22, .27, .32))  # VCMMR_estimation.R

#############################################################################
# 3. Inference on multilevel risk factors via bootstrap
#############################################################################

# NOTE: Obtaining the bootstrap confidence intervals may take several hours
#       depending on the size of the data set and processor speed. 

boot_Est <- VCMMR_bootstrap(nboot = 200, data = data, sigma2 = 1, hThetaBeta = .36, hEta = .8, hGamma = c(.22, .27, .32))  # VCMMR_bootstrap.R
 
#############################################################################
# 4. Hypothesis test of facility-specific fixed effects
#############################################################################

# NOTE: The hypothesis testing procedure will take approximately five minutes.  

# Produce p-value from hypothesis testing procedure for testing H0:gamma1(t)=0 for facility 1

pVal <- VCMMR_hypothesistest(Fid = 1, data = data, VCMMREst = VCMMREst, numBS = 500, hGamma = c(.22, .27, .32))  # VCMMR_hypothesistest.R

#############################################################################
# 5. Prediction of patient- and facility-level risk trajectories
#############################################################################  

gridPoints <- seq(0,1,1/19) # Grid points for theta(t) and beta(t) 
etaGrid <- seq(0,1,1/19) # Grid points for eta(c)

# True varying coefficient functions

# Facility fixed effects 
flatFxn <- function(x) {
  return(rep(-1,length(x)))
}
sqrtFxn <- function(x) {
  return(-sqrt(x) - 1)
}
quadraticFxn <- function(x) {
  return((-x-0.5)^2- 2)
}

# Effects of subject-level risk factors
beta1Fxn <- function(x) {
  return(cos(pi*x))
}
beta2Fxn <- function(x) {
  return(-cos(pi*x))
}

# Effects of facility-level risk factors
theta1Fxn <- function(x) {
  return(sin(pi*x))
}
theta2Fxn <- function(x) {
  return(-sin(pi*x))
}

# Calendar year effect
etaFxn <- function(x) {
  return(-0.2*x+0.1)
}

invLogit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

theta1F <- theta1Fxn(gridPoints)
theta2F <- theta2Fxn(gridPoints)
beta1F <- beta1Fxn(gridPoints)
beta2F <- beta2Fxn(gridPoints) 
etaF <- etaFxn(etaGrid)
ntheta <- 2
nbeta <- 2

# Facility shapes for 100 facilities
# 1 - Sqrt, 2 - Flat, 3 - Quad
facilityShapes <- c(rep(1,12),rep(2,11),rep(3,11),rep(1,11),rep(2,11),rep(3,11),rep(1,11),rep(2,11),rep(3,11))
# If you want 500 facilities comment the above line and uncomment the following line.
# facilityShapes <- c(rep(1,56),rep(2,56),rep(3,55),rep(1,56),rep(2,56),rep(3,55),rep(1,56),rep(2,55),rep(3,55))

df <- data
# Index for varying coefficient functions
df$r <- round(19*df$t)+1 
df$cr <- round(19*df$c)+1  

gammaVCMEst <- VCMMREst[[1]]
thetaVCMEst <- VCMMREst[[2]]
betaVCMEst <- VCMMREst[[3]]
etaVCMEst <- VCMMREst[[4]]
sigma2bEst <- VCMMREst[[5]]
bijEst <- VCMMREst[[6]]

# Patient-level prediction for the first patient from the first facility
Patid <- 1
ngrid <- length(gridPoints)
df_pat <- df[df$sid==Patid,]
Patfac <- df_pat$fid[1]
fac_covatiates <- df_pat[rep(1,ngrid), 5:(4+ntheta)]
pat_covariates <- df_pat[rep(1,ngrid), (5+ntheta):(4+ntheta+nbeta)]
pat_lpred <- gammaVCMEst[Patfac, ] +  rowSums(fac_covatiates * thetaVCMEst) + rowSums(pat_covariates * betaVCMEst) 
pred_pat <- invLogit(pat_lpred + etaVCMEst[df_pat$cr[1]]) # Predicted risk for the first patient

if(facilityShapes[Patfac]==1){
  gammaF <- sqrtFxn(gridPoints)
} else if(facilityShapes[Patfac]==2){
  gammaF <- flatFxn(gridPoints)
} else {
  gammaF <- quadraticFxn(gridPoints)
}

# True risk for the first patient
pat_ltrue <- gammaF + fac_covatiates[,1] * theta1F +  fac_covatiates[,2] * theta2F + pat_covariates[,1] * beta1F + pat_covariates[,2] * beta2F
true_pat <- invLogit(pat_ltrue + etaFxn(df_pat$c[1]))


# Facility-level prediction for the first facility
Facid <- 1
df_fac <- df[df$fid == Facid,] 
thetaPre <- thetaVCMEst[df_fac$r,]
betaPre <- betaVCMEst[df_fac$r,]
etaPre <- etaVCMEst[df_fac$cr]
bijPre <- bijEst[df_fac$sid]
thetaTrue <- cbind(theta1F[df_fac$r], theta2F[df_fac$r])
betaTrue <- cbind(beta1F[df_fac$r], beta2F[df_fac$r])
etaTrue <- etaF[df_fac$cr]
bijTrue <- df_fac$bi

if(facilityShapes[Facid]==1){
  gammaF <- sqrtFxn(gridPoints)
} else if(facilityShapes[Facid]==2){
  gammaF <- flatFxn(gridPoints)
} else {
  gammaF <- quadraticFxn(gridPoints)
} 
lfollow <- max(df_fac$r) # Maximum follow-up time of the first facility
  
pred_f <- c() # Predicted risk for the first facility
true_f <- c() # True risk for the first facility
for(j in 1:lfollow){
  spred <- rowSums(betaPre[df_fac$r==j,]*df_fac[df_fac$r==j,(5+ntheta):(4+ntheta+nbeta)])
  fpred <- rowSums(thetaPre[df_fac$r==j,]*df_fac[df_fac$r==j,5:(4+ntheta)])
  lpred_f <- spred + fpred + gammaVCMEst[Facid,j] + etaPre[df_fac$r==j] + bijPre[df_fac$r==j]
  pred_f <- c(pred_f,mean(invLogit(lpred_f)))
  strue <- rowSums(betaTrue[df_fac$r==j,]*df_fac[df_fac$r==j,(5+ntheta):(4+ntheta+nbeta)])
  ftrue <- rowSums(thetaTrue[df_fac$r==j,]*df_fac[df_fac$r==j,5:(4+ntheta)])
  ltrue_f <- strue + ftrue + gammaF[j] + etaTrue[df_fac$r==j] + bijTrue[df_fac$r==j]
  true_f <- c(true_f,mean(invLogit(ltrue_f)))
}  




#############################################################################
# 6. Visualization of VCMMR results
#############################################################################  

# Form 95% bootstrap confidence intervals for varying coefficient functions
theta1Boot = boot_Est[[1]]  
theta2Boot = boot_Est[[2]] 
beta1Boot = boot_Est[[3]]
beta2Boot = boot_Est[[4]]
etaBoot = boot_Est[[5]] 

theta1Boot_u <- apply(theta1Boot,2,quantile,0.975)
theta1Boot_l <- apply(theta1Boot,2,quantile,0.025)
theta2Boot_u <- apply(theta2Boot,2,quantile,0.975)
theta2Boot_l <- apply(theta2Boot,2,quantile,0.025)
beta1Boot_u <- apply(beta1Boot,2,quantile,0.975)
beta1Boot_l <- apply(beta1Boot,2,quantile,0.025)
beta2Boot_u <- apply(beta2Boot,2,quantile,0.975)
beta2Boot_l <- apply(beta2Boot,2,quantile,0.025)
etaBoot_u <- apply(etaBoot,2,quantile,0.975)
etaBoot_l <- apply(etaBoot,2,quantile,0.025)



# Plot estimates of varying coefficient functions and bootstrap confidence intervals
par(mfrow=c(3,2))

# Plot 95% confidence interval for theta1(t)
u <- max(theta1Boot_u)
l <- min(theta1Boot_l)
plot(gridPoints,thetaVCMEst[,1],'l',col="black",lwd=2,main="(a)", ylim = c(l,u),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(theta)[1](t)), line=2, cex.lab=1.5)
lines(gridPoints,theta1F,lwd=2,lty=1,col="grey")
lines(gridPoints,theta1Boot_l,col="black",lty=2,lwd=2) 
lines(gridPoints,theta1Boot_u,col="black",lty=2,lwd=2) 

# Plot 95% confidence interval for theta2(t)
u <- max(theta2Boot_u)
l <- min(theta2Boot_l)
plot(gridPoints,thetaVCMEst[,2],'l',col="black",lwd=2,main="(b)", ylim = c(l,u),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(theta)[2](t)), line=2, cex.lab=1.5)
lines(gridPoints,theta2F,lwd=2,lty=1,col="grey")
lines(gridPoints,theta2Boot_l,col="black",lty=2,lwd=2) 
lines(gridPoints,theta2Boot_u,col="black",lty=2,lwd=2) 

# Plot 95% confidence interval for beta1(t)
u <- max(beta1Boot_u)
l <- min(beta1Boot_l)
plot(gridPoints,betaVCMEst[,1],'l',col="black",lwd=2,main="(c)", ylim = c(l,u),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(beta)[1](t)), line=2, cex.lab=1.5)
lines(gridPoints,beta1F,lwd=2,lty=1,col="grey")
lines(gridPoints,beta1Boot_l,col="black",lty=2,lwd=2) 
lines(gridPoints,beta1Boot_u,col="black",lty=2,lwd=2) 

# Plot 95% confidence interval for beta2(t)
u <- max(beta2Boot_u)
l <- min(beta2Boot_l)
plot(gridPoints,betaVCMEst[,2],'l',col="black",lwd=2,main="(d)", ylim = c(l,u),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(beta)[2](t)), line=2, cex.lab=1.5)
lines(gridPoints,beta2F,lwd=2,lty=1,col="grey")
lines(gridPoints,beta2Boot_l,col="black",lty=2,lwd=2) 
lines(gridPoints,beta2Boot_u,col="black",lty=2,lwd=2) 

# Plot 95% confidence interval for eta(c)
u <- max(etaBoot_u)
l <- min(etaBoot_l)
plot(etaGrid,etaVCMEst,'l',col="black",lwd=2,main="(e)", ylim = c(l,u),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "c",ylab=expression(widehat(eta)(c)), line=2, cex.lab=1.5)
lines(etaGrid,etaF,lwd=2,lty=1,col="grey")
lines(etaGrid,etaBoot_l,col="black",lty=2,lwd=2) 
lines(etaGrid,etaBoot_u,col="black",lty=2,lwd=2) 

# Print the estimated variance of subject-specific random effects
cat(paste("Estimated variance of subject-specific random effects:", sigma2bEst), "\n", paste("True variance of subject-specific random effects:", 1.3))

dev.off()
# Visualization of the patient-level predicted risk trajectory for the first patient
plot(gridPoints,true_pat,"l",lwd=2,ylim = c(0,1),main = "Patient-level prediction",xlim=c(0,1),xaxs = "i",cex.main=1.5,xlab = "",ylab = "")
title(xlab = "t",ylab=expression(paste(widehat("p"),"'")[ij](t)), line=2, cex.lab=1.2)
lines(gridPoints,pred_pat,col = "grey")
legend("topleft",c("true risk","predicted risk"),lty = c(1,1),col = c("black","grey"),lwd=2,cex = 1,bty = "n",y.intersp=0.8)

# Visualization of the facility-level predicted risk trajectory for the first facility
plot(gridPoints[1:lfollow],true_f,"l",ylim = c(0,1),lwd=2,main = "Facility-level prediction",xlim=c(0,1),xaxs = "i", xlab="",ylab="", cex.main=1.5) 
title(xlab = "t",ylab=expression(paste(widehat("p"),"''")[i](t)), line=2, cex.lab=1.2)
lines(gridPoints[1:lfollow],pred_f,col = "grey")
legend("topleft",c("true risk","predicted risk"),lty = c(1,1),col = c("black","grey"),lwd=2,cex = 1,bty = "n",y.intersp=0.8)


 