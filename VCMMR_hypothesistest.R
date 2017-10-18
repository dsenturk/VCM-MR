VCMMR_hypothesistest <- function(Fid,        # id of the facility to be tested
                                 data,       # data.frame in long format
                                             # and row length equal to the length of the vectorized observations across all 
                                             # subjects and facilities
                                             # DATA.FRAME COLUMNS: 
                                             # fid: facility IDs (vector of length sum(Nij))
                                             # sid: subject IDs (vector of length sum(Nij))
                                             # y: hospitalization outcome data (vector of length sum(Nij))
                                             # t: follow-up time (vector of length sum(Nij)) 
                                             # z1: facility-level covariate (vector of length sum(Nij))
                                             # z2: facility-level covariate (vector of length sum(Nij))
                                             # x1: subject-level covariate (vector of length sum(Nij))
                                             # x2: subject-level covariate (vector of length sum(Nij))
                                             # c: calendar time at initiation of dialysis (vector of length sum(Nij))
                                             # size: facility sizes (vector of length sum(Nij))
                                             # bi: subject-specific random effects (vector of length sum(Nij))
                                 
                                 VCMMREst,   # VCMMR output from VCMMR_estimation.R
                                 numBS,      # number of resamples within facility (scalar) 
                                 hGamma      # bandwidth for estimating gammai(t)
                                             # of large, medium and small facilities (vector of length 3)
                                 ){
  
  #############################################################################
  ## Description: Function for performing bootstrap hypothesis testing for H0: gammai(t)=0 for the significance 
  ##              of the facility-specific fixed effects as described in Web Appendix B of "Modeling Time-Varying 
  ##              Effects of Multilevel Risk Factors of Hospitalizations in Patients on Dialysis". 
  ## Args:        see above
  ## Returns:     p-value: A p-value for the testing H0: gammai(t)=0 by resampling from within facility data. (scalar)
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
  
  # Number of subject-level risk factors
  nbeta <- 2
  
  # Number of facility-level risk factors
  ntheta <- 2
  
  # Define the grid points used to estimate gammai(t)
  gridPoints <- VCMMREst[[7]]
  ngrid <- length(gridPoints)
  
  # Maximum number of iterations allowed for convergence
  maxNRIter <- 50
  
  # Get the weights used in the Gauss-Hermite quadrature integral approximation
  ghrule <- gauss.quad(20,kind="hermite")
  numQuadPoints <- length(ghrule$nodes)
  
  # Read data from input
  df <- data
  
  # Bandwidth for estimating gammai(t)
  df$h <- (df$size == 1)*hGamma[1] + (df$size == 2)*hGamma[2] + (df$size == 3)*hGamma[3]
  
  # Index for updating varying coefficient functions
  df$r <- round(19*df$t)+1 
  df$cr <- round(19*df$c)+1 
  
  numF <- max(df$fid) # Number of facilities
  sumNi <- max(df$sid) # Number of subjects
  numOB <- as.numeric(table(df$sid)) # Number of observations per subject 
  Nij <- rep(numOB,numOB)   
  numDisPF <- as.numeric(table(df$fid)) 
  dft1 <- df[df$r==1,]
  numSubPF <- aggregate(dft1$sid,by=list(dft1$fid),FUN=length)[,2]
  df$numSubPF <- rep(numSubPF,numDisPF) 
  df$numOB <- rep(numOB,numOB)
  
  #######################################################################################################
  # Implement the hypothesis testing algorithm for testing H0: gammai(t)=0 as described in Web Appendix B
  #######################################################################################################
  
  # a. Read estimation results from input and calculate the test statistic riO for the facility specified for hypothesis testing in the input by Fid
  gammaVCMEst <- VCMMREst[[1]]
  thetaEstM <- VCMMREst[[2]][df$r,]
  betaEstM <- VCMMREst[[3]][df$r,]
  df$etaEst0 <- VCMMREst[[4]][df$cr]
  
  # Extract within facility data
  dfT <- df[df$fid==Fid,]
  covariates <- dfT[ diff(c(0,dfT$sid)) != 0, ]
  sumN1 <- covariates$numSubPF[1]
  numOBFac1 <- covariates$numOB
  
  riO <- (getArea(gridPoints,(gammaVCMEst[Fid,]-0)^2))^(1/2) # riO
  
  # b. Resample subject-specific random effects from the posterior distribution under H0
  ri <- rep(NA,numBS) # Test statistic under H0
  oldlpred <- rowSums(thetaEstM[df$fid==Fid,]*dfT[,5:(4+ntheta)]) + rowSums(betaEstM[df$fid==Fid,]*dfT[,(5+ntheta):(4+ntheta+nbeta)])
  # Calculate posterior means bij0 and variances vij0 under H0 
  # Gauss-Hermite quadrature
  temp <- matrix(NA, dim(dfT)[1], numQuadPoints)
  for(i in 1:numQuadPoints) {
    x <- ghrule$nodes[i]*sqrt(2*sigma2bEst)  
    temp[,i] <- exp(dfT$y*(x + oldlpred + dfT$etaEst0))/(1+exp(x + oldlpred + dfT$etaEst0))
  }
  ttt <- aggregate(temp,by=list(dfT$sid),FUN=prod)[-1]
  # Posterior means bij0 and variances vij0 under H0 
  meansLij <- rep(0,sumN1)
  means2Lij <- rep(0,sumN1)
  const <- rep(0,sumN1)
  for(i in 1:numQuadPoints) {
    p <- ttt[,i]*ghrule$weights[i]
    const <- const + p
    meansLij <- meansLij + ghrule$nodes[i]*p
    means2Lij <- means2Lij + ghrule$nodes[i]^2*p
  }
  meansLijH0 <- meansLij*sqrt(2*sigma2bEst)/const
  means2LijH0 <- means2Lij*2*sigma2bEst/const # bij0
  varsLijH0 <- means2LijH0 - meansLijH0^2 # vij0
  
  for(bs in 1:numBS) {
    subRandomEffectsBS <- rnorm(sumN1,meansLijH0,sqrt(varsLijH0)) # Resample subject-specific effects 
    
  # c. Generate outcome from a Bernoulli distribution under H0
    dfT$yb <- rbinom(dim(dfT)[1],1,invLogit(0+rep(subRandomEffectsBS,numOBFac1)+oldlpred + dfT$etaEst0)) 
    dfT$bisEst <- rep(0, dim(dfT)[1])
    dfT$gammaEst0 <- rep(0, dim(dfT)[1])
    dfT$gammaEst1 <- rep(0, dim(dfT)[1])
    gammaVCMEstTemp <- rep(0,ngrid)
    gammaVCMEst1Temp <- rep(0,ngrid)
    
  # d. Estimate bij0, vij0, gammai(t) and ri based on the resampled dataset 
    diff <- 1
    numIter <- 1
    while(diff > .0001 & numIter <= maxNRIter) {
      oldPijks <- invLogit(dfT$gammaEst0 + dfT$bisEst + oldlpred + dfT$etaEst0) # P0ijk from the previous iteration
      # Gauss-Hermite quadrature integral approximation
      temp <- matrix(NA, dim(dfT)[1], numQuadPoints)
      for(i in 1:numQuadPoints) {
        x <- ghrule$nodes[i]*sqrt(2*sigma2bEst)  
        temp[,i] <- exp(dfT$yb*(dfT$gammaEst0 + x + oldlpred + dfT$etaEst0))/(1+exp(dfT$gammaEst0 + x + oldlpred + dfT$etaEst0))
      }
      ttt <- aggregate(temp,by=list(dfT$sid),FUN=prod)[-1]
      
      # Update posterior means bij0 and variances vij0 based on the resampled dataset 
      meansLij2 <- rep(0,sumN1)
      means2Lij2 <- rep(0,sumN1)
      const2 <- rep(0,sumN1)
      for(i in 1:numQuadPoints) {
        p <- ttt[,i]*ghrule$weights[i]
        const2 <- const2 + p
        meansLij2 <- meansLij2 + ghrule$nodes[i]*p
        means2Lij2 <- means2Lij2 + ghrule$nodes[i]^2*p
      }
      meansLij2 <- meansLij2*sqrt(2*sigma2bEst)/const2
      means2Lij2 <- means2Lij2*2*sigma2bEst/const2
      varsLij2 <- means2Lij2 - meansLij2^2
      dfT$bisEst <- rep(meansLij2,numOBFac1) # bij0
      dfT$v <- rep(varsLij2,numOBFac1) # vij0
      
      # Update gammai(t) by maximizing the approximate expected local log-likelihood 
      # of the data from the facility based on the resampled dataset 
      for(j in 1:ngrid) {
        to <- gridPoints[j]
        indices <- (abs(dfT$t-to) < dfT$h)
        dfTemp <- dfT[indices,] # Data used for local MLE
        lpredTemp <- oldlpred[indices]
        tDiff <- dfTemp$t - to
        pijks <- invLogit(gammaVCMEstTemp[j] + gammaVCMEst1Temp[j]*tDiff + dfTemp$bisEst + lpredTemp + dfTemp$etaEst0)
        qijks <- 1 - pijks
        epan <- tDiff/dfTemp$h
        kernels <- 0.75*(1-epan^2)
        kernels <- (kernels > 0)*kernels/dfTemp$h # Kernel function
        a1 <- pijks*qijks + dfTemp$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
        a2 <- dfTemp$v/2*(pijks*qijks^2 - pijks^2*qijks)
        # Construct Hessian and Gradient 
        h11 <- sum(a1*kernels)
        h12 <- sum(a1*kernels*tDiff)
        h22 <- sum(a1*kernels*tDiff^2)
        v1 <- sum((dfTemp$yb - pijks - a2)*kernels)
        v2 <- sum((dfTemp$yb - pijks - a2)*kernels*tDiff)
        detVal <- 1/(h11*h22-h12^2)
        gammaVCMEstTemp[j] <- gammaVCMEstTemp[j] + detVal*(h22*v1-h12*v2) # Update gammai(t)
        gammaVCMEst1Temp[j] <- gammaVCMEst1Temp[j] + detVal*(-h12*v1+h11*v2)
        indices <- dfT$r==j
        dfT$gammaEst0[indices] <- gammaVCMEstTemp[j]
        dfT$gammaEst1[indices] <- gammaVCMEst1Temp[j]
      }
      newPijks <- invLogit(dfT$gammaEst0 + dfT$bisEst + oldlpred + dfT$etaEst0) # P0ijk from the current iteration
      diff <- max(newPijks - oldPijks)  
      numIter <- numIter + 1
    }
    # Calculate ri based on the resampled dataset
    ri[bs] <- (getArea(gridPoints,(gammaVCMEstTemp-0)^2))^(1/2) 
  }
 # e. Calculate the nominal p-value
  pVal <- mean(riO < ri) 
  return(pVal)
}