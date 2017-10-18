VCMMR_estimation <- function(data,        # data.frame in long format with eleven labeled columns (described below)
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
                             
                             sigma2,      # initial value for estimating variance of subject-specific random effects (scalar)
                             hThetaBeta,  # bandwidth for estimating theta(t) and beta(t) (scalar)
                             hEta,        # bandwidth for estimating eta(c) (scalar)
                             hGamma       # bandwidth for estimating gammai(t)
                                          # of large, medium and small facilities (vector of length 3)
                             ){
  
  #############################################################################
  ## Description: Function for estimation of VCM-MR model parameters described in "Modeling Time-Varying Effects 
  ##              of Multilevel Risk Factors of Hospitalizations in Patients on Dialysis", including estimation of facility-specific fixed effects,
  ##              time-varying effects of multilevel risk factors, calendar year effect and variance of subject specific random effects. 
  ## Args:        see above
  ## Returns:     list()
  ##              gamma: estimated facility-specific fixed effects (matrix of dimension numF*20)
  ##              theta: estimated facility-level risk factor effects (matrix of dimension 2*20)
  ##              beta: estimated subject-level risk factor effects (matrix of dimension 2*20)
  ##              eta: estimated calendar year effect (vector of length 20) 
  ##              sigma: estimated variance of subject specific random effects (scalar)
  ##              bijEst: estimated posterior mean of subject-level random effects (vector of length sum(Ni))
  ##              grid: grid points used to estimate the varying coefficient functions gammai(t), theta(t) and beta(t) (vector of length 20)
  ##              etagrid: grid points used to estimate eta(c) (vector of length 20)
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

  # Define functions 
  invLogit <- function(x) {
    return(exp(x)/(1+exp(x)))
  }
  getArea <- function(x,y) {
    if(length(x) != length(y)) {
      stop("The number of elements in the domain and range must equal.")
    }
    if(sum(sort(x) != x) > 0) {
      stop("The endpoints of the trapezoids must be ordered.")
    }
    k <- length(x)
    return(1/2*sum((y[-k]+y[-1])*(x[-1]-x[-k])))
  }
  
  # Create a list for returns as outlined above
  VCMMREst <- list(8)
  
  # Number of subject-level risk factors
  nbeta <- 2
  
  # Number of facility-level risk factors
  ntheta <- 2
  
  # Empty matrix for constructing Hessian matrix of theta(t) and beta(t) used for estimation step 5 in Web Appendix A
  h <- diag(nbeta+ntheta)   
  
  # Define the grid points used to estimate the varying coefficient functions gammai(t), theta(t) and beta(t)
  gridPoints <- seq(0,1,1/19)
  ngrid <- length(gridPoints)
  
  # Grid points used to estimate eta(c)
  etaGrid <- seq(0,1,1/19)
  nEtagrid <- length(etaGrid)
  
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

  

  
  ###########################################################################
  # Implement the approximate EM algorithm as described in Web Appendix A
  ###########################################################################
  
  # Find initial values for theta(t) and beta(t) from the non-time-varying risk effects
  # generalized linear model
  
  # Initial values for fitting non-time-varying risk effects generalized linear model
  df$bisEst <- rep(0, dim(df)[1])
  df$gammaEst0 <- rep(0, dim(df)[1])
  df$gammaEst1 <- rep(0, dim(df)[1])
  betaEsts <- rep(0,nbeta)
  thetaEsts <- rep(0,ntheta)
  df$etaEst0 <- rep(0, dim(df)[1])
  df$etaEst1 <- rep(0, dim(df)[1])
  gammaVCMEst <- matrix(0,numF,ngrid)
  gammaVCMEst1 <- matrix(0,numF,ngrid)
  sigma2bEst <- sigma2
  
  numIter <- 1
  Diff <- 1  
  while(Diff > 0.001 & numIter <= maxNRIter) { 
    oldlpred <- betaEsts[1]*df$x1 + betaEsts[2]*df$x2 + thetaEsts[1]*df$z1 + thetaEsts[2]*df$z2
    oldPijks <- invLogit(df$gammaEst0 + df$bisEst + oldlpred + df$etaEst0) # Pijk from the previous iteration
    
    # Gauss-Hermite quadrature
    temp <- matrix(NA, dim(df)[1], numQuadPoints)
    for(i in 1:numQuadPoints) {
      x <- ghrule$nodes[i]*sqrt(2*sigma2bEst)  
      temp[,i] <- exp(df$y*(df$gammaEst0 + x + oldlpred + df$etaEst0))/(1+exp(df$gammaEst0 + x + oldlpred + df$etaEst0))
    }
    ttt <- aggregate(temp,by=list(df$sid),FUN=prod)[-1]
   
    # Update sigma2 in the non-time-varying risk effects generalized linear model
    meansLij <- rep(0,sumNi)
    means2Lij <- rep(0,sumNi)
    const <- rep(0,sumNi)
    for(i in 1:numQuadPoints) {
      p <- ttt[,i]*ghrule$weights[i]
      const <- const + p
      meansLij <- meansLij + ghrule$nodes[i]*p
      means2Lij <- means2Lij + ghrule$nodes[i]^2*p
    }
    meansLij <- meansLij*sqrt(2*sigma2bEst)/const
    means2Lij <- means2Lij*2*sigma2bEst/const
    varsLij <- means2Lij - meansLij^2
    sigma2bEst <- 1/sumNi*sum(means2Lij)
    df$bisEst <- rep(meansLij,numOB) 
    df$v <- rep(varsLij,numOB)
    
    # Update gammai(t) in the non-time-varying risk effects generalized linear model
    for(j in 1:ngrid) {
      to <- gridPoints[j]
      indices <- (abs(df$t-to) < df$h)
      dfTemp <- df[indices,] 
      tDiff <- dfTemp$t - to
      lpredTemp <- oldlpred[indices]
      obj <- gammaVCMEst[dfTemp$fid,j] + gammaVCMEst1[dfTemp$fid,j]*tDiff + dfTemp$bisEst + lpredTemp + dfTemp$etaEst0
      pijks <- invLogit(obj)
      qijks <- 1 - pijks
      epan <- tDiff/dfTemp$h
      kernels <- 0.75*(1-epan^2)
      kernels <- (kernels > 0)*kernels/dfTemp$h
      # Construct Hessian and Gradient for gammai(t) in the non-time-varying risk effects generalized linear model
      a2 <- pijks*qijks + dfTemp$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
      a1 <- dfTemp$v/2*(pijks*qijks^2 - pijks^2*qijks)
      ttt <- cbind(a2*kernels, a2*kernels*tDiff, a2*kernels*tDiff^2,(dfTemp$y - pijks - a1)*kernels,(dfTemp$y - pijks - a1)*kernels*tDiff)
      aggTTT <- aggregate(ttt,by=list(dfTemp$fid),FUN=sum)
      h11 <- aggTTT[,2]
      h12 <- aggTTT[,3]
      h22 <- aggTTT[,4]
      v1 <- aggTTT[,5]
      v2 <- aggTTT[,6]
      detVal <- 1/(h11*h22-h12^2)
      gammaVCMEst[,j] <- gammaVCMEst[,j] + detVal*(h22*v1-h12*v2) # Update gammai(t)
      gammaVCMEst1[,j] <- gammaVCMEst1[,j] + detVal*(-h12*v1+h11*v2) 
      indices <- df$r==j
      df$gammaEst0[indices] <- gammaVCMEst[df$fid[indices],j]
      df$gammaEst1[indices] <- gammaVCMEst1[df$fid[indices],j]
    }   
    
    # Update theta and beta in the non-time-varying risk effects generalized linear model
    pijks <- invLogit(df$gammaEst0 + df$bisEst +  oldlpred + df$etaEst0)
    qijks <- 1 - pijks 
    a1 <- df$v/2*(pijks*qijks^2 - pijks^2*qijks)
    a2 <- pijks*qijks + df$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
    zxMatrix <- data.matrix(df[,5:(4+nbeta+ntheta)])
    zxVector <- t(zxMatrix)%*%(df$y - pijks - a1) # Gradient for theta and beta
    zzxxMatrix <- t(as.matrix(zxMatrix))
    for(mIndex in 1:(nbeta+ntheta)) {
      zzxxMatrix[mIndex,] <- zzxxMatrix[mIndex,]*a2
    }
    zzxxMatrix <- zzxxMatrix%*%zxMatrix # Hessian for theta and beta
    Phi_OneStep <- solve(zzxxMatrix) %*% zxVector
    thetaEsts <- matrix(thetaEsts,ntheta,1) + Phi_OneStep[1:ntheta]
    betaEsts <- matrix(betaEsts,nbeta,1) + Phi_OneStep[(ntheta+1):(ntheta+nbeta)]
    newlpred <- betaEsts[1]*df$x1 + betaEsts[2]*df$x2 + thetaEsts[1]*df$z1 + thetaEsts[2]*df$z2
    
    # Update eta(c) in the non-time-varying risk effects generalized linear model
    step0 <- rep(0,nEtagrid)
    step1 <- rep(0,nEtagrid)
    for(j in 1:nEtagrid) {
      co <- etaGrid[j]
      indices2 <- (abs(df$c-co) < hEta)
      dfTemp <- df[indices2,]
      cDiff <- dfTemp$c - co
      obj <- dfTemp$gammaEst0 + dfTemp$bisEst + newlpred[indices2] + dfTemp$etaEst0 + dfTemp$etaEst1*cDiff
      pijks <- invLogit(obj)
      qijks <- 1 - pijks
      epan <- cDiff/hEta
      kernels <- 0.75*(1-epan^2)
      kernels <- (kernels > 0)*kernels/hEta
      # Construct Hessian and Gradient for eta(c) in the non-time-varying risk effects generalized linear model
      a2 <- pijks*qijks + dfTemp$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
      a1 <- dfTemp$v/2*(pijks*qijks^2 - pijks^2*qijks)
      ttt <- cbind(a2*kernels, a2*kernels*cDiff, a2*kernels*cDiff^2,(dfTemp$y - pijks - a1)*kernels,(dfTemp$y - pijks - a1)*kernels*cDiff)
      aggTTT <- colSums(ttt)
      h11 <- aggTTT[1]
      h12 <- aggTTT[2]
      h22 <- aggTTT[3]
      v1 <- aggTTT[4]
      v2 <- aggTTT[5]
      detVal <- 1/(h11*h22-h12^2)
      step0[j] <- detVal*(h22*v1-h12*v2)
      step1[j] <- detVal*(-h12*v1+h11*v2)
    }
    
    # Begin line search for eta(c) in the non-time-varying risk effects generalized linear model
    Olddf <- df
    OldObj <-df$gammaEst0 + df$bisEst + newlpred + df$etaEst0
    OldPijks <- invLogit(OldObj)
    OldQijks <- 1-OldPijks
    OldlogLikelihood <- sum(OldObj*df$y+log(OldQijks)-df$v*OldPijks*OldQijks/2-(df$bisEst^2+df$v)/(2*sigma2bEst*Nij)-0.5*log(2*pi*sigma2bEst)/Nij)
    for(s in 1:10){
      smin <- (1 / 2) ^ (s - 1)
      for(j in 1:nEtagrid){
        indices <- df$cr==j
        df$etaEst0[indices] <- Olddf$etaEst0[indices] +  smin*step0[j]
        df$etaEst1[indices] <- Olddf$etaEst1[indices] +  smin*step1[j]
      }
      newObj <-df$gammaEst0 + df$bisEst + newlpred + df$etaEst0
      newPijks <- invLogit(newObj)
      newQijks <- 1-newPijks
      newlogLikelihood <- sum(newObj*df$y+log(newQijks)-df$v*newPijks*newQijks/2-(df$bisEst^2+df$v)/(2*sigma2bEst*Nij)-0.5*log(2*pi*sigma2bEst)/Nij)
      if(newlogLikelihood > OldlogLikelihood){
        break
      }
    } 
    # End line search
    etaVCMEst <- aggregate(df$etaEst0,by=list(df$cr),FUN=mean)[,2]
    df$etaEst0 <- df$etaEst0 - 1/2*sum((etaVCMEst[-1]+etaVCMEst[-nEtagrid]))/(nEtagrid-1) # Normalize eta(c)
    
    newObj <-df$gammaEst0 + df$bisEst + newlpred + df$etaEst0
    newPijks <- invLogit(newObj) # Pijk from the current iteration
    Diff <- max(newPijks - oldPijks)   
    numIter <- numIter + 1   
  }
  # End of finding initial values for theta(t) and beta(t)
  
  # Implement estimation steps as described in Web Appendix A
  # Step 1: set initial values for fitting VCM-MR
  df$bisEst <- rep(0, dim(df)[1])  
  df$etaEst0 <- rep(0, dim(df)[1])
  df$etaEst1 <- rep(0, dim(df)[1]) 
  df$gammaEst0 <- rep(0, dim(df)[1])
  df$gammaEst1 <- rep(0, dim(df)[1])
  gammaVCMEst <- matrix(0,numF,ngrid)
  gammaVCMEst1 <- matrix(0,numF,ngrid)
  thetaEstM <- t(matrix(thetaEsts,ntheta,dim(df)[1]))
  thetaEstM1 <- t(matrix(0,ntheta,dim(df)[1]))
  betaEstM <- t(matrix(betaEsts,nbeta,dim(df)[1]))
  betaEstM1 <- t(matrix(0,nbeta,dim(df)[1]))
  sigma2bEst <- sigma2
  
  numIter <- 0
  Diff <- 1
  while(Diff > 0.001 & numIter <= maxNRIter) {
    print(paste("Iter",numIter)) # Print number of iterations

    oldlpred <- rowSums(thetaEstM*df[,5:(4+ntheta)]) + rowSums(betaEstM*df[,(5+ntheta):(4+ntheta+nbeta)])
    oldObj <- df$gammaEst0 + df$bisEst + oldlpred + df$etaEst0
    oldPijks <- invLogit(oldObj) # P0ijk from the previous iteration
    
    # Step 2 (E-step): update posterior means bij0 and variances vij0 of subject-specific random effects
    # Gauss-Hermite quadrature integral approximation
    temp <- matrix(NA, dim(df)[1], numQuadPoints) 
    for(i in 1:numQuadPoints) {
      x <- ghrule$nodes[i]*sqrt(2*sigma2bEst)  
      temp[,i] <- exp(df$y*(df$gammaEst0 + x + oldlpred + df$etaEst0))/(1+exp(df$gammaEst0 + x + oldlpred + df$etaEst0))
    }
    ttt <- aggregate(temp,by=list(df$sid),FUN=prod)[-1]
    
    # Update posterior means bij0 and variances vij0 
    meansLij <- rep(0,sumNi)
    means2Lij <- rep(0,sumNi)
    const <- rep(0,sumNi)
    for(i in 1:numQuadPoints) {
      p <- ttt[,i]*ghrule$weights[i]
      const <- const + p
      meansLij <- meansLij + ghrule$nodes[i]*p
      means2Lij <- means2Lij + ghrule$nodes[i]^2*p
    }
    meansLij <- meansLij*sqrt(2*sigma2bEst)/const
    means2Lij <- means2Lij*2*sigma2bEst/const
    varsLij <- means2Lij - meansLij^2
    SumMean2 <- sum(means2Lij)
    df$bisEst <- rep(meansLij,numOB) # bij0
    df$v <- rep(varsLij,numOB)  # vij0
    
    # Step 3 (M-step): update sigma2b by maximizing the approximate expected log-likelihood
    sigma2bEst <- 1/sumNi*SumMean2

    # Step 4 (M-step): update gammai(t) by maximizing the approximate expected local log-likelihood 
    # of the data from the ith facility
    for(j in 1:ngrid) {
      to <- gridPoints[j]
      indices <- (abs(df$t-to) < df$h)
      dfTemp <- df[indices,] # Data used for local MLE
      lpredTemp <- oldlpred[indices]
      tDiff <- dfTemp$t - to
      obj <- gammaVCMEst[dfTemp$fid,j] + gammaVCMEst1[dfTemp$fid,j]*tDiff + dfTemp$bisEst + lpredTemp + dfTemp$etaEst0
      pijks <- invLogit(obj)
      qijks <- 1 - pijks
      epan <- tDiff/dfTemp$h
      kernels <- 0.75*(1-epan^2) 
      kernels <- (kernels > 0)*kernels/dfTemp$h # Kernel function
      a2 <- pijks*qijks + dfTemp$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
      a1 <- dfTemp$v/2*(pijks*qijks^2 - pijks^2*qijks)
      # Construct Hessian and Gradient 
      ttt <- cbind(a2*kernels, a2*kernels*tDiff, a2*kernels*tDiff^2,(dfTemp$y - pijks - a1)*kernels,(dfTemp$y - pijks - a1)*kernels*tDiff)
      aggTTT <- aggregate(ttt,by=list(dfTemp$fid),FUN=sum) # Use only within facility data
      h11 <- aggTTT[,2]
      h12 <- aggTTT[,3]
      h22 <- aggTTT[,4]
      v1 <- aggTTT[,5]
      v2 <- aggTTT[,6]
      detVal <- 1/(h11*h22-h12^2)
      gammaVCMEst[,j] <- gammaVCMEst[,j] + detVal*(h22*v1-h12*v2) # Update gammai(t)
      gammaVCMEst1[,j] <- gammaVCMEst1[,j] + detVal*(-h12*v1+h11*v2) 
      indices <- df$r==j
      df$gammaEst0[indices] <- gammaVCMEst[df$fid[indices],j]
      df$gammaEst1[indices] <- gammaVCMEst1[df$fid[indices],j]
    }  
    
    # Step 5 (M-step): update theta(t) and beta(t) by maximizing the approximate expected local log-likelihood
    Phi_OneStep <- matrix(0,ngrid,2*(ntheta+nbeta))
    for(j in 1:ngrid) {
      to <- gridPoints[j]
      indices1 <- (abs(df$t-to) < hThetaBeta)
      dfTemp <- df[indices1,] # Data used for local MLE
      tDiff <- dfTemp$t - to
      lpredTemp <- oldlpred[indices1] + rowSums(thetaEstM1[indices1,]*dfTemp[,5:(4+ntheta)]*tDiff) + rowSums(betaEstM1[indices1,]*dfTemp[,(5+ntheta):(4+ntheta+nbeta)]*tDiff)
      obj <- dfTemp$gammaEst0 + dfTemp$bisEst + lpredTemp + dfTemp$etaEst0
      pijks <- invLogit(obj)
      qijks <- 1 - pijks
      epan <- tDiff/hThetaBeta
      kernels <- 0.75*(1-epan^2)
      kernels <- (kernels > 0)*kernels/hThetaBeta # Kernel function
      a1 <- dfTemp$v/2*(pijks*qijks^2 - pijks^2*qijks)
      a2 <- pijks*qijks + dfTemp$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
      # Construct Gradient
      zxMatrix <- data.matrix(cbind(dfTemp[,5:(4+ntheta+nbeta)],dfTemp[,5:(4+ntheta+nbeta)]*tDiff))
      zxVector <- t(zxMatrix) %*% ((dfTemp$y - pijks - a1)*kernels)
      # Construct Hessian
      zzxxMatrix <- t(as.matrix(zxMatrix))      
      for(mIndex in 1:(2*(ntheta+nbeta))) {
        zzxxMatrix[mIndex,] <- zzxxMatrix[mIndex,]*a2*kernels
      }
      zzxxMatrix <- zzxxMatrix%*%zxMatrix
      # Newton-Raphson updating direction
      invHessian <- solve(zzxxMatrix)
      Phi_OneStep[j,] <- invHessian %*% zxVector
    } 
    # Begin line search for theta(t) and beta(t) as described in Web Appendix A
    # Calculate global log-likelihood
    OldbetaM <- betaEstM 
    OldthetaM <- thetaEstM 
    OldbetaM1 <- betaEstM1
    OldthetaM1 <- thetaEstM1
    OldObj <-df$gammaEst0 + df$bisEst + oldlpred + df$etaEst0
    OldPijks <- invLogit(OldObj)
    OldQijks <- 1-OldPijks
    OldlogLikelihood <- sum(OldObj*df$y+log(OldQijks)-df$v*OldPijks*OldQijks/2-(df$bisEst^2+df$v)/(2*sigma2bEst*Nij)-0.5*log(2*pi*sigma2bEst)/Nij)
    for(s in 1:10){
      # Step size
      ss <- (1 / 2) ^ (s - 1)
      # Same step size is used for all time points
      for(j in 1:ngrid) {
        indices <- df$r==j
        thetaEstM[indices,] <- t(t(OldthetaM[indices,]) + ss*Phi_OneStep[j,1:ntheta]) # Update theta(t)
        betaEstM[indices,] <- t(t(OldbetaM[indices,]) + ss*Phi_OneStep[j,(ntheta+1):(ntheta+nbeta)]) # Update beta(t)
        thetaEstM1[indices,] <- t(t(OldthetaM1[indices,]) + ss*Phi_OneStep[j,(ntheta+nbeta+1):(2*ntheta+nbeta)])
        betaEstM1[indices,] <- t(t(OldbetaM1[indices,]) + ss*Phi_OneStep[j,(2*ntheta+nbeta+1):(2*ntheta+2*nbeta)])
      } 
      # Calculate global log-likelihood with updated of theta(t) and beta(t)
      newlpred <- rowSums(thetaEstM*df[,5:(4+ntheta)]) + rowSums(betaEstM*df[,(5+ntheta):(4+ntheta+nbeta)])
      newObj <-df$gammaEst0 + df$bisEst + newlpred + df$etaEst0
      newPijks <- invLogit(newObj)
      newQijks <- 1-newPijks
      newlogLikelihood <- sum(newObj*df$y+log(newQijks)-df$v*newPijks*newQijks/2-(df$bisEst^2+df$v)/(2*sigma2bEst*Nij)-0.5*log(2*pi*sigma2bEst)/Nij)
      if(newlogLikelihood > OldlogLikelihood){ 
        break
      }
    }
    # End line search for theta(t) and beta(t)
    
    # Step 6 (M-step): update eta(c) by maximizing the approximate expected local likelihood
    step0 <- rep(0,nEtagrid)
    step1 <- rep(0,nEtagrid)
    for(j in 1:nEtagrid) {
      co <- etaGrid[j]
      indices2 <- (abs(df$c-co) < hEta)
      dfTemp <- df[indices2,]  # Data used for local MLE
      cDiff <- dfTemp$c - co
      obj <- dfTemp$gammaEst0 + dfTemp$bisEst + newlpred[indices2] + dfTemp$etaEst0 + dfTemp$etaEst1*cDiff
      pijks <- invLogit(obj)
      qijks <- 1 - pijks
      epan <- cDiff/hEta
      kernels <- 0.75*(1-epan^2)
      kernels <- (kernels > 0)*kernels/hEta # Kernel function
      a2 <- pijks*qijks + dfTemp$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
      a1 <- dfTemp$v/2*(pijks*qijks^2 - pijks^2*qijks)
      # Construct Hessian and Gradient
      ttt <- cbind(a2*kernels, a2*kernels*cDiff, a2*kernels*cDiff^2,(dfTemp$y - pijks - a1)*kernels,(dfTemp$y - pijks - a1)*kernels*cDiff)
      aggTTT <- colSums(ttt)
      h11 <- aggTTT[1]
      h12 <- aggTTT[2]
      h22 <- aggTTT[3]
      v1 <- aggTTT[4]
      v2 <- aggTTT[5]
      detVal <- 1/(h11*h22-h12^2)
      step0[j] <- detVal*(h22*v1-h12*v2)
      step1[j] <- detVal*(-h12*v1+h11*v2) 
    }
    # Begin line search for eta(c) as described in Web Appendix A
    # Calculate global log-likelihood
    Olddf <- df
    OldObj <-df$gammaEst0 + df$bisEst + newlpred + df$etaEst0
    OldPijks <- invLogit(OldObj)
    OldQijks <- 1-OldPijks
    OldlogLikelihood <- sum(OldObj*df$y+log(OldQijks)-df$v*OldPijks*OldQijks/2-(df$bisEst^2+df$v)/(2*sigma2bEst*Nij)-0.5*log(2*pi*sigma2bEst)/Nij)
    for(s in 1:10){
      ss <- (1 / 2) ^ (s - 1) 
      for(j in 1:nEtagrid){
        indices <- df$cr==j
        df$etaEst0[indices] <- Olddf$etaEst0[indices] +  ss*step0[j] # Update eta(c)
        df$etaEst1[indices] <- Olddf$etaEst1[indices] +  ss*step1[j]
      }
      # Calculate global log-likelihood with updated of eta(c)
      newObj <-df$gammaEst0 + df$bisEst + newlpred + df$etaEst0
      newPijks <- invLogit(newObj)
      newQijks <- 1-newPijks
      newlogLikelihood <- sum(newObj*df$y+log(newQijks)-df$v*newPijks*newQijks/2-(df$bisEst^2+df$v)/(2*sigma2bEst*Nij)-0.5*log(2*pi*sigma2bEst)/Nij)
      if(newlogLikelihood > OldlogLikelihood){
        break
      }
    } 
    # End line search for eta(c)
    etaVCMEst <- aggregate(df$etaEst0,by=list(df$cr),FUN=mean)[,2]
    df$etaEst0 <- df$etaEst0 - 1/2*sum((etaVCMEst[-1]+etaVCMEst[-nEtagrid]))/(nEtagrid-1) # Normalize eta(c)
    
    newObj <-df$gammaEst0 + df$bisEst + newlpred + df$etaEst0
    newPijks <- invLogit(newObj) # P0ijk from the current iteration
    Diff <- max(newPijks - oldPijks)   
    if(numIter == maxNRIter) { stop("WARNING: NO CONVERGENCE") } 
    print(paste("Diff=",Diff))
    numIter <- numIter + 1   
  }
  
  # Store the estimates of varying coefficient functions and variance of subject-specific random effects
  thetaVCMEst <- aggregate(thetaEstM,by=list(df$r),FUN=mean)[,2:3]
  betaVCMEst <- aggregate(betaEstM,by=list(df$r),FUN=mean)[,2:3]
  etaVCMEst <- aggregate(df$etaEst0,by=list(df$cr),FUN=mean)[,2]
  bijEst <- aggregate(df$bisEst,by=list(df$sid),FUN=mean)[,2]
  
  # Output estimation results
  VCMMREst[[1]] <- gammaVCMEst
  VCMMREst[[2]] <- thetaVCMEst
  VCMMREst[[3]] <- betaVCMEst
  VCMMREst[[4]] <- etaVCMEst
  VCMMREst[[5]] <- sigma2bEst
  VCMMREst[[6]] <- bijEst 
  VCMMREst[[7]] <- gridPoints
  VCMMREst[[8]] <- etaGrid
  names(VCMMREst) <- c("gamma", "theta", "beta", "eta", "sigma", "bijEst", "grid", "etagrid")
  return(VCMMREst)
}
 