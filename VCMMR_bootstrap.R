VCMMR_bootstrap <- function(nboot,       # number of bootstrap samples (scalar) 
                            data,        # data.frame in long format
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
  ## Description: Function for performing bootstrap inference for effects of multilevel factors 
  ##              and calendar year effects as described in "Modeling Time-Varying Effects 
  ##              of Multilevel Risk Factors of Hospitalizations in Patients on Dialysis". 
  ## Args:        see above
  ## Returns:     list()
  ##              theta1: estimated effect of facility-level covariate 1 from bootstrap samples (matrix of dimension nboot*20)
  ##              theta2: estimated effect of facility-level covariate 2 from bootstrap samples (matrix of dimension nboot*20)
  ##              beta1: estimated effect of subject-level covariate 1 from bootstrap samples (matrix of dimension nboot*20)
  ##              beta2: estimated effect of subject-level covariate 2 from bootstrap samples (matrix of dimension nboot*20)
  ##              eta: estimated calendar year effect from bootstrap samples (matrix of dimension nboot*20)
  ##              sigma: estimated variance of subject specific random effects from bootstrap samples (vector of length nboot)
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
  BootEst <- list(8)
  
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
  numDisPF <- as.numeric(table(df$fid)) 
  dft1 <- df[df$r==1,]
  numSubPF <- aggregate(dft1$sid,by=list(dft1$fid),FUN=length)[,2]
  df$numOB <- rep(numOB,numOB)

  # Start bootstrap 
  theta1Est_boot <- c()
  theta2Est_boot <- c()
  beta1Est_boot <- c()
  beta2Est_boot <- c() 
  etaEst_boot <- c()
  sigma2bEst_boot <- c()
  
  for(bt in 1:nboot){
    print(paste("Bootstrap Sample: ", bt))
    ###########################################################################
    # Generate a bootstrap sample as described in Section 2.3
    ###########################################################################
    boot_Fac <- sample(1:numF,numF,replace = TRUE) # Resample from facilities
    boot_numSubPF <- numSubPF[boot_Fac] 
    boot_sumNi <- sum(boot_numSubPF) # Number of subjects in the bootstrap sample
    boot_Ni <- rep(1:numF,boot_numSubPF)
    boot_index <- NULL
    boot_numOB <- NULL 
    for(i in 1:numF){
      if(boot_Fac[i] == 1){
        FacIndex <- 1:numDisPF[1]
      } else {
        FacIndex <- (sum(numDisPF[1:(boot_Fac[i]-1)]) + 1) : sum(numDisPF[1:(boot_Fac[i])])
      }
      boot_index <- c(boot_index, FacIndex)
      boot_numOB <- c(boot_numOB, numOB[unique(df[FacIndex,]$sid)])
    }
    df_boot <- df[boot_index,]
    df_boot$fid <- rep(rep(1:numF,boot_numSubPF),boot_numOB)
    df_boot$sid <- rep(1:boot_sumNi,boot_numOB)
    boot_Nij <- df_boot$numOB # Number of observations per subject in the bootstrap sample
    
    ###########################################################################
    # Implement the approximate EM algorithm as described in Web Appendix A
    ###########################################################################
    # Find initial values for theta(t) and beta(t) from the non-time-varying risk effects
    # generalized linear model
    
    # Initial values for fitting non-time-varying risk effects generalized linear model
    df_boot$bisEst <- rep(0, dim(df_boot)[1])
    df_boot$gammaEst0 <- rep(0, dim(df_boot)[1])
    df_boot$gammaEst1 <- rep(0, dim(df_boot)[1])
    betaEsts <- rep(0,nbeta)
    thetaEsts <- rep(0,ntheta)
    df_boot$etaEst0 <- rep(0, dim(df_boot)[1])
    df_boot$etaEst1 <- rep(0, dim(df_boot)[1])
    gammaVCMEst <- matrix(0,numF,ngrid)
    gammaVCMEst1 <- matrix(0,numF,ngrid)
    sigma2bEst <- sigma2
    
    numIter <- 1
    Diff <- 1  
    while(Diff > 0.001 & numIter <= maxNRIter) { 
      oldlpred <- betaEsts[1]*df_boot$x1 + betaEsts[2]*df_boot$x2 + thetaEsts[1]*df_boot$z1 + thetaEsts[2]*df_boot$z2
      oldPijks <- invLogit(df_boot$gammaEst0 + df_boot$bisEst + oldlpred + df_boot$etaEst0) # Pijk from the previous iteration
      
      # Gauss-Hermite quadrature
      temp <- matrix(NA, dim(df_boot)[1], numQuadPoints)
      for(i in 1:numQuadPoints) {
        x <- ghrule$nodes[i]*sqrt(2*sigma2bEst)  
        temp[,i] <- exp(df_boot$y*(df_boot$gammaEst0 + x + oldlpred + df_boot$etaEst0))/(1+exp(df_boot$gammaEst0 + x + oldlpred + df_boot$etaEst0))
      }
      ttt <- aggregate(temp,by=list(df_boot$sid),FUN=prod)[-1]
      
      # Update sigma2 in the non-time-varying risk effects generalized linear model
      meansLij <- rep(0,boot_sumNi)
      means2Lij <- rep(0,boot_sumNi)
      const <- rep(0,boot_sumNi)
      for(i in 1:numQuadPoints) {
        p <- ttt[,i]*ghrule$weights[i]
        const <- const + p
        meansLij <- meansLij + ghrule$nodes[i]*p
        means2Lij <- means2Lij + ghrule$nodes[i]^2*p
      }
      meansLij <- meansLij*sqrt(2*sigma2bEst)/const
      means2Lij <- means2Lij*2*sigma2bEst/const
      varsLij <- means2Lij - meansLij^2
      sigma2bEst <- 1/boot_sumNi*sum(means2Lij)
      df_boot$bisEst <- rep(meansLij,boot_numOB) 
      df_boot$v <- rep(varsLij,boot_numOB)
      
      # Update gammai(t) in the non-time-varying risk effects generalized linear model
      for(j in 1:ngrid) {
        to <- gridPoints[j]
        indices <- (abs(df_boot$t-to) < df_boot$h)
        df_bootTemp <- df_boot[indices,] 
        tDiff <- df_bootTemp$t - to
        lpredTemp <- oldlpred[indices]
        obj <- gammaVCMEst[df_bootTemp$fid,j] + gammaVCMEst1[df_bootTemp$fid,j]*tDiff + df_bootTemp$bisEst + lpredTemp + df_bootTemp$etaEst0
        pijks <- invLogit(obj)
        qijks <- 1 - pijks
        epan <- tDiff/df_bootTemp$h
        kernels <- 0.75*(1-epan^2)
        kernels <- (kernels > 0)*kernels/df_bootTemp$h
        # Construct Hessian and Gradient for gammai(t) in the non-time-varying risk effects generalized linear model
        a2 <- pijks*qijks + df_bootTemp$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
        a1 <- df_bootTemp$v/2*(pijks*qijks^2 - pijks^2*qijks)
        ttt <- cbind(a2*kernels, a2*kernels*tDiff, a2*kernels*tDiff^2,(df_bootTemp$y - pijks - a1)*kernels,(df_bootTemp$y - pijks - a1)*kernels*tDiff)
        aggTTT <- aggregate(ttt,by=list(df_bootTemp$fid),FUN=sum)
        h11 <- aggTTT[,2]
        h12 <- aggTTT[,3]
        h22 <- aggTTT[,4]
        v1 <- aggTTT[,5]
        v2 <- aggTTT[,6]
        detVal <- 1/(h11*h22-h12^2)
        gammaVCMEst[,j] <- gammaVCMEst[,j] + detVal*(h22*v1-h12*v2) # Update gammai(t)
        gammaVCMEst1[,j] <- gammaVCMEst1[,j] + detVal*(-h12*v1+h11*v2) 
        indices <- df_boot$r==j
        df_boot$gammaEst0[indices] <- gammaVCMEst[df_boot$fid[indices],j]
        df_boot$gammaEst1[indices] <- gammaVCMEst1[df_boot$fid[indices],j]
      }   
      
      # Update theta and beta in the non-time-varying risk effects generalized linear model
      pijks <- invLogit(df_boot$gammaEst0 + df_boot$bisEst +  oldlpred + df_boot$etaEst0)
      qijks <- 1 - pijks 
      a1 <- df_boot$v/2*(pijks*qijks^2 - pijks^2*qijks)
      a2 <- pijks*qijks + df_boot$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
      zxMatrix <- data.matrix(df_boot[,5:(4+nbeta+ntheta)])
      zxVector <- t(zxMatrix)%*%(df_boot$y - pijks - a1) # Gradient for theta and beta
      zzxxMatrix <- t(as.matrix(zxMatrix))
      for(mIndex in 1:(nbeta+ntheta)) {
        zzxxMatrix[mIndex,] <- zzxxMatrix[mIndex,]*a2
      }
      zzxxMatrix <- zzxxMatrix%*%zxMatrix # Hessian for theta and beta
      Phi_OneStep <- solve(zzxxMatrix) %*% zxVector
      thetaEsts <- matrix(thetaEsts,ntheta,1) + Phi_OneStep[1:ntheta]
      betaEsts <- matrix(betaEsts,nbeta,1) + Phi_OneStep[(ntheta+1):(ntheta+nbeta)]
      newlpred <- betaEsts[1]*df_boot$x1 + betaEsts[2]*df_boot$x2 + thetaEsts[1]*df_boot$z1 + thetaEsts[2]*df_boot$z2
      
      # Update eta(c) in the non-time-varying risk effects generalized linear model
      step0 <- rep(0,nEtagrid)
      step1 <- rep(0,nEtagrid)
      for(j in 1:nEtagrid) {
        co <- etaGrid[j]
        indices2 <- (abs(df_boot$c-co) < hEta)
        df_bootTemp <- df_boot[indices2,]
        cDiff <- df_bootTemp$c - co
        obj <- df_bootTemp$gammaEst0 + df_bootTemp$bisEst + newlpred[indices2] + df_bootTemp$etaEst0 + df_bootTemp$etaEst1*cDiff
        pijks <- invLogit(obj)
        qijks <- 1 - pijks
        epan <- cDiff/hEta
        kernels <- 0.75*(1-epan^2)
        kernels <- (kernels > 0)*kernels/hEta
        # Construct Hessian and Gradient for eta(c) in the non-time-varying risk effects generalized linear model
        a2 <- pijks*qijks + df_bootTemp$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
        a1 <- df_bootTemp$v/2*(pijks*qijks^2 - pijks^2*qijks)
        ttt <- cbind(a2*kernels, a2*kernels*cDiff, a2*kernels*cDiff^2,(df_bootTemp$y - pijks - a1)*kernels,(df_bootTemp$y - pijks - a1)*kernels*cDiff)
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
      Olddf_boot <- df_boot
      OldObj <-df_boot$gammaEst0 + df_boot$bisEst + newlpred + df_boot$etaEst0
      OldPijks <- invLogit(OldObj)
      OldQijks <- 1-OldPijks
      OldlogLikelihood <- sum(OldObj*df_boot$y+log(OldQijks)-df_boot$v*OldPijks*OldQijks/2-(df_boot$bisEst^2+df_boot$v)/(2*sigma2bEst*boot_Nij)-0.5*log(2*pi*sigma2bEst)/boot_Nij)
      for(s in 1:10){
        smin <- (1 / 2) ^ (s - 1)
        for(j in 1:nEtagrid){
          indices <- df_boot$cr==j
          df_boot$etaEst0[indices] <- Olddf_boot$etaEst0[indices] +  smin*step0[j]
          df_boot$etaEst1[indices] <- Olddf_boot$etaEst1[indices] +  smin*step1[j]
        }
        newObj <-df_boot$gammaEst0 + df_boot$bisEst + newlpred + df_boot$etaEst0
        newPijks <- invLogit(newObj)
        newQijks <- 1-newPijks
        newlogLikelihood <- sum(newObj*df_boot$y+log(newQijks)-df_boot$v*newPijks*newQijks/2-(df_boot$bisEst^2+df_boot$v)/(2*sigma2bEst*boot_Nij)-0.5*log(2*pi*sigma2bEst)/boot_Nij)
        if(newlogLikelihood > OldlogLikelihood){
          break
        }
      } 
      # End line search
      etaVCMEst <- aggregate(df_boot$etaEst0,by=list(df_boot$cr),FUN=mean)[,2]
      df_boot$etaEst0 <- df_boot$etaEst0 - 1/2*sum((etaVCMEst[-1]+etaVCMEst[-nEtagrid]))/(nEtagrid-1)  # Normalize eta(c)

      newObj <-df_boot$gammaEst0 + df_boot$bisEst + newlpred + df_boot$etaEst0
      newPijks <- invLogit(newObj)  # Pijk from the current iteration
      Diff <- max(newPijks - oldPijks)  
      numIter <- numIter + 1   
    } 
    # End of finding initial values for theta(t) and beta(t)
    
    # Implement estimation steps as described in Web Appendix A
    # Step 1: set initial values for fitting VCM-MR
    df_boot$bisEst <- rep(0, dim(df_boot)[1])  
    df_boot$etaEst0 <- rep(0, dim(df_boot)[1])
    df_boot$etaEst1 <- rep(0, dim(df_boot)[1]) 
    df_boot$gammaEst0 <- rep(0, dim(df_boot)[1])
    df_boot$gammaEst1 <- rep(0, dim(df_boot)[1]) 
    gammaVCMEst <- matrix(0,numF,ngrid)
    gammaVCMEst1 <- matrix(0,numF,ngrid)
    thetaEstM <- t(matrix(thetaEsts,ntheta,dim(df_boot)[1]))
    thetaEstM1 <- t(matrix(0,ntheta,dim(df_boot)[1]))
    betaEstM <- t(matrix(betaEsts,nbeta,dim(df_boot)[1]))
    betaEstM1 <- t(matrix(0,nbeta,dim(df_boot)[1]))
    sigma2bEst <- sigma2
    
    numIter <- 0
    Diff <- 1
    while(Diff > 0.001 & numIter <= maxNRIter) {
      print(paste("Iter",numIter)) # Print number of iterations
      
      oldlpred <- rowSums(thetaEstM*df_boot[,5:(4+ntheta)]) + rowSums(betaEstM*df_boot[,(5+ntheta):(4+ntheta+nbeta)])
      oldObj <- df_boot$gammaEst0 + df_boot$bisEst + oldlpred + df_boot$etaEst0
      oldPijks <- invLogit(oldObj) # P0ijk from the previous iteration
      
      # Step 2 (E-step): update posterior means bij0 and variances vij0 of subject-specific random effects
      # Gauss-Hermite quadrature integral approximation
      temp <- matrix(NA, dim(df_boot)[1], numQuadPoints)
      for(i in 1:numQuadPoints) {
        x <- ghrule$nodes[i]*sqrt(2*sigma2bEst)  
        temp[,i] <- exp(df_boot$y*(df_boot$gammaEst0 + x + oldlpred + df_boot$etaEst0))/(1+exp(df_boot$gammaEst0 + x + oldlpred + df_boot$etaEst0))
      }
      ttt <- aggregate(temp,by=list(df_boot$sid),FUN=prod)[-1]
      
      # Update posterior means bij0 and variances vij0 
      meansLij <- rep(0,boot_sumNi)
      means2Lij <- rep(0,boot_sumNi)
      const <- rep(0,boot_sumNi)
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
      df_boot$bisEst <- rep(meansLij,boot_numOB)  # bij0
      df_boot$v <- rep(varsLij,boot_numOB)  # vij0
      
      # Step 3 (M-step): update sigma2b by maximizing the approximate expected log-likelihood
      sigma2bEst <- 1/boot_sumNi*SumMean2
      
      # Step 4 (M-step): update gammai(t) by maximizing the approximate expected local log-likelihood 
      # of the data from the ith facility
      for(j in 1:ngrid) {
        to <- gridPoints[j]
        indices <- (abs(df_boot$t-to) < df_boot$h)
        df_bootTemp <- df_boot[indices,] # Data used for local MLE
        lpredTemp <- oldlpred[indices]
        tDiff <- df_bootTemp$t - to
        obj <- gammaVCMEst[df_bootTemp$fid,j] + gammaVCMEst1[df_bootTemp$fid,j]*tDiff + df_bootTemp$bisEst + lpredTemp + df_bootTemp$etaEst0
        pijks <- invLogit(obj)
        qijks <- 1 - pijks
        epan <- tDiff/df_bootTemp$h
        kernels <- 0.75*(1-epan^2)
        kernels <- (kernels > 0)*kernels/df_bootTemp$h # Kernel function
        a2 <- pijks*qijks + df_bootTemp$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
        a1 <- df_bootTemp$v/2*(pijks*qijks^2 - pijks^2*qijks)
        # Construct Hessian and Gradient 
        ttt <- cbind(a2*kernels, a2*kernels*tDiff, a2*kernels*tDiff^2,(df_bootTemp$y - pijks - a1)*kernels,(df_bootTemp$y - pijks - a1)*kernels*tDiff)
        aggTTT <- aggregate(ttt,by=list(df_bootTemp$fid),FUN=sum) # Use only within facility data
        h11 <- aggTTT[,2]
        h12 <- aggTTT[,3]
        h22 <- aggTTT[,4]
        v1 <- aggTTT[,5]
        v2 <- aggTTT[,6]
        detVal <- 1/(h11*h22-h12^2)
        gammaVCMEst[,j] <- gammaVCMEst[,j] + detVal*(h22*v1-h12*v2) # Update gammai(t)
        gammaVCMEst1[,j] <- gammaVCMEst1[,j] + detVal*(-h12*v1+h11*v2) 
        indices <- df_boot$r==j
        df_boot$gammaEst0[indices] <- gammaVCMEst[df_boot$fid[indices],j]
        df_boot$gammaEst1[indices] <- gammaVCMEst1[df_boot$fid[indices],j]
      }  
      
      # Step 5 (M-step): update theta(t) and beta(t) by maximizing the approximate expected local log-likelihood
      Phi_OneStep <- matrix(0,ngrid,2*(ntheta+nbeta))
      for(j in 1:ngrid) {
        to <- gridPoints[j]
        indices1 <- (abs(df_boot$t-to) < hThetaBeta)
        df_bootTemp <- df_boot[indices1,]  # Data used for local MLE
        tDiff <- df_bootTemp$t - to
        lpredTemp <- oldlpred[indices1] + rowSums(thetaEstM1[indices1,]*df_bootTemp[,5:(4+ntheta)]*tDiff) + rowSums(betaEstM1[indices1,]*df_bootTemp[,(5+ntheta):(4+ntheta+nbeta)]*tDiff)
        obj <- df_bootTemp$gammaEst0 + df_bootTemp$bisEst + lpredTemp + df_bootTemp$etaEst0
        pijks <- invLogit(obj)
        qijks <- 1 - pijks
        epan <- tDiff/hThetaBeta
        kernels <- 0.75*(1-epan^2)
        kernels <- (kernels > 0)*kernels/hThetaBeta # Kernel function
        a1 <- df_bootTemp$v/2*(pijks*qijks^2 - pijks^2*qijks)
        a2 <- pijks*qijks + df_bootTemp$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
        # Construct Gradient
        zxMatrix <- data.matrix(cbind(df_bootTemp[,5:(4+ntheta+nbeta)],df_bootTemp[,5:(4+ntheta+nbeta)]*tDiff))
        zxVector <- t(zxMatrix) %*% ((df_bootTemp$y - pijks - a1)*kernels)
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
      OldObj <-df_boot$gammaEst0 + df_boot$bisEst + oldlpred + df_boot$etaEst0
      OldPijks <- invLogit(OldObj)
      OldQijks <- 1-OldPijks
      OldlogLikelihood <- sum(OldObj*df_boot$y+log(OldQijks)-df_boot$v*OldPijks*OldQijks/2-(df_boot$bisEst^2+df_boot$v)/(2*sigma2bEst*boot_Nij)-0.5*log(2*pi*sigma2bEst)/boot_Nij)
      #begin line search
      for(s in 1:10){
        # Step size
        ss <- (1 / 2) ^ (s - 1)
        # Same step size is used for all time points
        for(j in 1:ngrid) {
          indices <- df_boot$r==j
          thetaEstM[indices,] <- t(t(OldthetaM[indices,]) + ss*Phi_OneStep[j,1:ntheta]) # Update theta(t)
          betaEstM[indices,] <- t(t(OldbetaM[indices,]) + ss*Phi_OneStep[j,(ntheta+1):(ntheta+nbeta)]) # Update beta(t)
          thetaEstM1[indices,] <- t(t(OldthetaM1[indices,]) + ss*Phi_OneStep[j,(ntheta+nbeta+1):(2*ntheta+nbeta)])
          betaEstM1[indices,] <- t(t(OldbetaM1[indices,]) + ss*Phi_OneStep[j,(2*ntheta+nbeta+1):(2*ntheta+2*nbeta)])
        } 
        # Calculate global log-likelihood with updated of theta(t) and beta(t)
        newlpred <- rowSums(thetaEstM*df_boot[,5:(4+ntheta)]) + rowSums(betaEstM*df_boot[,(5+ntheta):(4+ntheta+nbeta)])
        newObj <-df_boot$gammaEst0 + df_boot$bisEst + newlpred + df_boot$etaEst0
        newPijks <- invLogit(newObj)
        newQijks <- 1-newPijks
        newlogLikelihood <- sum(newObj*df_boot$y+log(newQijks)-df_boot$v*newPijks*newQijks/2-(df_boot$bisEst^2+df_boot$v)/(2*sigma2bEst*boot_Nij)-0.5*log(2*pi*sigma2bEst)/boot_Nij)
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
        indices2 <- (abs(df_boot$c-co) < hEta)
        df_bootTemp <- df_boot[indices2,]   # Data used for local MLE
        cDiff <- df_bootTemp$c - co
        obj <- df_bootTemp$gammaEst0 + df_bootTemp$bisEst + newlpred[indices2] + df_bootTemp$etaEst0 + df_bootTemp$etaEst1*cDiff
        pijks <- invLogit(obj)
        qijks <- 1 - pijks
        epan <- cDiff/hEta
        kernels <- 0.75*(1-epan^2)
        kernels <- (kernels > 0)*kernels/hEta  # Kernel function
        a2 <- pijks*qijks + df_bootTemp$v/2*(pijks*qijks^3 - 4*pijks^2*qijks^2 + pijks^3*qijks)
        a1 <- df_bootTemp$v/2*(pijks*qijks^2 - pijks^2*qijks)
        # Construct Hessian and Gradient
        ttt <- cbind(a2*kernels, a2*kernels*cDiff, a2*kernels*cDiff^2,(df_bootTemp$y - pijks - a1)*kernels,(df_bootTemp$y - pijks - a1)*kernels*cDiff)
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
      Olddf_boot <- df_boot
      OldObj <-df_boot$gammaEst0 + df_boot$bisEst + newlpred + df_boot$etaEst0
      OldPijks <- invLogit(OldObj)
      OldQijks <- 1-OldPijks
      OldlogLikelihood <- sum(OldObj*df_boot$y+log(OldQijks)-df_boot$v*OldPijks*OldQijks/2-(df_boot$bisEst^2+df_boot$v)/(2*sigma2bEst*boot_Nij)-0.5*log(2*pi*sigma2bEst)/boot_Nij)
      for(s in 1:10){
        ss <- (1 / 2) ^ (s - 1)
        for(j in 1:nEtagrid){
          indices <- df_boot$cr==j
          df_boot$etaEst0[indices] <- Olddf_boot$etaEst0[indices] +  ss*step0[j] # Update eta(c)
          df_boot$etaEst1[indices] <- Olddf_boot$etaEst1[indices] +  ss*step1[j]
        }
        # Calculate global log-likelihood with updated of eta(c)
        newObj <-df_boot$gammaEst0 + df_boot$bisEst + newlpred + df_boot$etaEst0
        newPijks <- invLogit(newObj)
        newQijks <- 1-newPijks
        newlogLikelihood <- sum(newObj*df_boot$y+log(newQijks)-df_boot$v*newPijks*newQijks/2-(df_boot$bisEst^2+df_boot$v)/(2*sigma2bEst*boot_Nij)-0.5*log(2*pi*sigma2bEst)/boot_Nij)
        if(newlogLikelihood > OldlogLikelihood){ 
          break
        }
      }
      # End line search for eta(c)
      etaVCMEst <- aggregate(df_boot$etaEst0,by=list(df_boot$cr),FUN=mean)[,2] 
      df_boot$etaEst0 <- df_boot$etaEst0 - 1/2*sum((etaVCMEst[-1]+etaVCMEst[-nEtagrid]))/(nEtagrid-1)  # Normalize eta(c)
      
      newObj <-df_boot$gammaEst0 + df_boot$bisEst + newlpred + df_boot$etaEst0
      newPijks <- invLogit(newObj) # P0ijk from the current iteration
      Diff <- max(newPijks - oldPijks) 
      print(paste("Diff=",Diff))
      numIter <- numIter + 1   
    }
    
    # Store the estimates of varying coefficient functions and variance of subject-specific random effects   
    # from the bootstrap sample
    thetaVCMEst <- aggregate(thetaEstM,by=list(df_boot$r),FUN=mean)[,2:3]
    betaVCMEst <- aggregate(betaEstM,by=list(df_boot$r),FUN=mean)[,2:3]   
    theta1Est_boot <- rbind(theta1Est_boot,t(thetaVCMEst[1]))
    theta2Est_boot <- rbind(theta2Est_boot,t(thetaVCMEst[2]))
    beta1Est_boot <- rbind(beta1Est_boot,t(betaVCMEst[1]))
    beta2Est_boot <- rbind(beta2Est_boot,t(betaVCMEst[2])) 
    etaEst_boot <- rbind(etaEst_boot,t(aggregate(df_boot$etaEst0,by=list(df_boot$cr),FUN=mean)[2]))
    sigma2bEst_boot <- c(sigma2bEst_boot, sigma2bEst)    
  } 
  
  # Output bootstrap estimation results
  BootEst[[1]] <- theta1Est_boot
  BootEst[[2]] <- theta2Est_boot
  BootEst[[3]] <- beta1Est_boot
  BootEst[[4]] <- beta2Est_boot
  BootEst[[5]] <- etaEst_boot
  BootEst[[6]] <- sigma2bEst_boot
  BootEst[[7]] <- gridPoints
  BootEst[[8]] <- etaGrid
  names(BootEst) <- c("theta1", "theta2", "beta1", "beta2", "eta", "sigma", "grid", "etagrid")
  return(BootEst)
}