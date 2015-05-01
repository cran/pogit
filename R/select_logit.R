#### Bayesian variable selection for (binomial) logit regression model 
#### with random intercept using Dirac spike and Student-t/ or normal slab
#### (invisible function)
#### last change: 2015/03/11

## -------------------------
# ... (list) objects needed:
#	- y ... number of binomial successes
#	- N ... number of binomial trials
# - W ... design matrix (including intercept)
# - H
# - compmix.bin ("logit": part I of data_augmentation for binomial logit models, "pogit": NULL)
# model [incl. deltafixf, gammafix, ri, clusterID, d, family, Zl]
# prior [incl. slab, psi.nu, V, M0, m0, aj0, w, pi, psi.Q, a0, invM0]
# mcmc  [incl. M, burnin, startsel, thin, verbose, msave, nmc]
# par=  [incl. parameters of binomial logit selection model]
# imc=	[= number of MCMC iteration]

# ... parameters needed (list object => par.logit):
# alpha, atilde, theta
# psi
# delta, pdelta, gamma, pgamma
# omega
# pi

select_logit <- function(y, N, W, H=NULL, compmix.bin=NULL, model, prior, mcmc, param, imc){   
      
    # linear predictior in binomial logit model
    muL <- W%*%param$alpha
    
    if (model$ri==1){
      ranEff <- param$atilde[model$Zl]
      linp <- muL + ranEff*param$theta 
    } else linp <- muL
          
          
    #### Step A --- data augmentation for the binomial logit model (binomial dRUM)
    
    ## DATA AUGMENTATION (part I) for binomial logit models taken from "binomlogit"
    if (model$family=="pogit"){
      compmix.bin <- dataug_binom_dRUM1(y, N)
      #compmix.bin <- do.call(dataug_binom_dRUM1, list(y, N))
    }
    
    ## DATA AUGMENTATION (part II) for binomial logit models taken from "binomlogit"
    augBinom <- dataug_binom_dRUM2(y, N, linp, compmix.bin)
    z        <- augBinom$ystar
    invSig   <- augBinom$invSig
    yS       <- z*invSig  # =sqrt(Sigma^-1)*z
     
    if (model$ri==0){
      Wall <- W*kronecker(matrix(1, 1, model$d + 1), invSig) 
    } else {
      Wall <- cbind(W,ranEff)*kronecker(matrix(1, 1, model$d + model$ri + 1), invSig)  
    }
    
    # inverse prior variance of regression effects (updated)
    invA0 <- diag(c(prior$invM0, 1/param$psi), nrow = model$d + model$ri + 1)    
       
    
    #### Step B --- starts Bayesian variable selection
    if (imc > mcmc$startsel && sum(sum(!model$deltafix) + sum(!model$gammafix))>0){
        
        ## (1) update mixture weights
        incfix <- sum(param$delta==1)  #==TODO: sum(param$delta)???
        omega  <- rbeta(1, prior$w['wa0'] + incfix, prior$w['wb0'] + model$d - incfix)
        
        if (model$ri==1){
          incran <- sum(param$gamma==1) #==TODO: sum(param$gamma)???
          pi     <- rbeta(1, prior$pi['pa0'] + incran, prior$pi['pb0'] + model$ri - incran)
        } else {
        	pi <- NULL
        } 

        ## (2) sample the indicators, the regression coefficients and the scale parameters
        ## --- (i) sample the indicators delta_{alpha,j}, gamma_{alpha} for the slab component
        
        ##==TODO generalize draw_indices for logit, pogit, poisson!!!
        indic <- draw_indicators(yS, Wall, param$delta, param$gamma, omega, pi, model, prior, invA0)
        delta  <- indic$deltanew 
        pdelta <- indic$pdeltanew
        gamma  <- indic$gammanew
        pgamma <- indic$pgammanew
    } else {
    	delta  <- param$delta
    	pdelta <- param$pdelta
    	gamma  <- param$gamma
    	pgamma <- param$pgamma
    	omega  <- param$omega
    	pi     <- param$pi
    }
       
    ## --- (ii_A) sample the (selected) regression effects 
    if (model$ri==0){
      index <- c(1, which(delta==1)+1) 
    } else {
      index <- c(1, which(c(delta, gamma)==1)+1)  
    }
    
    Zsel      <- Wall[, index, drop=FALSE]  # Z*=[1, W^delta,atilde]*sqrt(Sigma^-1)
    dsel      <- length(index)
    invA0_sel <- invA0[index, index, drop=FALSE]
    a0_sel    <- prior$a0[index,,drop=FALSE]
    
    AL    <- solve(invA0_sel + t(Zsel)%*%Zsel)  # A = (A0^-1 + (Z*)'Sigma^-1 Z*)
    aL    <- AL%*%(invA0_sel%*%a0_sel + t(Zsel)%*%yS) # a = A(A0^-1*a0 + (Z*)'Sigma^-1*y
    zetaL <- t(chol(AL))%*%matrix(rnorm(dsel), dsel, 1) + aL 
    
    v1        <- matrix(0, 1, model$d + model$ri + 1) 
    v1[index] <- t(zetaL)
    mu_alpha  <- v1[1]
  
    if (model$d > 0){
      alpha  <- v1[2:(model$d+1)]
    } else {
      alpha  <- matrix(0, 1, model$d)   
    }
    muL <- W%*%c(mu_alpha, alpha)
#    muL <- mu_alpha + W[,-1,drop=FALSE]%*%alpha #2015/04/14
              
    
    ## --- (ii_B) sample the random intercepts
    if (model$ri==1){
      theta <- v1[((model$d+1)+1):(model$d+model$ri+1)]   
       
      yh <- (z-muL)*invSig 
      Xh <- (H*theta)*kronecker(matrix(1, 1, max(model$Zl)), invSig) 
       
      BL <- solve(prior$invB0 + t(Xh)%*%Xh)   # B = (theta^2*H'*Sigma^-1 + B0)
      bL <- BL%*%t(Xh)%*%yh
      atilde <- t(chol(BL))%*%matrix(rnorm(max(model$Zl)), max(model$Zl),1) + bL
    
      # perform a random sign-switch (sign-switching step)
      sswitch <- sign(runif(1) - 0.5)
      theta   <- theta*sswitch
      atilde  <- atilde*sswitch 
    } else {  
        atilde <- NULL  
        theta  <- NULL #==TODO: or 0 
    }
    
    ## --- (iii) sample the variance parameter of Student-t/ or normal slab
    ind  <- c(delta, gamma)
    if (model$ri==1){
      effAlpha <- c(t(alpha), t(theta))
    } else {
      effAlpha <- t(alpha)
    }
    
    psi <- draw_psi(effAlpha, ind, prior)
    
    # returns updated par-list with parameters used for subsequent step
    return(list(alpha = c(mu_alpha, alpha), delta = delta, pdelta = pdelta, 
                omega = omega, psi = psi, atilde = atilde, theta = theta, 
                gamma = gamma, pgamma = pgamma, pi = pi))
}



