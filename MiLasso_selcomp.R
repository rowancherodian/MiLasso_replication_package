#######################################################
# MiLasso_selcomp_sims.R
# Replicate Simulation results for set-up A and B (manuscript-Table 1 and supplement-Table 2-13) in 
# "Moran's I Lasso for models with spatially correlated data"
######################################################

remove(list = ls()) #clear enviroment
#Install packages 
# install.packages("glmnet") #lasso package
library(glmnet)
#install.packages("penalized") #lasso package (CV)
library(penalized)

# Print version information about R, the OS and attached or loaded packages.
sessionInfo()


MiLasso_comp <- function(n1, seedzzz, n2, mu, Rhos, Beta1, psi, Alpha=0, 
                         tol=0.1, A=2, no_con=TRUE){
  
  #arguments
  #  rm(list = ls())
  # n1 <- 100 # sample size
  # seedzzz <- 999 # seed for loop
  # n2 <- 300 #number of replications
  # mu <- 8  #expected number of links
  # Rhos <- 0.9 #coef of spatial lags of y
  # Beta1 <- 1 # coef X_1
  # Alpha <- 0 #coef of constant
  # psi <- 0.9
  # tol <- 0.1
  #  A <- 2
  #  no_con=TRUE
  # # 
  
  
  
  GetMoranStat <- function(MSM, degfree) {
    #MSM    : M %*% S %*% M matrix
    #         M : projection matrix
    #         S : coded symmetric spatial link matrix
    #degfree: degrees of freedom
    
    MSM <- as.matrix(MSM)
    t1 <- sum(diag(MSM))
    t2 <- sum(diag(MSM %*% MSM))
    
    E <- t1 / degfree #equ 8 tg07
    Va <- 2 * (degfree * t2 - t1 * t1)/(degfree * degfree * (degfree + 2)) # equ 9 tg07
    return(list(Mean=E,Var=Va))     
  }
  
  grif_ESF <-function(y, X, W, V, vec, val, ExactEV = FALSE, tol = 0.1){# replacion of TG07 alg
    
    #function arguments
    # y=y_sp #dependent var
    # X=X # exog var
    # W=W #SWM
    # V=V #forced symeteric SWM
    # vec=E # eienvectors to use
    # ExactEV = FALSE #update E(m) and var(m) in search loop (Very computationally demanding)
    # tol = 0.1
    
    #Compute inital statistics 
    m_X <- diag(nrow(X)) - X %*% qr.solve(crossprod(X), t(X)) # projection matix for swm
    MStat <- GetMoranStat(MSM=m_X %*% V %*% m_X, degfree=nrow(X)-ncol(X)) 
    betagam <- solve(crossprod(X),crossprod(X,y)) #b=(X'X)^(-1)(X'Y)
    res <- y - X %*% betagam # residuals = Y - Xb
    IthisTime <- (crossprod(res, V) %*% res) / crossprod(res) #mi = (res'S)res/res'res                     
    zIthisTime <- (IthisTime - MStat$Mean) / sqrt(MStat$Var) #calcuating standarised morans I of residuals                     
    sigma2_hat <- crossprod(res) / (nrow(X) - ncol(X))
    se<-sqrt(diag(c(sigma2_hat) * solve(crossprod(X)))) #calculate se
    
    
    out <- c(	0,
              betagam[1,1],
              se[1],
              0,
              IthisTime,
              zIthisTime)
    Aout <- out
    
    
    
    
    #sel eigen vector number, eig, +1 -1 if +Ve or -ve, 
    sel <- cbind(val,0)
    #sel[,2] <- (val > abs(zerovalue)) - (val < -abs(zerovalue))
    sel[,2] <- (val > 0) - (val < 0)
    
    #Compute the Moran's I of the aspatial model (without any eigenvector)
    #i.e., the sign of autocorrelation if MI is positive, then acsign = 1 or if MI is negative, then acsign = -1 
    acsign <- 1
    if (((crossprod(res, V) %*% res) / crossprod(res)) < 0) acsign <- -1
    
    
    oldZMinMi <- Inf
    for (i in 1:ncol(vec)) { #Outer Loop
      z <- Inf
      idx <- 0
      for (j in 1:ncol(vec)) { #Inner Loop - Find next eigenvector
        if (sel[j,2] == acsign ) { #Use only feasible unselected evecs 
          
          xe <- cbind(X, vec[,j])  #Add test eigenvector
          res1 <- y - xe %*% solve(crossprod(xe), crossprod(xe, y)) #res = y - xe(ex'ex)^(-1)(ex'y)
          mi <- (crossprod(res1, V) %*% res1) / crossprod(res1) #mi = (res'S)res/res'res
          
          if (ExactEV) {
            M <- diag(nrow(X)) - xe %*% solve(crossprod(xe),t(xe))
            MStat <- GetMoranStat(MSM=M %*% V %*% M, degfree=nrow(X)-ncol(xe)) 
          }
          
          if (abs((mi - MStat$Mean) / sqrt(MStat$Var)) < z) { #Identify min z(Moran)
            MinMi = mi
            z <- (MinMi - MStat$Mean) / sqrt(MStat$Var)
            idx =j
          }
        }
      }  #End inner loop
      
      #Update design matrix permanently by selected eigenvector
      X <- cbind(X,vec[,idx])
      
      M <- diag(nrow(X)) - X %*% solve(crossprod(X),t(X))
      
      #Update Expectation and Variance
      MStat <- GetMoranStat(MSM= M %*% V %*% M, degfree=nrow(X) - ncol(X)) 
      ZMinMi <- ((MinMi - MStat$Mean) / sqrt(MStat$Var))
      
      
      betagam <- solve(crossprod(X),crossprod(X,y)) #b=(X'X)^(-1)(X'Y)
      res <- y - X %*% betagam # residuals = Y - Xb
      sigma2_hat <- crossprod(res) / (nrow(X) - ncol(X))
      se<-sqrt(diag(c(sigma2_hat) * solve(crossprod(X)))) #calculate se
      
      
      #Add results of i-th step
      out <- c(i,
               betagam[1,1],
               se[1],
               idx,
               MinMi,
               ZMinMi)
      
      Aout <- rbind(Aout, out)
      
      #To exclude the selected eigenvector in the next loop
      sel[idx,2] <- 0 
      
      
      if (is.na(ZMinMi)) {
        break 
      } else if (abs(ZMinMi) < tol){
        break
      }
      
    } # End Outer Loop
    colnames(Aout) <- c("Step","beta1_cof","beta1_se","SelEvec","MinMi","ZMinMi")
    return(Aout)
  }
  
  MiLasso <- function(y, X, V, E, A=A, no_con=no_con){

    #arguments 
    # y <-y_sp # dependant variable
    # X <- X1 # Matrix of Exonegnous variables
    # V <- V # forced symetric SWM
    # E <- E # Eigenvectors

    if (no_con==TRUE){
      Xcon <- X
    }else{
      Xcon <- cbind(X,rep(1,nrow(X)))
    }
    
    M_X <- diag(1, nrow(X)) - Xcon %*% solve(crossprod(Xcon),t(Xcon)) #projection matrix
    MStat <- GetMoranStat(MSM = M_X%*%V%*%M_X, degfree = nrow(X) - ncol(Xcon)) #Moran E and Var
    res <- y - Xcon %*% solve(crossprod(Xcon),crossprod(Xcon,y)) # residuals = y - xhat*b
    mI <- (crossprod(res, V) %*% res) / crossprod(res) #mi = (res'S)res/res'res
    zI <- (mI - MStat[["Mean"]]) / sqrt(MStat[["Var"]])
    
    theta <- abs(zI)^{-A} # absolute vlaue of the inverse of the morans I 
    pen <- c(rep(0,ncol(X)),rep(1,ncol(E)))
    XE <- cbind(X,E) # design matrix
    fit = glmnet(XE, y, lambda = theta, penalty.factor =pen) # Milasso
    ###
    #find selected evecs for Miplasso
    if((fit$df - ncol(X))==0){ # no eigenvectors selected
      
      if (no_con==TRUE){
        Miplasso <- lm(y ~ X-1) # simple ols model
      }else{
        Miplasso <- lm(y ~ X) # simple ols model
      }
      
      selV=NA
    } else {
      Temp <- as.matrix(cbind(fit[["beta"]],seq(1,nrow(fit[["beta"]])),0))
      Temp[,3] <- (Temp[,1] > 0 | Temp[,1] < 0)  # asign non-zero evecs value 1
      sel_no  <-  Temp[,2]*Temp[,3]
      Ext	= subset(sel_no,sel_no!=0)		#extracted eigenvectors
      Exx		= XE[,Ext]
      selV		= Exx[,-(1:ncol(X))]
      if (no_con==TRUE){
        Miplasso <- lm(y ~ X + selV - 1) # # post lasso
      }else{
        Miplasso <- lm(y ~ X + selV) # post lasso
      }
    }
    
    y_bar=y-E%*%as.vector(fit[["beta"]][-1])
    
    if (is.na(selV)[1]){ #no eigenvectors selected
      if(no_con==TRUE){ #exclude constant
        MiLasso <- lm(y_bar ~ X-1)
      } else{  #include constant
        MiLasso <- lm(y_bar ~ X)
      }
    }else{
      m_E <- diag(nrow(E)) - selV %*%  t(selV) # projection matix for swm
      if(no_con==TRUE){ #exclude constant
        MiLasso <- lm(y_bar ~ m_E%*%X-1)
      } else{  #include constant
        MiLasso <- lm(y_bar ~ m_E%*%X)
      }
    }
    return(list("p.Mi.Lasso"=Miplasso, "n.Mi.Lasso"=fit, "Mi.Lasso"=MiLasso, no_selected=(fit$df - ncol(X)), mI= mI, 
                zI=zI, E_s=as.matrix(selV)))

  }
  

  
  set.seed(seedzzz) # setseed
  p <- length(Rhos) #number of lags of y
  ones <- matrix(1,n1,1) #constant
  
  ################
  #1) SWM 
  ########
    r <- mu/n1 #probs for bernoulli draws of links of correct adj matrix
    
    G<- matrix(0,n1,n1) #initalise adj matrix
    
    for (i in 1:ncol(G)-1){ # upper triangular matrix  of  bernoulli draws of links
      G[i,(1+i):ncol(G)] <-rbinom((n1-i),1,r)
    }
    G <- G+t(G) # make symeteric 
    isSymmetric(G)

  
  W <- G/max(apply(G, 1, sum)) # max row standardised
  
  
  V <- (W + t(W)) / 2 # make SWM symeteric
  EigenV <- eigen(V, symmetric = TRUE)
  E <- as.matrix(EigenV$vectors) # eigenvectors
  vals <- EigenV$values
  
  
  #loops to make a sequence of lags for DGp of y
  L <- list()
  rhoL <-list(diag(n1))
  for(i in 1:(p+1)){
    Wlist <- lapply(seq(i), function(x){return(W)}) 
    L[[i]] <- Reduce("%*%", Wlist) #W^p
    rhoL[[i+1]] <-Rhos[i]*L[[i]] #rho_i*W^p for S
  }
  rhoL <- rhoL[-length(rhoL)]
  
  #form S
  S <- Reduce("-", rhoL) 
  
  
  results <- matrix(0, n2, 17) #colums OLS(x2), Niave-Mi-Lasso(x3), Post-Mi-Lasso(x2), Mi-Lasso(x2) ,Niave-CV-Lasso(x3),Post-CV-Lasso(x2), FstepZ(x3)
  outCI95_ho <- matrix(0, n2, 5) #coverage results
  outCI99_ho <- matrix(0, n2, 5) #coverage results
  
  for (i in 1:n2) { #rep loop
   # i=1
    
    X1 <- as.matrix(rnorm(n1,0,1)) # X1
    epsilon <- rnorm(n1,0,1) # generate error
    y_sp <- solve(S, (ones%*%Alpha+ X1%*%Beta1 + W%*%X1%*%psi  + epsilon))  # spatial endog varible DGP
    
    
    #ols
    
    if(no_con==TRUE){ #exclude constant
      simpOLS <- lm(y_sp ~ X1-1) #ignore all spatial process
    } else{  #include constant
      simpOLS <- lm(y_sp ~ X1) #ignore all spatial process
    }
    
    
    #####################  
    # Mi-Lasso
    
    res_mi <- MiLasso(y=y_sp, X=X1, V=V, E=E, A=A, no_con=no_con)
    
    
    #######################
    # # CV lasso seya
    step1 <-  tryCatch(optL1(response=y_sp, penalized=E, unpenalized=X1, model = "linear", fold=20,standardize=T, trace=FALSE),
                       error=function(e) NA)#Kfold CV
    if (length(step1)<5) step1 <- tryCatch(optL1(response=y_sp, penalized=E, unpenalized=X1, model = "linear", fold=20,standardize=T, trace=FALSE),
                                           error=function(e) NA)#Kfold CV
    if (length(step1)<5){
      CVlasso <- NA
    }else{
      CVlasso <- tryCatch(penalized(response=y_sp, penalized=E, unpenalized=X1, lambda1=step1$lambda, model = "linear",standardize=T, trace=FALSE),
                          error=function(e) NA)
    }
    
    
    # CV lasso
    ###
    # CV POST lasso seya
    #extract selected evects
    if (suppressWarnings(is.na(CVlasso))) {
      CVplasso <- simpOLS
    }else{
      if(sum(CVlasso@penalized!=0)==0){ # no eigenvectors selected
        CVplasso <- simpOLS
      } else {
        Temp <- as.matrix(cbind(as.matrix(CVlasso@penalized),seq(1,nrow(as.matrix(CVlasso@penalized))),0))
        Temp[,3] <- (Temp[,1] > 0 | Temp[,1] < 0) # asign non-zero evecs value 1
        
        sel_no  <-  Temp[,2]*Temp[,3]
        Ext	= subset(sel_no,sel_no!=0)		#extracted eigenvectors
        selV		= E[,Ext]
        
        if(no_con==TRUE){ #exclude constant
          CVplasso <- lm(y_sp ~ X1 + selV-1) # post lasso 
        } else{  #induce constant
          CVplasso <- lm(y_sp ~ X1 + selV) # post lasso 
        }
        
      }
    }
    
    
    
    #########################################
    #FstepZ (TG07 alg)
    grif <- grif_ESF(y=y_sp, X=X1, W=W, V=V, vec=E, val=vals, ExactEV = FALSE, tol = tol)
    grfE <- E[,grif[2:nrow(grif),4]]
    if(no_con==TRUE){ #exclude constant
      pgrif <- lm(y_sp ~ X1 + grfE-1) 
    } else{  #induce constant
      pgrif <- lm(y_sp ~ X1 + grfE)
    }
    
    if(no_con==TRUE){ #exclude constant
      g <- 1
    } else{  #induce constant
      g <- 2
    }
    
    #simple OLS (ignore spatial)
    results[i,1:2] <- coef(summary(simpOLS))[g,1:2] # save Beta and SE
    
    
    
    #MiLasso results  
    results[i,3] <- res_mi[["n.Mi.Lasso"]][["beta"]][1,] # lasso beta estimate
    results[i,4]<- res_mi[["no_selected"]] # Number of vec's selected Milasso
    results[i,5] <- res_mi[["n.Mi.Lasso"]][["lambda"]] # MiLasso tuning estiamte
    results[i,6:7] <- coef(summary(res_mi[["p.Mi.Lasso"]]))[g,1:2] #POST Milasso beta estiamte 
    results[i,8:9] <- coef(summary(res_mi[["Mi.Lasso"]]))[g,1:2] # Milasso beta estiamte (partial regression)
    
    #CV-Lasso results 
    if (suppressWarnings(is.na(CVlasso))) {
      results[i,10] <- NA # save lasso beta estiamte
      results[i,11] <- NA # save lasso theta estiamte
      results[i,12] <- NA #number of vec's selected lasso
    }else{
      results[i,10] <- CVlasso@unpenalized[1] # save lasso beta estiamte
      results[i,11] <- step1$lambda # save lasso theta estiamte
      results[i,12] <- sum(CVlasso@penalized!=0) #number of vec's selected lasso
    }
    results[i,13:14] <- coef(summary(CVplasso))[g,1:2] #POST CVlasso beta estiamte
    #FstepZ results    
    results[i,15:16] <- coef(summary(pgrif))[g,1:2] # save FstepZ beta estiamte
    results[i,17] <- grif[nrow(grif),1] #number of vec's FstepZ


    
    
    
    model_list <- list(simpOLS, res_mi[["p.Mi.Lasso"]],res_mi[["Mi.Lasso"]], CVplasso, pgrif)
    for(k in 1:length(model_list)){
      outCI95_ho[i,k] <- Beta1 <= confint(model_list[[k]])[g,1] |
        Beta1 >= confint(model_list[[k]])[g,2] 
      outCI99_ho[i,k] <- Beta1 <= confint(model_list[[k]], level=0.99)[g,1] |
        Beta1 >= confint(model_list[[k]], level=0.99)[g,2] 
    }
    
  }   
    
  
  res <- data.frame(matrix(NA, 7 , 18+p ))
  
  pars <- c(n2, n1, p, mu, Beta1 , Alpha, psi, Rhos)
  
  #ols
  res[1,] <- c(pars, NA, no_con, "OLS-simple",
               mean(results[,1]) - Beta1, mean((results[,1]-Beta1)^2), mean(results[,2]), #beta1 - bias/mse/ase
               sd(results[,1]), #beta1 sd (sims)
               (1-mean(outCI95_ho[,1])), (1-mean(outCI99_ho[,1])), NA, NA)
  #mi-lasso's
  res[2,] <- c(pars, A, NA, "Naive-Mi-Lasso",
               mean(results[,3]) - Beta1, mean((results[,3]-Beta1)^2), NA, #beta1 - bias/mse
               sd(results[,3]), #beta1 sd (sims)
               NA, NA, mean(results[,4]) , mean(results[,5]))
  res[3,] <- c(pars, A, no_con, "Post-Mi-Lasso",
               mean(results[,6]) - Beta1, mean((results[,6]-Beta1)^2), mean(results[,7]), #beta1 - bias/mse/ase
               sd(results[,6]), #beta1 sd (sims)
               (1-mean(outCI95_ho[,2])), (1-mean(outCI99_ho[,2])), mean(results[,4]) , NA)
  res[4,] <- c(pars, A, no_con, "Mi-Lasso",
               mean(results[,8]) - Beta1, mean((results[,8]-Beta1)^2), mean(results[,9]), #beta1 - bias/mse/ase
               sd(results[,8]), #beta1 sd (sims)
               (1-mean(outCI95_ho[,3])), (1-mean(outCI99_ho[,3])), mean(results[,4]) , NA)
  #CV-Lasso
  res[5,] <- c(pars, NA, NA, "Naive-CV-Lasso",
               mean(results[,10], na.rm=TRUE) - Beta1, mean((results[,10]-Beta1)^2, na.rm=TRUE), NA, #beta1 - bias/mse
               sd(results[,10], na.rm=TRUE), #beta1 sd (sims)
               NA, NA, mean(results[,12], na.rm=TRUE) , mean(results[,11], na.rm=TRUE))
  res[6,] <- c(pars, NA, no_con, "Post-CV-Lasso",
               mean(results[,13]) - Beta1, mean((results[,13]-Beta1)^2), mean(results[,14]), #beta1 - bias/mse/ase
               sd(results[,13]), #beta1 sd (sims)
               (1-mean(outCI95_ho[,4])), (1-mean(outCI99_ho[,4])), mean(results[,12], na.rm=TRUE) , NA)
  #FstepZ
  res[7,] <- c(pars, NA, no_con, "FstepZ",
               mean(results[,15]) - Beta1, mean((results[,15]-Beta1)^2), mean(results[,16]), #beta1 - bias/mse/ase
               sd(results[,15]), #beta1 sd (sims)
               (1-mean(outCI95_ho[,5])), (1-mean(outCI99_ho[,5])), mean(results[,17]), NA)


  
  
  if (length(Rhos)==1) rho_names= c("rho1")
  if (length(Rhos)==2) rho_names= c("rho1", "rho2")
  if (length(Rhos)==3) rho_names= c("rho1", "rho2", "rho3")
  if (length(Rhos)==4) rho_names= c("rho1", "rho2", "rho3", "rho4")
  
  colnames(res) <- c("reps", "sample_size", "p", "mu", "beta1", "alpha", "psi", rho_names, 
                     "A", "no_con",  "Estimator",
                     "beta1_Bias", "beta1_MSE", "beta1_AASE","beta1_SD", "CI95", "CI99", "No_Evecs", "Tuning")
  
  
  res$beta1_Bias<-as.numeric(res$beta1_Bias)
  res$beta1_SD<-as.numeric(res$beta1_SD)
  res$beta1_MSE<-as.numeric(res$beta1_MSE)
  res$beta1_AASE<-as.numeric(res$beta1_AASE)
  res$Tuning<-as.numeric(res$Tuning)
  res$No_Evecs<-as.numeric(res$No_Evecs)
  res$CI95 <- as.numeric(res$CI95)
  res$CI99 <- as.numeric(res$CI99)
  
  return(res)
}

n2 <- 1000
seedzz <- 999
Beta_val <- 1
rho_vals <- c(0.3,0.6,0.9)
n1 <- c(100,250,500)
psi_val <- 0.9


############################
#set up A
rho_vals <- c(0.3,0.6,0.9)
#binary SWM
res <- list()
mu <- c(4,8,12) 
for(i in 1: length(n1)){
  res[[i]] <- MiLasso_comp(n1=n1[i], seedzzz=seedzz, n2=n2, mu=mu[1], Rhos=rho_vals[1],
                           Beta1=Beta_val, psi=psi_val )
  for(p in 1:length(mu)){
     for(j in 1:length(rho_vals)){
      res[[i]] <- rbind(res[[i]],MiLasso_comp(n1=n1[i], seedzzz=seedzz, n2=n2, mu=mu[p], Rhos=rho_vals[j],
                                              Beta1=Beta_val, psi=psi_val ))
      
    }
    
  }
  res[[i]] <- res[[i]][-(1:7),]
}

save(res, file="SDM_setA.Rdata")

resA <- res

# produce tables for setup A in supplement (supplement tables 2-10)
for(i in 1:length(res)) res[[i]] <-split(res[[i]], res[[i]]$mu)

for(i in 1:length(res)){
  for(j in 1:length(res[[i]])){
    tab1 <- cbind("rho"=res[[i]][[j]]$rho1, 
                  "Estimator"=res[[i]][[j]]$Estimator,
                  "Bias"=format(round(res[[i]][[j]]$beta1_Bias,3), nsmall = 3),
                  "MSE"=format(round(res[[i]][[j]]$beta1_MSE,3), nsmall = 3), 
                  "SD"=format(round(res[[i]][[j]]$beta1_SD,3), nsmall = 3),
                  "AASE"=format(round(res[[i]][[j]]$beta1_AASE,3), nsmall = 3),
                  "CI95"=format(round(res[[i]][[j]]$CI95,3), nsmall = 3), 
                  "CI99"=format(round(res[[i]][[j]]$CI99,3), nsmall = 3),
                  "vecs#"=round(res[[i]][[j]]$No_Evecs,0))
    print(xtable::xtable(tab1, caption=paste0("Bias, MSE, SD, AASE, and 95/99 pre cent coverage and selected eigenvectors for mu=",
                                              res[[i]][[j]]$mu[1], " and n=",res[[i]][[j]]$sample_size[1], " (setup A)")),include.rownames=FALSE)
  }
}

res <- resA

# produce table for setup A in main paper (table 1)
for(i in 1:length(res)){
  res[[i]]<-res[[i]][!(res[[i]]$Estimator=="Naive-Mi-Lasso"|
                         res[[i]]$Estimator=="Naive-CV-Lasso"),] # drop raw lasso estimates
  res[[i]] <-split(res[[i]], res[[i]]$mu)
} 


tab1 <- list()
for(i in 1:length(res)){
  tab1[[i]] <- list()
  for(j in 1:length(res[[i]])){
    tab1[[i]][[j]] <- cbind("n"=res[[i]][[j]]$sample_size,
                            "rho"=res[[i]][[j]]$rho1, 
                            "Estimator"=res[[i]][[j]]$Estimator,
                            "Bias"=format(round(res[[i]][[j]]$beta1_Bias,3), nsmall = 3),
                            "MSE"=format(round(res[[i]][[j]]$beta1_MSE,3), nsmall = 3), 
                            "SD"=format(round(res[[i]][[j]]$beta1_SD,3), nsmall = 3),
                            "AASE"=format(round(res[[i]][[j]]$beta1_AASE,3), nsmall = 3),
                            "CI95"=format(round(res[[i]][[j]]$CI95,3), nsmall = 3), 
                            "CI99"=format(round(res[[i]][[j]]$CI99,3), nsmall = 3),
                            "vecs#"=round(res[[i]][[j]]$No_Evecs,0))
  }
}

mu4 <- rbind(tab1[[1]][[2]],tab1[[2]][[2]],tab1[[3]][[2]])
print(xtable::xtable(mu4),include.rownames=FALSE)



#########################
#set up B
rho_vals <- c(0.6,0.4,0.5)
#binary SWM
mu <- c(4,8,12) 
res <- list()
for(i in 1: length(n1)){
  res[[i]] <- list()
  for(p in 1:length(mu)){
    res[[i]][[p]] <- MiLasso_comp(n1=n1[i], seedzzz=seedzz, n2=n2, mu=mu[p], Rhos=rho_vals,
                             Beta1=Beta_val, psi=psi_val)
    
  }
}

save(res, file="SDM_setB.Rdata")

# produce tables for setup B in supplement (supplement tables 11-13)
for (i in 1:length(res)) {
  res[[i]] <-do.call("rbind", res[[i]])
  
  tab1 <- cbind("mu"=res[[i]]$mu,
                "Estimator"=res[[i]]$Estimator,
                "Bias"=format(round(as.numeric(res[[i]]$beta1_Bias),3), nsmall = 3),
                "MSE"=format(round(as.numeric(res[[i]]$beta1_MSE),3), nsmall = 3), 
                "SD"=format(round(as.numeric(res[[i]]$beta1_SD),3), nsmall = 3),
                "AASE"=format(round(as.numeric(res[[i]]$beta1_AASE),3), nsmall = 3),
                "CI95"=format(round(as.numeric(res[[i]]$CI95),3), nsmall = 3), 
                "CI99"=format(round(as.numeric(res[[i]]$CI99),3), nsmall = 3),
                "vecs#"=round(as.numeric(res[[i]]$No_Evecs),0) )
  print(xtable::xtable(tab1, caption=paste0("Bias, MSE, SD, AASE, and 95/99 pre cent coverage and selected eigenvectors for distance SWM (setup B) n=",
                                            res[[i]]$sample_size[1])),
        include.rownames=FALSE)
}

