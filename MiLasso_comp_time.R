############################################
# Milasso_comp_time.rdata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCC
# compares computational time of Mi-lasso, FstepZ and CV-lasso
# DGP SAR(1) with one co-variate 
# replicates (supplement) Table 14 in
# "Moran's I Lasso for models with spatially correlated data"
#################################################


remove(list = ls()) #clear enviroment

#Install packages 
# install.packages("glmnet") #lasso package
library(glmnet)
#install.packages("penalized") #lasso package
library(penalized)

# Print version information about R, the OS and attached or loaded packages.
sessionInfo()


set.seed(999)


#functions
SWMmu <- function(mu,n1){
  
  r <- mu/n1 #probs for bernoulli draws of links of correct adj matrix
  
  G<- matrix(0,n1,n1) #initalise adj matrix
  
  for (i in 1:ncol(G)-1){ # upper triangular matrix  of  bernoulli draws of links
    G[i,(1+i):ncol(G)] <-rbinom((n1-i),1,r)
  }
  G <- G+t(G) # make symeteric 
  isSymmetric(G)
  return(G)
}
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
MiLasso <- function(y, X, V, E, A=2, no_con=TRUE){
  
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
CVlasso <- function(y,  X , E, no_con=TRUE){
  step1 <-  optL1(response=y, penalized=E, unpenalized=X, model = "linear", fold=20, standardize=T, trace=FALSE) #Kfold CV
  fit <- penalized(response=y, penalized=E, unpenalized=X, lambda1=step1$lambda, model = "linear",standardize=T, trace=FALSE)
  Temp <- as.matrix(cbind(as.matrix(fit@penalized),seq(1,nrow(as.matrix(fit@penalized))),0))
  Temp[,3] <- (Temp[,1] > 0 | Temp[,1] < 0) # asign non-zero evecs value 1
  
  sel_no  <-  Temp[,2]*Temp[,3]
  Ext	= subset(sel_no,sel_no!=0)		#extracted eigenvectors
  selV		= E[,Ext]
  
  if(no_con==TRUE){ #exclude constant
    CVplasso <- lm(y ~ X + selV-1) # post lasso 
  } else{  #induce constant
    CVplasso <- lm(y ~ X + selV) # post lasso 
  }
  
  
  return(list(CVlasso=fit,CVplasso=CVplasso))
}
grif_ESF <-function(y, X, W, vec, val, ExactEV = FALSE, tol = 0.1){# replacion of TG07 alg
  
  #function arguments
  # y=y_sp #dependent var
  # X=X # exog var
  # W=W #SWM
  # V=V #forced symeteric SWM
  # vec=E # eienvectors to use
  # ExactEV = FALSE #update E(m) and var(m) in search loop (Very computationally demanding)
  # tol = 0.1
  V <- (W + t(W)) / 2 # make SWm symeteric
  
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

n1 <- c(250, 500, 1000, 10000)
n1names <- c("n=250", "n=500", "n=1000", "n=10000")

#swm
mu_vals <- 8
swm <- list()
  for (i in 1:length(n1)){
    swm[[i]] <- SWMmu(mu=mu_vals,n1=n1[i])
  }


#spectoral decomposition
eigenV	<- list()
E <- list()
vals <- list()
  for (i in 1:length(n1)){
   eigenV[[i]]	<- eigen((swm[[i]]), symmetric=TRUE)
   E[[i]] <- as.matrix(eigenV[[i]][["vectors"]])
   vals[[i]] <- eigenV[[i]][["values"]]
  }

#generate X variables
Xlist <-list()
for(i in 1:length(n1)){
  set.seed(999) # set seed
  Xlist[[i]] <- as.matrix(rnorm(n=n1[i],mean=0,sd=1))# generate 100 obs n(12,2^2)
}

#generate error
epsilon <-list()
for(i in 1:length(n1)){
  set.seed(1) # set seed
  epsilon[[i]] <- rnorm(n1[i],0,1) # generate error
}

#Beta true values 
Beta <-1
rho <- 0.3 #spatial lag y cofs
psi <- 0.9


#gen dep var y
y_sp <-list()
for(i in 1:length(n1)){
    y_sp[[i]] <- solve((diag(n1[i]) - rho*swm[[i]]), (Xlist[[i]]%*%Beta + swm[[i]]%*%Xlist[[i]]%*%psi
                                                                     + epsilon[[i]])) #spatial lag
  }


timecomp <- matrix(NA, nrow=length(n1),ncol=3)
rownames(timecomp)<-n1
colnames(timecomp)<-c("Mi-Lasso",  "CV-lasso", "FstepZ")


for(i in 1:length(n1)){
    timecomp[i,1] <- system.time(MiLasso(y=y_sp[[i]], X=Xlist[[i]], V=swm[[i]], E=E[[i]]))[3]
    timecomp[i,2] <- system.time(CVlasso(y=y_sp[[i]], X=Xlist[[i]], E=E[[i]]))[3]
    if(n1[i]>=10000){
      timecomp[i,3] <- NA # not computed as to computationally demanding
    }else{
      timecomp[i,3] <- system.time(grif_ESF(y=y_sp[[i]], X=Xlist[[i]], W=swm[[i]], 
                                                 vec = E[[i]], val=vals[[i]]))[3]
    }
}

# print results in (supplement) Table 14
timecomp # print times
timecomp/timecomp[,1] # prints times relative to Mi-Lasso



