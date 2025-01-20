###########################################
# MiLasso_bhousing_app.R
# Comparison of forward stepwise (FstepZ), CV, and Mi-Lasso using Boston housing dataset
# Replication of the application (manuscript-Table 3  and supplement-Table 15) in 
#"Moran's I Lasso for models with spatially correlated data"
#############################################


remove(list = ls()) #clear enviroment
#######################################################################################################
###Package installation
#install.packages("spdep") # package for data 
library(spdep)
# install.packages("penalized") #lasso package CV lasso
library(penalized)
# install.packages("glmnet") #lasso package
library(glmnet)
#install.packages("lmtest") 
library(lmtest)
#install.packages("sandwich") 
library(sandwich)

# Print version information about R, the OS and attached or loaded packages.
sessionInfo()

set.seed(1)	#Results depend on the seed because we use CV (may be a bit different from the result in the reference)

#######################################################################################################
#####Reading data
data(boston)
data = boston.c[,c("MEDV","LON","LAT","CRIM","ZN","INDUS","CHAS","NOX","RM","AGE","DIS","RAD","TAX","PTRATIO","B","LSTAT")]
data["MEDV"]	= log(data["MEDV"])
data["B2"]	= data["B"]
data["B"]	=  (0.63 - (data["B2"]/1000)^{1/2})*100
data["NOX"]		= (data["NOX"])^2
n = 506						#Sample size
k = 13						#Number of explanatory variables
attach(data)
Y = as.numeric(MEDV)				#Dependent variables vector [506 x 1]
X1 = as.matrix(data[,4:16])				
X = matrix(as.numeric(X1),n,13)		#Explanatory variables matrix [506 x 13] 
colnames(X) <- colnames(X1)

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

grif_ESF <-function(y, X, V, vec, val, ExactEV = FALSE, tol = 0.1){# replacion of TG07 alg
  
  #function arguments
  # y=Y #dependent var
  # X=X # exog var
  # V=W_mm #forced symeteric SWM
  # vec = egvec  # eienvectors to use
  # val = egval
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
  colnames(Aout) <- c("Step","SelEvec","MinMi","ZMinMi")
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
  
  y_bar=y-E%*%as.vector(fit[["beta"]][-(1:ncol(X))])
  
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


#######################################################################################################
##spatial weight matrix construction and spectoral decompostion

W 	= as.matrix(nb2mat(boston.soi,style="B")) # generate saptial weight matrix
W_mm = W/max(rowSums(W)) # max sum normalise

egW 	= eigen(W_mm)	#Extracting eigenvectors from modified SWM
egvec = egW$vectors
egval = egW$values
colnames(egvec)=rep(1:n)

######################################################################################
#OLS
ols <- lm(Y~X)


#########################################################################################################
##FstepZ
FstepZ_start <- Sys.time() # start timer
FstepZ <- grif_ESF(y=Y, X=X, V=W_mm, vec=egvec, val=egval, ExactEV = FALSE, tol = 0.1)
Fstepz_sel_no <- sort(FstepZ[-1,2]) #extracted eigenvectors number
Fstepz_sel		= egvec[,Fstepz_sel_no] 
FstepZ_reg <- lm(Y~X+Fstepz_sel)
FstepZ_end <- Sys.time() # end timer

#######################################################################################################
##CV-LASSO

cvlasso_start <- Sys.time() # start timer
LASSOecv_x 	= optL1(Y,penalized=egvec,unpenalized=X,fold=20, model="linear",standardize=T)	#Kfold CV
LASSOe_x 	= penalized(Y,penalized=egvec,unpenalized=X,model="linear",lambda1=LASSOecv_x$lambda,standardize=T)
coef_x		= as.data.frame(coefficients(LASSOe_x,"penalized"))
rownames(coef_x) = colnames(egvec)
#Eigenvectors which is not zero
xx			= as.numeric(rownames(coef_x))*ifelse(coef_x!=0,1,0)
CVlasso_sel_no	= subset(xx,xx!=0)		#extracted eigenvectors number
CVlasso_sel		= egvec[,CVlasso_sel_no]	 
CVplasso	= lm(Y~X+CVlasso_sel)	 # post lasso 
cvlasso_end <- Sys.time() # end timer

#######################################################################################################
##Mi-Lasso 

Milasso_start <- Sys.time() # start timer
milasso <-  MiLasso(y=Y, X=X, V=W_mm, E=egvec, A=2, no_con=FALSE)
Milasso_end <- Sys.time() # end timer



#computational times
FstepZ_end - FstepZ_start 
cvlasso_end - cvlasso_start
Milasso_end - Milasso_start


#print the regression results for the four models
models <- list("simple.OLS"=ols, "FstepZ"=FstepZ_reg,"post.CV.Lasso"= CVplasso,
               "Mi.Lasso"=milasso[["Mi.Lasso"]])

for(i in 1:length(models)){ # prints regression results
  print("********************************")
  print(paste0("Result in colum ", i, " of Table 3 in main paper and 15 in supplement"))
  print(names(models)[i])
  print(format(round(coeftest(models[[i]], vcov = vcovHC(models[[i]], type = "HC1")),3), nsmall = 3))
  # print(coeftest(models[[i]], vcov = vcovHC(models[[i]], type = "HC1"))) # uncomment for unrounded results
}


