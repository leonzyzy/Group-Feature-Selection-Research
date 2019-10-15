#library support
library(mc2d)
library(MASS)
library(mvtnorm)
library(gglasso)
library(grpreg)
library(ClustOfVar)
library(Matrix)
library(magic)
library(SGL)
library(mixmeta)
library(rlist)

#### ===================== define tuning para
beta = matrix(rep(c(1.2, -9, 0, 5.2, -4.1, 0, 1.2, 0, -0.95, 0),
                  c(1, 1, 6, 1, 1, 12, 1, 1, 1, 5)), nrow = 30)
block.group = c(3,4,4,3,4,3,4,3,2) 
block.num = 1:30
nrow = 30
ncol = 30
r = c(0:50)

#define true coefficient 
true.coef.matrix = matrix(0, rep(30), nrow = 1, ncol = 30)
true.coef.matrix[, which(beta != 0)] = 1
colnames(true.coef.matrix) = c(sprintf("X%d", seq(1,30)))


#### ===================== sim data function
sim_data = function(r, vec, group, nrow, ncol, beta, seed){
  ##define cov matrix
  cov.matrix = matrix(nrow=nrow, ncol=ncol)
  
  #fit values using r^|i-j|
  for(i in 1:30){
    for(j in 1:30){
      cov.matrix[i,j] = r^(abs(i-j))
    }
  }
  #create a list
  block.matrix.list = list()
  div.list = list()
  div.list[[1]] = vec[1:group[1]]
  
  for(i in 2:length(group)){
    left.indx = div.list[[i-1]][length(div.list[[i-1]])] + 1
    right.indx = sum(group[1:i])
    
    div.list[[i]] = vec[(left.indx):(right.indx)]
  }
  
  for(i in 1:length(group)){
    block.matrix.list[[i]] = cov.matrix[div.list[[i]], div.list[[i]]]
  }
  
  #create a block matrix
  block.matrix = bdiagMat(block.matrix.list)
  
  set.seed(seed)
  X = mvrnorm(100, rep(0, 30), block.matrix, empirical = TRUE)
  
  #add colnames
  features <- c(sprintf("X%d", seq(1,30)))
  colnames(X) = features
  
  #define signal-to-noise ratio/compute variance of error
  snr = 1.8
  var.error = (t(beta) %*% block.matrix %*% beta)/snr
  
  #simulate y = X^T * beta
  set.seed(-1)
  y = scale(as.numeric(X %*% beta + rnorm(100, 0, sqrt(var.error))))
  
  #return object
  sim.data = cbind(y,X)
  colnames(sim.data)[1] = "y"
  
  return(sim.data)
}





#### ===================== function for compute sensitivity/specificity
metrics = function(pred_vec, true_vec){
  TP = sum(pred_vec  == 1 & true_vec == 1)
  FP = sum(pred_vec  == 1 & true_vec == 0)
  FN = sum(pred_vec  == 0 & true_vec == 1)
  TN = sum(pred_vec  == 0 & true_vec == 0)
  
  sensitivity = TP/(TP+FN)
  specificity = TN/(FP+TN)
  
  return(c(sensitivity, specificity))
}

#### ===================== function for compute cve for group Lasso/MCP/SCAD
cve = function(x, y, group){
  data = list(x = x, y = y)
  
  set.seed(-2)
  cvglasso = cv.grpreg(x, y, group, penalty = "grLasso")
  
  set.seed(-2)
  cvgMCP = cv.grpreg(x, y, group, penalty = "grMCP")
  
  set.seed(-2)
  cvgSCAD = cv.grpreg(x, y, group, penalty = "grSCAD")
  
  set.seed(-2)
  cvSGL = cvSGL(data = data, index = group, type = "linear")
  
  cve.data = c(
    min(cvglasso$cve),
    min(cvgMCP$cve),
    min(cvgSCAD$cve),
    min(cvSGL$lldiff)*2/nrow(y)
  )
  return(cve.data)
}

#### ===================== function for compute coefficient for group Lasso/MCP/SCAD
coef = function(x, y, group){
  data = list(x = x, y = y)
  
  set.seed(-2)
  cvglasso = cv.grpreg(x, y, group, penalty = "grLasso")
  
  set.seed(-2)
  cvgMCP = cv.grpreg(x, y, group, penalty = "grMCP")
  
  set.seed(-2)
  cvgSCAD = cv.grpreg(x, y, group, penalty = "grSCAD")
  
  set.seed(-2)
  cvSGL = cvSGL(data = data, index = group, type = "linear")
  
  coef = list(
    ifelse(cvglasso$fit$beta[, which.min(cvglasso$cve)][-1] != 0 , 1, 0),
    ifelse(cvgMCP$fit$beta[, which.min(cvgMCP$cve)][-1] != 0 , 1, 0),
    ifelse(cvgSCAD$fit$beta[, which.min(cvgSCAD$cve)][-1] != 0 , 1, 0),
    ifelse(cvSGL$fit$beta[, which.min(cvSGL$lldiff)] != 0, 1, 0)
  )
  return(coef)
}








