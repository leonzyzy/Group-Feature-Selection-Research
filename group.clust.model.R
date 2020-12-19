#### ===================== define tuning para
#define beta matrix; 1, 2, 9, 10, 23, 25 not equals to 0
beta = matrix(rep(c(3.8, -2.8, 3.5, 0, -0.8, 2.2, -3.2, -4.2, 0, 1.8, 0.3, 1.1, -0.6, 0),
                  c(1, 1, 1, 4, 1, 1, 1, 1, 10, 1, 1, 1, 1, 5)), nrow = 30)

block.group = c(3,4,4,3,4,3,4,3,2) 
block.num = 1:30
nrow = 30
ncol = 30
r = seq(0.8, 0.95, by = 0.05)


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
  X = scale(mvrnorm(100, rep(0, 30), block.matrix, empirical = TRUE))
  
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

total.list = list()

for(r in r){
  sim.data.list = lapply(1:100, function(x) sim_data(r, block.num, block.group, nrow, ncol, beta, x))
  total.list = list.append(total.list, sim.data.list)
}

#define coefficent matrix
glasso.coef.matrix = matrix(NA, rep(30*100), nrow = 100, ncol = 30)
gMCP.coef.matrix = matrix(NA, rep(30*100), nrow = 100, ncol = 30)
gSCAD.coef.matrix = matrix(NA, rep(30*100), nrow = 100, ncol = 30)
SGL.coef.matrix = matrix(NA, rep(30*100), nrow = 100, ncol = 30)

colnames(glasso.coef.matrix) = c(sprintf("X%d", seq(1,30)))
colnames(gMCP.coef.matrix) = c(sprintf("X%d", seq(1,30)))
colnames(gSCAD.coef.matrix) = c(sprintf("X%d", seq(1,30)))
colnames(SGL.coef.matrix) = c(sprintf("X%d", seq(1,30)))

#define cve vectors
cve.glasso = c()
cve.gMCP = c()
cve.gSCAD = c()
cve.SGL = c()

#define coefficent matrix list
glasso.coef.list = list()
gMCP.coef.list = list()
gSCAD.coef.list= list()
SGL.coef.list = list()

#define model avg cve list
cve.avg.list = c()

for(i in 1:2){
  for(j in 1:100){
    model.data = total.list[[i]][[j]]
    X = model.data[,-1]
    y = model.data[, 1]
    data = list(x = X, y = y)
    
    #clustering method
    set.seed(j)
    tree = hclustvar(X.quanti = X) 
    stab = stability(tree, 50)
    
    if(max(stab$meanCR) <= 0.1){
      cluster = seq(1:ncol(X))
      
      #cross validation for group lasso/MCP/SCAD
      set.seed(j)
      cvglasso = cv.grpreg(X, y, cluster, penalty = "grLasso")
      cve.glasso[j] = min(cvglasso$cve)
      
      set.seed(j)
      cvgMCP = cv.grpreg(X, y, cluster, penalty = "grMCP")
      cve.gMCP[j] = min(cvgMCP$cve)
      
      set.seed(j)
      cvgSCAD = cv.grpreg(X, y, cluster, penalty = "grSCAD")
      cve.gSCAD[j] = min(cvgSCAD$cve)
      
      set.seed(j)
      cvSGL = cvSGL(data = data, index = cluster, type = "linear")
      cve.SGL[j] = 2*min(cvSGL$lldiff)/100
      
      # #coefficient matirx
      glasso.coef.matrix[j,] = ifelse(cvglasso$fit$beta[, which.min(cvglasso$cve)][-1] != 0 , 1, 0)
      gMCP.coef.matrix[j,] = ifelse(cvgMCP$fit$beta[, which.min(cvgMCP$cve)][-1] != 0 , 1, 0)
      gSCAD.coef.matrix[j,] = ifelse(cvgSCAD$fit$beta[, which.min(cvgSCAD$cve)][-1] != 0 , 1, 0)
      SGL.coef.matrix[j,] = ifelse(cvSGL$fit$beta[, which.min(cvSGL$lldiff)] != 0, 1, 0)
    }
    else{
      #cut tree using the optimal adjusted rand index 
      P <- cutreevar(tree, as.numeric(which.max(stab$meanCR) + 1))
      var.name = names(P$cluster)
      cluster = as.integer(P$cluster)
      cluster.data = data.frame(var.name, cluster)
      cluster.data = cluster.data[order(cluster.data$cluster),]
      columns.order = as.integer(rownames(cluster.data))
      
      #cross validation for group lasso/MCP/SCAD
      set.seed(j)
      cvglasso = cv.grpreg(X, y, cluster.data$cluster, penalty = "grLasso")
      cve.glasso[j] = min(cvglasso$cve)
      
      set.seed(j)
      cvgMCP = cv.grpreg(X, y, cluster.data$cluster, penalty = "grMCP")
      cve.gMCP[j] = min(cvgMCP$cve)
      
      set.seed(j)
      cvgSCAD = cv.grpreg(X, y, cluster.data$cluster, penalty = "grSCAD")
      cve.gSCAD[j] = min(cvgSCAD$cve)
      
      set.seed(j)
      cvSGL = cvSGL(data = data, index = cluster.data$cluster, type = "linear")
      cve.SGL[j] = min(2*cvSGL$lldiff)/100
      
      # #coefficient matirx
      glasso.coef.matrix[j,] = ifelse(cvglasso$fit$beta[, which.min(cvglasso$cve)][-1] != 0 , 1, 0)
      gMCP.coef.matrix[j,] = ifelse(cvgMCP$fit$beta[, which.min(cvgMCP$cve)][-1] != 0 , 1, 0)
      gSCAD.coef.matrix[j,] = ifelse(cvgSCAD$fit$beta[, which.min(cvgSCAD$cve)][-1] != 0 , 1, 0)
      SGL.coef.matrix[j,] = ifelse(cvSGL$fit$beta[, which.min(cvSGL$lldiff)] != 0, 1, 0)
    }
  }
  
  #data frame for cve
  cve.data.avg = data.frame("GLasso" = mean(cve.glasso), "GLasso.std"= sd(cve.glasso),
                            "GMCP" = mean(cve.gMCP), "GMCP.std"= sd(cve.gMCP),
                            "GSCAD" = mean(cve.gSCAD), "GSCAD.std"= sd(cve.gSCAD),
                            "SparseGLasso" = mean(cve.SGL), "SparseGLasso.std" = sd(cve.SGL))
  
  cve.avg.list = list.append(cve.avg.list, cve.data.avg)
  
  #add the coef into list
  glasso.coef.list = list.append(glasso.coef.list, glasso.coef.matrix)
  gMCP.coef.list = list.append(gMCP.coef.list, gMCP.coef.matrix)
  gSCAD.coef.list= list.append(gSCAD.coef.list, gSCAD.coef.matrix)
  SGL.coef.list = list.append(SGL.coef.list, SGL.coef.matrix)
}

#data frame to load the cve.avg
glasso.cve = data.frame("GLasso" = cve.avg.list[[1]], "GLasso.std"= cve.avg.list[[2]],
                        "GMCP" = cve.avg.list[[3]], "GMCP.std"= cve.avg.list[[4]],
                        "GSCAD" = cve.avg.list[[5]], "GSCAD.std"= cve.avg.list[[6]],
                        "SparseGLasso" = cve.avg.list[[7]], "SparseGLasso.std" = cve.avg.list[[8]])

group.model.cve = data.frame("GLasso" = c(cve.avg.list[[1]], cve.avg.list[[9]][[1]]),
                             "GMCP" = c(cve.avg.list[[3]], cve.avg.list[[9]][[3]]),
                             "GSCAD" = c(cve.avg.list[[5]], cve.avg.list[[9]][[5]]),
                             "SparseGLasso" = c(cve.avg.list[[7]], cve.avg.list[[9]][[7]]))

write.csv(group.model.cve, "group.model.cve3.csv", row.names = F)



#define the matrix as true coefficient for testing 
true.coef.matrix = matrix(0, rep(30), nrow = 1, ncol = 30)
true.coef.matrix[, which(rep(c(3.8, -2.8, 3.5, 0, -0.8, 2.2, -3.2, -4.2, 0, 1.8, 0.3, 1.1, -0.6, 0),
                             c(1, 1, 1, 4, 1, 1, 1, 1, 10, 1, 1, 1, 1, 5)) != 0)] = 1
colnames(true.coef.matrix) = c(sprintf("X%d", seq(1,30)))


#function for compute sensitivity/specificity
metrics = function(pred_vec, true_vec){
  TP = sum(pred_vec  == 1 & true_vec == 1)
  FP = sum(pred_vec  == 1 & true_vec == 0)
  FN = sum(pred_vec  == 0 & true_vec == 1)
  TN = sum(pred_vec  == 0 & true_vec == 0)
  
  sensitivity = TP/(TP+FN)
  specificity = TN/(FP+TN)
  
  return(c(sensitivity, specificity))
}

#loop for group lasso sensitivity/specificity
sensitivity.glasso = matrix(NA, nrow = 100, ncol = 8)
specificity.glasso = matrix(NA, nrow = 100, ncol = 8)

sensitivity.gMCP = matrix(NA, nrow = 100, ncol = 8)
specificity.gMCP = matrix(NA, nrow = 100, ncol = 8)

sensitivity.gSCAD = matrix(NA, nrow = 100, ncol = 8)
specificity.gSCAD = matrix(NA, nrow = 100, ncol = 8)

sensitivity.SGL = matrix(NA, nrow = 100, ncol = 8)
specificity.SGL = matrix(NA, nrow = 100, ncol = 8)

for(i in 1:2){ 
  sensitivity.glasso[,i] = sapply(1:100, function(x) metrics(glasso.coef.list[[i]][x,], true.coef.matrix[1,])[1])
  specificity.glasso[,i] = sapply(1:100, function(x) metrics(glasso.coef.list[[i]][x,], true.coef.matrix[1,])[2])
  
  sensitivity.gMCP[,i] = sapply(1:100, function(x) metrics(gMCP.coef.list[[i]][x,], true.coef.matrix[1,])[1])
  specificity.gMCP[,i] = sapply(1:100, function(x) metrics(gMCP.coef.list[[i]][x,], true.coef.matrix[1,])[2])
  
  sensitivity.gSCAD[,i] = sapply(1:100, function(x) metrics(gSCAD.coef.list[[i]][x,], true.coef.matrix[1,])[1])
  specificity.gSCAD[,i] = sapply(1:100, function(x) metrics(gSCAD.coef.list[[i]][x,], true.coef.matrix[1,])[2])
  
  sensitivity.SGL[,i] = sapply(1:100, function(x) metrics(SGL.coef.list[[i]][x,], true.coef.matrix[1,])[1])
  specificity.SGL[,i] = sapply(1:100, function(x) metrics(SGL.coef.list[[i]][x,], true.coef.matrix[1,])[2])
}
#group lasso
sensitivity.glasso.avg = sapply(1:2, function(x) mean(sensitivity.glasso[,x]))
specificity.glasso.avg = sapply(1:2, function(x) mean(specificity.glasso[,x]))
sensitivity.glasso.std = sapply(1:2, function(x) sd(sensitivity.glasso[,x]))
specificity.glasso.std = sapply(1:2, function(x) sd(specificity.glasso[,x]))

#group MCP
sensitivity.gMCP.avg = sapply(1:2, function(x) mean(sensitivity.gMCP[,x]))
specificity.gMCP.avg = sapply(1:2, function(x) mean(specificity.gMCP[,x]))
sensitivity.gMCP.std = sapply(1:2, function(x) sd(sensitivity.gMCP[,x]))
specificity.gMCP.std = sapply(1:2, function(x) sd(specificity.gMCP[,x]))

#group SCAD
sensitivity.gSCAD.avg = sapply(1:2, function(x) mean(sensitivity.gSCAD[,x]))
specificity.gSCAD.avg = sapply(1:2, function(x) mean(specificity.gSCAD[,x]))
sensitivity.gSCAD.std = sapply(1:2, function(x) sd(sensitivity.gSCAD[,x]))
specificity.gSCAD.std = sapply(1:2, function(x) sd(specificity.gSCAD[,x]))

#SGL
sensitivity.SGL.avg = sapply(1:2, function(x) mean(sensitivity.SGL[,x]))
specificity.SGL.avg = sapply(1:2, function(x) mean(specificity.SGL[,x]))
sensitivity.SGL.std = sapply(1:2, function(x) sd(sensitivity.SGL[,x]))
specificity.SGL.std = sapply(1:2, function(x) sd(specificity.SGL[,x]))


#conclude a data frame
performance.coef.TPR = data.frame("Correlation" = seq(0.8, 0.85, by = 0.05)[1:2],
                                  "GLasso.TPR" = sensitivity.glasso.avg,
                                  "GMCP.TPR" = sensitivity.gMCP.avg,
                                  "GSCAD.TPR" = sensitivity.gSCAD.avg,
                                  "SGL.TPR" = sensitivity.SGL.avg)

performance.coef.TNR = data.frame("Correlation" = seq(0.8, 0.85, by = 0.05)[1:2],
                                  "GLasso.TNR" = specificity.glasso.avg,
                                  "GMCP.TNR" = specificity.gMCP.avg,
                                  "GSCAD.TNR" = specificity.gSCAD.avg,
                                  "SGL.TNR" = specificity.SGL.avg)

write.csv(performance.coef.TPR, "performance.coef.TPR3.csv", row.names = F)
write.csv(performance.coef.TNR, "performance.coef.TNR3.csv", row.names = F)

###combine csv
performance_coef_TPR1 <- read_csv("performance.coef.TPR1.csv")
performance_coef_TNR1 <- read_csv("performance.coef.TNR1.csv")

performance.coef.TPR.overall = rbind(performance_coef_TPR1, performance.coef.TPR)
performance.coef.TNR.overall = rbind(performance_coef_TNR1, performance.coef.TNR)
group_model_cve = rbind(group_model_cve1, group.model.cve[-9,])


write.csv(performance.coef.TPR.overall, "performance.coef.TPR.overall.csv", row.names = F)
write.csv(performance.coef.TNR.overall, "performance.coef.TNR.overall.csv", row.names = F)
write.csv(group_model_cve, "group_model_cve.csv", row.names = F)


##compute porportion for 0.50
true.index = which(rep(c(3.8, -2.8, 3.5, 0, -0.8, 2.2, -3.2, -4.2, 0, 1.8, 0.3, 1.1, -0.6, 0),
                       c(1, 1, 1, 4, 1, 1, 1, 1, 10, 1, 1, 1, 1, 5)) != 0)

glasso.coef.percent = (colSums(glasso.coef.list[[4]])/100)[true.index]
gMCP.coef.percent = (colSums(gMCP.coef.list[[4]])/100)[true.index]
gSCAD.coef.percent = (colSums(gSCAD.coef.list[[4]])/100)[true.index]
SGL.coef.percent = (colSums(SGL.coef.list[[4]])/100)[true.index]

coef.percent3 = data.frame("Glasso" = glasso.coef.percent,
                           "GMCP" = gMCP.coef.percent,
                           "GSCAD" = gSCAD.coef.percent,
                           "SGL" = SGL.coef.percent)
write.csv(coef.percent3, "coef.percent3.csv")
