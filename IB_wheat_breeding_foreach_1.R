library("parallel")
library("foreach")
library("doParallel")
library("regress")

res_num = 1
# G = read.table(paste0('num_impute.txt'), header = T)
data = read.table("./pheno_data2011.txt", header = T)

myY = data[, c(2, 3, 5, 7 + res_num)]
# myY=cbind(myY,G[data$trt_id,])
Names = colnames(myY)
resp_num = 4
Names[resp_num] = "response"

# Names[7:(7+102323)]=seq(1,102324)

colnames(myY) = Names

# for (cv_num in c(1, 3, 5, 10)) {

cv_num = 3

myKI = read.csv(paste0("./Kin_VanRaden_", res_num, ".csv"), header = FALSE)
myCV = read.csv(paste0("./PCA_10.csv"), header = TRUE)[, 1:(1 + cv_num)]

myY = cbind(myY, myCV[, 2:(1 + cv_num)])

myKI = myKI[, -1]

myKI = myKI[complete.cases(myY[, resp_num]), complete.cases(myY[, resp_num])]

myY = myY[complete.cases(myY[, resp_num]), ]

n = dim(myY)[1]

myY = as.data.frame(myY)

myY$rep = factor(myY$rep)

myY$bl = factor(myY$bl)

myY$trt_id = factor(myY$trt_id)

IB_wheat_breeding = foreach(sim_num = 1:1000) %dopar% {
    
    library("regress")
    
    s = 1/7  #sample of inference, e.g. set it to 1/7 for seven fold cross validation
    n.missing = round(n * s)
    set.seed(sim_num)
    
    sample.missing = sample(1:n, n.missing)
    sample.obs = c(1:n)[-sample.missing]
    
    myY_obs = myY[sample.obs, ]
    myKI_obs = myKI[sample.obs, sample.obs]
    
    A_obs = model.matrix(~rep, data = myY_obs)
    V1_obs = tcrossprod(A_obs)  # fast way to compute AA'
    B_obs = model.matrix(~bl:rep, data = myY_obs)
    V2_obs = tcrossprod(B_obs)
    K_obs = as.matrix(myKI_obs)
    
    fit = regress(response ~ PC1 + PC2 + PC3, Vformula = ~K_obs + V1_obs + V2_obs, data = myY_obs)
    
    A = model.matrix(~rep, data = myY)
    V1 = tcrossprod(A)  # fast way to compute AA'
    B = model.matrix(~bl:rep, data = myY)
    V2 = tcrossprod(B)
    K = as.matrix(myKI)
    
    G = fit$sigma[1] * K + fit$sigma[2] * V1 + fit$sigma[3] * V2
    
    Simga = G + fit$sigma[4] * diag(dim(myY)[1])
    
    Simgaobs = Simga[sample.obs, sample.obs]
    
    pred_IB = as.matrix(cbind(1, myY[sample.missing, 5:7])) %*% matrix(fit$beta, ncol = 1) + 
        G[sample.missing, sample.obs] %*% solve(Simgaobs) %*% (myY_obs$response - as.matrix(cbind(1, 
            myY[sample.obs, 5:7])) %*% matrix(fit$beta, ncol = 1))
    
    result = cor(pred_IB, myY$response[sample.missing])
    
    saveRDS(result, paste0("./IB/2011IB_res_", res_num, "_cv_", cv_num, "_sim_", sim_num, 
        ".rds"))
    
}
stopCluster(cl)
# }
