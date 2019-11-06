library("parallel")
library("foreach")
library("doParallel")

data = read.table("pheno_data2011.txt", header = T)
res_num = 1

data = data[complete.cases(data[, (7 + res_num)]), ]

getUV_info = function(U, V) {
    
    n = length(U)
    
    # m1 for row numbers, m2 for column numbers
    m1 = length(unique(V))
    m2 = length(unique(U))
    m1p = m1 + 4
    m2p = m2 + 4
    
    return(list(U = U, V = V, n = n, m1 = m1, m2 = m2, m1p = m1p, m2p = m2p))
}

get_cov_Field_no_sub_no_spatial = function(DistanceX, tot_num, n, para) {
    SigmaX = para[1] * as.matrix(exp(-DistanceX/(para[3]^2 * tot_num)))
    
    SigmaE = diag(rep(para[2]^2, n))
    
    Sigma = SigmaX + SigmaE
    
    return(list(SigmaX = SigmaX, SigmaE = SigmaE, Sigma = Sigma))
}


getmuhat = function(n, cholSigma, Y) {
    
    Sigmay = backsolve(cholSigma, Y, transpose = TRUE)
    
    Sigma1 = backsolve(cholSigma, rep(1, n), transpose = TRUE)
    
    muhat = sum(Sigmay * Sigma1)/sum(Sigma1^2)
    
    return(muhat)
}

UV_info = getUV_info(data$col, data$row)

myY = data[, (7 + res_num)]

Ymu = (myY - mean(myY))/sd(myY)

DistanceX = as.matrix(readRDS(paste0("./dist.rds")))
DistanceX_all = DistanceX[data$trt_id, data$trt_id]

tot_num = 102324

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl, cores = detectCores() - 1)

GRF_wheat_breeding = foreach(sim_num = 1:1000) %dopar% {
    
    library("Matrix")
    
    optimpara = rep(0, 3)
    result = rep(0, 2)
    
    s = 1/7  #sample of inference, e.g. set it to 1/7 for seven fold cross validation
    n.missing = round(UV_info$n * s)
    
    set.seed(sim_num)
    
    sample.missing = sample(1:UV_info$n, n.missing)
    sample.obs = c(1:UV_info$n)[-sample.missing]
    
    Ymu_obs = Ymu[sample.obs]
    
    sep_loglikelihood_obs = function(para) {
        
        Sigma = get_cov_Field_no_sub_no_spatial(DistanceX_all[sample.obs, sample.obs], 
            tot_num, length(sample.obs), para)$Sigma
        
        cholSigma = chol(forceSymmetric(Sigma))
        
        muhat = getmuhat(length(sample.obs), cholSigma, Ymu_obs)
        
        xmu = backsolve(cholSigma, Ymu_obs - muhat * rep(1, length(sample.obs)), transpose = TRUE)
        
        return(2 * sum(log(diag(cholSigma))) + sum(xmu^2))
    }
    
    t1 = proc.time()
    optimY_obs = optim(c(0.5, 0.6, 0.6), sep_loglikelihood_obs, method = c("L-BFGS-B"), 
        lower = c(0.01, 0.001, 0.001), upper = c(2, 5, 2))
    tused = (proc.time() - t1)[3]
    
    optimpara = optimY_obs$par
    print(optimpara)
    
    pred_GRF = function(para) {
        
        
        SigFarmCPU = get_cov_Field_no_sub_no_spatial(DistanceX_all, tot_num, UV_info$n, 
            para)$Sigma
        
        cholSigmaFarmCPU = chol(forceSymmetric(SigFarmCPU[sample.obs, sample.obs]))
        
        muhat = getmuhat(length(sample.obs), cholSigmaFarmCPU, Ymu_obs)
        
        pred_GRF = muhat + as.vector(SigFarmCPU[sample.missing, sample.obs] %*% solve(SigFarmCPU[sample.obs, 
            sample.obs]) %*% as.matrix(Ymu_obs - muhat))
        
        return(pred_Yhat = pred_GRF)
    }
    
    pred_Y_GRF = pred_GRF(optimpara)
    
    result[1] = sum((Ymu[sample.missing] - pred_Y_GRF)^2 * sd(myY)^2)/length(sample.missing)
    result[2] = cor(pred_Y_GRF, Ymu[sample.missing])
    
    saveRDS(optimpara, paste0("./GRF/2011GRF_no_spatial_optimpara_res_", res_num, "_sim_", 
        sim_num, ".rds"))
    saveRDS(result, paste0("./GRF/2011GRF_no_spatial_res_", res_num, "_sim_", sim_num, 
        ".rds"))
}

stopCluster(cl)
