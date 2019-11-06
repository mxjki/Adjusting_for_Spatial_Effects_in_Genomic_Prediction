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

WeightMat = function(k) {
    Weight = matrix(0, k, k)
    Weight[1, 1:2] = c(1, -1)
    Weight[k, (k - 1):k] = c(-1, 1)
    for (i in 2:(k - 1)) {
        Weight[i, (i - 1):(i + 1)] = c(-1, 2, -1)
    }
    return(Weight)
}

spationinfo = function(U, V, m1p, m2p) {
    Nu = kronecker(diag(1, m2p), WeightMat(m1p))
    Nv = kronecker(WeightMat(m2p), diag(1, m1p))
    
    eigenWm1p = eigen(WeightMat(m1p))
    lambda = eigenWm1p$values
    lambdaz = eigenWm1p$vectors
    
    eigenWm2p = eigen(WeightMat(m2p))
    nu = eigenWm2p$values
    nuz = eigenWm2p$vectors
    
    return(list(Nu = Nu, Nv = Nv, lambda = lambda, nu = nu, lambdaz = lambdaz, nuz = nuz))
}

WW = function(m1p, m2p, lambda00, beta, Nu, Nv) {
    return(lambda00 * diag(1, m1p * m2p) + beta * Nu + (0.5 - beta) * Nv)
}

WWinverse = function(m1p, m2p, lambda00, beta, lambda, nu, combinez) {
    return(combinez %*% diag(1/(lambda00 + beta * rep(lambda, m2p) + (0.5 - beta) * rep(nu, 
        each = m1p))) %*% t(combinez))
}

getH = function(n, U, V, m1p, m2p) {
    H = matrix(0, n, m1p * m2p)
    Vorder = V - min(V) + 1
    Uorder = U - min(U) + 1
    
    for (i in 1:n) {
        # expand along columns
        H[i, (1 + Uorder[i]) * m1p + Vorder[i] + 2] = 1
    }
    
    return(H)
}

get_cov_Field_no_sub = function(DistanceX, SigmaS_std, tot_num, n, para) {
    SigmaX = para[1] * as.matrix(exp(-DistanceX/(para[4]^2 * tot_num)))
    
    SigmaS = para[3] * SigmaS_std
    
    SigmaE = diag(rep(para[2]^2, n))
    
    Sigma = SigmaX + SigmaS + SigmaE
    
    return(list(SigmaX = SigmaX, SigmaS = SigmaS, SigmaE = SigmaE, Sigma = Sigma))
}


getmuhat = function(n, cholSigma, Y) {
    
    Sigmay = backsolve(cholSigma, Y, transpose = TRUE)
    
    Sigma1 = backsolve(cholSigma, rep(1, n), transpose = TRUE)
    
    muhat = sum(Sigmay * Sigma1)/sum(Sigma1^2)
    
    return(muhat)
}

UV_info = getUV_info(data$col, data$row)
spatio = spationinfo(UV_info$U, UV_info$V, UV_info$m1p, UV_info$m2p)
combinez = kronecker(spatio$nuz, spatio$lambdaz)
H = getH(UV_info$n, UV_info$U, UV_info$V, UV_info$m1p, UV_info$m2p)

myY = data[, (7 + res_num)]

Ymu = (myY - mean(myY))/sd(myY)

DistanceX = as.matrix(readRDS(paste0("./dist.rds")))
DistanceX_all = DistanceX[data$trt_id, data$trt_id]

tot_num = 102324

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl, cores = detectCores() - 1)

GRF_wheat_breeding = foreach(sim_num = 1:1000) %dopar% {
    
    library("Matrix")
    
    optimpara = rep(0, 5)
    result = rep(0, 2)
    
    s = 1/7  #sample of inference, e.g. set it to 1/7 for seven fold cross validation
    n.missing = round(UV_info$n * s)
    
    set.seed(sim_num)
    
    sample.missing = sample(1:UV_info$n, n.missing)
    sample.obs = c(1:UV_info$n)[-sample.missing]
    
    Ymu_obs = Ymu[sample.obs]
    
    sep_loglikelihood_obs = function(para) {
        SigmaS = as.matrix(H %*% WWinverse(UV_info$m1p, UV_info$m2p, 0.001, para[1], 
            spatio$lambda, spatio$nu, combinez) %*% t(H))
        diagS = diag(SigmaS)
        SigmaS_std = 1/sqrt(diagS) * SigmaS/sqrt(diagS)
        
        Sigma = get_cov_Field_no_sub(DistanceX_all[sample.obs, sample.obs], SigmaS_std[sample.obs, 
            sample.obs], tot_num, length(sample.obs), para[-1])$Sigma
        
        cholSigma = chol(forceSymmetric(Sigma))
        
        muhat = getmuhat(length(sample.obs), cholSigma, Ymu_obs)
        
        xmu = backsolve(cholSigma, Ymu_obs - muhat * rep(1, length(sample.obs)), transpose = TRUE)
        
        return(2 * sum(log(diag(cholSigma))) + sum(xmu^2))
    }
    
    t1 = proc.time()
    optimY_obs = optim(c(0.25, 0.5, 0.6, 0.3, 0.6), sep_loglikelihood_obs, method = c("L-BFGS-B"), 
        lower = c(0, 0.01, 0.001, 0, 0.001), upper = c(0.5, 2, 5, 2, 2))
    tused = (proc.time() - t1)[3]
    
    optimpara = optimY_obs$par
    print(optimpara)
    
    pred_GRF = function(para) {
        
        SigmaS = as.matrix(H %*% WWinverse(UV_info$m1p, UV_info$m2p, 0.001, para[1], 
            spatio$lambda, spatio$nu, combinez) %*% t(H))
        diagS = diag(SigmaS)
        SigmaS_std = 1/sqrt(diagS) * SigmaS/sqrt(diagS)
        
        SigFarmCPU = get_cov_Field_no_sub(DistanceX_all, SigmaS_std, tot_num, UV_info$n, 
            para[-1])$Sigma
        
        cholSigmaFarmCPU = chol(forceSymmetric(SigFarmCPU[sample.obs, sample.obs]))
        
        muhat = getmuhat(length(sample.obs), cholSigmaFarmCPU, Ymu_obs)
        
        pred_GRF = muhat + as.vector(SigFarmCPU[sample.missing, sample.obs] %*% solve(SigFarmCPU[sample.obs, 
            sample.obs]) %*% as.matrix(Ymu_obs - muhat))
        
        return(pred_Yhat = pred_GRF)
    }
    
    pred_Y_GRF = pred_GRF(optimpara)
    
    result[1] = sum((Ymu[sample.missing] - pred_Y_GRF)^2 * sd(myY)^2)/length(sample.missing)
    result[2] = cor(pred_Y_GRF, Ymu[sample.missing])
    
    saveRDS(optimpara, paste0("./GRF/2011GRF_optimpara_res_", res_num, "_sim_", sim_num, 
        ".rds"))
    saveRDS(result, paste0("./GRF/2011GRF_res_", res_num, "_sim_", sim_num, ".rds"))
}

stopCluster(cl)
