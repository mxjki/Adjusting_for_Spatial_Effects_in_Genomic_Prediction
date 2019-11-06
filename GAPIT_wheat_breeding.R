for (cv_num in c(1, 3, 5, 10)) {
    
    library(gplots)
    library(genetics)
    library(EMMREML)
    library(compiler)  #this library is already installed in R
    source("http://zzlab.net/GAPIT/gapit_functions.txt")
    source("http://zzlab.net/GAPIT/emma.txt")
    
    library("Matrix")
    data = read.table("pheno_data2011.txt", header = T)
    result = matrix(0, 8, 2)
    
    for (res_num in 1:8) {
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
            
            return(list(Nu = Nu, Nv = Nv, lambda = lambda, nu = nu, lambdaz = lambdaz, 
                nuz = nuz))
        }
        
        WW = function(m1p, m2p, lambda00, beta, Nu, Nv) {
            return(lambda00 * diag(1, m1p * m2p) + beta * Nu + (0.5 - beta) * Nv)
        }
        
        WWinverse = function(m1p, m2p, lambda00, beta, lambda, nu, combinez) {
            return(combinez %*% diag(1/(lambda00 + beta * rep(lambda, m2p) + (0.5 - beta) * 
                rep(nu, each = m1p))) %*% t(combinez))
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
        
        # GAPIT functions
        pred_GAPIT = function(n, sample.missing, myKI, myCV, myY) {
            
            Y.raw = myY  #choos a trait
            Y.raw = Y.raw[!is.na(Y.raw[, 2]), ]  #Remove missing data
            
            if (length(sample.missing) > 0) {
                Y0 = Y.raw[-sample.missing, ]
            } else {
                Y0 = Y.raw
            }
            
            myGAPIT <- GAPIT(Y = Y0, KI = myKI, CV = myCV, group.from = n, group.to = n, 
                kinship.cluster = c("average"), kinship.group = c("Mean"))
            prediction = myGAPIT$Pred
            
            return(list(pred_yhat = prediction))
            
        }
        
        # GAPIT
        
        UV_info = getUV_info(data$col, data$row)
        n = UV_info$n
        spatio = spationinfo(UV_info$U, UV_info$V, UV_info$m1p, UV_info$m2p)
        combinez = kronecker(spatio$nuz, spatio$lambdaz)
        H = getH(UV_info$n, UV_info$U, UV_info$V, UV_info$m1p, UV_info$m2p)
        
        myY = data.frame(Taxa = paste0("X", data$trt_id), Y = data[, (7 + res_num)])
        
        s = 1/7  #sample of inference, e.g. set it to 1/7 for seven fold cross validation
        n.missing = round(UV_info$n * s)
        
        set.seed(sim_num)
        
        sample.missing = sample(1:UV_info$n, n.missing)
        sample.obs = c(1:UV_info$n)[-sample.missing]
        
        myKI = read.csv(paste0("./Kin_VanRaden_", res_num, ".csv"), header = FALSE)
        myCV = read.csv(paste0("./PCA_10.csv"), header = TRUE)[, 1:(1 + cv_num)]
        
        pred_GAPIT_Y = pred_GAPIT(UV_info$n, sample.missing, myKI, myCV, myY)$pred_yhat
        
        result[res_num, 1] = sum((myY[sample.missing, 2] - pred_GAPIT_Y$Prediction[data$trt_id[sample.missing]])^2)/length(sample.missing)
        result[res_num, 2] = cor(pred_GAPIT_Y$Prediction[data$trt_id[sample.missing]], 
            myY[sample.missing, 2])
    }
    saveRDS(result, paste0("./GAPIT/2011GAPIT_cv_", cv_num, "_sim_", sim_num, ".rds"))
}
