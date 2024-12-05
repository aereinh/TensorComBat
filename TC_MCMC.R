##### Functions to run longitudinal Tensor Combat harmonization model
## Model (Voxel-level): Y_ij(v) = M(v) + \Beta_j(v) + \Gamma_i(v) + \sum_{q=1}^K \Theta_q x_jq + \epsilon_ij(v), i indexes scanner and j indexes subject
## M ~ Population Intercept
## Beta_j ~ Subject-specific intercept
## Gamma_i ~ Scanner-specific intercept
## Theta_q ~ Population-level covariate effects
## epsilon_ij(v) ~ Residual; for given scanner (i), follows finite mixture across voxels (v)

library(Rcpp)
library(e1071)
library(parallel)
library(glmnet)
# set directory where 'BTRR_harm_long'
#setwd('/Users/alecreinhardt/Dropbox/MD Anderson/lBTRR/')

sourceCpp('TC_MCMC_c.cpp')


# Tensor Combat Functions -------------------------------------------------

mdev_cred_int <- function(Gamma_k_mcmc,alpha=0.05,missing_vox=NULL) {
  if (!is.null(missing_vox)) {
    Gamma_k_mcmc[,missing_vox] <- NA
  }
  nvox <- dim(Gamma_k_mcmc)[2]
  mean_vox <- apply(Gamma_k_mcmc,2,mean)
  s_alpha_vox <- apply(Gamma_k_mcmc,2,function(x) quantile(x,c(alpha/2,1-alpha/2),na.rm = T))
  slow_alpha <- max(mean_vox-s_alpha_vox[1,],na.rm=T)
  shigh_alpha <- max(s_alpha_vox[2,]-mean_vox,na.rm=T)
  credint_vox <- rbind(mean_vox-slow_alpha,
                       mean_vox+shigh_alpha)
  return(credint_vox[1,]*credint_vox[2,]>0)
}

getCoef_mcmc_R = function(gamma_store) {
  niter = length(gamma_store)
  gamma_iter1 = gamma_store[[1]]
  isMulti = is.list(gamma_iter1[[1]])
  if (isMulti==FALSE) {
    p = unlist(lapply(gamma_iter1,function(x) dim(x)[1]))
    Gamma_mcmc = array(dim=c(niter,prod(p)))
    for (iter in 1:niter) {
      Gamma_mcmc[iter,] = getGamma_cpp(gamma_store[[iter]])
    }
  } else {
    Q = length(gamma_iter1)
    p = unlist(lapply(gamma_iter1[[1]], function(x) dim(x)[1]))
    Gamma_mcmc = array(dim=c(niter,Q,prod(p)))
    for (iter in 1:niter) {
      for (q in 1:Q) {
        Gamma_mcmc[iter,q,] = getGamma_cpp(gamma_store[[iter]][[q]])
      }
    }
  }
  return(Gamma_mcmc)
}

getSigSq_mcmc_R = function(Scanner_ID, ssq_store, Zeta_store) {
  niter = length(ssq_store)
  nobs = length(Scanner_ID)
  V = dim(Zeta_store[[1]])[2]
  SigSq_mcmc = array(dim=c(niter,nobs,V))
  for (iter in 1:niter) {
    SigSq_mcmc[iter,,] = getSigSq_harm_cpp(Scanner_ID,ssq_store[[iter]],Zeta_store[[iter]])
  }
  return(SigSq_mcmc)
}

# compacting/uncompacting Zeta and Z_subj
uncomp_Zeta = function(Zeta_comp, nScan, V, H) {
  vec = rep(0,nScan*V)
  for (h in 1:H) {
    vec[Zeta_comp[[h]]] = h
  }
  arr = array(vec,dim=c(nScan,V))
  return(arr)
}

comp_Zeta = function(Zeta, H) {
  Zeta_vec = c(Zeta)
  Zeta_comp = lapply(1:H,function(h) which(Zeta_vec==h))
  return(Zeta_comp)
}

uncomp_Z_subj = function(Z_subj_comp, nSubj, V, na_vox) {
  vec = rep(0,nSubj*V)
  vec[Z_subj_comp] = 1
  arr = array(vec,dim=c(nSubj,V))
  arr[,na_vox==T] = -1
  return(arr)
}

comp_Z_subj = function(Z_subj) {
  Z_subj_vec = c(Z_subj)
  Z_subj_comp = which(Z_subj_vec==1)
  return(Z_subj_comp)
}

# Initialize MCMC samples
initStorageVars_harm_long_R = function(Yvec, p, X_multi, Scanner_ID, Subj_ID, niter, R, H, subj_coef_type) {
  # subj_coef_type = 1: Lowrank
  # subj_coef_type = 2: Equal voxels
  # subj_coef_type = 3: Spike and slab
  
  #All parameters are lists with niter observations. The iter^th iteration of each has the following sizes:
  
  # population intercept margin:
  # gamma_pop ~ List of length D, each element is matrix of size p(d) x R(0)
  # w_pop / lambda_pop / alpha_pop ~ Matrix of size D x R(0)
  # tau_pop ~ Scalar
  
  #scanner intercept margins:
  # gamma_scan ~ List of length nScanners, each element is list of length D, each element is matrix of size p(d) x R(1)
  # w_scan / lambda_scan / alpha_scan ~ Matrix of size D x R(1)
  # tau_scan ~ Vector of length nScanners
  
  # covariate effects margins:
  # gamma_cov ~ List of length Q, each element is list of length D, each element is matrix of size p(d) x R(2)
  # w_cov / lambda_cov / alpha_cov ~ Matrix of size D x R(2)
  # tau_cov ~ Vector of length Q
  
  # subject intercept margins:
  # if subj_coef_type = 1:
  # gamma_subj ~ List of length nSubj, each element is list of length D, each element is matrix of size p(d) x R(3)
  # w_subj / lambda_subj / alpha_subj ~ Matrix of size D x R(3)
  # tau_subj ~ Vector of length nSubj
  # Z_subj ~ NULL
  # rho_subj ~ NULL
  # if subj_coef_type = 2:
  # gamma_subj ~ Vector of length nSubj
  # w_subj / lambda_subj / alpha_subj ~ NULL
  # tau_subj ~ Scalar
  # Z_subj ~ NULL
  # rho_subj ~ NULL
  # if subj_coef_type = 3:
  # gamma_subj ~ Vector of length nSubj
  # w_subj / lambda_subj / alpha_subj ~ NULL
  # tau_subj ~ Scalar
  # Z_subj ~ Matrix of size nSubj x V
  # rho_subj ~ Vector of length nSubj
  
  # residual variance:
  # Zeta ~ Matrix of size nobs x V
  # rho / ssq ~ Matrix of size nobs x H
  
  probs = rep(.5,2)
  nScanners = length(unique(Scanner_ID))
  nSubj = length(unique(Subj_ID))
  Q = ncol(X_multi)
  D = length(p)
  V = ncol(Yvec)
  probs_H = rep(1/H, H)
  if (length(R)!=4) {
    R = rep(R[1],4)
  }
  
  # Initialize lists
  gamma_pop_store = w_pop_store = lambda_pop_store = alpha_pop_store = tau_pop_store = vector(mode="list",niter+1)
  gamma_scan_store = w_scan_store = lambda_scan_store = alpha_scan_store = tau_scan_store = vector(mode="list",niter+1)
  gamma_cov_store = w_cov_store = lambda_cov_store = alpha_cov_store = tau_cov_store = vector(mode="list",niter+1)
  gamma_subj_store = w_subj_store = lambda_subj_store = alpha_subj_store = tau_subj_store = Z_subj_comp_store = rho_subj_store = vector(mode="list",niter+1)
  Zeta_comp_store = rho_store = ssq_store = vector(mode="list",niter+1)
  
  # Fill lists
  for (i in 1:(niter+1)) {
    
    # gamma's
    gamma_pop_store[[i]] = vector(mode="list",D)
    for (d in 1:D) {
      gamma_pop_store[[i]][[d]] = array(0,dim=c(p[d],R[1]))
      if (i==1) {
        for (r in 1:R[1]) {
          gamma_pop_store[[i]][[d]][,r] = rnorm(p[d],0,1)
        }
      }
    }
    
    gamma_scan_store[[i]] = vector(mode="list",nScanners)
    for (s in 1:nScanners) {
      gamma_scan_store[[i]][[s]] = vector(mode="list",D)
      for (d in 1:D) {
        gamma_scan_store[[i]][[s]][[d]] = array(0,dim=c(p[d],R[2]))
        if (i==1) {
          for (r in 1:R[2]) {
            gamma_scan_store[[i]][[s]][[d]][,r] = rnorm(p[d],0,1)
          }
        }
      }
    }
    
    gamma_cov_store[[i]] = vector(mode="list",Q)
    for (q in 1:Q) {
      gamma_cov_store[[i]][[q]] = vector(mode="list",D)
      for (d in 1:D) {
        gamma_cov_store[[i]][[q]][[d]] = array(0,dim=c(p[d],R[3]))
        if (i==1) {
          for (r in 1:R[3]) {
            gamma_cov_store[[i]][[q]][[d]][,r] = rnorm(p[d],0,1)
          }
        }
      }
    }
    
    if (subj_coef_type == 1) {
      gamma_subj_store[[i]] = vector(mode="list",nSubj)
      for (u in 1:nSubj) {
        gamma_subj_store[[i]][[u]] = vector(mode="list",D)
        for (d in 1:D) {
          gamma_subj_store[[i]][[u]][[d]] = array(0,dim=c(p[d],R[4]))
          if (i==1) {
            for (r in 1:R[4]) {
              gamma_subj_store[[i]][[u]][[d]][,r] = rnorm(p[d],0,1)
            }
          }
        }
      }
    } else {
      gamma_subj_store[[i]] = rep(0,nSubj)
      if (i==1) {
        gamma_subj_store[[i]] = rnorm(nSubj,0,1)
      }
    }
    
    # w/lambda/alpha's
    w_pop_store[[i]] = array(0,dim=c(D,R[1]))
    w_scan_store[[i]] = array(0,dim=c(D,R[2]))
    w_cov_store[[i]] = array(0,dim=c(D,R[3]))
    w_subj_store[[i]] = array(0,dim=c(D,R[4]))
    
    lambda_pop_store[[i]] = array(0,dim=c(D,R[1]))
    lambda_scan_store[[i]] = array(0,dim=c(D,R[2]))
    lambda_cov_store[[i]] = array(0,dim=c(D,R[3]))
    lambda_subj_store[[i]] = array(0,dim=c(D,R[4]))
    
    alpha_pop_store[[i]] = array(0,dim=c(D,R[1]))
    alpha_scan_store[[i]] = array(0,dim=c(D,R[2]))
    alpha_cov_store[[i]] = array(0,dim=c(D,R[3]))
    alpha_subj_store[[i]] = array(0,dim=c(D,R[4]))
    
    if (i==1) {
      w_pop_store[[i]] = w_pop_store[[i]] + 1
      w_scan_store[[i]] = w_scan_store[[i]] + 1
      w_cov_store[[i]] = w_cov_store[[i]] + 1
      w_subj_store[[i]] = w_subj_store[[i]] + 1
      
      lambda_pop_store[[i]] = lambda_pop_store[[i]] + 1
      lambda_scan_store[[i]] = lambda_scan_store[[i]] + 1
      lambda_cov_store[[i]] = lambda_cov_store[[i]] + 1
      lambda_subj_store[[i]] = lambda_subj_store[[i]] + 1
      
      alpha_pop_store[[i]] = alpha_pop_store[[i]] + 1
      alpha_scan_store[[i]] = alpha_scan_store[[i]] + 1
      alpha_cov_store[[i]] = alpha_cov_store[[i]] + 1
      alpha_subj_store[[i]] = alpha_subj_store[[i]] + 1
    }
    
    # tau's
    tau_pop_store[[i]] = 0
    tau_scan_store[[i]] = rep(0,nScanners)
    tau_cov_store[[i]] = rep(0,Q)
    if (subj_coef_type==1) {
      tau_subj_store[[i]] = rep(0,nSubj)
    } else {
      tau_subj_store[[i]] = 0
    }
    if (i==1) {
      tau_pop_store[[i]] = tau_pop_store[[i]]+1
      tau_scan_store[[i]] = tau_scan_store[[i]]+1
      tau_cov_store[[i]] = tau_cov_store[[i]]+1
      tau_subj_store[[i]] = tau_subj_store[[i]]+1
    }
    
    # Z_subj_comp & rho (spike-and-slab)
    if (subj_coef_type == 3) {
      rho_subj_store[[i]] = rep(0,nSubj)
      if (i==1) {
        rho_subj_store[[i]] = rep(0.5,nSubj)
        Z_subj_i = array(0,dim=c(nSubj,V))
        for (u in 1:nSubj) {
          for (v in 1:V) {
            Z_subj_i[u,v] = sample_1int_rcpp(probs)-1
          }
        }
        Z_subj_comp_store[[i]] = comp_Z_subj(Z_subj_i)
      }
    }
    
    # Zeta
    if (i==1) {
      Zeta_i = array(0,dim=c(nScanners,V))
      for (n in 1:nScanners) {
        for (v in 1:V) {
          Zeta_i[n,v] = sample_1int_rcpp(probs_H)
        }
      }
      Zeta_comp_store[[i]] = comp_Zeta(Zeta_i,H)
    }
    
    # rho & ssq
    rho_store[[i]] = array(0,dim=c(nScanners,H))
    ssq_store[[i]] = array(0,dim=c(nScanners,H))
    if (i==1) {
      for (n in 1:nScanners) {
        rho_store[[i]][n,] = probs_H
        ssq_store[[i]][n,] = 1
      }
    }
  }
  
  # save all
  all_gammas = list(gamma_pop_store, gamma_scan_store, gamma_cov_store, gamma_subj_store)
  all_ws = list(w_pop_store, w_scan_store, w_cov_store, w_subj_store)
  all_lambdas = list(lambda_pop_store, lambda_scan_store, lambda_cov_store, lambda_subj_store)
  all_alphas = list(alpha_pop_store, alpha_scan_store, alpha_cov_store, alpha_subj_store)
  all_taus = list(tau_pop_store, tau_scan_store, tau_cov_store, tau_subj_store)
  allparams = list(all_gammas, all_ws, all_lambdas, all_alphas, all_taus,
                   Zeta_comp_store, rho_store, ssq_store, Z_subj_comp_store, rho_subj_store)
  
  return(allparams)
}


# MCMC sampler
runMCMC_harm_long_R <- function(Yvec, p, X_cov, Scanner_ID, Subj_ID=NULL, niter, R, H, subj_coef_type=3,
                                a_lam=1, b_lam=1, a_tau=1, b_tau=1, a_alpha=1, b_alpha=1, sigma_log_alpha=0.01, beta_rho=1, a_s=1, b_s=1,
                                prog_count=1, show_all_steps=F, null_pop=F, null_scan=F, null_subj=F, null_cov=F,
                                initParams=NULL, condense_output=F, burn=.3) {
  
  # a_lam = b_lam = a_tau = b_tau = a_alpha = b_alpha = 1
  # sigma_log_alpha = .01
  # beta_rho = a_s = b_s = 1
  
  
  ## If voxel index (v1, v2, v3) ~ NA, then do not need to have an element \gamma_{v1}, \gamma_{v2}, \gamma_{v3}
  nobs = nrow(Yvec)
  na_vox = apply(Yvec,2,anyNA)
  Yarray = array(Yvec,dim=c(nobs,p))
  na_vox_all = apply(Yvec,2,function(x) all(is.na(x)))
  
  #na_vox_array_ind = which(na_vox_all==T)
  
  ## for d = 1, ..., D, have a array: elements of d-th margin x matrix of remaining 
  
  ## Presave indices of each tensor slice
  sliceInds_alldims = lapply(1:length(p), function(dd) t(sapply(1:p[dd],function(jj) getSliceIndices_iota_cpp(p,dd-1,jj-1))))
  
  ## For each element, find which remaining elements are missing
  sliceInds_alldim_exclude_na = lapply(1:length(p), function(d) {
    Dlist_allj = lapply(1:p[d], function(j) {
      na_dj = slice_tensor_cpp(na_vox_all, p, d-1, j-1);
      return(sliceInds_alldims[[d]][j,na_dj==F])
    })
    return(Dlist_allj)
  })
  
  a_rho = b_rho = 1
  
  if (length(R)!=4) {
    R = rep(R[1],4)
  }
  nobs = nrow(Yvec)
  if (is.null(Subj_ID)) Subj_ID = 1:nobs
  if (!is.null(X_cov) && !is.matrix(X_cov)) X_cov = as.matrix(X_cov)
  Q = ncol(X_cov)
  nScanners = length(unique(Scanner_ID)) #max(Scanner_ID)
  nSubj = length(unique(Subj_ID)) #max(Subj_ID)
  D = length(p)
  X_intercept = rep(1,nobs)
  allparams_store = initStorageVars_harm_long_R(Yvec,p,X_cov,Scanner_ID,Subj_ID,niter,R,H,subj_coef_type)
  if (!is.null(initParams)) {
    # gamma's
    allparams_store[[1]][[1]][[1]] = initParams$gamma_pop
    allparams_store[[1]][[2]][[1]] = initParams$gamma_scan
    allparams_store[[1]][[3]][[1]] = initParams$gamma_cov
    allparams_store[[1]][[4]][[1]] = initParams$gamma_subj
    allparams_store[[2]][[1]][[1]] = initParams$w_pop
    allparams_store[[2]][[2]][[1]] = initParams$w_scan
    allparams_store[[2]][[3]][[1]] = initParams$w_cov
    if (null_subj==F) {allparams_store[[2]][[4]][[1]] = initParams$w_subj}
    allparams_store[[3]][[1]][[1]] = initParams$lambda_pop
    allparams_store[[3]][[2]][[1]] = initParams$lambda_scan
    allparams_store[[3]][[3]][[1]] = initParams$lambda_cov
    if (null_subj==F) {allparams_store[[3]][[4]][[1]] = initParams$lambda_subj}
    allparams_store[[4]][[1]][[1]] = initParams$alpha_pop
    allparams_store[[4]][[2]][[1]] = initParams$alpha_scan
    allparams_store[[4]][[3]][[1]] = initParams$alpha_cov
    if (null_subj==F) {allparams_store[[4]][[4]][[1]] = initParams$alpha_subj}
    allparams_store[[5]][[1]][[1]] = initParams$tau_pop
    allparams_store[[5]][[2]][[1]] = initParams$tau_scan
    allparams_store[[5]][[3]][[1]] = initParams$tau_cov
    allparams_store[[5]][[4]][[1]] = initParams$tau_subj
    #allparams_store[[6]][[1]] = initParams$Zeta
    allparams_store[[6]][[1]] = initParams$Zeta_comp
    allparams_store[[7]][[1]] = initParams$rho
    allparams_store[[8]][[1]] = initParams$ssq
    #allparams_store[[9]][[1]] = initParams$Z_subj
    allparams_store[[9]][[1]] = initParams$Z_subj_comp
    allparams_store[[10]][[1]] = initParams$rho_subj
  }
  
  gamma_pop_store = allparams_store[[1]][[1]]
  w_pop_store = allparams_store[[2]][[1]]
  lambda_pop_store = allparams_store[[3]][[1]]
  alpha_pop_store = allparams_store[[4]][[1]]
  tau_pop_store = allparams_store[[5]][[1]]
  
  gamma_scan_store = allparams_store[[1]][[2]]
  w_scan_store = allparams_store[[2]][[2]]
  lambda_scan_store = allparams_store[[3]][[2]]
  alpha_scan_store = allparams_store[[4]][[2]]
  tau_scan_store = allparams_store[[5]][[2]]
  
  gamma_cov_store = allparams_store[[1]][[3]]
  w_cov_store = allparams_store[[2]][[3]]
  lambda_cov_store = allparams_store[[3]][[3]]
  alpha_cov_store = allparams_store[[4]][[3]]
  tau_cov_store = allparams_store[[5]][[3]]
  
  gamma_subj_store = allparams_store[[1]][[4]]
  w_subj_store = allparams_store[[2]][[4]]
  lambda_subj_store = allparams_store[[3]][[4]]
  alpha_subj_store = allparams_store[[4]][[4]]
  tau_subj_store = allparams_store[[5]][[4]]
  #Z_subj_store = allparams_store[[9]]
  Z_subj_comp_store = allparams_store[[9]]
  rho_subj_store = allparams_store[[10]]
  
  #Zeta_store = allparams_store[[6]]
  Zeta_comp_store = allparams_store[[6]]
  rho_store = allparams_store[[7]]
  ssq_store = allparams_store[[8]]
  
  # Get initial coefficients
  if (null_pop) {
    Gamma_pop_iter = rep(0,prod(p))
  } else {
    Gamma_pop_iter = getGamma_cpp(gamma_pop_store[[1]])
  }
  
  if (null_scan) {
    Gamma_scan_iter = vector(mode="list",nScanners)
    for (s in 1:nScanners) {
      Gamma_scan_iter[[s]] = rep(0,prod(p))
    }
  } else {
    Gamma_scan_iter = getGamma_multi_cpp(gamma_scan_store[[1]])
  }
  
  if (Q==0) {
    Q = 1
    null_cov = T
  }
  
  if (null_cov) {
    Gamma_cov_iter = vector(mode="list",Q)
    for (q in 1:Q) {
      Gamma_cov_iter[[q]] = rep(0,prod(p))
    }
  } else {
    Gamma_cov_iter = getGamma_multi_cpp(gamma_cov_store[[1]])
  }
  
  if (null_subj) {
    Gamma_subj_iter = vector(mode="list",nSubj)
    for (u in 1:nSubj) {
      Gamma_subj_iter[[u]] = rep(0,prod(p))
    }
  } else {
    if (subj_coef_type==1) {
      Gamma_subj_iter = getGamma_multi_cpp(gamma_subj_store[[1]])
    } else if (subj_coef_type==2) {
      Gamma_subj_iter = lapply(1:nSubj, function(u) rep(gamma_subj_store[[1]][u],prod(p)))
    } else if (subj_coef_type==3) {
      Z_subj_iter = uncomp_Z_subj(Z_subj_comp_store[[1]], nSubj, prod(p), na_vox_all)
      Gamma_subj_iter = lapply(1:nSubj, function(j) getCoef_j_spsl_cpp(gamma_subj_store[[1]],Z_subj_iter,j-1))
    }
  }
  
  for (iter in 1:niter) {
    if (iter %% prog_count == 0) {
      print(paste0('Iteration: ', iter))
    }
    
    gamma_pop_iter = gamma_pop_store[[iter]]
    w_pop_iter = w_pop_store[[iter]]
    lambda_pop_iter = lambda_pop_store[[iter]]
    alpha_pop_iter = alpha_pop_store[[iter]]
    tau_pop_iter = tau_pop_store[[iter]]
    
    gamma_scan_iter = gamma_scan_store[[iter]]
    w_scan_iter = w_scan_store[[iter]]
    lambda_scan_iter = lambda_scan_store[[iter]]
    alpha_scan_iter = alpha_scan_store[[iter]]
    tau_scan_iter = tau_scan_store[[iter]]
    
    gamma_cov_iter = gamma_cov_store[[iter]]
    w_cov_iter = w_cov_store[[iter]]
    lambda_cov_iter = lambda_cov_store[[iter]]
    alpha_cov_iter = alpha_cov_store[[iter]]
    tau_cov_iter = tau_cov_store[[iter]]
    
    gamma_subj_iter = gamma_subj_store[[iter]]
    tau_subj_iter = tau_subj_store[[iter]]
    if (subj_coef_type==1) {
      w_subj_iter = w_subj_store[[iter]]
      lambda_subj_iter = lambda_subj_store[[iter]]
      alpha_subj_iter = alpha_subj_store[[iter]]
    }
    if (subj_coef_type==3) {
      #Z_subj_iter = Z_subj_store[[iter]]
      Z_subj_comp_iter = Z_subj_comp_store[[iter]]
      Z_subj_iter = uncomp_Z_subj(Z_subj_comp_iter, nSubj, prod(p), na_vox_all)
      rho_subj_iter = rho_subj_store[[iter]]
    }
    
    #Zeta_iter = Zeta_store[[iter]]
    Zeta_comp_iter = Zeta_comp_store[[iter]]
    Zeta_iter = uncomp_Zeta(Zeta_comp_iter, nScanners, prod(p), H)
    rho_iter = rho_store[[iter]]
    ssq_iter = ssq_store[[iter]]
    SigSq_iter = getSigSq_cpp(ssq_iter, Zeta_iter)
    
    # population intercept
    if (!null_pop) {
      if (show_all_steps) {
        print(paste0(iter, ': Population intercept'))
      }
      Y_intercept = getY_intercept_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_cov_iter, Gamma_scan_iter, Gamma_subj_iter)
      for (r in 1:R[1]) {
        Yrvec_pop = getYr_cpp(Y_intercept, X_intercept, gamma_pop_iter, r-1)
        for (d in 1:D) {
          gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_withInds2_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1, sliceInds_alldim_exclude_na[[d]])
          #gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1, sliceInds_alldims[[d]])
          #gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_parallel_R(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d, r, sliceInds_alldims[[d]])
          #gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1)
          w_pop_iter[d,r] = sample_w_dr_cpp(gamma_pop_iter, tau_pop_iter, alpha_pop_iter, lambda_pop_iter, p, d-1, r-1)
          lambda_pop_iter[d,r] = sample_lambda_dr_cpp(w_pop_iter, p, a_lam, b_lam, d-1, r-1)
          alpha_pop_iter[d,r] = sample_alpha_dr_cpp(alpha_pop_iter, gamma_pop_iter, tau_pop_iter, w_pop_iter, p, a_alpha, b_alpha, sigma_log_alpha, d-1, r-1)
        }
      }
      tau_pop_iter = sample_tau_cpp(gamma_pop_iter, w_pop_iter, alpha_pop_iter, a_tau, b_tau)
      Gamma_pop_iter = getGamma_cpp(gamma_pop_iter)
    }
    
    # scanner intercepts
    if (!null_scan) {
      if (show_all_steps) {
        print(paste0(iter, ': Scanner intercept'))
      }
      for (s in 1:nScanners) {
        Y_scan_s = getY_scan_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_pop_iter, Gamma_cov_iter, Gamma_subj_iter, s)
        X_s = rep(1,nrow(Y_scan_s))
        Scanner_s = Scanner_ID[Scanner_ID==s]
        for (r in 1:R[2]) {
          Yrsvec = getYr_cpp(Y_scan_s, X_s, gamma_scan_iter[[s]], r-1)
          for (d in 1:D) {
            gamma_scan_iter[[s]][[d]][,r] = sample_gamma_dr_harm_withInds2_cpp(Yrsvec, p, Scanner_s, X_s, gamma_scan_iter[[s]], tau_scan_iter[s], w_scan_iter, alpha_scan_iter, SigSq_iter, d-1, r-1,sliceInds_alldim_exclude_na[[d]])
            #gamma_scan_iter[[s]][[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yrsvec, p, Scanner_s, X_s, gamma_scan_iter[[s]], tau_scan_iter[s], w_scan_iter, alpha_scan_iter, SigSq_iter, d-1, r-1,sliceInds_alldims[[d]])
            #gamma_scan_iter[[s]][[d]][,r] = sample_gamma_dr_harm_parallel_R(Yrsvec, p, Scanner_s, X_s, gamma_scan_iter[[s]], tau_scan_iter[s], w_scan_iter, alpha_scan_iter, SigSq_iter, d, r,sliceInds_alldims[[d]])
            #gamma_scan_iter[[s]][[d]][,r] = sample_gamma_dr_harm_cpp(Yrsvec, p, Scanner_s, X_s, gamma_scan_iter[[s]], tau_scan_iter[s], w_scan_iter, alpha_scan_iter, SigSq_iter, d-1, r-1)
          }
        }
      }
      for (r in 1:R[2]) {
        for (d in 1:D) {
          w_scan_iter[d,r] = sample_w_dr_multi_cpp(gamma_scan_iter, tau_scan_iter, alpha_scan_iter, lambda_scan_iter, p, d-1, r-1)
          lambda_scan_iter[d,r] = sample_lambda_dr_cpp(w_scan_iter, p, a_lam, b_lam, d-1, r-1)
          alpha_scan_iter[d,r] = sample_alpha_dr_multi_cpp(alpha_scan_iter, gamma_scan_iter, tau_scan_iter, w_scan_iter, p, a_alpha,b_alpha, sigma_log_alpha, d-1, r-1)
        }
      }
      tau_scan_iter = sample_tau_multi_cpp(gamma_scan_iter, w_scan_iter, alpha_scan_iter, a_tau, b_tau)
      Gamma_scan_iter = getGamma_multi_cpp(gamma_scan_iter)
    }
    
    # covariate effects
    if (!null_cov) {
      if (show_all_steps) {
        print(paste0(iter,': Covariate effects'))
      }
      Y_cov = getY_cov_long_cpp(Yvec, Scanner_ID, Subj_ID, Gamma_pop_iter, Gamma_scan_iter, Gamma_subj_iter)
      for (q in 1:Q) {
        if (q > 1) {
          Gamma_cov_iter[[q-1]] = getGamma_cpp(gamma_cov_iter[[q-1]])
        }
        Y_cov_q = getYq_multi_cpp(Y_cov, X_cov, Gamma_cov_iter, q-1)
        for (r in 1:R[3]) {
          Yrqvec = getYr_cpp(Y_cov_q, X_cov[,q], gamma_cov_iter[[q]], r-1)
          for (d in 1:D) {
            gamma_cov_iter[[q]][[d]][,r] = sample_gamma_dr_harm_withInds2_cpp(Yrqvec, p, Scanner_ID, X_cov[,q], gamma_cov_iter[[q]], tau_cov_iter[q], w_cov_iter, alpha_cov_iter, SigSq_iter, d-1, r-1,sliceInds_alldim_exclude_na[[d]])
            #gamma_cov_iter[[q]][[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yrqvec, p, Scanner_ID, X_cov[,q], gamma_cov_iter[[q]], tau_cov_iter[q], w_cov_iter, alpha_cov_iter, SigSq_iter, d-1, r-1,sliceInds_alldims[[d]])
            #gamma_cov_iter[[q]][[d]][,r] = sample_gamma_dr_harm_parallel_R(Yrqvec, p, Scanner_ID, X_cov[,q], gamma_cov_iter[[q]], tau_cov_iter[q], w_cov_iter, alpha_cov_iter, SigSq_iter, d, r,sliceInds_alldims[[d]])
            #gamma_cov_iter[[q]][[d]][,r] = sample_gamma_dr_harm_cpp(Yrqvec, p, Scanner_ID, X_cov[,q], gamma_cov_iter[[q]], tau_cov_iter[q], w_cov_iter, alpha_cov_iter, SigSq_iter, d-1, r-1)
          }
        }
      }
      for (r in 1:R[3]) {
        for (d in 1:D) {
          w_cov_iter[d,r] = sample_w_dr_multi_cpp(gamma_cov_iter, tau_cov_iter, alpha_cov_iter, lambda_cov_iter, p, d-1, r-1)
          lambda_cov_iter[d,r] = sample_lambda_dr_cpp(w_cov_iter, p, a_lam, b_lam, d-1, r-1)
          alpha_cov_iter[d,r] = sample_alpha_dr_multi_cpp(alpha_cov_iter, gamma_cov_iter, tau_cov_iter, w_cov_iter, p, a_alpha, b_alpha, sigma_log_alpha, d-1, r-1)
        }
      }
      tau_cov_iter = sample_tau_multi_cpp(gamma_cov_iter, w_cov_iter, alpha_cov_iter, a_tau, b_tau)
      Gamma_cov_iter = getGamma_multi_cpp(gamma_cov_iter)
    }
    
    # subject intercept
    if (!null_subj) {
      if (show_all_steps) {
        print(paste0(iter, ': Subject intercept'))
      }
      if (subj_coef_type==1) {
        for (u in 1:nSubj) {
          Y_subj_u = getY_subj_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_pop_iter, Gamma_cov_iter, Gamma_scan_iter, u)
          X_u = rep(1,nrow(Y_subj_u))
          Scanner_u = Scanner_ID[Subj_ID==u]
          for (r in 1:R[4]) {
            Yruvec = getYr_cpp(Y_subj_u, X_u, gamma_subj_iter[[u]], r-1)
            for (d in 1:D) {
              gamma_subj_iter[[u]][[d]][,r] = sample_gamma_dr_harm_withInds2_cpp(Yruvec, p, Scanner_u, X_u, gamma_subj_iter[[u]], tau_subj_iter[u], w_subj_iter, alpha_subj_iter, SigSq_iter, d-1, r-1,sliceInds_alldim_exclude_na[[d]])
              #gamma_subj_iter[[u]][[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yruvec, p, Scanner_u, X_u, gamma_subj_iter[[u]], tau_subj_iter[u], w_subj_iter, alpha_subj_iter, SigSq_iter, d-1, r-1,sliceInds_alldims[[d]])
              #gamma_subj_iter[[u]][[d]][,r] = sample_gamma_dr_harm_parallel_R(Yruvec, p, Scanner_u, X_u, gamma_subj_iter[[u]], tau_subj_iter[u], w_subj_iter, alpha_subj_iter, SigSq_iter, d, r,sliceInds_alldims[[d]])
              #gamma_subj_iter[[u]][[d]][,r] = sample_gamma_dr_harm_cpp(Yruvec, p, Scanner_u, X_u, gamma_subj_iter[[u]], tau_subj_iter[u], w_subj_iter, alpha_subj_iter, SigSq_iter, d-1, r-1)
            }
          }
        }
        for (r in 1:R[4]) {
          for (d in 1:D) {
            w_subj_iter[d,r] = sample_w_dr_multi_cpp(gamma_subj_iter, tau_subj_iter, alpha_subj_iter, lambda_subj_iter, p, d-1, r-1)
            lambda_subj_iter[d,r] = sample_lambda_dr_cpp(w_subj_iter, p, a_lam, b_lam, d-1, r-1)
            alpha_subj_iter[d,r] = sample_alpha_dr_multi_cpp(alpha_subj_iter, gamma_subj_iter, tau_subj_iter, w_subj_iter, p, a_alpha, b_alpha, sigma_log_alpha, d-1, r-1)
          }
        }
        tau_subj_iter = sample_tau_multi_cpp(gamma_subj_iter, w_subj_iter, alpha_subj_iter, a_tau, b_tau)
        Gamma_subj_iter = getGamma_multi_cpp(gamma_subj_iter)
        #### *** Equal subject intercept across voxels
      } else if (subj_coef_type==2) {
        for (u in 1:nSubj) {
          Y_subj_u = getY_subj_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_pop_iter, Gamma_cov_iter, Gamma_scan_iter, u)
          Scanner_u = Scanner_ID[Subj_ID==u]
          gamma_subj_iter[u] = sample_gamma_u_eqvox_harm_cpp(Y_subj_u,Scanner_u,tau_subj_iter,SigSq_iter)
          Gamma_subj_iter[[u]] = rep(gamma_subj_iter[u],prod(p))
        }
        tau_subj_iter = sample_tau_eqvox_cpp(gamma_subj_iter, a_tau, b_tau)
      } else if (subj_coef_type==3) {
        for (u in 1:nSubj) {
          Y_subj_u = getY_subj_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_pop_iter, Gamma_cov_iter, Gamma_scan_iter, u)
          Scanner_u = Scanner_ID[Subj_ID==u]
          Z_subj_iter[u,] = sample_Zj_cpp(Y_subj_u, Scanner_u, gamma_subj_iter, SigSq_iter, rho_subj_iter, na_vox, u-1)
          rho_subj_iter[u] = sample_rhoj_cpp(Z_subj_iter,a_rho,b_rho,u-1)
          gamma_subj_iter[u] = sample_betaj_cpp(Y_subj_u, Scanner_u, Z_subj_iter, SigSq_iter, tau_subj_iter, u-1)
          Gamma_subj_iter[[u]] = getCoef_j_spsl_cpp(gamma_subj_iter,Z_subj_iter,u-1)
        }
        tau_subj_iter = sample_tau_spsl_cpp(gamma_subj_iter,a_tau,b_tau)
      }
    }
    
    # residual noise
    if (show_all_steps) {
      print(paste0(iter, ': Residual noise'))
    }
    Residvec_iter = getResid_harm_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_pop_iter, Gamma_cov_iter, Gamma_scan_iter, Gamma_subj_iter)
    for (scan in 1:nScanners) {
      Zeta_iter[scan,] = sample_Zeta_scan_cpp(Residvec_iter, Scanner_ID, ssq_iter, rho_iter, scan)
      rho_iter[scan,] = sample_rho_scan_cpp(Zeta_iter, beta_rho, H, scan)
      ssq_iter[scan,] = sample_ssq_scan_cpp(Residvec_iter, Scanner_ID, Zeta_iter, a_s, b_s, H, scan)
    }
    
    # store sampled parameters
    gamma_pop_store[[iter+1]] = gamma_pop_iter
    w_pop_store[[iter+1]] = w_pop_iter
    lambda_pop_store[[iter+1]] = lambda_pop_iter
    alpha_pop_store[[iter+1]] = alpha_pop_iter
    tau_pop_store[[iter+1]] = tau_pop_iter
    
    gamma_scan_store[[iter+1]] = gamma_scan_iter
    w_scan_store[[iter+1]] = w_scan_iter
    lambda_scan_store[[iter+1]] = lambda_scan_iter
    alpha_scan_store[[iter+1]] = alpha_scan_iter
    tau_scan_store[[iter+1]] = tau_scan_iter
    
    gamma_cov_store[[iter+1]] = gamma_cov_iter
    w_cov_store[[iter+1]] = w_cov_iter
    lambda_cov_store[[iter+1]] = lambda_cov_iter
    alpha_cov_store[[iter+1]] = alpha_cov_iter
    tau_cov_store[[iter+1]] = tau_cov_iter
    
    gamma_subj_store[[iter+1]] = gamma_subj_iter
    tau_subj_store[[iter+1]] = tau_subj_iter
    if (subj_coef_type==1) {
      w_subj_store[[iter+1]] = w_subj_iter
      lambda_subj_store[[iter+1]] = lambda_subj_iter
      alpha_subj_store[[iter+1]] = alpha_subj_iter
    }
    if (subj_coef_type==3) {
      #Z_subj_store[[iter+1]] = Z_subj_iter
      Z_subj_comp_store[[iter+1]] = comp_Z_subj(Z_subj_iter)
      rho_subj_store[[iter+1]] = rho_subj_iter
    }
    
    #Zeta_store[[iter+1]] = Zeta_iter
    Zeta_comp_store[[iter+1]] = comp_Zeta(Zeta_iter, H)
    rho_store[[iter+1]] = rho_iter
    ssq_store[[iter+1]] = ssq_iter
  }
  
  if (subj_coef_type==1) {
    output_list = list(gamma_pop=gamma_pop_store,w_pop=w_pop_store,lambda_pop=lambda_pop_store,alpha_pop=alpha_pop_store,tau_pop=tau_pop_store,
                       gamma_scan=gamma_scan_store,w_scan=w_scan_store,lambda_scan=lambda_scan_store,alpha_scan=alpha_scan_store,tau_scan=tau_scan_store,
                       gamma_cov=gamma_cov_store,w_cov=w_cov_store,lambda_cov=lambda_cov_store,alpha_cov=alpha_cov_store,tau_cov=tau_cov_store,
                       gamma_subj=gamma_subj_store,w_subj=w_subj_store,lambda_subj=lambda_subj_store,alpha_subj=alpha_subj_store,tau_subj=tau_subj_store,
                       Zeta_comp=Zeta_comp_store,rho=rho_store,ssq=ssq_store)
  } else if (subj_coef_type==2) {
    output_list = list(gamma_pop=gamma_pop_store,w_pop=w_pop_store,lambda_pop=lambda_pop_store,alpha_pop=alpha_pop_store,tau_pop=tau_pop_store,
                       gamma_scan=gamma_scan_store,w_scan=w_scan_store,lambda_scan=lambda_scan_store,alpha_scan=alpha_scan_store,tau_scan=tau_scan_store,
                       gamma_cov=gamma_cov_store,w_cov=w_cov_store,lambda_cov=lambda_cov_store,alpha_cov=alpha_cov_store,tau_cov=tau_cov_store,
                       gamma_subj=gamma_subj_store,tau_subj=tau_subj_store,
                       Zeta_comp=Zeta_comp_store,rho=rho_store,ssq=ssq_store)
  } else if (subj_coef_type==3) {
    output_list = list(gamma_pop=gamma_pop_store,w_pop=w_pop_store,lambda_pop=lambda_pop_store,alpha_pop=alpha_pop_store,tau_pop=tau_pop_store,
                       gamma_scan=gamma_scan_store,w_scan=w_scan_store,lambda_scan=lambda_scan_store,alpha_scan=alpha_scan_store,tau_scan=tau_scan_store,
                       gamma_cov=gamma_cov_store,w_cov=w_cov_store,lambda_cov=lambda_cov_store,alpha_cov=alpha_cov_store,tau_cov=tau_cov_store,
                       gamma_subj=gamma_subj_store,tau_subj=tau_subj_store,Z_subj_comp=Z_subj_comp_store,rho_subj=rho_subj_store,
                       Zeta_comp=Zeta_comp_store,rho=rho_store,ssq=ssq_store)
  }
  return(output_list)
}

# more efficient coefficient extraction when many missing voxels are present
getCoef_mcmc_missingvox_R = function(gamma_store, p, missing_vox_vec, burn.in=.5,show.prog=T) {
  
  obs_vox_vec = which(missing_vox_vec==F,arr.ind=T)
  missing_vox_tens = array(missing_vox_vec,dim=p)
  obs_vox_tens = which(missing_vox_tens==F,arr.ind=T)
  
  niter = length(gamma_store)
  mcmc_ss = round(burn.in*niter):niter
  niter2 = length(mcmc_ss)
  gamma_iter1 = gamma_store[[1]]
  isMulti = is.list(gamma_iter1[[1]])
  
  if (isMulti==FALSE) {
    Gamma_mcmc = array(dim=c(niter2,prod(p)))
    for (iter in 1:niter2) {
      if (show.prog==T) print(mcmc_ss[iter])
      for (vv in 1:length(obs_vox_vec)) {
        vox_vec = obs_vox_vec[vv]
        vox_tens = obs_vox_tens[vv,]
        Gamma_mcmc[iter,vox_vec] = getLowRankCoefAtVox_cpp(gamma_store[[mcmc_ss[iter]]],vox_tens)
      }
    }
  } else {
    Q = length(gamma_iter1)
    Gamma_mcmc = array(dim=c(niter2,Q,prod(p)))
    for (iter in 1:niter2) {
      if (show.prog==T) print(mcmc_ss[iter])
      for (q in 1:Q) {
        for (vv in 1:length(obs_vox_vec)) {
          vox_vec = obs_vox_vec[vv]
          vox_tens = obs_vox_tens[vv,]
          Gamma_mcmc[iter,q,vox_vec] = getLowRankCoefAtVox_cpp(gamma_store[[mcmc_ss[iter]]][[q]],vox_tens)
        }
      }
    }
  }
  return(Gamma_mcmc)
}

getSigSq_harm_mcmc_R = function(Scanner_ID, ssq_store, Zeta_store, burn.in = .5, show.prog = T) {
  niter = length(ssq_store)
  H = ncol(ssq_store[[1]])
  nScanners = length(unique(Scanner_ID))
  V = dim(Zeta_store[[1]])[2]
  mcmc_ss = round(burn.in*niter):niter
  niter2 = length(mcmc_ss)
  SigSq_mcmc = array(dim=c(niter2,nScanners,V))
  
  for (iter in 1:niter2) {
    if (show.prog==T) print(mcmc_ss[iter])
    SigSq_mcmc[iter,,] = getSigSq_cpp(ssq_store[[mcmc_ss[iter]]], Zeta_store[[mcmc_ss[iter]]])
  }
  SigSq_mcmc[SigSq_mcmc<.000001] = NA
  return(SigSq_mcmc)
}

# more efficient residual noise extraction when many missing voxels are present
getSigSq_harm_mcmc_missingvox_R = function(Scanner_ID, ssq_store, Zeta_store, missing_vox_vec, burn.in=.5, show.prog=T) {
  niter = length(ssq_store)
  H = ncol(ssq_store[[1]])
  nScanners = length(unique(Scanner_ID))
  V = dim(Zeta_store[[1]])[2]
  mcmc_ss = round(burn.in*niter):niter
  niter2 = length(mcmc_ss)
  SigSq_mcmc = array(dim=c(niter2,nScanners,V))
  
  for (iter in 1:niter2) {
    if (show.prog==T) print(mcmc_ss[iter])
    for (v in 1:V) {
      if (missing_vox_vec[v]==F) {
        if (H>1) {
          SigSq_mcmc[iter,,v] = getSigSq_harm_atVox_cpp(Scanner_ID,ssq_store[[mcmc_ss[iter]]],Zeta_store[[mcmc_ss[iter]]],v)
        } else if (H==1) {
          SigSq_mcmc[iter,,v] = ssq_store[[mcmc_ss[iter]]]
        }
      }
    }
  }
  
  
  return(SigSq_mcmc)
}


# Competing Harmonization Methods & Applying Harmonization -----------------------------------------

# Competing harmonization methods -- Adjusted Residuals
adj_resid = function(Yvec, X, Scanner_ID, maxZeros=Inf) {
  V = ncol(Yvec)
  nScanners = length(unique(Scanner_ID))
  Q = ncol(X)
  Gamma_pop_est = rep(0,V)
  Gamma_scan_est = array(0,dim=c(nScanners,V))
  Gamma_cov_est = array(0,dim=c(Q,V))
  SigSq_est = array(1,dim=c(nScanners,V))
  allmissing = apply(Yvec,2,function(x) all(is.na(x)))
  nzeros = apply(Yvec,2,function(x) sum(x==0,na.rm=T))
  for (v in 1:V) {
    if (allmissing[v]==F && nzeros[v]<maxZeros) {
      lm_model = lmer(Yvec[,v]~X+(1|factor(Scanner_ID)))
      lm_model_summ = summary(lm_model)
      Gamma_pop_est[v] = lm_model_summ$coefficients[1,1]
      Gamma_scan_est[,v] = ranef(lm_model)$`factor(Scanner_ID)`[[1]]
      Gamma_cov_est[,v] = lm_model_summ$coefficients[-1,1]
      SigSq_est[,v] = lm_model_summ$sigma^2
    }
  }
  return(list(Gamma_pop_est=Gamma_pop_est,Gamma_scan_est=Gamma_scan_est, Gamma_cov_est=Gamma_cov_est, SigSq_est=SigSq_est))
}

adj_resid_long = function(Ysl, X, Scanner_ID, Subj_ID, maxZeros=Inf) {
  V = ncol(Ysl)
  nScanners = length(unique(Scanner_ID))
  nSubj = length(unique(Subj_ID))
  Q = ncol(X)
  Gamma_pop_est = rep(0,V)
  Gamma_scan_est = array(0,dim=c(nScanners,V))
  Gamma_subj_est = array(0,dim=c(nSubj,V))
  Gamma_cov_est = array(0,dim=c(Q,V))
  SigSq_est = array(1,dim=c(nScanners,V))
  allmissing = apply(Ysl,2,function(x) all(is.na(x)))
  nzeros = apply(Ysl,2,function(x) sum(x==0,na.rm=T))
  for (v in 1:V) {
    if (allmissing[v]==F && nzeros[v]<maxZeros) {
      lm_model = lmer(Ysl[,v]~X+(1|factor(Scanner_ID))+(1|factor(Subj_ID)))
      lm_model_summ = summary(lm_model)
      Gamma_pop_est[v] = lm_model_summ$coefficients[1,1]
      Gamma_scan_est[,v] = ranef(lm_model)$`factor(Scanner_ID)`[[1]]
      Gamma_subj_est[,v] = ranef(lm_model)$`factor(Subj_ID)`[[1]]
      Gamma_cov_est[,v] = lm_model_summ$coefficients[-1,1]
      SigSq_est[,v] = lm_model_summ$sigma^2
    }
  }
  return(list(Gamma_pop_est=Gamma_pop_est,Gamma_scan_est=Gamma_scan_est, Gamma_subj_est = Gamma_subj_est, Gamma_cov_est=Gamma_cov_est, SigSq_est=SigSq_est))
}

# harmonization functions
harmonize_cs = function(Yvec, Scanner_ID, Xq, M, Gamma_i, Theta_q=NULL, SigSq, cov_adj=T, thresh=.001) {
  N = nrow(Yvec)
  V = ncol(Yvec)
  Q = ncol(Xq)
  if (cov_adj == F) {
    Theta_q = array(0,dim=c(Q,V))
  }
  Yharm = 0*Yvec
  for (n in 1:N) {
    Yharm[n,] = (Yvec[n,]-M-Gamma_i[Scanner_ID[n],]-Xq[n,]%*%Theta_q)/sqrt(SigSq[Scanner_ID[n],]) + M + Xq[n,]%*%Theta_q
  }
  return(Yharm)
}


harmonize_long = function(Yvec, Scanner_ID, Xq, Subj_ID, M, Gamma_i, Theta_q, Beta_u, SigSq, cov_adj=T, thresh=.001, apply_missing=T) {
  if (apply_missing==T) {
    missing_vox = apply(SigSq,2,function(x) anyNA(x) || sd(x,na.rm=T)<thresh)
    M[missing_vox] = NA
    Gamma_i[,missing_vox] = NA
    Theta_q[,missing_vox] = NA
    Beta_u[,missing_vox] = NA
    SigSq[,missing_vox] = NA
  }
  N = nrow(Yvec)
  V = ncol(Yvec)
  Q = ncol(Xq)
  if (cov_adj == F) {
    Theta_q = array(0,dim=c(Q,V))
  }
  Yharm = 0*Yvec
  for (n in 1:N) {
    Yharm[n,] = (Yvec[n,]-M-Gamma_i[Scanner_ID[n],]-Xq[n,]%*%Theta_q-Beta_u[Subj_ID[n],])/sqrt(SigSq[Scanner_ID[n],]) + M + Xq[n,]%*%Theta_q + Beta_u[Subj_ID[n],]
  }
  return(Yharm)
}

harmonize_long_eqvox = function(Yvec, Scanner_ID, Xq, Subj_ID, M, Gamma_i, Theta_q, Beta_u, SigSq, cov_adj=T, thresh=.001, apply_missing=T) {
  if (apply_missing==T) {
    missing_vox = apply(SigSq,2,function(x) anyNA(x) || sd(x,na.rm=T)<thresh)
    M[missing_vox] = NA
    Gamma_i[,missing_vox] = NA
    Theta_q[,missing_vox] = NA
    #Beta_u[,missing_vox] = NA
    SigSq[,missing_vox] = NA
  }
  N = nrow(Yvec)
  V = ncol(Yvec)
  Q = ncol(Xq)
  if (cov_adj == F) {
    Theta_q = array(0,dim=c(Q,V))
  }
  Yharm = 0*Yvec
  for (n in 1:N) {
    Yharm[n,] = (Yvec[n,]-M-Gamma_i[Scanner_ID[n],]-Xq[n,]%*%Theta_q-Beta_u[Subj_ID[n]])/sqrt(SigSq[Scanner_ID[n],]) + M + Xq[n,]%*%Theta_q + Beta_u[Subj_ID[n]]
  }
  return(Yharm)
}