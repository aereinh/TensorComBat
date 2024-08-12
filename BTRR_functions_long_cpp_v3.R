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

sourceCpp('BTRR_harm_long_Rcpp_v3.cpp')


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

getCoef_mcmc_R2 = function(gamma_store) {
  niter = length(gamma_store)
  gamma_iter1 = gamma_store[[1]]
  isMulti = is.list(gamma_iter1[[1]])
  if (isMulti==FALSE) {
    Gamma_mcmc_cpp = getGammaUniv_mcmc_cpp(gamma_store)
    Gamma_mcmc = t(array(unlist(Gamma_mcmc_cpp),dim=c(length(Gamma_mcmc_cpp[[1]]),niter)))
  } else {
    Q = length(gamma_iter1)
    Gamma_mcmc_cpp = getGammaMulti_mcmc_cpp(gamma_store)
    Gamma_mcmc = aperm(array(unlist(Gamma_mcmc_cpp),dim=c(length(Gamma_mcmc_cpp[[1]][[1]]),Q,niter)),c(3,2,1))
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

sample_gamma_dr_harm_parallel_R = function(Yrvec, p, Scanner_ID, X, gamma, tau, w, alpha, SigSq, d, r, slice_inds_d) {
  nobs = nrow(Yrvec)
  D = length(p)
  gamma_d = gamma[[d]]
  gamma_dr = gamma_d[,r]
  #gamma_dr_new = rep(0,p[d])
  #c_dr = rep(0,p[d])
  #d_dr = rep(0,p[d])
  w_dr = w[d,r];
  alpha_dr = alpha[d,r];
  
  if (D==3) {
    gamma_x = gamma[[1]]
    gamma_xr = gamma_x[,r]
    gamma_y = gamma[[2]]
    gamma_yr = gamma_y[,r]
    gamma_z = gamma[[3]]
    gamma_zr = gamma_z[,r]
    if (d==1) {
      outer_dr = TP_2D_cpp(gamma_yr, gamma_zr)
    } else if (d==2) {
      outer_dr = TP_2D_cpp(gamma_xr, gamma_zr)
    } else if (d==3) {
      outer_dr = TP_2D_cpp(gamma_xr, gamma_yr)
    }
  } else if (D==2) {
    gamma_x = gamma[[1]]
    gamma_xr = gamma_x[,r]
    gamma_y = gamma[[2]]
    gamma_yr = gamma_y[,r]
    if (d==1) {
      outer_dr = gamma_yr
    } else if (d==2) {
      outer_dr = gamma_xr
    }
  }
  
  slice_inds_d = sliceInds_alldims[[d]]

  # parallel computation over margin elements
  gamma_dr_new = mclapply(1:p[d], function(vd) {
    slice_inds_dvd = slice_inds_d[vd,]+1
    c_drv = d_drv = 0
    
    # parallel computation over samples
    mn_all_vd = mclapply(1:nobs, function(i) {
      SigSq_sl_div = SigSq[Scanner_ID[i],slice_inds_dvd];
      Y_sl_driv = Yrvec[i,slice_inds_dvd];
      Ysl_isNA = is.na(Y_sl_driv);
    
      m_drv_i = sum(((X[i]*outer_dr)^2/SigSq_sl_div)[Ysl_isNA==F]);
      n_drv_i = sum(((X[i]*Y_sl_driv*outer_dr)/SigSq_sl_div)[Ysl_isNA==F]);
      return(c(m_drv_i,n_drv_i))
    }) %>% unlist %>% array(dim=c(2,nobs)) %>% t()
    
    m_drv = sum(mn_all_vd[,1])
    n_drv = sum(mn_all_vd[,2])
  
    if (vd==1) {
      c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)))
      d_drv = n_drv+(exp(-alpha_dr)*gamma_dr[vd+1])/(tau*w_dr*(1-exp(-2*alpha_dr)));
    } else if (vd==p[d]) {
      c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*gamma_dr[vd-1])/(tau*w_dr*(1-exp(-2*alpha_dr)));
    } else {
      c_drv = m_drv+(1+exp(-2*alpha_dr))/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*(gamma_dr[vd-1]+gamma_dr[vd+1]))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    }
  
    if (c_drv==0) {
      return(0);
    } else {
      return(rnorm(n = 1, mean = d_drv/c_drv, sd = sqrt(1.0/c_drv)))
    }
  }) %>% unlist()
  
  return(gamma_dr_new)
}

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
  Q = ncol(X_cov)
  nScanners = max(Scanner_ID)
  nSubj = max(Subj_ID)
  D = length(p)
  X_intercept = rep(1,nobs)
  allparams_store = initStorageVars_harm_long_cpp(Yvec,p,X_cov,Scanner_ID,Subj_ID,niter,R,H,subj_coef_type)
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
    allparams_store[[6]][[1]] = initParams$Zeta
    allparams_store[[7]][[1]] = initParams$rho
    allparams_store[[8]][[1]] = initParams$ssq
    allparams_store[[9]][[1]] = initParams$Z_subj
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
  Z_subj_store = allparams_store[[9]]
  rho_subj_store = allparams_store[[10]]
  
  Zeta_store = allparams_store[[6]]
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
      Gamma_subj_iter = lapply(1:nSubj, function(j) getCoef_j_spsl_cpp(gamma_subj_store[[1]],Z_subj_store[[1]],j-1))
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
      Z_subj_iter = Z_subj_store[[iter]]
      rho_subj_iter = rho_subj_store[[iter]]
    }
    
    Zeta_iter = Zeta_store[[iter]]
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
      Z_subj_store[[iter+1]] = Z_subj_iter
      rho_subj_store[[iter+1]] = rho_subj_iter
    }
    
    Zeta_store[[iter+1]] = Zeta_iter
    rho_store[[iter+1]] = rho_iter
    ssq_store[[iter+1]] = ssq_iter
  }
  
  if (subj_coef_type==1) {
    output_list = list(gamma_pop=gamma_pop_store,w_pop=w_pop_store,lambda_pop=lambda_pop_store,alpha_pop=alpha_pop_store,tau_pop=tau_pop_store,
                       gamma_scan=gamma_scan_store,w_scan=w_scan_store,lambda_scan=lambda_scan_store,alpha_scan=alpha_scan_store,tau_scan=tau_scan_store,
                       gamma_cov=gamma_cov_store,w_cov=w_cov_store,lambda_cov=lambda_cov_store,alpha_cov=alpha_cov_store,tau_cov=tau_cov_store,
                       gamma_subj=gamma_subj_store,w_subj=w_subj_store,lambda_subj=lambda_subj_store,alpha_subj=alpha_subj_store,tau_subj=tau_subj_store,
                       Zeta=Zeta_store,rho=rho_store,ssq=ssq_store)
  } else if (subj_coef_type==2) {
    output_list = list(gamma_pop=gamma_pop_store,w_pop=w_pop_store,lambda_pop=lambda_pop_store,alpha_pop=alpha_pop_store,tau_pop=tau_pop_store,
                       gamma_scan=gamma_scan_store,w_scan=w_scan_store,lambda_scan=lambda_scan_store,alpha_scan=alpha_scan_store,tau_scan=tau_scan_store,
                       gamma_cov=gamma_cov_store,w_cov=w_cov_store,lambda_cov=lambda_cov_store,alpha_cov=alpha_cov_store,tau_cov=tau_cov_store,
                       gamma_subj=gamma_subj_store,tau_subj=tau_subj_store,
                       Zeta=Zeta_store,rho=rho_store,ssq=ssq_store)
  } else if (subj_coef_type==3) {
    output_list = list(gamma_pop=gamma_pop_store,w_pop=w_pop_store,lambda_pop=lambda_pop_store,alpha_pop=alpha_pop_store,tau_pop=tau_pop_store,
                       gamma_scan=gamma_scan_store,w_scan=w_scan_store,lambda_scan=lambda_scan_store,alpha_scan=alpha_scan_store,tau_scan=tau_scan_store,
                       gamma_cov=gamma_cov_store,w_cov=w_cov_store,lambda_cov=lambda_cov_store,alpha_cov=alpha_cov_store,tau_cov=tau_cov_store,
                       gamma_subj=gamma_subj_store,tau_subj=tau_subj_store,Z_subj=Z_subj_store,rho_subj=rho_subj_store,
                       Zeta=Zeta_store,rho=rho_store,ssq=ssq_store)
  }
  return(output_list)
}

runMCMC_harm_long_R2 <- function(Yvec, p, X_cov, Scanner_ID, Subj_ID=NULL, niter, R, H, subj_coef_type=3,
                                a_lam=1, b_lam=1, a_tau=1, b_tau=1, a_alpha=1, b_alpha=1, sigma_log_alpha=0.01, beta_rho=1, a_s=1, b_s=1,
                                prog_count=1, show_all_steps=F, null_pop=F, null_scan=F, null_subj=F, null_cov=F,
                                initParams=NULL, condense_output=F, burn=.3) {
  
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
  Q = ncol(X_cov)
  nScanners = max(Scanner_ID)
  nSubj = max(Subj_ID)
  D = length(p)
  X_intercept = rep(1,nobs)
  allparams_store = initStorageVars_harm_long_cpp(Yvec,p,X_cov,Scanner_ID,Subj_ID,niter,R,H,subj_coef_type)
  
  #### ** FIX THIS TO ALLOW FOR OTHER INITIALIZATION (e.g. if we have partially run the MCMC)
  # if (!is.null(initParams)) {
  #   for (l in 1:length(allparams_store)) {
  #     allparams_store[[l]][[1]] = initParams[[l]]
  #   }
  # }
  
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
  Z_subj_store = allparams_store[[9]]
  rho_subj_store = allparams_store[[10]]
  
  Zeta_store = allparams_store[[6]]
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
      Gamma_subj_iter = lapply(1:nSubj, function(j) getCoef_j_spsl_cpp(gamma_subj_store[[1]],Z_subj_store[[1]],j-1))
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
      Z_subj_iter = Z_subj_store[[iter]]
      rho_subj_iter = rho_subj_store[[iter]]
    }
    
    Zeta_iter = Zeta_store[[iter]]
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
          #gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_withInds2_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1, sliceInds_alldim_exclude_na[[d]])
          gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1, sliceInds_alldims[[d]])
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
            #gamma_scan_iter[[s]][[d]][,r] = sample_gamma_dr_harm_withInds2_cpp(Yrsvec, p, Scanner_s, X_s, gamma_scan_iter[[s]], tau_scan_iter[s], w_scan_iter, alpha_scan_iter, SigSq_iter, d-1, r-1,sliceInds_alldim_exclude_na[[d]])
            gamma_scan_iter[[s]][[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yrsvec, p, Scanner_s, X_s, gamma_scan_iter[[s]], tau_scan_iter[s], w_scan_iter, alpha_scan_iter, SigSq_iter, d-1, r-1,sliceInds_alldims[[d]])
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
            #gamma_cov_iter[[q]][[d]][,r] = sample_gamma_dr_harm_withInds2_cpp(Yrqvec, p, Scanner_ID, X_cov[,q], gamma_cov_iter[[q]], tau_cov_iter[q], w_cov_iter, alpha_cov_iter, SigSq_iter, d-1, r-1,sliceInds_alldim_exclude_na[[d]])
            gamma_cov_iter[[q]][[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yrqvec, p, Scanner_ID, X_cov[,q], gamma_cov_iter[[q]], tau_cov_iter[q], w_cov_iter, alpha_cov_iter, SigSq_iter, d-1, r-1,sliceInds_alldims[[d]])
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
              #gamma_subj_iter[[u]][[d]][,r] = sample_gamma_dr_harm_withInds2_cpp(Yruvec, p, Scanner_u, X_u, gamma_subj_iter[[u]], tau_subj_iter[u], w_subj_iter, alpha_subj_iter, SigSq_iter, d-1, r-1,sliceInds_alldim_exclude_na[[d]])
              gamma_subj_iter[[u]][[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yruvec, p, Scanner_u, X_u, gamma_subj_iter[[u]], tau_subj_iter[u], w_subj_iter, alpha_subj_iter, SigSq_iter, d-1, r-1,sliceInds_alldims[[d]])
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
      Z_subj_store[[iter+1]] = Z_subj_iter
      rho_subj_store[[iter+1]] = rho_subj_iter
    }
    
    Zeta_store[[iter+1]] = Zeta_iter
    rho_store[[iter+1]] = rho_iter
    ssq_store[[iter+1]] = ssq_iter
  }
  
  
  if (subj_coef_type==1) {
    output_list = list(gamma_pop=gamma_pop_store,w_pop=w_pop_store,lambda_pop=lambda_pop_store,alpha_pop=alpha_pop_store,tau_pop=tau_pop_store,
                       gamma_scan=gamma_scan_store,w_scan=w_scan_store,lambda_scan=lambda_scan_store,alpha_scan=alpha_scan_store,tau_scan=tau_scan_store,
                       gamma_cov=gamma_cov_store,w_cov=w_cov_store,lambda_cov=lambda_cov_store,alpha_cov=alpha_cov_store,tau_cov=tau_cov_store,
                       gamma_subj=gamma_subj_store,w_subj=w_subj_store,lambda_subj=lambda_subj_store,alpha_subj=alpha_subj_store,tau_subj=tau_subj_store,
                       Zeta=Zeta_store,rho=rho_store,ssq=ssq_store)
  } else if (subj_coef_type==2) {
    output_list = list(gamma_pop=gamma_pop_store,w_pop=w_pop_store,lambda_pop=lambda_pop_store,alpha_pop=alpha_pop_store,tau_pop=tau_pop_store,
                       gamma_scan=gamma_scan_store,w_scan=w_scan_store,lambda_scan=lambda_scan_store,alpha_scan=alpha_scan_store,tau_scan=tau_scan_store,
                       gamma_cov=gamma_cov_store,w_cov=w_cov_store,lambda_cov=lambda_cov_store,alpha_cov=alpha_cov_store,tau_cov=tau_cov_store,
                       gamma_subj=gamma_subj_store,tau_subj=tau_subj_store,
                       Zeta=Zeta_store,rho=rho_store,ssq=ssq_store)
  } else if (subj_coef_type==3) {
    output_list = list(gamma_pop=gamma_pop_store,w_pop=w_pop_store,lambda_pop=lambda_pop_store,alpha_pop=alpha_pop_store,tau_pop=tau_pop_store,
                       gamma_scan=gamma_scan_store,w_scan=w_scan_store,lambda_scan=lambda_scan_store,alpha_scan=alpha_scan_store,tau_scan=tau_scan_store,
                       gamma_cov=gamma_cov_store,w_cov=w_cov_store,lambda_cov=lambda_cov_store,alpha_cov=alpha_cov_store,tau_cov=tau_cov_store,
                       gamma_subj=gamma_subj_store,tau_subj=tau_subj_store,Z_subj=Z_subj_store,rho_subj=rho_subj_store,
                       Zeta=Zeta_store,rho=rho_store,ssq=ssq_store)
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
adj_resid = function(Yvec, X, Scanner_ID, maxZeros=Inf, show.prog = F) {
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
    if (show.prog == T) print(v)
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

adj_resid_long = function(Ysl, X, Scanner_ID, Subj_ID, maxZeros=Inf, show.prog=F) {
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
    if (show.prog==T) print(v)
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


# Examining Site Effects --------------------------------------------------

getScannerMeans = function(Yharm_nona, Scanner_ID) {
  if (class(Scanner_ID)!="numeric") Scanner_ID = as.numeric(Scanner_ID)
  nScanners = length(unique(Scanner_ID))
  Yharm_nona_scanmean = array(dim=c(nScanners,ncol(Yharm_nona)))
  for (i in 1:nScanners) {
    print(i)
    Yharm_nona_scanmean[i,] = apply(Yharm_nona[Scanner_ID==i,],2,function(x) mean(x,na.rm=T))
  }
  return(Yharm_nona_scanmean)
}


partialCor = function(x,y) {
  ind = which(!is.na(x) & !is.na(y))
  return(cor(x[ind],y[ind]))
}

getPairwiseScannerMetrics = function(Yharm_nona_scanmean,Scanner_ID) {
  if (class(Scanner_ID)!="numeric") Scanner_ID = as.numeric(Scanner_ID)
  nScanners = length(unique(Scanner_ID))
  cormat = rmsemat = array(dim=c(nScanners,nScanners))
  for (i in 2:nScanners-1) {
    print(i)
    for (j in (i+1):nScanners) {
      rmsemat[j,i] = sqrt(mean((Yharm_nona_scanmean[i,]-Yharm_nona_scanmean[j,])^2,na.rm=T))
      cormat[j,i] = partialCor(Yharm_nona_scanmean[i,],Yharm_nona_scanmean[j,])
    }
  }
  return(list(cors=cormat,rmses=rmsemat))
}


getResidualsVxlWise = function(Ylist, X_cov, prog_count = 100) {
  if (!is.list(Ylist)) Ylist = list(Ylist)
  L = length(Ylist)
  N = nrow(Ylist[[1]])
  #V = ncol(Ylist[[1]])
  V_allmethods = unlist(lapply(Ylist,ncol))
  V = max(V_allmethods)
  Residlist = lapply(1:L,function(x) array(dim=c(N,V_allmethods[x])))
  
  for (v in 1:V) {
    if (v %% prog_count==0) print(paste0(signif(100*v/V,2),'% Complete'))
    for (l in 1:L) {
      if (V_allmethods[l]>=v) {
        yv = Ylist[[l]][,v]
        mod_v = lm(yv ~ X_cov)
        Residlist[[l]][,v] = mod_v$residuals
      }
    }
  }
  return(Residlist)
}

makeResidPlotsBySite = function(Residlist, Scanner_ID, Methods, Lq = 1000, type="boxplot") {
  L = length(Residlist)
  nScanners = length(unique(Scanner_ID))
  quants = seq(from=.001, to=.999, length.out=Lq)
  Residlist_ds = lapply(1:L, function(x) array(dim=c(nScanners,Lq)))
  for (i in 1:nScanners) {
    for (l in 1:L) {
      Residlist_ds[[l]][i,] = quantile(c(Residlist[[l]][Scanner_ID==i,]),quants)
    }
  }
  #scanner_order = order(apply(Residlist_ds[[1]],1,median))
  #Residlist_ds_ordered = lapply(Residlist_ds, function(x) x[scanner_order,])
  
  n_tot = prod(dim(Residlist_ds[[1]]))
  Method_reps = c()
  Resids_reps = c()
  Sites_reps = c()
  for (l in 1:L) {
    Method_reps = c(Method_reps, rep(Methods[l],n_tot))
    for (i in 1:nScanners) {
      Resids_reps = c(Resids_reps, Residlist_ds[[l]][i,])
      Sites_reps = c(Sites_reps, rep(i,Lq))
    }
  }
  
  plot_data = data.frame(Method = factor(Method_reps), Site = factor(Sites_reps), Resid = Resids_reps)
  if (type=="boxplot") {
    plot_data %>% 
      ggplot(aes(x=Site, y=Resid, fill=Method))+
      geom_boxplot()+
      facet_wrap(vars(Method))+
      theme_bw()+
      ylab('Residuals')
  } else {
    plot_data %>%
      ggplot(aes(x=Site, y=Resid, color=Method))+
      geom_jitter()+
      facet_wrap(vars(Method))+
      theme_bw()+
      ylab('Residuals')
  }
}

testSiteMeanAcrossVoxels = function(Ylist, Scanner_ID, prog_count=100) {
  if (!is.list(Ylist)) Ylist = list(Ylist)
  if (!is.factor(Scanner_ID)) Scanner_ID = factor(Scanner_ID)
  L = length(Ylist)
  V = ncol(Ylist[[1]])
  pvals_list = lapply(1:L, function(x) rep(NA,V))
  for (v in 1:V) {
    if (v %% prog_count==0) print(paste0(signif(100*v/V,3),'% Complete'))
    for (l in 1:L) {
      anova_lv = aov(Ylist[[l]][,v] ~ Scanner_ID)
      pvals_list[[l]][v] = summary(anova_lv)[[1]][1,5]
    }
  }
  return(pvals_list)
}

testSiteVarianceAcrossVoxels = function(Ylist, Scanner_ID, prog_count=100) {
  if (!is.list(Ylist)) Ylist = list(Ylist)
  if (!is.factor(Scanner_ID)) Scanner_ID = factor(Scanner_ID)
  L = length(Ylist)
  V = ncol(Ylist[[1]])
  pvals_list = lapply(1:L, function(x) rep(NA,V))
  for (v in 1:V) {
    if (v %% prog_count==0) print(paste0(signif(100*v/V,3),'% Complete'))
    for (l in 1:L) {
      barlett_lv = bartlett.test(formula = Ylist[[l]][,v] ~ Scanner_ID)
      pvals_list[[l]][v] = barlett_lv$p.value
    }
  }
  return(pvals_list)
}

makePValuePlot = function(pvals_list, Methods, type="boxplot") {
  L = length(pvals_list)
  n_tot = length(pvals_list[[1]])
  Method_reps = c()
  pval_reps = unlist(pvals_list)
  for (l in 1:L) {
    Method_reps = c(Method_reps, rep(Methods[l],n_tot))
  }
  plot_data = data.frame(Method = factor(Method_reps), p = pval_reps)
  
  if (type=="boxplot") {
    plot_data %>%
      ggplot(aes(x = Method, y = p, fill=Method))+
      geom_boxplot()+
      geom_hline(yintercept = 0.05, linetype=2)+
      theme_bw()+
      ylab('p-value')
  } else {
    plot_data %>%
      ggplot(aes(x=Method, y=p, color=Method))+
      geom_jitter()+
      geom_hline(yintercept = 0.05, linetype=2)+
      theme_bw()+
      ylab('p-value')
  }
}

# Biological Prediction ---------------------------------------------------

get_kfold_inds = function(Scanner_ID,kfold) {
  kfold_inds = rep(NA,length(Scanner_ID))
  for (scan in unique(Scanner_ID)) {
    scan_inds = which(Scanner_ID==scan)
    kfold_inds[scan_inds] = sample(1:kfold,length(scan_inds),replace=T)
  }
  return(kfold_inds)
}

get_kfold_inds2 = function(Scanner_ID,Subj_ID,kfold) {
  kfold_inds = rep(NA,length(Scanner_ID))
  for (scan in unique(Scanner_ID)) {
    scan_inds = which(Scanner_ID==scan)
    for (subj in unique(Subj_ID[scan_inds])) {
      subj_scan_inds = which(Scanner_ID==scan & Subj_ID==subj)
      kfold_inds[subj_scan_inds] = sample(1:kfold,length(subj_scan_inds),replace=T)
    }
  }
  return(kfold_inds)
}

SVMgridSearchRMSE = function(xtrain, ytrain, xtest, ytest, kernel="radial", ep, C) {
  rmse = array(dim=c(length(ep),length(C)))
  for (e in ep) {
    for (cc in C) {
      mod = svm(x=xtrain,y=ytrain,kernel=kernel,cost=cc,epsilon=e)
      rmse[which(ep==e),which(C==cc)] = sqrt(mean((as.numeric(ytest) - as.numeric(predict(mod, xtest)))^2))
    }
  }
  return(min(rmse))
}

catPrediction <- function(imgdata, cat_cov, Scanner_ID, Subj_ID=NULL, nreps = 10, kfolds = 3, seed = 123, show.prog=T) {
  catcov_acc_reps = array(dim=c(nreps,kfolds))
  cat_cov = as.factor(cat_cov)
  set.seed(seed)
  for (i in 1:nreps) {
    if (show.prog==T) print(i)
    if (is.null(Subj_ID)) {
      k_inds = get_kfold_inds(Scanner_ID,kfolds)
    } else {
      k_inds = get_kfold_inds2(Scanner_ID,Subj_ID,kfolds)
    }
    for (k in 1:kfolds) {
      if (sum(k_inds==k)>=2) {
        svm_catcov = svm(x = imgdata[k_inds!=k,], y=cat_cov[k_inds!=k])
        catcov_acc_reps[i,k] = sum(as.numeric(predict(svm_catcov,imgdata[k_inds==k,]))==as.numeric(cat_cov[k_inds==k]))/length(cat_cov[k_inds==k])
      }
    }
  }
  return(catcov_acc_reps)
}

contPrediction <- function(imgdata, cont_cov, Scanner_ID, Subj_ID=NULL, nreps = 10, kfolds = 3, ep = .1, C = 1, seed = 123, show.prog = T) {
  contcov_rmse_reps = array(dim=c(nreps,kfolds))
  set.seed(seed)
  for (i in 1:nreps) {
    if (show.prog==T) print(i)
    if (is.null(Subj_ID)) {
      k_inds = get_kfold_inds(Scanner_ID,kfolds)
    } else {
      k_inds = get_kfold_inds2(Scanner_ID,Subj_ID,kfolds)
    }
    for (k in 1:kfolds) {
      if (sum(k_inds==k)>=2) {
        contcov_rmse_reps[i,k] = SVMgridSearchRMSE(imgdata[k_inds!=k,],cont_cov[k_inds!=k],imgdata[k_inds==k,],cont_cov[k_inds==k],kernel="radial",ep,C)
      }
    }
  }
  return(contcov_rmse_reps)
}

lassoPrediction_cont <- function(imgdata, cont_cov, Scan, Subj=NULL, nreps = 10, kfolds = 3, seed=123, show.prog = T) {
  # remove any NAs
  if (anyNA(cont_cov)) {
    imgdata = imgdata[!is.na(cont_cov),]
    Scan = Scan[!is.na(cont_cov)]
    if (!is.null(Subj)) Subj = Subj[!is.na(cont_cov)]
    cont_cov = cont_cov[!is.na(cont_cov)]
  }
  # repeat over replicates
  rmses_allreps_allfolds = array(dim=c(nreps,kfolds))
  for (i in 1:nreps) {
    if (show.prog==T) print(i)
    if (is.null(Subj)) {
      k_inds = get_kfold_inds(Scan,kfolds)
    } else {
      k_inds = get_kfold_inds2(Scan,Subj,kfolds)
    }
    for (k in 1:kfolds) {
      training_data = data.frame(y = cont_cov[k_inds!=k], X = imgdata[k_inds!=k,])
      testing_data = data.frame(y = cont_cov[k_inds==k], X = imgdata[k_inds==k,])
      model_k_cv = cv.glmnet(x = as.matrix(training_data[,-1]), y = training_data[,1], alpha = 1)
      model_k = glmnet(x = as.matrix(training_data[,-1]), y = training_data[,1], alpha = 1, lambda = model_k_cv$lambda.min)
      ypred_rk = predict(model_k, newx = as.matrix(testing_data[,-1]), type = "response", s = model_k_cv$lambda.min)
      rmses_allreps_allfolds[i,k] = sqrt(mean((ypred_rk - testing_data[,1])^2,na.rm=T))
    }
  }
  return(rmses_allreps_allfolds)
}

lassoPrediction_cat <- function(imgdata, cat_cov, Scan, Subj=NULL, nreps = 10, kfolds = 3, seed=123, show.prog = T) {
  if (!is.factor(cat_cov)) cat_cov = factor(cat_cov)
  
  # remove any NAs
  if (anyNA(cat_cov)) {
    imgdata = imgdata[!is.na(cat_cov),]
    Scan = Scan[!is.na(cat_cov)]
    if (!is.null(Subj)) Subj = Subj[!is.na(cat_cov)]
    cat_cov = cat_cov[!is.na(cat_cov)]
  }
  
  # repeat over replicates
  accs_allreps_allfolds = array(dim=c(nreps,kfolds))
  for (i in 1:nreps) {
    if (show.prog==T) print(i)
    if (is.null(Subj)) {
      k_inds = get_kfold_inds(Scan,kfolds)
    } else {
      k_inds = get_kfold_inds2(Scan,Subj,kfolds)
    }
    for (k in 1:kfolds) {
      training_data = data.frame(y = cat_cov[k_inds!=k], X = imgdata[k_inds!=k,])
      testing_data = data.frame(y = cat_cov[k_inds==k], X = imgdata[k_inds==k,])
      model_k_cv = cv.glmnet(x = as.matrix(training_data[,-1]), y = training_data[,1], alpha = 1, family="binomial")
      model_k = glmnet(x = as.matrix(training_data[,-1]), y = training_data[,1], alpha = 1, lambda = model_k_cv$lambda.min, family="binomial")
      ypred_rk = predict(model_k, newx = as.matrix(testing_data[,-1]), type = "response", s = model_k_cv$lambda.min, family="binomial")
      ypred_rk_bin = ifelse(ypred_rk >= .5, 1, 0)
      accs_allreps_allfolds[i,k] = sum(diag(table(ypred_rk_bin, testing_data[,1]))) / nrow(testing_data)
    }
  }
  return(accs_allreps_allfolds)
}

pairwiseTTests = function(xlist,Methods) {
  L = length(xlist)
  pairwise_pvals = array(dim=c(L,L))
  for (i in 1:(L-1)) {
    xi = c(xlist[[i]]) %>% na.omit()
    for (j in (i+1):L) {
      xj = c(xlist[[j]]) %>% na.omit()
      pairwise_pvals[i,j] = t.test(xi,xj)$p.value
    }
  }
  colnames(pairwise_pvals) = Methods
  rownames(pairwise_pvals) = Methods
  return(pairwise_pvals)
}

showPredictionBoxplots = function(xlist, Methods, type="RMSE") {
  L = length(xlist)
  pred_reps = c()
  Method_reps = c()
  for (l in 1:L) {
    len = sum(!is.na(xlist[[l]]))
    xdata = xlist[[l]] %>% na.omit()
    pred_reps = c(pred_reps, xdata)
    Method_reps = c(Method_reps, rep(Methods[l],length(xdata)))
  }
  plot_data = data.frame(Method = factor(Method_reps), predval = pred_reps)
  plot_data %>%
    ggplot(aes(x=Method, y=predval, fill=Method))+
    geom_boxplot()+
    theme_bw()+
    ylab(type)
}




#TEST --------------------------------------------------------------------
# subj_coef_type = 3
# 
# #simulate data
# set.seed(12345)
# N=100
# nScanners=5
# nSubj = 50
# p=rep(10,3)
# D=length(p)
# R=2
# Q=5
# H=10
# Scanner_ID = sort(rep(c(1:nScanners),N/nScanners)) #sample(1:nScanners,N,replace=T)
# Subj_ID = sort(rep(c(1:nSubj),N/nSubj))
# #Subj_ID = sort(rep(c(1:nSubj),N/nSubj))[sample(N,N,replace=F)] #sample(1:nSubj,N,replace=T)
# table(Scanner_ID)
# table(Subj_ID)
# X_cov=array(rnorm(N*Q),dim=c(N,Q))
# hist(X_cov[,1])
# hist(X_cov[,2])
# 
# # generate random coefficients from low-rank decomp
# set.seed(12345)
# 
# # population intercept
# marg_pop = vector(mode="list",D)
# for (d in 1:D) {marg_pop[[d]] = array(rnorm(p[d],R),dim=c(p[d],R))}
# Gamma_pop_true = getGamma_cpp(marg_pop)
# 
# # covariate effects
# marg_cov = vector(mode="list",Q)
# for (q in 1:Q) {
#   marg_cov[[q]] = vector(mode="list",D)
#   for (d in 1:D) {
#     marg_cov[[q]][[d]] = array(rnorm(p[d]*R),dim=c(p[d],R))
#   }
# }
# Gamma_cov_true = getGamma_multi_cpp(marg_cov)
# 
# # scanner intercept
# marg_scan = vector(mode="list",nScanners)
# for (s in 1:nScanners) {
#   marg_scan[[s]] = vector(mode="list",D)
#   for (d in 1:D) {
#     marg_scan[[s]][[d]] = array(rnorm(p[d]*R),dim=c(p[d],R))
#   }
# }
# Gamma_scan_true = getGamma_multi_cpp(marg_scan)
# 
# # subject intercept
# mean_subj_devs = rnorm(nSubj,mean=0,sd=1)
# Gamma_subj_true = lapply(1:nSubj, function(u) {z=rep(0,prod(p)); z[sample(1:prod(p),round(prod(p)/2))] = rnorm(round(prod(p)/2),mean=mean_subj_devs[u],sd=0); return(z)})#lapply(1:nSubj,function(u) rnorm(prod(p),mean=mean_subj_devs[u],sd=.1))
# 
# SigSq_true=array(sample(c(1:H)^2,N*prod(p),replace=T),dim=c(N,prod(p)))
# 
# Yvec=array(dim=c(N,prod(p)))
# for (n in 1:N) {
#   Yvec[n,] = Gamma_pop_true + Gamma_scan_true[[Scanner_ID[n]]] + Gamma_subj_true[[Subj_ID[n]]]
#   for (q in 1:Q) {
#     Yvec[n,] = Yvec[n,] + X_cov[n,q]*Gamma_cov_true[[q]]
#   }
#   Yvec[n,] = Yvec[n,] + rnorm(prod(p),mean=0,sd=sqrt(SigSq_true[n,]))
# }


# Fit Model ---------------------------------------------------------------
# library(tictoc)
# library(profvis)
# set.seed(123)
# 
# niter = 100
# Rfitted = 2
# res_harm_A = runMCMC_harm_long_R(Yvec = Yvec, p = p, X_cov = X_cov, Scanner_ID = Scanner_ID, Subj_ID = Subj_ID,
#                                  niter = niter, R = Rfitted, H = H, prog_count = 1, show_all_steps = T, subj_coef_type = 3)
# 
# res_harm_B = runMCMC_harm_long_parallel_R(Yvec = Yvec, p = p, X_cov = X_cov, Scanner_ID = Scanner_ID, Subj_ID = Subj_ID,
#                                  niter = niter, R = Rfitted, H = H, prog_count = 1, show_all_steps = T, subj_coef_type = 3, null_cov = T, null_scan=T, null_subj=T)


# niter = 1
# ranks = c(1:5)
# reps = 5
# times_noparallel = times_parallel = array(dim=c(reps,length(ranks)))
# 
# for (rr in 1:reps) {
#   for (Rfitted in ranks) {
# 
#   tic()
#   res_harm_A = runMCMC_harm_long_R(Yvec = Yvec, p = p, X_cov = X_cov, Scanner_ID = Scanner_ID, Subj_ID = Subj_ID,
#                                 niter = niter, R = Rfitted, H = H, prog_count = 1, show_all_steps = T, subj_coef_type = 3)
#   t1 = toc()
#   times_noparallel[rr,Rfitted] = t1$toc-t1$tic
# 
#   tic()
#   res_harm_B = runMCMC_harm_long_parallel_R(Yvec = Yvec, p = p, X_cov = X_cov, Scanner_ID = Scanner_ID, Subj_ID = Subj_ID,
#                                niter = niter, R = Rfitted, H = H, prog_count = 1, show_all_steps = T, subj_coef_type = 3)
#   t2 = toc()
#   times_parallel[rr,Rfitted] = t2$toc-t2$tic
# 
#   }
# }
# 
# times_noparallel_avg = apply(times_noparallel, 2, mean)
# times_parallel_avg = apply(times_parallel, 2, mean)
# plot(ranks, times_noparallel_avg, type="l")
# plot(ranks, times_parallel_avg, type="l")
# 
# plot(ranks, times_noparallel_avg/times_parallel_avg, type="l" )




# profvis({res_harm = runMCMC_harm_long_R(Yvec = Yvec, p = p, X_cov = X_cov, Scanner_ID = Scanner_ID, Subj_ID = Subj_ID,
#                                niter = niter, R = Rfitted, H = H, prog_count = 1, show_all_steps = T, subj_coef_type = 3)
# })
# 
# 
# 
# # Rcpp (fastest version)
# profvis({gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1, sliceInds_alldims[[d]])})
# 
# 
# profvis({
# nobs = nrow(Yrvec_pop)
# D = length(p)
# gamma_d = gamma_pop_iter[[d]]
# gamma_dr = gamma_d[,r]
# gamma_dr_new = rep(0,p[d])
# c_dr = rep(0,p[d])
# d_dr = rep(0,p[d])
# w_dr = w_pop_iter[d,r];
# alpha_dr = alpha_pop_iter[d,r];
# if (D==3) {
#   gamma_x = gamma_pop_iter[[1]];
#   gamma_xr = gamma_x[,r];
#   gamma_y = gamma_pop_iter[[2]];
#   gamma_yr = gamma_y[,r];
#   gamma_z = gamma_pop_iter[[3]];
#   gamma_zr = gamma_z[,r];
#   if (d==1) {
#     outer_dr = TP_2D_cpp(gamma_yr, gamma_zr);
#   } else if (d==2) {
#     outer_dr = TP_2D_cpp(gamma_xr, gamma_zr);
#   } else if (d==3) {
#     outer_dr = TP_2D_cpp(gamma_xr, gamma_yr);
#   }
# } else if (D==2) {
#   gamma_x = gamma_pop_iter[[1]];
#   gamma_xr = gamma_x[,r];
#   gamma_y = gamma_pop_iter[[2]];
#   gamma_yr = gamma_y[,r];
#   if (d==1) {
#     outer_dr = gamma_yr;
#   } else if (d==2) {
#     outer_dr = gamma_xr;
#   }
# }
# 
# slice_inds_d = sliceInds_alldims[[d]]
# for (vd in 1:p[d]) {
#   slice_inds_dvd = slice_inds_d[vd,]+1
#   c_drv = 0;
#   d_drv = 0;
#   m_drv = 0;
#   n_drv = 0;
#   
#   for (i in 1:nobs) {
#     #SigSq_i = SigSq_iter[Scanner_ID[i],];
#     #Yrvec_i = Yrvec_pop[i,];
#     #SigSq_sl_div = SigSq_i[slice_inds_dvd]
#     #Y_sl_driv = Yrvec_i[slice_inds_dvd]
#     SigSq_sl_div = SigSq_iter[Scanner_ID[i],slice_inds_dvd];
#     Y_sl_driv = Yrvec_pop[i,slice_inds_dvd];
#     Ysl_isNA = is.na(Y_sl_driv);
# 
#     m_drv = m_drv + sum(((X[i]*outer_dr)^2/SigSq_sl_div)[Ysl_isNA==F])
#     n_drv = n_drv + sum(((X[i]*Y_sl_driv*outer_dr)/SigSq_sl_div)[Ysl_isNA==F])
# 
#     # for (j in 1:length(outer_dr)) {
#     #   if (Ysl_isNA[j]==FALSE) {
#     #     m_drv = m_drv+(X[i]*outer_dr[j])^2/SigSq_sl_div[j];
#     #     n_drv = n_drv+(X[i]*Y_sl_driv[j]*outer_dr[j])/SigSq_sl_div[j];
#     #   }
#     # }
# 
#   }
#   
# # mn_all = mclapply(1:nobs, function(i) {
# #     SigSq_sl_div = SigSq_iter[Scanner_ID[i],slice_inds_dvd];
# #     Y_sl_driv = Yrvec_pop[i,slice_inds_dvd];
# #     Ysl_isNA = is.na(Y_sl_driv);
# #     
# #     m_drv_i = sum(((X[i]*outer_dr)^2/SigSq_sl_div)[Ysl_isNA==F]);
# #     n_drv_i = sum(((X[i]*Y_sl_driv*outer_dr)/SigSq_sl_div)[Ysl_isNA==F]);
# #     return(c(m_drv_i,n_drv_i))
# #     }) %>% unlist %>% array(dim=c(2,nobs)) %>% t()
# #   m_drv = sum(mn_all[,1])
# #   n_drv = sum(mn_all[,2])
#   
#   
#   
#   if (vd==1) {
#     c_drv = m_drv+1/(tau_pop_iter*w_dr*(1-exp(-2*alpha_dr)))
#     d_drv = n_drv+(exp(-alpha_dr)*gamma_dr[vd+1])/(tau_pop_iter*w_dr*(1-exp(-2*alpha_dr)));
#   } else if (vd==p[d]) {
#     c_drv = m_drv+1/(tau_pop_iter*w_dr*(1-exp(-2*alpha_dr)));
#     d_drv = n_drv+(exp(-alpha_dr)*gamma_dr[vd-1])/(tau_pop_iter*w_dr*(1-exp(-2*alpha_dr)));
#   } else {
#     c_drv = m_drv+(1+exp(-2*alpha_dr))/(tau_pop_iter*w_dr*(1-exp(-2*alpha_dr)));
#     d_drv = n_drv+(exp(-alpha_dr)*(gamma_dr[vd-1]+gamma_dr[vd+1]))/(tau_pop_iter*w_dr*(1-exp(-2*alpha_dr)));
#   }
#   
#   if (c_drv==0) {
#     gamma_dr_new[vd] = 0;
#   } else {
#     gamma_dr_new[vd] = rnorm(n = 1, mean = d_drv/c_drv, sd = sqrt(1.0/c_drv))
#   }
# }
# })
# 
# 
# # R (with parallel computing)
# profvis({gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_parallel_R(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d, r, sliceInds_alldims[[d]])})
# 
# 
# 
# profvis({
#   nobs = nrow(Yrvec_pop)
#   D = length(p)
#   gamma_d = gamma_pop_iter[[d]]
#   gamma_dr = gamma_d[,r]
#   gamma_dr_new = rep(0,p[d])
#   c_dr = rep(0,p[d])
#   d_dr = rep(0,p[d])
#   w_dr = w_pop_iter[d,r];
#   alpha_dr = alpha_pop_iter[d,r];
#   if (D==3) {
#     gamma_x = gamma_pop_iter[[1]];
#     gamma_xr = gamma_x[,r];
#     gamma_y = gamma_pop_iter[[2]];
#     gamma_yr = gamma_y[,r];
#     gamma_z = gamma_pop_iter[[3]];
#     gamma_zr = gamma_z[,r];
#     if (d==1) {
#       outer_dr = TP_2D_cpp(gamma_yr, gamma_zr);
#     } else if (d==2) {
#       outer_dr = TP_2D_cpp(gamma_xr, gamma_zr);
#     } else if (d==3) {
#       outer_dr = TP_2D_cpp(gamma_xr, gamma_yr);
#     }
#   } else if (D==2) {
#     gamma_x = gamma_pop_iter[[1]];
#     gamma_xr = gamma_x[,r];
#     gamma_y = gamma_pop_iter[[2]];
#     gamma_yr = gamma_y[,r];
#     if (d==1) {
#       outer_dr = gamma_yr;
#     } else if (d==2) {
#       outer_dr = gamma_xr;
#     }
#   }
#   
#   slice_inds_d = sliceInds_alldims[[d]]
#   
#   gamma_dr_new = mclapply(1:p[d], function(vd) {
#     slice_inds_dvd = slice_inds_d[vd,]+1
#     c_drv = d_drv = 0
#     mn_all = mclapply(1:nobs, function(i) {
#       SigSq_sl_div = SigSq_iter[Scanner_ID[i],slice_inds_dvd];
#       Y_sl_driv = Yrvec_pop[i,slice_inds_dvd];
#       Ysl_isNA = is.na(Y_sl_driv);
#       
#       m_drv_i = sum(((X[i]*outer_dr)^2/SigSq_sl_div)[Ysl_isNA==F]);
#       n_drv_i = sum(((X[i]*Y_sl_driv*outer_dr)/SigSq_sl_div)[Ysl_isNA==F]);
#       return(c(m_drv_i,n_drv_i))
#     }) %>% unlist %>% array(dim=c(2,nobs)) %>% t()
#     m_drv = sum(mn_all[,1])
#     n_drv = sum(mn_all[,2])
#     
#     if (vd==1) {
#       c_drv = m_drv+1/(tau_pop_iter*w_dr*(1-exp(-2*alpha_dr)))
#       d_drv = n_drv+(exp(-alpha_dr)*gamma_dr[vd+1])/(tau_pop_iter*w_dr*(1-exp(-2*alpha_dr)));
#     } else if (vd==p[d]) {
#       c_drv = m_drv+1/(tau_pop_iter*w_dr*(1-exp(-2*alpha_dr)));
#       d_drv = n_drv+(exp(-alpha_dr)*gamma_dr[vd-1])/(tau_pop_iter*w_dr*(1-exp(-2*alpha_dr)));
#     } else {
#       c_drv = m_drv+(1+exp(-2*alpha_dr))/(tau_pop_iter*w_dr*(1-exp(-2*alpha_dr)));
#       d_drv = n_drv+(exp(-alpha_dr)*(gamma_dr[vd-1]+gamma_dr[vd+1]))/(tau_pop_iter*w_dr*(1-exp(-2*alpha_dr)));
#     }
#     
#     if (c_drv==0) {
#       return(0);
#     } else {
#       return(rnorm(n = 1, mean = d_drv/c_drv, sd = sqrt(1.0/c_drv)))
#     }
#     
#   }) %>% unlist()
# })
# 
# 
# 
# ### Test sampling step for one term
# tic()
# Y_intercept = getY_intercept_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_cov_iter, Gamma_scan_iter, Gamma_subj_iter)
# for (r in 1:R[1]) {
#   Yrvec_pop = getYr_cpp(Y_intercept, X_intercept, gamma_pop_iter, r-1)
#   for (d in 1:D) {
#     gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1, sliceInds_alldims[[d]])
#     #gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_parallel_R(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d, r, sliceInds_alldims[[d]])
#     #gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1)
#     w_pop_iter[d,r] = sample_w_dr_cpp(gamma_pop_iter, tau_pop_iter, alpha_pop_iter, lambda_pop_iter, p, d-1, r-1)
#     lambda_pop_iter[d,r] = sample_lambda_dr_cpp(w_pop_iter, p, a_lam, b_lam, d-1, r-1)
#     alpha_pop_iter[d,r] = sample_alpha_dr_cpp(alpha_pop_iter, gamma_pop_iter, tau_pop_iter, w_pop_iter, p, a_alpha, b_alpha, sigma_log_alpha, d-1, r-1)
#   }
# }
# toc()
# 
# tic()
# Y_intercept = getY_intercept_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_cov_iter, Gamma_scan_iter, Gamma_subj_iter)
# for (r in 1:R[1]) {
#   Yrvec_pop = getYr_cpp(Y_intercept, X_intercept, gamma_pop_iter, r-1)
#   for (d in 1:D) {
#     #gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1, sliceInds_alldims[[d]])
#     gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_parallel_R(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d, r, sliceInds_alldims[[d]])
#     #gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1)
#     w_pop_iter[d,r] = sample_w_dr_cpp(gamma_pop_iter, tau_pop_iter, alpha_pop_iter, lambda_pop_iter, p, d-1, r-1)
#     lambda_pop_iter[d,r] = sample_lambda_dr_cpp(w_pop_iter, p, a_lam, b_lam, d-1, r-1)
#     alpha_pop_iter[d,r] = sample_alpha_dr_cpp(alpha_pop_iter, gamma_pop_iter, tau_pop_iter, w_pop_iter, p, a_alpha, b_alpha, sigma_log_alpha, d-1, r-1)
#   }
# }
# toc()
# 
# tic()
# Y_intercept = getY_intercept_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_cov_iter, Gamma_scan_iter, Gamma_subj_iter)
# mclapply(1:R[1], function(r) {
#   Yrvec_pop = getYr_cpp(Y_intercept, X_intercept, gamma_pop_iter, r-1);
#   mclapply(1:D, function(d) {
#     gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1, sliceInds_alldims[[d]]);
#     #gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_parallel_R(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d, r, sliceInds_alldims[[d]])
#     #gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1)
#     w_pop_iter[d,r] = sample_w_dr_cpp(gamma_pop_iter, tau_pop_iter, alpha_pop_iter, lambda_pop_iter, p, d-1, r-1);
#     lambda_pop_iter[d,r] = sample_lambda_dr_cpp(w_pop_iter, p, a_lam, b_lam, d-1, r-1);
#     alpha_pop_iter[d,r] = sample_alpha_dr_cpp(alpha_pop_iter, gamma_pop_iter, tau_pop_iter, w_pop_iter, p, a_alpha, b_alpha, sigma_log_alpha, d-1, r-1);
#   });
#   # for (d in 1:D) {
#   #   gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_withInds_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1, sliceInds_alldims[[d]])
#   #   #gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_parallel_R(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d, r, sliceInds_alldims[[d]])
#   #   #gamma_pop_iter[[d]][,r] = sample_gamma_dr_harm_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d-1, r-1)
#   #   w_pop_iter[d,r] = sample_w_dr_cpp(gamma_pop_iter, tau_pop_iter, alpha_pop_iter, lambda_pop_iter, p, d-1, r-1)
#   #   lambda_pop_iter[d,r] = sample_lambda_dr_cpp(w_pop_iter, p, a_lam, b_lam, d-1, r-1)
#   #   alpha_pop_iter[d,r] = sample_alpha_dr_cpp(alpha_pop_iter, gamma_pop_iter, tau_pop_iter, w_pop_iter, p, a_alpha, b_alpha, sigma_log_alpha, d-1, r-1)
#   # }
# })
# toc()
# 






# # Get Results -------------------------------------------------------------

# res_harm = res_harm_B # res_harm_B
# 
# ### Population intercept
# Gamma_pop_mcmc = getCoef_mcmc_R(res_harm$gamma_pop)
# burn.in = .5
# plot(Gamma_pop_mcmc[round(burn.in*(niter+1)):(niter+1),1])
# acf(Gamma_pop_mcmc[round(burn.in*(niter+1)):(niter+1),1],main='Gamma(v) - Harm model')
# Gamma_pop_est = apply(Gamma_pop_mcmc[round(burn.in*(niter+1)):(niter+1),],2,mean)
# plot(Gamma_pop_est,Gamma_pop_true)
# cor(Gamma_pop_est,Gamma_pop_true)
# 
# ### Scanner intercept
# Gamma_scan_mcmc = getCoef_mcmc_R(res_harm$gamma_scan)
# burn.in = .5
# plot(Gamma_scan_mcmc[round(burn.in*(niter+1)):(niter+1),1,1])
# acf(Gamma_scan_mcmc[round(burn.in*(niter+1)):(niter+1),1,1])
# Gamma_scan_est = apply(Gamma_scan_mcmc[round(burn.in*(niter+1)):(niter+1),,],c(2,3),mean)
# plot(Gamma_scan_est[1,],Gamma_scan_true[[1]])
# cor(Gamma_scan_est[1,],Gamma_scan_true[[1]])
# plot(Gamma_scan_est[2,],Gamma_scan_true[[2]])
# cor(Gamma_scan_est[2,],Gamma_scan_true[[2]])
# sapply(1:nScanners, function(i) cor(Gamma_scan_est[i,],Gamma_scan_true[[i]]))
# 
# ### Covariate effects
# Gamma_cov_mcmc = getCoef_mcmc_R(res_harm$gamma_cov)
# burn.in=.5
# plot(Gamma_cov_mcmc[round(burn.in*(niter+1)):(niter+1),1,1])
# acf(Gamma_cov_mcmc[round(burn.in*(niter+1)):(niter+1),1,1],main='Gamma(v) - Harm model')
# Gamma_cov_est = apply(Gamma_cov_mcmc[round(burn.in*(niter+1)):(niter+1),,],c(2,3),mean)
# plot(Gamma_cov_est[1,],Gamma_cov_true[[1]])
# cor(Gamma_cov_est[1,],Gamma_cov_true[[1]])
# plot(Gamma_cov_est[2,],Gamma_cov_true[[2]])
# cor(Gamma_cov_est[2,],Gamma_cov_true[[2]])
# sapply(1:Q,function(q) cor(Gamma_cov_est[q,], Gamma_cov_true[[q]]))
# 
# ### Subject intercept
# if (subj_coef_type==1) {
#   Gamma_subj_mcmc = getCoef_mcmc_R(res_harm$gamma_subj)
# } else {
#   Gamma_subj_mcmc = array(0,dim=c(niter+1,nSubj,prod(p)))
#   for (iter in 1:(niter+1)) {
#     for (u in 1:nSubj) {
#       if (subj_coef_type==2) {
#         Gamma_subj_mcmc[iter,u,] = res_harm$gamma_subj[[iter]][u]
#       } else if (subj_coef_type==3) {
#         Gamma_subj_mcmc[iter,u,res_harm$Z_subj[[iter]][u,]==1] = res_harm$gamma_subj[[iter]][u]
#       }
#     }
#   }
# }
# burn.in = .5
# plot(Gamma_subj_mcmc[round(burn.in*(niter+1)):(niter+1),1,1])
# Gamma_subj_est = apply(Gamma_subj_mcmc[round(burn.in*(niter+1)):(niter+1),,],c(2,3),mean)
# plot(Gamma_subj_est[1,],Gamma_subj_true[[1]])
# mean(Gamma_subj_est[1,])
# cor(Gamma_subj_est[1,],Gamma_subj_true[[1]])
# plot(Gamma_subj_est[2,],Gamma_subj_true[[2]])
# mean(Gamma_subj_est[2,])
# cor(Gamma_subj_est[2,],Gamma_subj_true[[2]])
# sapply(1:nSubj,function(u) cor(Gamma_subj_est[u,], Gamma_subj_true[[u]]))
# 
# 
# # get residual variance
# SigSq_mcmc = getSigSq_mcmc_R(Scanner_ID, res_harm$ssq, res_harm$Zeta)
# burn.in=.5
# plot(SigSq_mcmc[round(burn.in*(niter+1)):(niter+1),1,1])
# acf(SigSq_mcmc[round(burn.in*(niter+1)):(niter+1),1,1],main='SigSq(v)')
# SigSq_est = apply(SigSq_mcmc[round(burn.in*(niter+1)):(niter+1),,],c(2,3),mean)
# hist(SigSq_est[1,])
# hist(SigSq_true[1,])
# boxplot(SigSq_est[1,]~SigSq_true[1,])
# cor(SigSq_true[1,],SigSq_est[1,])
# 
