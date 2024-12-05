source('TC_MCMC.R')

tcombat_harmonize = function(Y, Scan = Scan, X = NULL, Subj = NULL, R = 5, H = 5, mcmc_iter = 1000, burn.in = .5, long = F, covar = T, verbose = T) {
  p = dim(Y)[-1]
  if (verbose) print('Running MCMC...')
  MCMCres = runMCMC_harm_long_R(Yvec = Yvec, p = p, X_cov = X, Scanner_ID = Scan, Subj_ID = Subj, 
                                niter = mcmc_iter, R = R, H = H, prog_count = 10, show_all_steps = F, null_subj = is.null(Subj), null_cov = !covar)
  
  nobs = nrow(Y)
  V = ncol(Y)
  na_vox = apply(Y,2,anyNA)
  uncomp = any(grepl("comp",names(MCMCres)))
  niter = length(MCMCres$gamma_pop)
  ss = round(burn.in*niter):niter
  nss = length(ss)
  if (is.null(na_vox)) na_vox = rep(F,V)
  if (is.null(X)) covar = F
  nCov = ncol(X)
  nScan = length(unique(Scan))
  nSubj = length(unique(Subj))
  
  PopInt = rep(0,V)
  ScanInt = array(0,dim=c(nScan,V))
  CovEff = SubjInt = 0
  if (covar == F) {
    X = as.matrix(rep(0,nobs))
  }
  if (covar == T) CovEff = array(0,dim=c(nCov,V))
  if (long == T) SubjInt = array(0,dim=c(nSubj,V))
  ScanResidVar = array(0,dim=c(nScan,V))
  for (i in 1:nss) {
    if (i %% 10 == 0) print(paste0('Iteration = ',ss[i]+1,'/',niter-1))
    Gamma_pop_i = TP_rankR_cpp(MCMCres$gamma_pop[[ss[i]]])
    PopInt = PopInt + Gamma_pop_i/nss
    for (j in 1:nScan) {
      Gamma_scan_ij = TP_rankR_cpp(MCMCres$gamma_scan[[ss[i]]][[j]])
      ScanInt[j,] =  ScanInt[j,] + Gamma_scan_ij/nss
    }
    if (covar == T) {
      for (k in 1:nCov) {
        Gamma_cov_ik = TP_rankR_cpp(MCMCres$gamma_cov[[ss[i]]][[k]])
        CovEff[k,] = CovEff[k,] + Gamma_cov_ik/nss
      }
    }
    if (long == T) {
      if (uncomp == T) {
        Zi = uncomp_Z_subj(MCMCres$Z_subj_comp[[ss[i]]], nSubj = nSubj, V = V, na_vox = na_vox)
      } else {
        Zi = MCMCres$Z_subj[[ss[i]]]
      }
      gamma_subj_i = MCMCres$gamma_subj[[ss[i]]]
      SubjInt_i = 0*SubjInt
      for (j in 1:nSubj) {
        SubjInt_i[j,Zi[j,]==1] = gamma_subj_i[j]
      }
      SubjInt = SubjInt + SubjInt_i/nss
    }
    if (uncomp == T) {
      Zeta_i = uncomp_Zeta(MCMCres$Zeta_comp[[ss[i]]],nScan,V,H = H)
    } else {
      Zeta_i = MCMCres$Zeta[[ss[i]]]
    }
    ssq_i = MCMCres$ssq[[ss[i]]]
    SigSq_i = getSigSq_cpp(ssq = ssq_i, Zeta = Zeta_i)
    ScanResidVar = ScanResidVar + SigSq_i/nss
  }
  
  ## get additive (~centered around 0) and multiplicative deviations (~centered around 1)
  weights = as.numeric(table(Scan))
  ScanIntAvg = apply(ScanInt,2,function(x) sum(weights*x)/sum(weights))
  ScanDev = t(apply(ScanInt,1,function(x) x-ScanIntAvg))
  PopInt_agg = PopInt+ScanIntAvg
  ScanResidVarAvg = apply(ScanResidVar,2,function(x) sum(weights*x)/sum(weights))
  ScanResidVarDev = t(apply(ScanResidVar,1,function(x) sqrt(x/ScanResidVarAvg)))
  
  ## Set terms to NA where na_vox==T
  PopInt_agg[na_vox==T] = NA
  ScanDev[,na_vox==T] = NA
  if (covar==T) CovEff[,na_vox==T] = NA
  if (long == T) SubjInt[,na_vox==T] = NA
  ScanResidVarDev[,na_vox==T] = NA
  
  ## apply harmonization step
  if (long == F) {
    Harm = t(sapply(1:nobs, function(n) (Y[n,] - PopInt_agg - ScanDev[Scan[n],] - t(X[n,]%*%CovEff))/sqrt(ScanResidVarDev[Scan[n],]) + PopInt_agg + t(X[n,]%*%CovEff)))
  } else {
    Harm = t(sapply(1:nobs, function(n) (Y[n,] - PopInt_agg - ScanDev[Scan[n],] - SubjInt[Subj[n],] - t(X[n,]%*%CovEff))/sqrt(ScanResidVarDev[Scan[n],]) + PopInt_agg + SubjInt[Subj[n],] + t(X[n,]%*%CovEff)))
  }
  return(Harm)
}