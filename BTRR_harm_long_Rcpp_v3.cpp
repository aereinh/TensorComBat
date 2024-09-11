//
//  BTRR_harm_Rcpp.cpp
//  
//
//  Created by Alec Reinhardt on 10/27/22.
//

//#include <stdio.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double gig_mode(double lambda, double omega) {
    if (lambda >= 1.)
        return (sqrt((lambda-1.)*(lambda-1.) + omega*omega)+(lambda-1.))/omega;
    else
        return omega / (sqrt((1.-lambda)*(1.-lambda) + omega*omega)+(1.-lambda));
}

// [[Rcpp::export]]
double rgig_ROU_shift_alt (double lambda, double lambda_old, double omega, double alpha) {
    double xm, nc;     /* location of mode; c=log(f(xm)) normalization constant */
    double s, t;       /* auxiliary variables */
    double U, V, X;    /* random variables */

    int count = 0;     /* counter for total number of iterations */

    double a, b, c;    /* coefficent of cubic */
    double p, q;       /* coefficents of depressed cubic */
    double fi, fak;    /* auxiliary results for Cardano's rule */

    double y1, y2;     /* roots of (1/x)*sqrt(f((1/x)+m)) */

    double uplus, uminus;  /* maximum and minimum of x*sqrt(f(x+m)) */
    
    double res;

    /* -- Setup -------------------------------------------------------------- */

    /* shortcuts */
    t = 0.5 * (lambda-1.);
    s = 0.25 * omega;

    /* mode = location of maximum of sqrt(f(x)) */
    xm = gig_mode(lambda, omega);

    /* normalization constant: c = log(sqrt(f(xm))) */
    nc = t*log(xm) - s*(xm + 1./xm);

    /* location of minimum and maximum of (1/x)*sqrt(f(1/x+m)):  */

    /* compute coeffients of cubic equation y^3+a*y^2+b*y+c=0 */
    a = -(2.*(lambda+1.)/omega + xm);       /* < 0 */
    b = (2.*(lambda-1.)*xm/omega - 1.);
    c = xm;

    /* we need the roots in (0,xm) and (xm,inf) */

    /* substitute y=z-a/3 for depressed cubic equation z^3+p*z+q=0 */
    p = b - a*a/3.;
    q = (2.*a*a*a)/27. - (a*b)/3. + c;

    /* use Cardano's rule */
    fi = acos(-q/(2.*sqrt(-(p*p*p)/27.)));
    fak = 2.*sqrt(-p/3.);
    y1 = fak * cos(fi/3.) - a/3.;
    y2 = fak * cos(fi/3. + 4./3.*M_PI) - a/3.;

    /* boundaries of minmal bounding rectangle:                  */
    /* we us the "normalized" density f(x) / f(xm). hence        */
    /* upper boundary: vmax = 1.                                 */
    /* left hand boundary: uminus = (y2-xm) * sqrt(f(y2)) / sqrt(f(xm)) */
    /* right hand boundary: uplus = (y1-xm) * sqrt(f(y1)) / sqrt(f(xm)) */
    uplus  = (y1-xm) * exp(t*log(y1) - s*(y1 + 1./y1) - nc);
    uminus = (y2-xm) * exp(t*log(y2) - s*(y2 + 1./y2) - nc);

    /* -- Generate sample ---------------------------------------------------- */

    do {
        ++count;
        U = uminus + unif_rand() * (uplus - uminus);    /* U(u-,u+)  */
        V = unif_rand();                                /* U(0,vmax) */
        X = U/V + xm;
    }                                         /* Acceptance/Rejection */
    while ((X <= 0.) || ((log(V)) > (t*log(X) - s*(X + 1./X) - nc)));

    /* store random point */
    res = (lambda_old < 0.) ? (alpha / X) : (alpha * X);

    /* -- End ---------------------------------------------------------------- */
    return res;
}

// [[Rcpp::export]]
double rgig_ROU_noshift (double lambda, double lambda_old, double omega, double alpha) {
    double xm, nc;     /* location of mode; c=log(f(xm)) normalization constant */
    double ym, um;     /* location of maximum of x*sqrt(f(x)); umax of MBR */
    double s, t;       /* auxiliary variables */
    double U, V, X;    /* random variables */

    int count = 0;     /* counter for total number of iterations */
    
    double res;

    /* -- Setup -------------------------------------------------------------- */

    /* shortcuts */
    t = 0.5 * (lambda-1.);
    s = 0.25 * omega;
  
    /* mode = location of maximum of sqrt(f(x)) */
    xm = gig_mode(lambda, omega);

    /* normalization constant: c = log(sqrt(f(xm))) */
    nc = t*log(xm) - s*(xm + 1./xm);

    /* location of maximum of x*sqrt(f(x)):           */
    /* we need the positive root of                   */
    /*    omega/2*y^2 - (lambda+1)*y - omega/2 = 0    */
    ym = ((lambda+1.) + sqrt((lambda+1.)*(lambda+1.) + omega*omega))/omega;

    /* boundaries of minmal bounding rectangle:                   */
    /* we us the "normalized" density f(x) / f(xm). hence         */
    /* upper boundary: vmax = 1.                                  */
    /* left hand boundary: umin = 0.                              */
    /* right hand boundary: umax = ym * sqrt(f(ym)) / sqrt(f(xm)) */
    um = exp(0.5*(lambda+1.)*log(ym) - s*(ym + 1./ym) - nc);

    /* -- Generate sample ---------------------------------------------------- */

    do {
      ++count;
      U = um * unif_rand();        /* U(0,umax) */
      V = unif_rand();             /* U(0,vmax) */
      X = U/V;
    }                              /* Acceptance/Rejection */
    while (((log(V)) > (t*log(X) - s*(X + 1./X) - nc)));

    /* store random point */
    res = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
  
  /* -- End ---------------------------------------------------------------- */

  return res;
}

// [[Rcpp::export]]
double rgig_newapproach1 (double lambda, double lambda_old, double omega, double alpha) {
    /* parameters for hat function */
    double A[3], Atot;  /* area below hat */
    double k0;          /* maximum of PDF */
    double k1, k2;      /* multiplicative constant */

    double xm;          /* location of mode */
    double x0;          /* splitting point T-concave / T-convex */
    double a;           /* auxiliary variable */

    double U, V, X;     /* random numbers */
    double hx;          /* hat at X */

    int count = 0;      /* counter for total number of iterations */
    
    double res;

    /* -- Check arguments ---------------------------------------------------- */

    if (lambda >= 1. || omega >1.)
        res = 0;
        //error ("invalid parameters");

    /* -- Setup -------------------------------------------------------------- */

    /* mode = location of maximum of sqrt(f(x)) */
    xm = gig_mode(lambda, omega);

    /* splitting point */
    x0 = omega/(1.-lambda);

    /* domain [0, x_0] */
    k0 = exp((lambda-1.)*log(xm) - 0.5*omega*(xm + 1./xm));     /* = f(xm) */
    A[0] = k0 * x0;

    /* domain [x_0, Infinity] */
    if (x0 >= 2./omega) {
        k1 = 0.;
        A[1] = 0.;
        k2 = pow(x0, lambda-1.);
        A[2] = k2 * 2. * exp(-omega*x0/2.)/omega;
    }
  
    else {
        /* domain [x_0, 2/omega] */
        k1 = exp(-omega);
        A[1] = (lambda == 0.)
        ? k1 * log(2./(omega*omega))
        : k1 / lambda * ( pow(2./omega, lambda) - pow(x0, lambda) );

        /* domain [2/omega, Infinity] */
        k2 = pow(2/omega, lambda-1.);
        A[2] = k2 * 2 * exp(-1.)/omega;
    }

    /* total area */
    Atot = A[0] + A[1] + A[2];

    /* -- Generate sample ---------------------------------------------------- */

    do {
      ++count;

      /* get uniform random number */
      V = Atot * unif_rand();
      
      do {
    
    /* domain [0, x_0] */
    if (V <= A[0]) {
      X = x0 * V / A[0];
      hx = k0;
      break;
    }
    
    /* domain [x_0, 2/omega] */
    V -= A[0];
    if (V <= A[1]) {
      if (lambda == 0.) {
        X = omega * exp(exp(omega)*V);
        hx = k1 / X;
      }
      else {
        X = pow(pow(x0, lambda) + (lambda / k1 * V), 1./lambda);
        hx = k1 * pow(X, lambda-1.);
      }
      break;
    }
    
    /* domain [max(x0,2/omega), Infinity] */
    V -= A[1];
    a = (x0 > 2./omega) ? x0 : 2./omega;
    X = -2./omega * log(exp(-omega/2. * a) - omega/(2.*k2) * V);
    hx = k2 * exp(-omega/2. * X);
    break;
    
      } while(0);
      
      /* accept or reject */
      U = unif_rand() * hx;

      if (log(U) <= (lambda-1.) * log(X) - omega/2. * (X+1./X)) {
    /* store random point */
    res = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
    break;
      }
    } while(1);
    
    /* -- End ---------------------------------------------------------------- */

    return res;
}

// [[Rcpp::export]]
double do_rgig(double lambda, double chi, double psi) {
    double res = 0.0;
    double ZTOL=.8*DBL_EPSILON;

    /* alternative parametrization */
    double omega = sqrt(psi*chi);
    double alpha = sqrt(chi/psi);
    double abs_lambda = (lambda >= 0.) ? lambda : -lambda;
    
    /* run generator */
    if (omega < ZTOL) {
        /* for very small values of omega we have serious round-off errors */
        if (lambda > 0.0) {
            /* special cases which are basically Gamma distribution */
            //NumericVector gvec = rgamma(1,lambda,2.0/psi);
            //res = gvec(0);
            res = R::rgamma(lambda,2.0/psi);
        }
        else if (lambda < 0.0) {
            /* special cases which are basically Gamma and Inverse Gamma distribution */
            //NumericVector gvec = rgamma(1,-lambda,2.0/chi);
            //res = 1.0/gvec(0);
            res = 1.0/R::rgamma(-lambda,2.0/chi);
        }
    } else if (abs_lambda > 2. || omega > 3.) {
        /* Ratio-of-uniforms with shift by 'mode', alternative implementation */
        res = rgig_ROU_shift_alt(abs_lambda, lambda, omega, alpha);
    } else if (abs_lambda >= 1.-2.25*omega*omega || omega > 0.2) {
        /* Ratio-of-uniforms without shift */
        res = rgig_ROU_noshift(abs_lambda, lambda, omega, alpha);
    } else if (omega > 0.) { /* remaining case */
        /* New approach, constant hat in log-concave part. */
        res = rgig_newapproach1(abs_lambda, lambda, omega, alpha);
    }

    return res;
}

//Functions for sampling MCMC

// [[Rcpp::export]]
NumericVector rdirichlet_cpp(NumericVector alpha_vec) {
    int K = alpha_vec.length();
    NumericVector y_vec(K);
    NumericVector x_vec(K);
    for (int k=0; k<K; ++k) {
        //NumericVector x = rgamma(1,alpha_vec(k),1);
        //y_vec(k) = x(0);
        y_vec(k) = R::rgamma(alpha_vec(k),1);
    }
    x_vec = y_vec/sum(y_vec);
    return x_vec;
}

// [[Rcpp::export]]
int sample_1int_rcpp(NumericVector probs) {
    int ncomp = probs.length();
    NumericVector cumprobs = cumsum(probs);
    //NumericVector rand_num_vec = runif(1,0,1);
    //double rand_num = rand_num_vec(0);
    double rand_num = R::runif(0,1);
    int k = 0;
    while (rand_num >= cumprobs(k) && k<ncomp) {
        k++;
    }
    return k+1;
}

// [[Rcpp::export]]
IntegerVector modVec_cpp(IntegerVector x, int n) {
  IntegerVector y = x - (x/n)*n;
  return y;
}

// [[Rcpp::export]]
IntegerVector getSliceIndices_iota_cpp(IntegerVector p, int d, int j){
  int D = p.length();
  IntegerVector sliceIndices;
  if (D==2) {
    if (d==0) {
      sliceIndices = IntegerVector(p(1));
      std::iota(sliceIndices.begin(),sliceIndices.end(),0);
      sliceIndices = j+sliceIndices*p(0);
    } else if (d==1) {
      sliceIndices = IntegerVector(p(0));
      std::iota(sliceIndices.begin(),sliceIndices.end(),0);
      sliceIndices = j*p(0)+modVec_cpp(sliceIndices,p(0));
    }
  } else if (D==3) {
    if (d==0) {
      sliceIndices = IntegerVector(p(1)*p(2));
      std::iota(sliceIndices.begin(),sliceIndices.end(),0);
      sliceIndices = j+p(0)*modVec_cpp(sliceIndices,p(1))+(sliceIndices/p(1))*p(0)*p(1);
    } else if (d==1) {
      sliceIndices = IntegerVector(p(0)*p(2));
      std::iota(sliceIndices.begin(),sliceIndices.end(),0);
      sliceIndices = j*p(0)+modVec_cpp(sliceIndices,p(0))+(sliceIndices/p(0))*p(0)*p(1);
    } else if (d==2) {
      sliceIndices = IntegerVector(p(0)*p(1));
      std::iota(sliceIndices.begin(),sliceIndices.end(),0);
      sliceIndices = j*p(0)*p(1)+modVec_cpp(sliceIndices,p(0)*p(1));
    }
  }
  return sliceIndices;
}

// [[Rcpp::export]]
IntegerVector getSliceIndices_loop_cpp(IntegerVector p, int d, int j){
  int D = p.length();
  IntegerVector sliceIndices;
  if (D==2) {
    if (d==0) {
      sliceIndices = IntegerVector(p(1));
      for (int vv=0; vv<p(1); ++vv) {
        sliceIndices(vv) = j+vv*p(0);
      }
    } else if (d==1) {
      sliceIndices = IntegerVector(p(0));
      for (int vv=0; vv<p(0); ++vv) {
        sliceIndices(vv) = j*p(0)+(vv%p(0));
      }
    }
  } else if (D==3) {
    if (d==0) {
      sliceIndices = IntegerVector(p(1)*p(2));
      for (int vv=0; vv<(p(1)*p(2)); ++vv) {
        sliceIndices(vv) = j+(vv%p(1))*p(0)+(vv/p(1))*p(0)*p(1);
      }
    } else if (d==1) {
      sliceIndices = IntegerVector(p(0)*p(2));
      for (int vv=0; vv<(p(0)*p(2)); ++vv) {
        sliceIndices(vv) = j*p(0)+vv%p(0)+(vv/p(0))*p(0)*p(1);
      }
    } else if (d==2) {
      sliceIndices = IntegerVector(p(0)*p(1));
      for (int vv=0; vv<(p(0)*p(1)); ++vv) {
        sliceIndices(vv) = j*p(0)*p(1)+(vv%(p(0)*p(1)));
      }
    }
  }
  return sliceIndices;
}

// [[Rcpp::export]]
NumericVector slice_tensor_cpp(NumericVector X, IntegerVector p, int d, int j) {
    int D = p.length();
    NumericVector tensor_sl;
    if (D==2) { // 2D case
        if (d==0) {
            tensor_sl = NumericVector(p(1));
            for (int vv=0; vv<p(1); ++vv) {
                int tensor_ind = j+vv*p(0);
                tensor_sl(vv) = X(tensor_ind);
            }
        } else if (d==1) {
            tensor_sl = NumericVector(p(0));
            for (int vv=0; vv<p(0); ++vv) {
                int tensor_ind = j*p(0)+(vv%p(0));
                tensor_sl(vv) = X(tensor_ind);
            }
        }
    } else if (D==3) {  // 3D case
        if (d==0) {
            tensor_sl = NumericVector(p(1)*p(2));
            for (int vv=0; vv<(p(1)*p(2)); ++vv) {
                int tensor_ind = j+(vv%p(1))*p(0)+(vv/p(1))*p(0)*p(1);
                tensor_sl(vv) = X(tensor_ind);
            }
        } else if (d==1) {
            tensor_sl = NumericVector(p(0)*p(2));
            for (int vv=0; vv<(p(0)*p(2)); ++vv) {
                int tensor_ind = j*p(0)+vv%p(0)+(vv/p(0))*p(0)*p(1);
                tensor_sl(vv) = X(tensor_ind);
            }
        } else if (d==2) {
            tensor_sl = NumericVector(p(0)*p(1));
            for (int vv=0; vv<(p(0)*p(1)); ++vv) {
                int tensor_ind = j*p(0)*p(1)+(vv%(p(0)*p(1)));
                tensor_sl(vv) = X(tensor_ind);
            }
        }
    }
    return tensor_sl;
}

// [[Rcpp::export]]
NumericVector slice_tensor_cpp2(NumericVector X, IntegerVector p, int d, int j) {
  IntegerVector slice_inds = getSliceIndices_iota_cpp(p, d, j);
  NumericVector y = X[slice_inds];
  return y;
}

// [[Rcpp::export]]
NumericVector slice_tensor_cpp3(NumericVector X, IntegerVector p, int d, int j) {
  IntegerVector slice_inds = getSliceIndices_loop_cpp(p, d, j);
  NumericVector y = X[slice_inds];
  return y;
}

// [[Rcpp::export]]
NumericVector TP_2D_cpp(NumericVector x, NumericVector y) {
    int px = x.length();
    int py = y.length();
    NumericVector TP_vec(px*py);
    int vx; int vy;
    for (int v=0; v<(px*py); ++v) {
        vx = v%px;
        vy = (v/px)%py;
        TP_vec(v) = x(vx)*y(vy);
    }
    return TP_vec;
}

// [[Rcpp::export]]
NumericVector TP_3D_cpp(NumericVector x, NumericVector y, NumericVector z) {
    int px = x.length();
    int py = y.length();
    int pz = z.length();
    NumericVector TP_vec(px*py*pz);
    int vx; int vy; int vz;
    for (int v=0; v<(px*py*pz); ++v) {
        vx = v%px;
        vy = (v/px)%py;
        vz = (v/(px*py))%pz;
        TP_vec(v) = x(vx)*y(vy)*z(vz);
    }
    return TP_vec;
}

// [[Rcpp::export]]
NumericVector TP_rankR_cpp(List marg) {
    int D = marg.length();
    NumericMatrix x = marg[0];
    int R = x.ncol();
    int V = 1;
    IntegerVector p(D);
    for (int d=0; d<D; ++d) {
        NumericMatrix z = marg[d];
        p(d) = z.nrow();
        V = V*z.nrow();
    }
    NumericVector TP_vec(V);
    int vx; int vy; int vz;
    for (int r=0; r<R; ++r) {
        for (int v=0; v<V; ++v) {
            if (D==2) { // 2D case
                vx = v%p(0);
                vy = (v/p(0))%p(1);
                NumericMatrix marg_x = marg[0];
                NumericMatrix marg_y = marg[1];
                TP_vec(v)=TP_vec(v)+marg_x(vx,r)*marg_y(vy,r);
            } else if (D==3) { // 3D case
                vx = v%p(0);
                vy = (v/p(0))%p(1);
                vz = (v/(p(0)*p(1)))%p(2);
                NumericMatrix marg_x = marg[0];
                NumericMatrix marg_y = marg[1];
                NumericMatrix marg_z = marg[2];
                TP_vec(v)=TP_vec(v)+marg_x(vx,r)*marg_y(vy,r)*marg_z(vz,r);
            }
        }
    }
    return TP_vec;
}

// [[Rcpp::export]]
NumericVector getGamma_cpp(List gamma) {
    NumericVector Gamma = TP_rankR_cpp(gamma);
    return Gamma;
}

// [[Rcpp::export]]
List getGamma_multi_cpp(List gamma_multi) {
    int Q = gamma_multi.length();
    List Gamma_multi(Q);
    for (int q=0; q<Q; ++q) {
        List gamma_q = gamma_multi[q];
        NumericVector Gamma_q = getGamma_cpp(gamma_q);
        Gamma_multi[q] = Gamma_q;
    }
    return Gamma_multi;
}

// [[Rcpp::export]]
NumericMatrix getSigSq_cpp(NumericMatrix ssq, IntegerMatrix Zeta) {
    int nobs = Zeta.nrow();
    int V = Zeta.ncol();
    NumericMatrix SigSq(nobs,V);
    for (int i=0; i<nobs; ++i) {
        for (int v=0; v<V; ++v) {
            SigSq(i,v) = ssq(i,Zeta(i,v)-1);
        }
    }
    return SigSq;
}

// [[Rcpp::export]]
NumericMatrix getSigSq_harm_cpp(IntegerVector Scanner_ID, NumericMatrix ssq, IntegerMatrix Zeta) {
    //int nScanners = max(Scanner_ID);
    int nobs = Scanner_ID.length();
    int V = Zeta.ncol();
    NumericMatrix SigSq(nobs,V);
    for (int i=0; i<nobs; ++i) {
        int scan = Scanner_ID(i);
        for (int v=0; v<V; ++v) {
            SigSq(i,v) = ssq(scan-1, Zeta(scan-1,v)-1);
        }
    }
    return SigSq;
}

// [[Rcpp::export]]
NumericMatrix getSigSq_harm_cpp2(IntegerVector Scanner_ID, NumericMatrix ssq, IntegerMatrix Zeta) {
  int nScanners = max(Scanner_ID);
  //int nobs = Scanner_ID.length();
  int V = Zeta.ncol();
  NumericMatrix SigSq(nScanners,V);
  for (int i=0; i<nScanners; ++i) {
    for (int v=0; v<V; ++v) {
      SigSq(i,v) = ssq(i, Zeta(i,v)-1);
    }
  }
  return SigSq;
}

// [[Rcpp::export]]
List getSigSq_harm_mcmc_cpp(IntegerVector Scanner_ID, List ssq_store, List Zeta_store) {
  int niter = ssq_store.length();
  List SigSq_mcmc(niter);
  for (int iter=0; iter<niter; ++iter) {
    NumericMatrix ssq_iter = ssq_store[iter];
    IntegerMatrix Zeta_iter = Zeta_store[iter];
    SigSq_mcmc[iter] = getSigSq_harm_cpp(Scanner_ID, ssq_iter, Zeta_iter);
  }
  return SigSq_mcmc;
}

// [[Rcpp::export]]
NumericMatrix getYr_cpp(NumericMatrix Yvec, NumericVector X, List gamma, int r) {
    NumericMatrix Yrvec = clone(Yvec);
    int nobs = Yvec.nrow();
    int D = gamma.length();
    NumericMatrix gamma0 = gamma[0];
    int R = gamma0.ncol();
    List gamma_nor = clone(gamma);
    if (R>1) {
        for (int d=0; d<D; ++d) {
            NumericMatrix gamma_nor_d = gamma_nor[d];
            for (int vd=0; vd<gamma_nor_d.nrow(); ++vd) {
                gamma_nor_d(vd,r) = 0;
            }
            gamma_nor[d] = gamma_nor_d;
        }
        NumericVector Gamma_nor = getGamma_cpp(gamma_nor);
        for (int i=0; i<nobs; ++i) {
            Yrvec(i,_) = Yvec(i,_)-Gamma_nor*X(i);
        }
    }
    return Yrvec;
}
    
// [[Rcpp::export]]
NumericMatrix getYq_multi_cpp(NumericMatrix Yvec, NumericMatrix X_multi, List Gamma_multi, int q) {
    NumericMatrix Yqvec = clone(Yvec);
    int nobs = Yvec.nrow();
    int Q = X_multi.ncol();
    for (int qq=0; qq<Q; ++qq) {
        if (qq!=q) {
            NumericVector Gamma_qq = Gamma_multi[qq];
            for (int i=0; i<nobs; ++i) {
                Yqvec(i,_) = Yqvec(i,_)-Gamma_qq*X_multi(i,qq);
            }
        }
    }
    return Yqvec;
}

// [[Rcpp::export]]
NumericVector sample_gamma_dr_cpp(NumericMatrix Yrvec, IntegerVector p, NumericVector X, List gamma, double tau, NumericMatrix w, NumericMatrix alpha, NumericMatrix SigSq, int d, int r) {
    int nobs = Yrvec.nrow();
    int D = p.length();
    NumericMatrix gamma_d = gamma[d];
    NumericVector gamma_dr = gamma_d(_,r);
    NumericVector gamma_dr_new(p(d));
    //NumericVector c_dr(p(d));
    //NumericVector d_dr(p(d));
    double w_dr = w(d,r);
    double alpha_dr = alpha(d,r);
    NumericVector outer_dr;
    if (D==3) {
        NumericMatrix gamma_x = gamma[0];
        NumericVector gamma_xr = gamma_x(_,r);
        NumericMatrix gamma_y = gamma[1];
        NumericVector gamma_yr = gamma_y(_,r);
        NumericMatrix gamma_z = gamma[2];
        NumericVector gamma_zr = gamma_z(_,r);
        if (d==0) {
            outer_dr = TP_2D_cpp(gamma_yr, gamma_zr);
        } else if (d==1) {
            outer_dr = TP_2D_cpp(gamma_xr, gamma_zr);
        } else if (d==2) {
            outer_dr = TP_2D_cpp(gamma_xr, gamma_yr);
        }
    } else if (D==2) {
        NumericMatrix gamma_x = gamma[0];
        NumericVector gamma_xr = gamma_x(_,r);
        NumericMatrix gamma_y = gamma[1];
        NumericVector gamma_yr = gamma_y(_,r);
        if (d==0) {
            outer_dr = gamma_yr;
        } else if (d==1) {
            outer_dr = gamma_xr;
        }
    }
    for (int vd=0; vd<p(d); ++vd) {
        double c_drv = 0;
        double d_drv = 0;
        double m_drv = 0;
        double n_drv = 0;
        for (int i=0; i<nobs; ++i) {
            NumericVector SigSq_i = SigSq(i,_);
            NumericVector Yrvec_i = Yrvec(i,_);
            NumericVector SigSq_sl_div = slice_tensor_cpp(SigSq_i,p,d,vd);
            NumericVector Y_sl_driv = slice_tensor_cpp(Yrvec_i,p,d,vd);
            LogicalVector Ysl_isNA = is_na(Y_sl_driv);
            for (int j=0; j<outer_dr.length(); ++j) {
                if (Ysl_isNA(j)==FALSE) {
                    m_drv = m_drv+pow(X(i)*outer_dr(j),2)/SigSq_sl_div(j);
                    n_drv = n_drv+(X(i)*Y_sl_driv(j)*outer_dr(j))/SigSq_sl_div(j);
                }
            }
        }
        if (vd==0) {
            c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
            d_drv = n_drv+(exp(-alpha_dr)*gamma_dr(vd+1))/(tau*w_dr*(1-exp(-2*alpha_dr)));
        } else if (vd==(p(d)-1)) {
            c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
            d_drv = n_drv+(exp(-alpha_dr)*gamma_dr(vd-1))/(tau*w_dr*(1-exp(-2*alpha_dr)));
        } else {
            c_drv = m_drv+(1+exp(-2*alpha_dr))/(tau*w_dr*(1-exp(-2*alpha_dr)));
            d_drv = n_drv+(exp(-alpha_dr)*(gamma_dr(vd-1)+gamma_dr(vd+1)))/(tau*w_dr*(1-exp(-2*alpha_dr)));
        }
        
        if (c_drv==0) {
            gamma_dr_new(vd) = 0;
        } else {
            gamma_dr_new(vd) = R::rnorm(d_drv/c_drv, pow(1.0/c_drv,0.5));
        }
        //c_dr(vd) = c_drv;
        //d_dr(vd) = d_drv;
    }
    return gamma_dr_new;
    //List res = List::create(c_dr,d_dr);
    //return res;
}

// [[Rcpp::export]]
NumericVector sample_gamma_dr_harm_cpp(NumericMatrix Yrvec, IntegerVector p, IntegerVector Scanner_ID, NumericVector X, List gamma, double tau, NumericMatrix w, NumericMatrix alpha, NumericMatrix SigSq, int d, int r) {
  int nobs = Yrvec.nrow();
  int D = p.length();
  NumericMatrix gamma_d = gamma[d];
  NumericVector gamma_dr = gamma_d(_,r);
  NumericVector gamma_dr_new(p(d));
  //NumericVector c_dr(p(d));
  //NumericVector d_dr(p(d));
  double w_dr = w(d,r);
  double alpha_dr = alpha(d,r);
  NumericVector outer_dr;
  if (D==3) {
    NumericMatrix gamma_x = gamma[0];
    NumericVector gamma_xr = gamma_x(_,r);
    NumericMatrix gamma_y = gamma[1];
    NumericVector gamma_yr = gamma_y(_,r);
    NumericMatrix gamma_z = gamma[2];
    NumericVector gamma_zr = gamma_z(_,r);
    if (d==0) {
      outer_dr = TP_2D_cpp(gamma_yr, gamma_zr);
    } else if (d==1) {
      outer_dr = TP_2D_cpp(gamma_xr, gamma_zr);
    } else if (d==2) {
      outer_dr = TP_2D_cpp(gamma_xr, gamma_yr);
    }
  } else if (D==2) {
    NumericMatrix gamma_x = gamma[0];
    NumericVector gamma_xr = gamma_x(_,r);
    NumericMatrix gamma_y = gamma[1];
    NumericVector gamma_yr = gamma_y(_,r);
    if (d==0) {
      outer_dr = gamma_yr;
    } else if (d==1) {
      outer_dr = gamma_xr;
    }
  }
  for (int vd=0; vd<p(d); ++vd) {
    double c_drv = 0;
    double d_drv = 0;
    double m_drv = 0;
    double n_drv = 0;
    for (int i=0; i<nobs; ++i) {
      NumericVector SigSq_i = SigSq(Scanner_ID(i)-1,_);
      NumericVector Yrvec_i = Yrvec(i,_);
      NumericVector SigSq_sl_div = slice_tensor_cpp(SigSq_i,p,d,vd);
      NumericVector Y_sl_driv = slice_tensor_cpp(Yrvec_i,p,d,vd);
      LogicalVector Ysl_isNA = is_na(Y_sl_driv);
      for (int j=0; j<outer_dr.length(); ++j) {
        if (Ysl_isNA(j)==FALSE) {
          m_drv = m_drv+pow(X(i)*outer_dr(j),2)/SigSq_sl_div(j);
          n_drv = n_drv+(X(i)*Y_sl_driv(j)*outer_dr(j))/SigSq_sl_div(j);
        }
      }
    }
    if (vd==0) {
      c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*gamma_dr(vd+1))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    } else if (vd==(p(d)-1)) {
      c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*gamma_dr(vd-1))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    } else {
      c_drv = m_drv+(1+exp(-2*alpha_dr))/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*(gamma_dr(vd-1)+gamma_dr(vd+1)))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    }
    
    if (c_drv==0) {
      gamma_dr_new(vd) = 0;
    } else {
      gamma_dr_new(vd) = R::rnorm(d_drv/c_drv, pow(1.0/c_drv,0.5));
    }
    //c_dr(vd) = c_drv;
    //d_dr(vd) = d_drv;
  }
  return gamma_dr_new;
  //List res = List::create(c_dr,d_dr);
  //return res;
}

// [[Rcpp::export]]
NumericVector sample_gamma_dr_harm_cpp2(NumericMatrix Yrvec, IntegerVector p, IntegerVector Scanner_ID, NumericVector X, List gamma, double tau, NumericMatrix w, NumericMatrix alpha, NumericMatrix SigSq, int d, int r) {
  int nobs = Yrvec.nrow();
  int D = p.length();
  NumericMatrix gamma_d = gamma[d];
  NumericVector gamma_dr = gamma_d(_,r);
  NumericVector gamma_dr_new(p(d));
  //NumericVector c_dr(p(d));
  //NumericVector d_dr(p(d));
  double w_dr = w(d,r);
  double alpha_dr = alpha(d,r);
  NumericVector outer_dr;
  if (D==3) {
    NumericMatrix gamma_x = gamma[0];
    NumericVector gamma_xr = gamma_x(_,r);
    NumericMatrix gamma_y = gamma[1];
    NumericVector gamma_yr = gamma_y(_,r);
    NumericMatrix gamma_z = gamma[2];
    NumericVector gamma_zr = gamma_z(_,r);
    if (d==0) {
      outer_dr = TP_2D_cpp(gamma_yr, gamma_zr);
    } else if (d==1) {
      outer_dr = TP_2D_cpp(gamma_xr, gamma_zr);
    } else if (d==2) {
      outer_dr = TP_2D_cpp(gamma_xr, gamma_yr);
    }
  } else if (D==2) {
    NumericMatrix gamma_x = gamma[0];
    NumericVector gamma_xr = gamma_x(_,r);
    NumericMatrix gamma_y = gamma[1];
    NumericVector gamma_yr = gamma_y(_,r);
    if (d==0) {
      outer_dr = gamma_yr;
    } else if (d==1) {
      outer_dr = gamma_xr;
    }
  }
  for (int vd=0; vd<p(d); ++vd) {
    double c_drv = 0;
    double d_drv = 0;
    double m_drv = 0;
    double n_drv = 0;
    for (int i=0; i<nobs; ++i) {
      NumericVector SigSq_i = SigSq(Scanner_ID(i)-1,_);
      NumericVector Yrvec_i = Yrvec(i,_);
      //NumericVector SigSq_sl_div = slice_tensor_cpp(SigSq_i,p,d,vd);
      //NumericVector Y_sl_driv = slice_tensor_cpp(Yrvec_i,p,d,vd);
      IntegerVector slice_inds_dvd = getSliceIndices_iota_cpp(p, d, vd);
      NumericVector SigSq_sl_div = SigSq_i[slice_inds_dvd];
      NumericVector Y_sl_driv = Yrvec_i[slice_inds_dvd];
      LogicalVector Ysl_isNA = is_na(Y_sl_driv);
      for (int j=0; j<outer_dr.length(); ++j) {
        if (Ysl_isNA(j)==FALSE) {
          m_drv = m_drv+pow(X(i)*outer_dr(j),2)/SigSq_sl_div(j);
          n_drv = n_drv+(X(i)*Y_sl_driv(j)*outer_dr(j))/SigSq_sl_div(j);
        }
      }
    }
    if (vd==0) {
      c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*gamma_dr(vd+1))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    } else if (vd==(p(d)-1)) {
      c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*gamma_dr(vd-1))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    } else {
      c_drv = m_drv+(1+exp(-2*alpha_dr))/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*(gamma_dr(vd-1)+gamma_dr(vd+1)))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    }
    
    if (c_drv==0) {
      gamma_dr_new(vd) = 0;
    } else {
      gamma_dr_new(vd) = R::rnorm(d_drv/c_drv, pow(1.0/c_drv,0.5));
    }
    //c_dr(vd) = c_drv;
    //d_dr(vd) = d_drv;
  }
  return gamma_dr_new;
  //List res = List::create(c_dr,d_dr);
  //return res;
}

// [[Rcpp::export]]
NumericVector sample_gamma_dr_harm_cpp3(NumericMatrix Yrvec, IntegerVector p, IntegerVector Scanner_ID, NumericVector X, List gamma, double tau, NumericMatrix w, NumericMatrix alpha, NumericMatrix SigSq, int d, int r) {
  int nobs = Yrvec.nrow();
  int D = p.length();
  NumericMatrix gamma_d = gamma[d];
  NumericVector gamma_dr = gamma_d(_,r);
  NumericVector gamma_dr_new(p(d));
  //NumericVector c_dr(p(d));
  //NumericVector d_dr(p(d));
  double w_dr = w(d,r);
  double alpha_dr = alpha(d,r);
  NumericVector outer_dr;
  if (D==3) {
    NumericMatrix gamma_x = gamma[0];
    NumericVector gamma_xr = gamma_x(_,r);
    NumericMatrix gamma_y = gamma[1];
    NumericVector gamma_yr = gamma_y(_,r);
    NumericMatrix gamma_z = gamma[2];
    NumericVector gamma_zr = gamma_z(_,r);
    if (d==0) {
      outer_dr = TP_2D_cpp(gamma_yr, gamma_zr);
    } else if (d==1) {
      outer_dr = TP_2D_cpp(gamma_xr, gamma_zr);
    } else if (d==2) {
      outer_dr = TP_2D_cpp(gamma_xr, gamma_yr);
    }
  } else if (D==2) {
    NumericMatrix gamma_x = gamma[0];
    NumericVector gamma_xr = gamma_x(_,r);
    NumericMatrix gamma_y = gamma[1];
    NumericVector gamma_yr = gamma_y(_,r);
    if (d==0) {
      outer_dr = gamma_yr;
    } else if (d==1) {
      outer_dr = gamma_xr;
    }
  }
  for (int vd=0; vd<p(d); ++vd) {
    double c_drv = 0;
    double d_drv = 0;
    double m_drv = 0;
    double n_drv = 0;
    for (int i=0; i<nobs; ++i) {
      NumericVector SigSq_i = SigSq(Scanner_ID(i)-1,_);
      NumericVector Yrvec_i = Yrvec(i,_);
      //NumericVector SigSq_sl_div = slice_tensor_cpp(SigSq_i,p,d,vd);
      //NumericVector Y_sl_driv = slice_tensor_cpp(Yrvec_i,p,d,vd);
      IntegerVector slice_inds_dvd = getSliceIndices_loop_cpp(p, d, vd);
      NumericVector SigSq_sl_div = SigSq_i[slice_inds_dvd];
      NumericVector Y_sl_driv = Yrvec_i[slice_inds_dvd];
      LogicalVector Ysl_isNA = is_na(Y_sl_driv);
      for (int j=0; j<outer_dr.length(); ++j) {
        if (Ysl_isNA(j)==FALSE) {
          m_drv = m_drv+pow(X(i)*outer_dr(j),2)/SigSq_sl_div(j);
          n_drv = n_drv+(X(i)*Y_sl_driv(j)*outer_dr(j))/SigSq_sl_div(j);
        }
      }
    }
    if (vd==0) {
      c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*gamma_dr(vd+1))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    } else if (vd==(p(d)-1)) {
      c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*gamma_dr(vd-1))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    } else {
      c_drv = m_drv+(1+exp(-2*alpha_dr))/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*(gamma_dr(vd-1)+gamma_dr(vd+1)))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    }
    
    if (c_drv==0) {
      gamma_dr_new(vd) = 0;
    } else {
      gamma_dr_new(vd) = R::rnorm(d_drv/c_drv, pow(1.0/c_drv,0.5));
    }
    //c_dr(vd) = c_drv;
    //d_dr(vd) = d_drv;
  }
  return gamma_dr_new;
  //List res = List::create(c_dr,d_dr);
  //return res;
}

// [[Rcpp::export]]
NumericVector sample_gamma_dr_harm_withInds_cpp(NumericMatrix Yrvec, IntegerVector p, IntegerVector Scanner_ID, NumericVector X, List gamma, double tau, NumericMatrix w, NumericMatrix alpha, NumericMatrix SigSq, int d, int r, IntegerMatrix slice_inds_d) {
  int nobs = Yrvec.nrow();
  int D = p.length();
  NumericMatrix gamma_d = gamma[d];
  NumericVector gamma_dr = gamma_d(_,r);
  NumericVector gamma_dr_new(p(d));
  //NumericVector c_dr(p(d));
  //NumericVector d_dr(p(d));
  double w_dr = w(d,r);
  double alpha_dr = alpha(d,r);
  NumericVector outer_dr;
  if (D==3) {
    NumericMatrix gamma_x = gamma[0];
    NumericVector gamma_xr = gamma_x(_,r);
    NumericMatrix gamma_y = gamma[1];
    NumericVector gamma_yr = gamma_y(_,r);
    NumericMatrix gamma_z = gamma[2];
    NumericVector gamma_zr = gamma_z(_,r);
    if (d==0) {
      outer_dr = TP_2D_cpp(gamma_yr, gamma_zr);
    } else if (d==1) {
      outer_dr = TP_2D_cpp(gamma_xr, gamma_zr);
    } else if (d==2) {
      outer_dr = TP_2D_cpp(gamma_xr, gamma_yr);
    }
  } else if (D==2) {
    NumericMatrix gamma_x = gamma[0];
    NumericVector gamma_xr = gamma_x(_,r);
    NumericMatrix gamma_y = gamma[1];
    NumericVector gamma_yr = gamma_y(_,r);
    if (d==0) {
      outer_dr = gamma_yr;
    } else if (d==1) {
      outer_dr = gamma_xr;
    }
  }
  for (int vd=0; vd<p(d); ++vd) {
    IntegerVector slice_inds_dvd = slice_inds_d(vd,_);
    double c_drv = 0;
    double d_drv = 0;
    double m_drv = 0;
    double n_drv = 0;
    for (int i=0; i<nobs; ++i) {
      NumericVector SigSq_i = SigSq(Scanner_ID(i)-1,_);
      NumericVector Yrvec_i = Yrvec(i,_);
      //NumericVector SigSq_sl_div = slice_tensor_cpp(SigSq_i,p,d,vd);
      //NumericVector Y_sl_driv = slice_tensor_cpp(Yrvec_i,p,d,vd);
      //IntegerVector slice_inds_dvd = getSliceIndices_loop_cpp(p, d, vd);
      NumericVector SigSq_sl_div = SigSq_i[slice_inds_dvd];
      NumericVector Y_sl_driv = Yrvec_i[slice_inds_dvd];
      LogicalVector Ysl_isNA = is_na(Y_sl_driv);
      for (int j=0; j<outer_dr.length(); ++j) {
        if (Ysl_isNA(j)==FALSE) {
          m_drv = m_drv+pow(X(i)*outer_dr(j),2)/SigSq_sl_div(j);
          n_drv = n_drv+(X(i)*Y_sl_driv(j)*outer_dr(j))/SigSq_sl_div(j);
        }
      }
    }
    if (vd==0) {
      c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*gamma_dr(vd+1))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    } else if (vd==(p(d)-1)) {
      c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*gamma_dr(vd-1))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    } else {
      c_drv = m_drv+(1+exp(-2*alpha_dr))/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*(gamma_dr(vd-1)+gamma_dr(vd+1)))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    }
    
    if (c_drv==0) {
      gamma_dr_new(vd) = 0;
    } else {
      gamma_dr_new(vd) = R::rnorm(d_drv/c_drv, pow(1.0/c_drv,0.5));
    }
    //c_dr(vd) = c_drv;
    //d_dr(vd) = d_drv;
  }
  return gamma_dr_new;
  //List res = List::create(c_dr,d_dr);
  //return res;
}

// [[Rcpp::export]]
NumericVector sample_gamma_dr_harm_withInds2_cpp(NumericMatrix Yrvec, IntegerVector p, IntegerVector Scanner_ID, NumericVector X, List gamma, double tau, NumericMatrix w, NumericMatrix alpha, NumericMatrix SigSq, int d, int r, List slice_inds_d) {
  int nobs = Yrvec.nrow();
  int D = p.length();
  //int V = Yrvec.ncol();
  NumericMatrix gamma_d = gamma[d];
  NumericVector gamma_dr = gamma_d(_,r);
  NumericVector gamma_dr_new(p(d));
  //NumericVector mean_dr(p(d));
  //NumericVector sd_dr(p(d));
  
  double w_dr = w(d,r);
  double alpha_dr = alpha(d,r);
  
  NumericVector Gamma_r;
  if (D==3) {
    NumericMatrix gamma_x = gamma[0];
    NumericVector gamma_xr = gamma_x(_,r);
    NumericMatrix gamma_y = gamma[1];
    NumericVector gamma_yr = gamma_y(_,r);
    NumericMatrix gamma_z = gamma[2];
    NumericVector gamma_zr = gamma_z(_,r);
    Gamma_r = TP_3D_cpp(gamma_xr, gamma_yr, gamma_zr);
  } else if (D==2) {
    NumericMatrix gamma_x = gamma[0];
    NumericVector gamma_xr = gamma_x(_,r);
    NumericMatrix gamma_y = gamma[1];
    NumericVector gamma_yr = gamma_y(_,r);
    Gamma_r = TP_2D_cpp(gamma_xr, gamma_yr);
  }

  for (int vd=0; vd<p(d); ++vd) {
    IntegerVector slice_inds_dvd = slice_inds_d[vd];
    int L_vd = slice_inds_dvd.length();
    
    double c_drv = 0;
    double d_drv = 0;
    double m_drv = 0;
    double n_drv = 0;
    
    if (L_vd > 0) {
      NumericVector Gamma_rvd = Gamma_r[slice_inds_dvd];
      Gamma_rvd = Gamma_rvd/gamma_dr(vd);
      
      for (int i=0; i<nobs; ++i) {
        NumericVector SigSq_i = SigSq(Scanner_ID(i)-1,_);
        NumericVector Yrvec_i = Yrvec(i,_);
        NumericVector SigSq_sl_div = SigSq_i[slice_inds_dvd];
        NumericVector Y_sl_driv = Yrvec_i[slice_inds_dvd];
        for (int j=0; j<L_vd; ++j) {
          m_drv = m_drv + pow(X(i)*Gamma_rvd(j),2)/SigSq_sl_div(j);
          n_drv = n_drv + (X(i)*Y_sl_driv(j)*Gamma_rvd(j))/SigSq_sl_div(j);
        }
      }
    }
      
    if (vd==0) {
      c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*gamma_dr(vd+1))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    } else if (vd==(p(d)-1)) {
      c_drv = m_drv+1/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*gamma_dr(vd-1))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    } else {
      c_drv = m_drv+(1+exp(-2*alpha_dr))/(tau*w_dr*(1-exp(-2*alpha_dr)));
      d_drv = n_drv+(exp(-alpha_dr)*(gamma_dr(vd-1)+gamma_dr(vd+1)))/(tau*w_dr*(1-exp(-2*alpha_dr)));
    }
    
    if (c_drv==0) {
     gamma_dr_new(vd) = 0;
    } else {
     gamma_dr_new(vd) = R::rnorm(d_drv/c_drv, pow(1.0/c_drv,0.5));
    }
    // if (c_drv==0) {
    //   mean_dr(vd) = 0;
    //   sd_dr(vd) = 0;
    // } else {
    //   mean_dr(vd) = d_drv/c_drv;
    //   sd_dr(vd) = pow(1.0/c_drv,.5);
    // }

  }
  
  //List results = List::create(mean_dr, sd_dr);
  //return results;
  return gamma_dr_new;
}


// [[Rcpp::export]]
double sample_gamma_u_eqvox_harm_cpp(NumericMatrix Yvec, IntegerVector Scanner_ID, double tau, NumericMatrix SigSq) {
  int nobs = Yvec.nrow();
  int V = Yvec.ncol();
  double gamma;
  double mu;
  double rho;
  double numer = 0;
  double denom = 1.0/tau;
  for (int i=0; i<nobs; ++i) {
    NumericVector SigSq_i = SigSq(Scanner_ID(i)-1,_);
    NumericVector Yvec_i = Yvec(i,_);
    LogicalVector Ysl_isNA = is_na(Yvec_i);
    for (int v=0; v<V; ++v) {
      if (Ysl_isNA(v)==FALSE) {
        numer = numer + Yvec_i(v)/SigSq_i(v);
        denom = denom + 1.0/SigSq_i(v);
      }
    }
  }
  mu = numer/denom;
  if (denom==0) {
    gamma = 0;
  } else {
    rho = pow(1.0/denom,0.5);
    gamma = R::rnorm(mu, rho);
  }
  return gamma;
}

// [[Rcpp::export]]
double get_gamma_sum(NumericVector gamma_dr, double alpha_dr) {
    int pd = gamma_dr.length();
    double sum = 0;
    for (int vv=0; vv<pd; ++vv) {
        if (vv==0 || vv==(pd-1)) {
            sum = sum + pow(gamma_dr(vv),2);
        } else {
            sum = sum + (1+exp(-2*alpha_dr))*pow(gamma_dr(vv),2);
        }
        if (vv<(pd-1)) {
            sum = sum - 2*exp(-alpha_dr)*gamma_dr(vv)*gamma_dr(vv+1);
        }
    }
    return sum;
}

// [[Rcpp::export]]
double sample_w_dr_cpp(List gamma, double tau, NumericMatrix alpha, NumericMatrix lambda, IntegerVector p, int d, int r) {
    double w_dr;
    double alpha_dr = alpha(d,r);
    double lambda_dr = lambda(d,r);
    double c_dr = 0;
    NumericMatrix gamma_d = gamma[d];
    NumericVector gamma_dr = gamma_d(_,r);
    c_dr = get_gamma_sum(gamma_dr,alpha_dr)/(tau*(1-exp(-2*alpha_dr)));
    w_dr = do_rgig(1-.5*p(d), c_dr, lambda_dr);
    return w_dr;
    //return c_dr;
}

// [[Rcpp::export]]
double sample_w_dr_multi_cpp(List gamma_multi, NumericVector tau_multi, NumericMatrix alpha, NumericMatrix lambda, IntegerVector p, int d, int r) {
    int Q = tau_multi.length();
    double w_dr;
    double alpha_dr = alpha(d,r);
    double lambda_dr = lambda(d,r);
    double c_dr = 0;
    for (int q=0; q<Q; ++q) {
        List gamma_q = gamma_multi[q];
        NumericMatrix gamma_dq = gamma_q[d];
        NumericVector gamma_drq = gamma_dq(_,r);
        c_dr = c_dr + get_gamma_sum(gamma_drq,alpha_dr)/(tau_multi(q)*(1-exp(-2*alpha_dr)));
    }
    w_dr = do_rgig(1-.5*p(d)*Q,c_dr,lambda_dr);
    return w_dr;
}
    
// [[Rcpp::export]]
double sample_lambda_dr_cpp(NumericMatrix w, IntegerVector p, double a_lam, double b_lam, int d, int r) {
    double lambda_dr = R::rgamma(a_lam+p(d), 1.0/(b_lam+.5*p(d)*w(d,r)));
    return lambda_dr;
}

// [[Rcpp::export]]
double get_log_Palpha_dr_cpp(double alpha_dr, List gamma, double tau, NumericMatrix w, IntegerVector p, double a_alpha, double b_alpha, int d, int r) {
    double w_dr = w(d,r);
    NumericMatrix gamma_d = gamma[d];
    NumericVector gamma_dr = gamma_d(_,r);
    double logy = (a_alpha-1)*log(alpha_dr)-.5*(p(d)-1)*log(1-exp(-2*alpha_dr))-b_alpha*alpha_dr-get_gamma_sum(gamma_dr, alpha_dr)/(2*tau*w_dr*(1-exp(-2*alpha_dr)));
    return logy;
}

// [[Rcpp::export]]
double get_log_Palpha_dr_multi_cpp(double alpha_dr, List gamma_multi, NumericVector tau_multi, NumericMatrix w, IntegerVector p, double a_alpha, double b_alpha, int d, int r) {
    int Q = tau_multi.length();
    double w_dr = w(d,r);
    double logy = (a_alpha-1)*log(alpha_dr)-.5*Q*(p(d)-1)*log(1-exp(-2*alpha_dr))-b_alpha*alpha_dr;
    for (int q=0; q<Q; ++q) {
        List gamma_q = gamma_multi[q];
        NumericMatrix gamma_dq = gamma_q[d];
        NumericVector gamma_drq = gamma_dq(_,r);
        logy = logy - get_gamma_sum(gamma_drq,alpha_dr)/(2*tau_multi(q)*w_dr*(1-exp(-2*alpha_dr)));
    }
    return logy;
}

// [[Rcpp::export]]
double sample_alpha_dr_cpp(NumericMatrix alpha, List gamma, double tau, NumericMatrix w, IntegerVector p, double a_alpha, double b_alpha, double sigma_log_alpha, int d, int r) {
    double alpha_dr = alpha(d,r);
    double alpha_dr_cand = exp(R::rnorm(log(alpha_dr),sigma_log_alpha));
    double alpha_dr_new;
    double acc_ratio = exp(get_log_Palpha_dr_cpp(alpha_dr_cand, gamma, tau, w, p, a_alpha, b_alpha, d, r)-get_log_Palpha_dr_cpp(alpha_dr, gamma, tau, w, p, a_alpha, b_alpha, d, r));
    double xunif = R::runif(0,1);
    if (xunif > acc_ratio) {
        alpha_dr_new = alpha_dr;
    } else {
        alpha_dr_new = alpha_dr_cand;
    }
    return alpha_dr_new;
}

// [[Rcpp::export]]
double sample_alpha_dr_multi_cpp(NumericMatrix alpha, List gamma_multi, NumericVector tau_multi, NumericMatrix w, IntegerVector p, double a_alpha, double b_alpha, double sigma_log_alpha, int d, int r) {
    double alpha_dr = alpha(d,r);
    double alpha_dr_cand = exp(R::rnorm(log(alpha_dr),sigma_log_alpha));
    double alpha_dr_new;
    double acc_ratio = exp(get_log_Palpha_dr_multi_cpp(alpha_dr_cand, gamma_multi, tau_multi, w, p, a_alpha, b_alpha, d, r)-get_log_Palpha_dr_multi_cpp(alpha_dr, gamma_multi, tau_multi, w, p, a_alpha, b_alpha, d, r));
    double xunif = R::runif(0,1);
    if (xunif > acc_ratio) {
        alpha_dr_new = alpha_dr;
    } else {
        alpha_dr_new = alpha_dr_cand;
    }
    return alpha_dr_new;
}

// [[Rcpp::export]]
double sample_tau_cpp(List gamma, NumericMatrix w, NumericMatrix alpha, double a_tau, double b_tau) {
    int D = w.nrow();
    int R = w.ncol();
    double sum_tau = 0;
    double sum_p = 0;
    for (int d=0; d<D; ++d) {
        NumericMatrix gamma_d = gamma[d];
        sum_p = sum_p + gamma_d.nrow();
        for (int r=0; r<R; ++r) {
            NumericVector gamma_dr = gamma_d(_,r);
            sum_tau = sum_tau + get_gamma_sum(gamma_dr, alpha(d,r))/(w(d,r)*(1-exp(-2*alpha(d,r))));
        }
    }
    double tau = do_rgig(a_tau-.5*R*sum_p,sum_tau,2*b_tau);
    return tau;
    //return sum_tau;
}

// [[Rcpp::export]]
NumericVector sample_tau_multi_cpp(List gamma_multi, NumericMatrix w, NumericMatrix alpha, double a_tau, double b_tau) {
    int Q = gamma_multi.length();
    int D = w.nrow();
    int R = w.ncol();
    NumericVector tau_multi(Q);
    for (int q=0; q<Q; ++q) {
        double sum_q, sum_p;
        List gamma_q = gamma_multi[q];
        for (int d=0; d<D; ++d) {
            NumericMatrix gamma_dq = gamma_q[d];
            sum_p = sum_p + gamma_dq.nrow();
            for (int r=0; r<R; ++r) {
                NumericVector gamma_drq = gamma_dq(_,r);
                sum_q = sum_q + get_gamma_sum(gamma_drq, alpha(d,r))/(w(d,r)*(1-exp(-2*alpha(d,r))));
            }
        }
        tau_multi(q) = do_rgig(a_tau-.5*R*sum_p,sum_q,2*b_tau);
    }
    return tau_multi;
}

// [[Rcpp::export]]
double sample_tau_eqvox_cpp(NumericVector gamma, double a_tau, double b_tau) {
  double tau;
  int nSubj = gamma.length();
  double aprime = a_tau + 0.5*nSubj;
  double sum_gamma_sq = 0;
  for (int u=0; u<nSubj; ++u) {
    sum_gamma_sq = sum_gamma_sq + pow(gamma(u),2);
  }
  double bprime = b_tau + .5*sum_gamma_sq;
  tau = 1.0/R::rgamma(aprime, 1.0/bprime);
  return tau;
}

// [[Rcpp::export]]
NumericMatrix getResid_cpp(NumericMatrix Yvec, NumericVector X, NumericVector Gamma) {
    NumericMatrix Residvec = clone(Yvec);
    int nobs = Yvec.nrow();
    for (int i=0; i<nobs; ++i) {
        Residvec(i,_) = Yvec(i,_)-Gamma*X(i);
    }
    return Residvec;
}

// [[Rcpp::export]]
NumericMatrix getResid_multi_cpp(NumericMatrix Yvec, NumericMatrix X_multi, List Gamma_multi) {
    NumericMatrix Residvec = clone(Yvec);
    int Q = X_multi.ncol();
    int nobs = Yvec.nrow();
    for (int q=0; q<Q; ++q) {
        NumericVector Gamma_q = Gamma_multi[q];
        for (int i=0; i<nobs; ++i) {
            Residvec(i,_) = Residvec(i,_) - Gamma_q*X_multi(i,q);
        }
    }
    return Residvec;
}

// [[Rcpp::export]]
NumericMatrix getResid_harm_long_cpp(NumericMatrix Yvec, NumericMatrix X_cov, IntegerVector Scanner_ID, IntegerVector Subj_ID, NumericVector Gamma_pop, List Gamma_cov, List Gamma_scan, List Gamma_subj) {
  NumericMatrix Residvec = getResid_multi_cpp(Yvec, X_cov, Gamma_cov);
  int nobs = Yvec.nrow();
  for (int n=0; n<nobs; ++n) {
    NumericVector Gamma_scan_n = Gamma_scan[Scanner_ID(n)-1];
    NumericVector Gamma_subj_n = Gamma_subj[Subj_ID(n)-1];
    Residvec(n,_) = Residvec(n,_) - Gamma_scan_n - Gamma_subj_n - Gamma_pop;
  }
  return Residvec;
}

// [[Rcpp::export]]
NumericMatrix getY_intercept_long_cpp(NumericMatrix Yvec, NumericMatrix X_cov, IntegerVector Scanner_ID, IntegerVector Subj_ID, List Gamma_cov, List Gamma_scan, List Gamma_subj) {
    NumericMatrix Y_intercept = getResid_multi_cpp(Yvec, X_cov, Gamma_cov);
    int nobs = Yvec.nrow();
    for (int n=0; n<nobs; ++n) {
        NumericVector Gamma_scan_n = Gamma_scan[Scanner_ID(n)-1];
        NumericVector Gamma_subj_n = Gamma_subj[Subj_ID(n)-1];
        Y_intercept(n,_) = Y_intercept(n,_) - Gamma_scan_n - Gamma_subj_n;
    }
    return Y_intercept;
}

// [[Rcpp::export]]
NumericMatrix getY_cov_long_cpp(NumericMatrix Yvec, IntegerVector Scanner_ID, IntegerVector Subj_ID, NumericVector Gamma_pop, List Gamma_scan, List Gamma_subj) {
    NumericMatrix Y_cov = clone(Yvec);
    int nobs = Yvec.nrow();
    for (int n=0; n<nobs; ++n) {
        NumericVector Gamma_scan_n = Gamma_scan[Scanner_ID(n)-1];
        NumericVector Gamma_subj_n = Gamma_subj[Subj_ID(n)-1];
        Y_cov(n,_) = Y_cov(n,_) - Gamma_pop - Gamma_scan_n - Gamma_subj_n;
    }
    return Y_cov;
}

// [[Rcpp::export]]
NumericMatrix getY_scan_long_cpp(NumericMatrix Yvec, NumericMatrix X_cov, IntegerVector Scanner_ID, IntegerVector Subj_ID, NumericVector Gamma_pop, List Gamma_cov, List Gamma_subj, int scan) {
    int nobs = Yvec.nrow();
    int V = Yvec.ncol();
    int Q = X_cov.ncol();
    int num_scan_obs = 0;
    for (int n=0; n<nobs; ++n) {
        if (Scanner_ID(n)==scan) {
            num_scan_obs++;
        }
    }
    NumericMatrix Yscan(num_scan_obs,V);
    NumericMatrix X_cov_scan(num_scan_obs,Q);
    IntegerVector Subj_ID_scan(num_scan_obs);
    int ind_scan = 0;
    for (int n=0; n<nobs; ++n) {
        if (Scanner_ID(n)==scan) {
            Yscan(ind_scan,_) = Yvec(n,_);
            X_cov_scan(ind_scan,_) = X_cov(n,_);
            Subj_ID_scan(ind_scan) = Subj_ID(n);
            ind_scan++;
        }
    }
    NumericMatrix Residscan = getResid_multi_cpp(Yscan,X_cov_scan,Gamma_cov);
    for (int nn=0; nn<num_scan_obs; ++nn) {
        NumericVector Gamma_subj_nn = Gamma_subj[Subj_ID_scan(nn)-1];
        Residscan(nn,_) = Residscan(nn,_) - Gamma_pop - Gamma_subj_nn;
    }
    return Residscan;
}

// [[Rcpp::export]]
NumericMatrix getY_subj_long_cpp(NumericMatrix Yvec, NumericMatrix X_cov, IntegerVector Scanner_ID, IntegerVector Subj_ID, NumericVector Gamma_pop, List Gamma_cov, List Gamma_scan, int subj) {
  int nobs = Yvec.nrow();
  int V = Yvec.ncol();
  int Q = X_cov.ncol();
  int num_subj_obs = 0;
  for (int n=0; n<nobs; ++n) {
    if (Subj_ID(n)==subj) {
      num_subj_obs++;
    }
  }
  NumericMatrix Ysubj(num_subj_obs,V);
  NumericMatrix X_cov_subj(num_subj_obs,Q);
  IntegerVector Scanner_ID_subj(num_subj_obs);
  int ind_subj = 0;
  for (int n=0; n<nobs; ++n) {
    if (Subj_ID(n)==subj) {
      Ysubj(ind_subj,_) = Yvec(n,_);
      X_cov_subj(ind_subj,_) = X_cov(n,_);
      Scanner_ID_subj(ind_subj) = Scanner_ID(n);
      ind_subj++;
    }
  }
  NumericMatrix Residsubj = getResid_multi_cpp(Ysubj,X_cov_subj,Gamma_cov);
  for (int nn=0; nn<num_subj_obs; ++nn) {
    NumericVector Gamma_scan_nn = Gamma_scan[Scanner_ID_subj(nn)-1];
    Residsubj(nn,_) = Residsubj(nn,_) - Gamma_pop - Gamma_scan_nn;
  }
  return Residsubj;
}

// [[Rcpp::export]]
NumericMatrix getY_scan_noadj_cpp(NumericMatrix Yvec, IntegerVector Scanner_ID, int scan) {
    int nobs = Yvec.nrow();
    int V = Yvec.ncol();
    int num_scan_obs = 0;
    for (int n=0; n<nobs; ++n) {
        if (Scanner_ID(n)==scan) {
            num_scan_obs++;
        }
    }
    NumericMatrix Yscan(num_scan_obs, V);
    int ind_scan = 0;
    for (int n=0; n<nobs; ++n) {
        if (Scanner_ID(n)==scan) {
            Yscan(ind_scan,_) = Yvec(n,_);
            ind_scan++;
        }
    }
    return Yscan;
}

// [[Rcpp::export]]
IntegerVector sample_Zeta_i_cpp(NumericMatrix Residvec, NumericMatrix ssq, NumericMatrix rho, int i) {
    int V = Residvec.ncol();
    int H = ssq.ncol();
    NumericVector Residvec_i = Residvec(i,_);
    LogicalVector Residvec_i_isNA = is_na(Residvec_i);
    IntegerVector Zeta_i(V);
    if (H==1) {
        for (int v=0; v<V; ++v) {
            if (Residvec_i_isNA(v)==FALSE) {
                Zeta_i(v) = 1;
            }
        }
    } else {
        for (int v=0; v<V; ++v) {
            if (Residvec_i_isNA(v)==FALSE) {
                NumericVector qiv(H);
                NumericVector qprime_iv(H);
                NumericVector probs_iv(H);
                bool allneg = TRUE;
                for (int h=0; h<H; ++h) {
                    qiv(h) = log(rho(i,h))-.5*log(ssq(i,h))-.5*pow(Residvec_i(v),2)/ssq(i,h);
                    if (qiv(h)>=0) {
                        allneg = FALSE;
                    }
                }
                if (allneg == TRUE) {
                    qprime_iv = qiv;
                } else {
                    qprime_iv = qiv-min(qiv);
                }
                probs_iv = exp(qprime_iv)/sum(exp(qprime_iv));
                LogicalVector isnan_probs = is_na(probs_iv);
                for (int hh=0; hh<H; ++hh) {
                    if (isnan_probs(hh)==TRUE) {
                        probs_iv(hh)=1;
                    }
                }
                Zeta_i(v) = sample_1int_rcpp(probs_iv);
            }
        }
    }
    return Zeta_i;
}

// [[Rcpp::export]]
bool anyNA_cpp(NumericVector x) {
    int L = x.length();
    LogicalVector x_isNA = is_na(x);
    bool anyNA = FALSE;
    for (int l=0; l<L; ++l) {
        if (x_isNA(l)==TRUE) {
            anyNA = TRUE;
        }
    }
    return anyNA;
}

// [[Rcpp::export]]
bool allZero_cpp(NumericVector x) {
    int L = x.length();
    bool allZero = TRUE;
    for (int l=0; l<L; ++l) {
        if (x(l)!=0) {
            allZero = FALSE;
        }
    }
    return allZero;
}

// [[Rcpp::export]]
NumericVector stl_sort(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

// [[Rcpp::export]]
IntegerVector sample_Zeta_scan_cpp(NumericMatrix Residvec, IntegerVector Scanner_ID, NumericMatrix ssq, NumericMatrix rho, int scan) {
    int V = Residvec.ncol();
    int H = ssq.ncol();
    IntegerVector Zeta_scan(V);
    NumericMatrix Resid_scan = getY_scan_noadj_cpp(Residvec, Scanner_ID, scan);
    int n_subj = Resid_scan.nrow();
    for (int v=0; v<V; ++v) {
        NumericVector Resid_scan_v = Resid_scan(_,v);
        LogicalVector Rsv_isNA = is_na(Resid_scan_v);
        int Rvcount = 0;
        double sumRvSq = 0;
        for (int j=0; j<n_subj; ++j) {
            if (Rsv_isNA(j)==FALSE) {
                Rvcount++;
                sumRvSq = sumRvSq + pow(Resid_scan_v(j),2);
            }
        }
        if (Rvcount>0) {
            if (H==1) {
                Zeta_scan(v) = 1;
            } else {
                NumericVector q_scanv(H);
                NumericVector qprime_scanv(H);
                NumericVector probs_scanv(H);
                for (int h=0; h<H; ++h) {
                    q_scanv(h) = log(rho(scan-1,h))-.5*Rvcount*log(ssq(scan-1,h))-.5*sumRvSq/ssq(scan-1,h);
                }
                NumericVector q_scanv_sort = stl_sort(q_scanv);
                int index = 0;
                double min_q_scanv = q_scanv_sort(index);
                qprime_scanv = q_scanv-min_q_scanv;
                probs_scanv = exp(qprime_scanv)/sum(exp(qprime_scanv));
                bool anyNaN = anyNA_cpp(probs_scanv);
                bool allZero = allZero_cpp(probs_scanv);
                while (anyNaN == TRUE || allZero == TRUE) {
                    index++;
                    min_q_scanv = q_scanv_sort(index);
                    qprime_scanv = q_scanv-min_q_scanv;
                    probs_scanv = exp(qprime_scanv)/sum(exp(qprime_scanv));
                    anyNaN = anyNA_cpp(probs_scanv);
                    allZero = allZero_cpp(probs_scanv);
                }
                Zeta_scan(v) = sample_1int_rcpp(probs_scanv);
            }
        }
    }
    return Zeta_scan;
}

// [[Rcpp::export]]
NumericVector sample_rho_i_cpp(IntegerMatrix Zeta, double beta_rho, int H, int i) {
    int V = Zeta.ncol();
    IntegerVector Zeta_i = Zeta(i,_);
    NumericVector rho_i(H);
    NumericVector diri_params(H);
    IntegerVector N_i(H);
    for (int v=0; v<V; ++v) {
        if (Zeta_i(v)>0) {
            N_i(Zeta_i(v)-1)++;
        }
    }
    for (int h=0; h<H; ++h) {
        diri_params(h) = N_i(h) + beta_rho/H;
    }
    rho_i = rdirichlet_cpp(diri_params);
    return rho_i;
}

// [[Rcpp::export]]
NumericVector sample_rho_scan_cpp(IntegerMatrix Zeta, double beta_rho, int H, int scan) {
    int V = Zeta.ncol();
    IntegerVector Zeta_scan = Zeta(scan-1,_);
    NumericVector rho_scan(H);
    NumericVector diri_params(H);
    IntegerVector N_scan(H);
    for (int v=0; v<V; ++v) {
        if (Zeta_scan(v)>0) {
            N_scan(Zeta_scan(v)-1)++;
        }
    }
    for (int h=0; h<H; ++h) {
        diri_params(h) = N_scan(h) + beta_rho/H;
    }
    rho_scan = rdirichlet_cpp(diri_params);
    return rho_scan;
}

// [[Rcpp::export]]
NumericVector sample_ssq_i_cpp(NumericMatrix Residvec, IntegerMatrix Zeta, double a_s, double b_s, int H, int i) {
    int V = Residvec.ncol();
    NumericVector ssq_i(H);
    NumericVector Resid_i = Residvec(i,_);
    for (int h=0; h<H; ++h) {
        int N_ih = 0;
        double sumsq_R_h = 0;
        for (int v=0; v<V; ++v) {
            if (Zeta(i,v)==(h+1)) {
                sumsq_R_h = sumsq_R_h + pow(Resid_i(v),2);
                N_ih++;
            }
        }
        ssq_i(h) = 1.0/R::rgamma(a_s+.5*N_ih, 1.0/(b_s+.5*sumsq_R_h));
    }
    return ssq_i;
}

// [[Rcpp::export]]
NumericVector sample_ssq_scan_cpp(NumericMatrix Residvec, IntegerVector Scanner_ID, IntegerMatrix Zeta, double a_s, double b_s, int H, int scan) {
    int V = Residvec.ncol();
    NumericVector ssq_scan(H);
    NumericMatrix Resid_scan = getY_scan_noadj_cpp(Residvec, Scanner_ID, scan);
    int n_subj = Resid_scan.nrow();
    for (int h=0; h<H; ++h) {
        int N_scanh = 0;
        double sumsq_R_scanh = 0;
        for (int v=0; v<V; ++v) {
            if (Zeta(scan-1,v)==(h+1)) {
                NumericVector Resid_scan_v = Resid_scan(_,v);
                LogicalVector R_scanv_isNA = is_na(Resid_scan_v);
                for (int j=0; j<n_subj; ++j) {
                    if (R_scanv_isNA(j)==FALSE) {
                        N_scanh++;
                        sumsq_R_scanh = sumsq_R_scanh + pow(Resid_scan_v(j),2);
                    }
                }
            }
        }
        ssq_scan(h) = 1.0/R::rgamma(a_s+.5*N_scanh, 1.0/(b_s+.5*sumsq_R_scanh));
    }
    return ssq_scan;
}

// [[Rcpp::export]]
List initStorageVars_harm_long_cpp(NumericMatrix Yvec, IntegerVector p, NumericMatrix X_multi, IntegerVector Scanner_ID, IntegerVector Subj_ID, int niter, IntegerVector R, int H, int subj_coef_type) {

// subj_coef_type = 1: Lowrank
// subj_coef_type = 2: Equal voxels
// subj_coef_type = 3: Spike and slab
  
// All parameters are lists with niter observations. The iter^th iteration of each has the following sizes:

// population intercept margin:
// gamma_pop ~ List of length D, each element is matrix of size p(d) x R(0)
// w_pop / lambda_pop / alpha_pop ~ Matrix of size D x R(0)
// tau_pop ~ Scalar

// scanner intercept margins:
// gamma_scan ~ List of length nScanners, each element is list of length D, each element is matrix of size p(d) x R(1)
// w_scan / lambda_scan / alpha_scan ~ Matrix of size D x R(1)
// tau_scan ~ Vector of length nScanners
    
// covariate effects margins:
// gamma_cov ~ List of length Q, each element is list of length D, each element is matrix of size p(d) x R(2)
// w_cov / lambda_cov / alpha_cov ~ Matrix of size D x R(2)
// tau_cov ~ Vector of length Q

// subject intercept margins:
// if subj_coef_type = 1:
//  gamma_subj ~ List of length nSubj, each element is list of length D, each element is matrix of size p(d) x R(3)
//  w_subj / lambda_subj / alpha_subj ~ Matrix of size D x R(3)
//  tau_subj ~ Vector of length nSubj
// if subj_coef_type = 2:
//  gamma_subj ~ Vector of length nSubj
//  w_subj / lambda_subj / alpha_subj ~ NULL
//  tau_subj ~ Scalar
// if subj_coef_type = 3:
//  gamma_subj ~ Vector of length nSubj
// w_subj / lambda_subj / alpha_subj ~ NULL
// tau_subj ~ Scalar
// Z_subj ~ Matrix of size nSubj x V
// rho_subj ~ Vector of length nSubj

// residual variance:
// Zeta ~ Matrix of size nobs x V
// rho / ssq ~ Matrix of size nobs x H
    
    NumericVector probs(2);
    probs(0) = 0.5;
    probs(1) = 0.5;
    //int nobs = Yvec.nrow();
    int nScanners = max(Scanner_ID);
    int nSubj = max(Subj_ID);
    int Q = X_multi.ncol();
    int D = p.length();
    int V = Yvec.ncol();
    NumericVector probs_H(H);
    for (int h=0; h<H; ++h) {
      probs_H(h) = 1.0/H;
    }
    if (R.length()!=4) {
      R = rep(R[1],4);
    }
    
    // initialize lists
    List gamma_pop_store(niter+1);
    List w_pop_store(niter+1);
    List lambda_pop_store(niter+1);
    List alpha_pop_store(niter+1);
    List tau_pop_store(niter+1);
    
    List gamma_scan_store(niter+1);
    List w_scan_store(niter+1);
    List lambda_scan_store(niter+1);
    List alpha_scan_store(niter+1);
    List tau_scan_store(niter+1);
    
    List gamma_cov_store(niter+1);
    List w_cov_store(niter+1);
    List lambda_cov_store(niter+1);
    List alpha_cov_store(niter+1);
    List tau_cov_store(niter+1);
    
    List gamma_subj_store(niter+1);
    List w_subj_store(niter+1);
    List lambda_subj_store(niter+1);
    List alpha_subj_store(niter+1);
    List tau_subj_store(niter+1);
    List Z_subj_store(niter+1);
    List rho_subj_store(niter+1);
    
    List Zeta_store(niter+1);
    List rho_store(niter+1);
    List ssq_store(niter+1);
    
    // fill lists
    for (int i=0; i<(niter+1); ++i) {

        // gamma's
        List gamma_pop_i(D);
        for (int d=0; d<D; ++d) {
            NumericMatrix gamma_pop_id(p(d),R(0));
            if (i==0) {
                for (int r=0; r<R(0); ++r) {
                    gamma_pop_id(_,r) = rnorm(p(d),0,1);
                }
            }
            gamma_pop_i[d] = gamma_pop_id;
        }
        gamma_pop_store[i] = gamma_pop_i;
        
        List gamma_scan_i(nScanners);
        for (int s=0; s<nScanners; ++s) {
            List gamma_scan_is(D);
            for (int d=0; d<D; ++d) {
                NumericMatrix gamma_scan_isd(p(d),R(1));
                if (i==0) {
                    for (int r=0; r<R(1); ++r) {
                        gamma_scan_isd(_,r) = rnorm(p(d),0,1);
                    }
                }
                gamma_scan_is[d] = gamma_scan_isd;
            }
            gamma_scan_i[s] = gamma_scan_is;
        }
        gamma_scan_store[i] = gamma_scan_i;
        
        List gamma_cov_i(Q);
        for (int q=0; q<Q; ++q) {
            List gamma_cov_iq(D);
            for (int d=0; d<D; ++d) {
                NumericMatrix gamma_cov_iqd(p(d),R(2));
                if (i==0) {
                    for (int r=0; r<R(2); ++r) {
                        gamma_cov_iqd(_,r) = rnorm(p(d),0,1);
                    }
                }
                gamma_cov_iq[d] = gamma_cov_iqd;
            }
            gamma_cov_i[q] = gamma_cov_iq;
        }
        gamma_cov_store[i] = gamma_cov_i;
        
        if (subj_coef_type == 1) {
          List gamma_subj_i(nSubj);
          for (int u=0; u<nSubj; ++u) {
            List gamma_subj_iu(D);
            for (int d=0; d<D; ++d) {
              NumericMatrix gamma_subj_iud(p(d),R(3));
              if (i==0) {
                for (int r=0; r<R(3); ++r) {
                  gamma_subj_iud(_,r) = rnorm(p(d),0,1);
                }
              }
              gamma_subj_iu[d] = gamma_subj_iud;
            }
            gamma_subj_i[u] = gamma_subj_iu;
          }
          gamma_subj_store[i] = gamma_subj_i;
        } else {
          NumericVector gamma_subj_i(nSubj);
          if (i==0) {
            gamma_subj_i = rnorm(nSubj,0,1);
          }
          gamma_subj_store[i] = gamma_subj_i;
        }

        // w/lambda/alpha's
        NumericMatrix w_pop_i(D,R(0));
        NumericMatrix w_scan_i(D,R(1));
        NumericMatrix w_cov_i(D,R(2));
        NumericMatrix w_subj_i(D,R(3));
        
        NumericMatrix lambda_pop_i(D,R(0));
        NumericMatrix lambda_scan_i(D,R(1));
        NumericMatrix lambda_cov_i(D,R(2));
        NumericMatrix lambda_subj_i(D,R(3));
        
        NumericMatrix alpha_pop_i(D,R(0));
        NumericMatrix alpha_scan_i(D,R(1));
        NumericMatrix alpha_cov_i(D,R(2));
        NumericMatrix alpha_subj_i(D,R(3));
        
        if (i==0) {
            w_pop_i = w_pop_i+1;
            w_scan_i = w_scan_i+1;
            w_cov_i = w_cov_i+1;
            w_subj_i = w_subj_i+1;
            
            lambda_pop_i = lambda_pop_i+1;
            lambda_scan_i = lambda_scan_i+1;
            lambda_cov_i = lambda_cov_i+1;
            lambda_subj_i = lambda_subj_i+1;
            
            alpha_pop_i = alpha_pop_i+1;
            alpha_scan_i = alpha_scan_i+1;
            alpha_cov_i = alpha_cov_i+1;
            alpha_subj_i = alpha_subj_i+1;
        }
        
        w_pop_store[i] = w_pop_i;
        w_scan_store[i] = w_scan_i;
        w_cov_store[i] = w_cov_i;
        w_subj_store[i] = w_subj_i;
        
        lambda_pop_store[i] = lambda_pop_i;
        lambda_scan_store[i] = lambda_scan_i;
        lambda_cov_store[i] = lambda_cov_i;
        lambda_subj_store[i] = lambda_subj_i;
        
        alpha_pop_store[i] = alpha_pop_i;
        alpha_scan_store[i] = alpha_scan_i;
        alpha_cov_store[i] = alpha_cov_i;
        alpha_subj_store[i] = alpha_subj_i;

        // tau's
        double tau_pop_i;
        NumericVector tau_scan_i(nScanners);
        NumericVector tau_cov_i(Q);
        //double tau_subj_i = 0;
        NumericVector tau_subj_i;
        if (subj_coef_type == 1) {
          tau_subj_i = NumericVector(nSubj);
        } else {
          tau_subj_i = NumericVector(1);
        }
        if (i==0) {
            tau_pop_i = 1;
            tau_scan_i = tau_scan_i+1;
            tau_cov_i = tau_cov_i+1;
            tau_subj_i = tau_subj_i+1;
        }
        tau_pop_store[i] = tau_pop_i;
        tau_scan_store[i] = tau_scan_i;
        tau_cov_store[i] = tau_cov_i;
        tau_subj_store[i] = tau_subj_i;
        
        // spike and slab terms (Z, rho)
        if (subj_coef_type == 3) {
          IntegerMatrix Z_subj_i(nSubj,V);
          NumericVector rho_subj_i(nSubj);
          if (i==0) {
            for (int i=0; i<nSubj; ++i) {
              for (int v=0; v<V; ++v) {
                Z_subj_i(i,v) = sample_1int_rcpp(probs)-1;
              }
            }
            rho_subj_i = rho_subj_i + 0.5;
          }
          Z_subj_store[i] = Z_subj_i;
          rho_subj_store[i] = rho_subj_i;
        }

        // Zeta
        IntegerMatrix Zeta_i(nScanners,V);
        if (i==0) {
            for (int n=0; n<nScanners; ++n) {
                for (int v=0; v<V; ++v) {
                    Zeta_i(n,v) = sample_1int_rcpp(probs_H);
                }
            }
        }
        Zeta_store[i] = Zeta_i;

        // rho and ssq
        NumericMatrix rho_i(nScanners,H);
        NumericMatrix ssq_i(nScanners,H);
        if (i==0) {
            for (int n=0; n<nScanners; ++n) {
                rho_i(n,_) = probs_H;
                ssq_i(n,_) = ssq_i(n,_)+1;
            }
        }
        rho_store[i] = rho_i;
        ssq_store[i] = ssq_i;
    }
    
    List all_gammas = List::create(gamma_pop_store, gamma_scan_store, gamma_cov_store,gamma_subj_store);
    List all_ws = List::create(w_pop_store,w_scan_store,w_cov_store,w_subj_store);
    List all_lambdas = List::create(lambda_pop_store,lambda_scan_store,lambda_cov_store,lambda_subj_store);
    List all_alphas = List::create(alpha_pop_store,alpha_scan_store,alpha_cov_store,alpha_subj_store);
    List all_taus = List::create(tau_pop_store,tau_scan_store,tau_cov_store,tau_subj_store);
  
    List allparams = List::create(all_gammas, all_ws, all_lambdas, all_alphas, all_taus,
                                  Zeta_store, rho_store, ssq_store, Z_subj_store, rho_subj_store);
    
    return allparams;
}

// [[Rcpp::export]]
List set_gamma_dr(List gamma, NumericVector gamma_dr, int d, int r) {
    List gamma_new = clone(gamma);
    NumericMatrix gamma_d = gamma_new[d];
    gamma_d(_,r) = gamma_dr;
    gamma_new[d] = gamma_d;
    return gamma_new;
}

// [[Rcpp::export]]
NumericMatrix set_w_dr(NumericMatrix w, double w_dr, int d, int r) {
    NumericMatrix w_new = clone(w);
    w_new(d,r) = w_dr;
    return w_new;
}

// [[Rcpp::export]]
NumericMatrix set_lambda_dr(NumericMatrix lambda, double lambda_dr, int d, int r) {
    NumericMatrix lambda_new = clone(lambda);
    lambda_new(d,r) = lambda_dr;
    return lambda_new;
}

// [[Rcpp::export]]
NumericMatrix set_alpha_dr(NumericMatrix alpha, double alpha_dr, int d, int r) {
    NumericMatrix alpha_new = clone(alpha);
    alpha_new(d,r) = alpha_dr;
    return alpha_new;
}

// [[Rcpp::export]]
IntegerMatrix set_Zeta_i(IntegerMatrix Zeta, IntegerVector Zeta_i, int i) {
    IntegerMatrix Zeta_new = clone(Zeta);
    Zeta_new(i,_) = Zeta_i;
    return Zeta_new;
}

// [[Rcpp::export]]
NumericMatrix set_rho_i(NumericMatrix rho, NumericVector rho_i, int i) {
    NumericMatrix rho_new = clone(rho);
    rho_new(i,_) = rho_i;
    return rho_new;
}

// [[Rcpp::export]]
NumericMatrix set_ssq_i(NumericMatrix ssq, NumericVector ssq_i, int i) {
    NumericMatrix ssq_new = clone(ssq);
    ssq_new(i,_) = ssq_i;
    return ssq_new;
}

// [[Rcpp::export]]
List set_gamma_q(List gamma, List gamma_q, int q) {
    List gamma_new = clone(gamma);
    gamma_new[q] = gamma_q;
    return gamma_new;
}

// [[Rcpp::export]]
List getGammaUniv_mcmc_cpp(List gamma_mcmc) {
  int niter = gamma_mcmc.length();
  List Gamma_mcmc(niter);
  for (int iter=0; iter<niter; ++iter) {
    List gamma_iter = gamma_mcmc[iter];
    Gamma_mcmc[iter] = getGamma_cpp(gamma_iter);
  }
  return Gamma_mcmc;
}

// [[Rcpp::export]]
List getGammaMulti_mcmc_cpp(List gamma_mcmc) {
  int niter = gamma_mcmc.length();
  List Gamma_mcmc(niter);
  for (int iter=0; iter<niter; ++iter) {
    List gamma_iter = gamma_mcmc[iter];
    int Q = gamma_iter.length();
    List Gamma_iter(Q);
    for (int q=0; q<Q; ++q) {
      List gamma_q_iter = gamma_iter[q];
      Gamma_iter[q] = getGamma_cpp(gamma_q_iter);
    }
    Gamma_mcmc[iter] = Gamma_iter;
  }
  return Gamma_mcmc;
}

// [[Rcpp::export]]
double getDev_harm_cpp(NumericMatrix Yvec, NumericMatrix X_cov, IntegerVector Scanner_ID, IntegerVector Subj_ID, NumericVector Gamma_pop, List Gamma_cov, List Gamma_scan, List Gamma_subj, NumericMatrix SigSq) {
    int N = Yvec.nrow();
    int V = Yvec.ncol();
    NumericMatrix Resid = getResid_harm_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_pop, Gamma_cov, Gamma_scan, Gamma_subj);
    double L = 0;
    for (int n=0; n<N; ++n) {
        LogicalVector isNA_vec1 = is_na(Resid(n,_));
        LogicalVector isNA_vec2 = is_na(SigSq(Scanner_ID(n)-1,_));
        for (int v=0; v<V; ++v) {
            if (isNA_vec1(v)==FALSE && isNA_vec2(v)==FALSE) {
                L = L - .5*(log(SigSq(Scanner_ID(n)-1,v))+pow(Resid(n,v),2)/SigSq(Scanner_ID(n)-1,v));
            }
        }
    }
    double D = -2*L;
    return D;
}

// [[Rcpp::export]]
NumericVector getDev_allmcmc_harm_cpp(NumericMatrix Yvec, NumericMatrix X_cov, IntegerVector Scanner_ID, IntegerVector Subj_ID, List Gamma_pop_mcmc, List Gamma_cov_mcmc, List Gamma_scan_mcmc, List Gamma_subj_mcmc, List SigSq_mcmc, int prog_count=10) {
    int niter = Gamma_pop_mcmc.length();
    NumericVector D_mcmc(niter);
    for (int iter=0; iter<niter; ++iter) {
        if ((iter % prog_count) == 0) {
            std::cout << "Iteration: " << iter << "\n";
        }
        D_mcmc(iter) = getDev_harm_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_pop_mcmc[iter], Gamma_cov_mcmc[iter], Gamma_scan_mcmc[iter], Gamma_subj_mcmc[iter], SigSq_mcmc[iter]);
    }
    return D_mcmc;
}

// [[Rcpp::export]]
void change_gamma_list(List gamma, NumericVector gamma_dr, int d, int r) {
  NumericMatrix gamma_d = gamma(d);
  gamma_d(_,r) = gamma_dr;
  gamma(d) = gamma_d;
}

// [[Rcpp::export]]
void change_gamma_multi_list(List gamma, NumericVector gamma_kdr, int d, int r, int k) {
  List gamma_k = gamma(k);
  NumericMatrix gamma_kd = gamma_k(d);
  gamma_kd(_,r) = gamma_kdr;
  gamma_k(d) = gamma_kd;
  gamma(k) = gamma_k;
}

//FUNCTIONS FOR SPIKE-AND-SLAB COEFs

// [[Rcpp::export]]
NumericVector getCoef_j_spsl_cpp(NumericVector beta, IntegerMatrix Z, int j) {
  int V = Z.ncol();
  NumericVector Coef_j(V);
  Coef_j[Z(j,_)==1] = beta(j);
  return Coef_j;
}

// [[Rcpp::export]]
IntegerVector sample_Zj_cpp(NumericMatrix Ytilde, IntegerVector Scanner_ID, NumericVector beta, NumericMatrix SigSq, NumericVector rho, LogicalVector na_vox, int j) {
  double beta_j = beta(j);
  double rho_j = rho(j);
  int V = Ytilde.ncol();
  int nobs = Ytilde.nrow();
  double log_prob0_un = log(1-rho_j);
  double log_prob1_un = log(rho_j);
  IntegerVector Zj(V);
  for (int v=0; v<V; ++v) {
    double log_prob0_un_v = log_prob0_un;
    double log_prob1_un_v = log_prob1_un;
    if (na_vox(v) == TRUE) {
      Zj(v) = -1;
    } else {
      for (int i=0; i<nobs; ++i) {
        log_prob1_un_v = log_prob1_un_v - .5*pow(beta_j,2)/SigSq(Scanner_ID(i)-1,v) + Ytilde(i,v)*beta_j/SigSq(Scanner_ID(i)-1,v);
      }
      double log_prob0_c_v = 0;
      double log_prob1_c_v = 0;
      if (log_prob1_un_v >= log_prob0_un_v) {
        log_prob0_c_v = log_prob0_un_v - log_prob1_un_v;
      } else {
        log_prob1_c_v = log_prob1_un_v - log_prob0_un_v;
      }
      double prob1_n_v = exp(log_prob1_c_v)/(exp(log_prob1_c_v)+exp(log_prob0_c_v));
      double prob0_n_v = 1-prob1_n_v;
      NumericVector probs_v(2);
      probs_v(0) = prob0_n_v;
      probs_v(1) = prob1_n_v;
      Zj(v) = sample_1int_rcpp(probs_v)-1;
    }
  }
  return Zj;
}

// [[Rcpp::export]]
double sample_betaj_cpp(NumericMatrix Ytilde, IntegerVector Scanner_ID, IntegerMatrix Z, NumericMatrix SigSq, double tausq, int j) {
  // need to include X in both numerator and denominator if sampling random slope
  int nobs = Ytilde.nrow();
  int V = Ytilde.ncol();
  IntegerVector Zj = Z(j,_);
  double numer;
  double denom = 1/tausq;
  for (int i=0; i<nobs; ++i) {
    LogicalVector isna_Y = is_na(Ytilde(i,_));
    for (int v=0; v<V; ++v) {
      if (isna_Y(v)==FALSE && Zj(v)==1) {
        numer = numer + Ytilde(i,v)/SigSq(Scanner_ID(i)-1,v);
        denom = denom + 1/SigSq(Scanner_ID(i)-1,v);
      }
    }
  }
  double mean = numer/denom;
  double var = 1/denom;
  double betaj = R::rnorm(mean, pow(var,0.5));
  return betaj;
}

// [[Rcpp::export]]
double sample_tau_spsl_cpp(NumericVector gamma, double a_tau, double b_tau) {
  double tau;
  int nSubj = gamma.length();
  double aprime = a_tau + 0.5*nSubj;
  double sum_gamma_sq = 0;
  for (int u=0; u<nSubj; ++u) {
    sum_gamma_sq = sum_gamma_sq + pow(gamma(u),2);
  }
  double bprime = b_tau + .5*sum_gamma_sq;
  tau = 1.0/R::rgamma(aprime, 1.0/bprime);
  return tau;
}

// [[Rcpp::export]]
double sample_rhoj_cpp(IntegerMatrix Z, double a_rho, double b_rho, int j) {
  IntegerVector Zj = Z(j,_);
  double aprime_rho = a_rho + sum(Zj==1);
  double bprime_rho = b_rho + sum(Zj==0);
  double rhoj = R::rbeta(aprime_rho, bprime_rho);
  return rhoj;
}

// [[Rcpp::export]]
double getLowRankCoefAtVox_cpp(List gamma, IntegerVector vox) {
  int D = vox.length();
  NumericMatrix gamma0 = gamma(0);
  int R = gamma0.ncol();
  double lowrank_sum = 0;
  for (int r=0; r<R; ++r) {
    double rank_r_prod = 1;
    for (int d=0; d<D; ++d) {
      NumericMatrix gamma_d = gamma(d);
      rank_r_prod = rank_r_prod * gamma_d(vox(d)-1,r);
    }
    lowrank_sum = lowrank_sum + rank_r_prod;
  }
  return lowrank_sum;
}

// [[Rcpp::export]]
NumericVector getSigSq_harm_atVox_cpp(IntegerVector Scanner_ID, NumericMatrix ssq, IntegerMatrix Zeta, int v) {
  int nScanners = max(Scanner_ID);
  NumericVector SigSq_vox(nScanners);
  for (int i=0; i<nScanners; ++i) {
    SigSq_vox(i) = ssq(i, Zeta(i,v)-1);
  }
  return SigSq_vox;
}

// MISC / OBSOLETE FUNCTIONS

// List runMCMC_harm_long_cpp(NumericMatrix Yvec, IntegerVector p, NumericMatrix X_cov, IntegerVector Scanner_ID, IntegerVector Subj_ID, int niter, IntegerVector R, int H, int subj_coef_type,
//                            double a_lam, double b_lam, double a_tau, double b_tau, double a_alpha, double b_alpha, double sigma_log_alpha, double beta_rho, double a_s, double b_s,
//                            int prog_count, bool show_all_steps, bool null_pop, bool null_scan, bool null_cov, bool null_subj, List sliceInds_alldims, LogicalVector na_vox) {
//   // set hyperparameters
//   // double a_lam = 1;
//   // double b_lam = 1;
//   // double a_tau = 1;
//   // double b_tau = 1;
//   // double a_alpha = 1;
//   // double b_alpha = 1;
//   // double sigma_log_alpha = 0.01;
//   // double beta_rho = 1;
//   // double a_s = 1;
//   // double b_s = 1;
// 
//   // establish rank if equal across coefficients
//   IntegerVector Rtemp(4);
//   if (R.length()!=4) {
//     for (int ii=0; ii<4; ++ii) {
//       Rtemp(ii) = R(0);
//     }
//     R = clone(Rtemp);
//   }
// 
//   // establish length terms
//   int D = p.length();
//   int Q = X_cov.ncol();
//   int nScanners = max(Scanner_ID);
//   int nSubj = max(Subj_ID);
//   int nobs = Yvec.nrow();
//   int V = Yvec.ncol();
//   NumericVector X_intercept(nobs);
//   X_intercept = X_intercept + 1;
// 
//   // initialize model parameters
//   List allparams_store = initStorageVars_harm_long_cpp(Yvec, p, X_cov, Scanner_ID, Subj_ID, niter, R, H, subj_coef_type);
// 
//   List gammas_store = allparams_store(0);
//   List gamma_pop_store = gammas_store(0);
//   List gamma_scan_store = gammas_store(1);
//   List gamma_cov_store = gammas_store(2);
//   List gamma_subj_store = gammas_store(3);
// 
//   List ws_store = allparams_store(1);
//   List w_pop_store = ws_store(0);
//   List w_scan_store = ws_store(1);
//   List w_cov_store = ws_store(2);
//   List w_subj_store = ws_store(3);
// 
//   List lambdas_store = allparams_store(2);
//   List lambda_pop_store = lambdas_store(0);
//   List lambda_scan_store = lambdas_store(1);
//   List lambda_cov_store = lambdas_store(2);
//   List lambda_subj_store = lambdas_store(3);
// 
//   List alphas_store = allparams_store(3);
//   List alpha_pop_store = alphas_store(0);
//   List alpha_scan_store = alphas_store(1);
//   List alpha_cov_store = alphas_store(2);
//   List alpha_subj_store = alphas_store(3);
// 
//   List taus_store = allparams_store(4);
//   List tau_pop_store = taus_store(0);
//   List tau_scan_store = taus_store(1);
//   List tau_cov_store = taus_store(2);
//   List tau_subj_store = taus_store(3);
// 
//   List Zeta_store = allparams_store(5);
//   List rho_store = allparams_store(6);
//   List ssq_store = allparams_store(7);
// 
//   // Get initial coefficients
//   NumericVector zeros_allvox(V);
//   NumericVector Gamma_pop_iter(V);
//   if (null_pop == FALSE) {
//     List gamma_pop_init = gamma_pop_store(0);
//     Gamma_pop_iter = getGamma_cpp(gamma_pop_init);
//   }
// 
//   List Gamma_scan_iter(nScanners);
//   for (int ss=0; ss<nScanners; ++ss) {
//     Gamma_scan_iter(ss) = clone(zeros_allvox);
//   }
//   if (null_scan == FALSE) {
//     List gamma_scan_init = gamma_scan_store(0);
//     Gamma_scan_iter = getGamma_multi_cpp(gamma_scan_init);
//   }
// 
//   if (Q==0) {
//     Q = 1;
//     null_cov = TRUE;
//   }
//   List Gamma_cov_iter(Q);
//   for (int qq=0; qq<Q; ++qq) {
//     Gamma_cov_iter(qq) = clone(zeros_allvox);
//   }
//   if (null_cov == FALSE) {
//     List gamma_cov_init = gamma_cov_store(0);
//     Gamma_cov_iter = getGamma_multi_cpp(gamma_cov_init);
//   }
// 
//   List Gamma_subj_iter(nSubj);
//   for (int uu=0; uu<nSubj; ++uu) {
//     Gamma_subj_iter(uu) = clone(zeros_allvox);
//   }
//   if (null_subj == FALSE) {
//     if (subj_coef_type == 1) {
//       List gamma_subj_init = gamma_subj_store(0);
//       Gamma_subj_iter = getGamma_multi_cpp(gamma_subj_init);
//     } else {
//       NumericVector gamma_subj_init = gamma_subj_store(0);
//       for (int uu=0; uu<nSubj; ++uu) {
//         NumericVector Gamma_subj_iter_uu(V);
//         Gamma_subj_iter_uu = Gamma_subj_iter_uu + gamma_subj_init(uu);
//         Gamma_subj_iter(uu) = Gamma_subj_iter_uu;
//       }
//     }
//   }
// 
//   // loop over MCMC iterations
//   for (int iter=0; iter<niter; ++iter) {
//     // print current iteration
//     if (((iter+1) % prog_count)==0) {
//       printf("Iteration: %d\n",iter+1);
//     }
// 
//     // initialize current iteration for model parameters
//     List gamma_pop_iter = clone(gamma_pop_store(iter));
//     NumericMatrix w_pop_iter = clone(w_pop_store(iter));
//     NumericMatrix lambda_pop_iter = clone(lambda_pop_store(iter));
//     NumericMatrix alpha_pop_iter = clone(alpha_pop_store(iter));
//     double tau_pop_iter = clone(tau_pop_store(iter));
// 
//     List gamma_scan_iter = clone(gamma_scan_store(iter));
//     NumericMatrix w_scan_iter = clone(w_scan_store(iter));
//     NumericMatrix lambda_scan_iter = clone(lambda_scan_store(iter));
//     NumericMatrix alpha_scan_iter = clone(alpha_scan_store(iter));
//     NumericVector tau_scan_iter = clone(tau_scan_store(iter));
// 
//     List gamma_cov_iter = clone(gamma_cov_store(iter));
//     NumericMatrix w_cov_iter = clone(w_cov_store(iter));
//     NumericMatrix lambda_cov_iter = clone(lambda_cov_store(iter));
//     NumericMatrix alpha_cov_iter = clone(alpha_cov_store(iter));
//     NumericVector tau_cov_iter = clone(tau_cov_store(iter));
// 
//     // if (subj_int_lowrank==TRUE) {
//     //   List gamma_subj_iter = gamma_subj_store(iter);
//     //   NumericVector tau_subj_iter = tau_subj_store(iter);
//     // } else {
//     //   NumericVector gamma_subj_iter = gamma_subj_store(iter);
//     //   double tau_subj_iter = tau_subj_store(iter);
//     // }
//     NumericMatrix w_subj_iter = clone(w_subj_store(iter));
//     NumericMatrix lambda_subj_iter = clone(lambda_subj_store(iter));
//     NumericMatrix alpha_subj_iter = clone(alpha_subj_store(iter));
// 
//     IntegerMatrix Zeta_iter = clone(Zeta_store(iter));
//     NumericMatrix rho_iter = clone(rho_store(iter));
//     NumericMatrix ssq_iter = clone(ssq_store(iter));
//     NumericMatrix SigSq_iter = getSigSq_cpp(ssq_iter, Zeta_iter);
// 
// 
//     // sample population intercept term
//     if (null_pop == FALSE) {
//       if (show_all_steps==TRUE) {
//         printf("%d: Population intercept\n",iter+1);
//       }
//       NumericMatrix Y_intercept = getY_intercept_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_cov_iter, Gamma_scan_iter, Gamma_subj_iter);
//       for (int r=0; r<R(0); ++r) {
//         NumericMatrix Yrvec_pop = getYr_cpp(Y_intercept, X_intercept, gamma_pop_iter, r);
//         for (int d=0; d<D; ++d) {
//           NumericVector gamma_pop_iter_dr = sample_gamma_dr_harm_withInds_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d, r, sliceInds_alldims(d));
//           //NumericVector gamma_pop_iter_dr = sample_gamma_dr_harm_cpp(Yrvec_pop, p, Scanner_ID, X_intercept, gamma_pop_iter, tau_pop_iter, w_pop_iter, alpha_pop_iter, SigSq_iter, d, r);
//           change_gamma_list(gamma_pop_iter,gamma_pop_iter_dr,d,r);
//           w_pop_iter(d,r) = sample_w_dr_cpp(gamma_pop_iter, tau_pop_iter, alpha_pop_iter, lambda_pop_iter, p, d, r);
//           lambda_pop_iter(d,r) = sample_lambda_dr_cpp(w_pop_iter, p, a_lam, b_lam, d, r);
//           alpha_pop_iter(d,r) = sample_alpha_dr_cpp(alpha_pop_iter, gamma_pop_iter, tau_pop_iter, w_pop_iter, p, a_alpha, b_alpha, sigma_log_alpha, d, r);
//         }
//       }
//       tau_pop_iter = sample_tau_cpp(gamma_pop_iter, w_pop_iter, alpha_pop_iter, a_tau, b_tau);
//       Gamma_pop_iter = getGamma_cpp(gamma_pop_iter);
//     }
// 
//     //sample scanner intercept terms
//     if (null_scan == FALSE) {
//       if (show_all_steps==TRUE) {
//         printf("%d: Scanner intercept\n",iter+1);
//       }
//       for (int s=0; s<nScanners; ++s) {
//         NumericMatrix Y_scan_s = getY_scan_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_pop_iter, Gamma_cov_iter, Gamma_subj_iter, s+1);
//         int Ns = Y_scan_s.nrow();
//         NumericVector X_s(Ns);
//         X_s = X_s + 1.0;
//         IntegerVector Scanner_s(Ns);
//         Scanner_s = Scanner_s + s + 1;
//         for (int r=0; r<R(1); ++r) {
//           NumericMatrix Yrsvec = getYr_cpp(Y_scan_s, X_s, gamma_scan_iter(s), r);
//           for (int d=0; d<D; ++d) {
//             NumericVector gamma_scan_iter_drs = sample_gamma_dr_harm_withInds_cpp(Yrsvec, p, Scanner_s, X_s, gamma_scan_iter(s), tau_scan_iter(s), w_scan_iter, alpha_scan_iter, SigSq_iter, d, r,sliceInds_alldims(d));
//             //NumericVector gamma_scan_iter_drs = sample_gamma_dr_harm_cpp(Yrsvec, p, Scanner_s, X_s, gamma_scan_iter(s), tau_scan_iter(s), w_scan_iter, alpha_scan_iter, SigSq_iter, d, r);
//             change_gamma_multi_list(gamma_scan_iter, gamma_scan_iter_drs, d, r, s);
//           }
//         }
//       }
//       for (int r=0; r<R(1); ++r) {
//         for (int d=0; d<D; ++d) {
//           w_scan_iter(d,r) = sample_w_dr_multi_cpp(gamma_scan_iter, tau_scan_iter, alpha_scan_iter, lambda_scan_iter, p, d, r);
//           lambda_scan_iter(d,r) = sample_lambda_dr_cpp(w_scan_iter, p, a_lam, b_lam, d, r);
//           alpha_scan_iter(d,r) = sample_alpha_dr_multi_cpp(alpha_scan_iter, gamma_scan_iter, tau_scan_iter, w_scan_iter, p, a_alpha,b_alpha, sigma_log_alpha, d, r);
//         }
//       }
//       tau_scan_iter = sample_tau_multi_cpp(gamma_scan_iter, w_scan_iter, alpha_scan_iter, a_tau, b_tau);
//       Gamma_scan_iter = getGamma_multi_cpp(gamma_scan_iter);
//     }
// 
//     // sample covariate effects terms
//     if (null_cov == FALSE) {
//       if (show_all_steps==TRUE) {
//         printf("%d: Covariate effects\n",iter+1);
//       }
//       NumericMatrix Y_cov = getY_cov_long_cpp(Yvec, Scanner_ID, Subj_ID, Gamma_pop_iter, Gamma_scan_iter, Gamma_subj_iter);
//       for (int q=0; q<Q; ++q) {
//         if (q > 0) {
//           Gamma_cov_iter(q-1) = getGamma_cpp(gamma_cov_iter(q-1));
//         }
//         NumericMatrix Y_cov_q = getYq_multi_cpp(Y_cov, X_cov, Gamma_cov_iter, q);
//         for (int r=0; r<R(2); ++r) {
//           NumericMatrix Yrqvec = getYr_cpp(Y_cov_q, X_cov(_,q), gamma_cov_iter(q), r);
//           for (int d=0; d<D; ++d) {
//             NumericVector gamma_cov_iter_drq = sample_gamma_dr_harm_withInds_cpp(Yrqvec, p, Scanner_ID, X_cov(_,q), gamma_cov_iter(q), tau_cov_iter(q), w_cov_iter, alpha_cov_iter, SigSq_iter, d, r,sliceInds_alldims(d));
//             //NumericVector gamma_cov_iter_drq = sample_gamma_dr_harm_cpp(Yrqvec, p, Scanner_ID, X_cov(_,q), gamma_cov_iter(q), tau_cov_iter(q), w_cov_iter, alpha_cov_iter, SigSq_iter, d, r);
//             change_gamma_multi_list(gamma_cov_iter, gamma_cov_iter_drq, d, r, q);
//           }
//         }
//       }
//       for (int r=0; r<R(2); ++r) {
//         for (int d=0; d<D; ++d) {
//           w_cov_iter(d,r) = sample_w_dr_multi_cpp(gamma_cov_iter, tau_cov_iter, alpha_cov_iter, lambda_cov_iter, p, d, r);
//           lambda_cov_iter(d,r) = sample_lambda_dr_cpp(w_cov_iter, p, a_lam, b_lam, d, r);
//           alpha_cov_iter(d,r) = sample_alpha_dr_multi_cpp(alpha_cov_iter, gamma_cov_iter, tau_cov_iter, w_cov_iter, p, a_alpha, b_alpha, sigma_log_alpha, d, r);
//         }
//       }
//       tau_cov_iter = sample_tau_multi_cpp(gamma_cov_iter, w_cov_iter, alpha_cov_iter, a_tau, b_tau);
//       Gamma_cov_iter = getGamma_multi_cpp(gamma_cov_iter);
//     }
// 
//     // sample subject intercept terms
//     if (null_subj == FALSE) {
//       if (show_all_steps==TRUE) {
//         printf("%d: Subject intercept\n",iter+1);
//       }
//       if (subj_coef_type == 1) {
//         List gamma_subj_iter = clone(gamma_subj_store(iter));
//         NumericVector tau_subj_iter = clone(tau_subj_store(iter));
//         for (int u=0; u<nSubj; ++u) {
//           NumericMatrix Y_subj_u = getY_subj_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_pop_iter, Gamma_cov_iter, Gamma_scan_iter, u+1);
//           int Nu = Y_subj_u.nrow();
//           NumericVector X_u(Nu);
//           X_u = X_u + 1.0;
//           LogicalVector Subj_u = (Subj_ID==(u+1));
//           IntegerVector Scanner_u = Scanner_ID[Subj_u];
//           for (int r=0; r<R(3); ++r) {
//             NumericMatrix Yruvec = getYr_cpp(Y_subj_u, X_u, gamma_subj_iter(u), r);
//             for (int d=0; d<D; ++d) {
//               NumericVector gamma_subj_iter_dru = sample_gamma_dr_harm_withInds_cpp(Yruvec, p, Scanner_u, X_u, gamma_subj_iter(u), tau_subj_iter(u), w_subj_iter, alpha_subj_iter, SigSq_iter, d, r,sliceInds_alldims(d));
//               //NumericVector gamma_subj_iter_dru = sample_gamma_dr_harm_cpp(Yruvec, p, Scanner_ID, X_u, gamma_subj_iter(u), tau_subj_iter(u), w_subj_iter, alpha_subj_iter, SigSq_iter, d, r);
//               change_gamma_multi_list(gamma_subj_iter, gamma_subj_iter_dru, d, r, u);
//             }
//           }
//         }
//         for (int r=0; r<R(3); ++r) {
//           for (int d=0; d<D; ++d) {
//             w_subj_iter(d,r) = sample_w_dr_multi_cpp(gamma_subj_iter, tau_subj_iter, alpha_subj_iter, lambda_subj_iter, p, d, r);
//             lambda_subj_iter(d,r) = sample_lambda_dr_cpp(w_subj_iter, p, a_lam, b_lam, d, r);
//             alpha_subj_iter(d,r) = sample_alpha_dr_multi_cpp(alpha_subj_iter, gamma_subj_iter, tau_subj_iter, w_subj_iter, p, a_alpha, b_alpha, sigma_log_alpha, d, r);
//           }
//         }
//         tau_subj_iter = sample_tau_multi_cpp(gamma_subj_iter, w_subj_iter, alpha_subj_iter, a_tau, b_tau);
//         Gamma_subj_iter = getGamma_multi_cpp(gamma_subj_iter);
//         //save storage parameters
//         gamma_subj_store(iter+1) = gamma_subj_iter;
//         tau_subj_store(iter+1) = tau_subj_iter;
//       } else {
//         NumericVector gamma_subj_iter = clone(gamma_subj_store(iter));
//         double tau_subj_iter = clone(tau_subj_store(iter));
//         for (int u=0; u<nSubj; ++u) {
//           NumericMatrix Y_subj_u = getY_subj_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_pop_iter, Gamma_cov_iter, Gamma_scan_iter, u+1);
//           LogicalVector Subj_u = (Subj_ID==(u+1));
//           IntegerVector Scanner_u = Scanner_ID[Subj_u];
//           gamma_subj_iter(u) = sample_gamma_u_eqvox_harm_cpp(Y_subj_u,Scanner_u,tau_subj_iter,SigSq_iter);
//           NumericVector Gamma_subj_iter_u(V);
//           Gamma_subj_iter_u = Gamma_subj_iter_u + gamma_subj_iter(u);
//           Gamma_subj_iter(u) = Gamma_subj_iter_u;
//         }
//         tau_subj_iter = sample_tau_eqvox_cpp(gamma_subj_iter, a_tau, b_tau);
//         // save storage parameters
//         gamma_subj_store(iter+1) = gamma_subj_iter;
//         tau_subj_store(iter+1) = tau_subj_iter;
//       }
//     }
// 
//     // sample residual terms
//     if (show_all_steps==TRUE) {
//       printf("%d: Residual noise\n",iter+1);
//     }
//     NumericMatrix Residvec_iter = getResid_harm_long_cpp(Yvec, X_cov, Scanner_ID, Subj_ID, Gamma_pop_iter, Gamma_cov_iter, Gamma_scan_iter, Gamma_subj_iter);
//     for (int scan = 0; scan<nScanners; ++scan) {
//       Zeta_iter(scan,_) = sample_Zeta_scan_cpp(Residvec_iter, Scanner_ID, ssq_iter, rho_iter, scan+1);
//       rho_iter(scan,_) = sample_rho_scan_cpp(Zeta_iter, beta_rho, H, scan+1);
//       ssq_iter(scan,_) = sample_ssq_scan_cpp(Residvec_iter, Scanner_ID, Zeta_iter, a_s, b_s, H, scan+1);
//     }
// 
//     // store sampled parameters
//     gamma_pop_store(iter+1) = gamma_pop_iter;
//     w_pop_store(iter+1) = w_pop_iter;
//     lambda_pop_store(iter+1) = lambda_pop_iter;
//     alpha_pop_store(iter+1) = alpha_pop_iter;
//     tau_pop_store(iter+1) = tau_pop_iter;
// 
//     gamma_scan_store(iter+1) = gamma_scan_iter;
//     w_scan_store(iter+1) = w_scan_iter;
//     lambda_scan_store(iter+1) = lambda_scan_iter;
//     alpha_scan_store(iter+1) = alpha_scan_iter;
//     tau_scan_store(iter+1) = tau_scan_iter;
// 
//     gamma_cov_store(iter+1) = gamma_cov_iter;
//     w_cov_store(iter+1) = w_cov_iter;
//     lambda_cov_store(iter+1) = lambda_cov_iter;
//     alpha_cov_store(iter+1) = alpha_cov_iter;
//     tau_cov_store(iter+1) = tau_cov_iter;
// 
//     //gamma_subj_store(iter+1) = gamma_subj_iter;
//     w_subj_store(iter+1) = w_subj_iter;
//     lambda_subj_store(iter+1) = lambda_subj_iter;
//     alpha_subj_store(iter+1) = alpha_subj_iter;
//     //tau_subj_store(iter+1) = tau_subj_iter;
// 
//     Zeta_store(iter+1) = Zeta_iter;
//     rho_store(iter+1) = rho_iter;
//     ssq_store(iter+1) = ssq_iter;
// 
//   }
// 
//   // save all MCMC samples
//   List all_gammas_samp = List::create(gamma_pop_store, gamma_scan_store, gamma_cov_store,gamma_subj_store);
//   List all_ws_samp = List::create(w_pop_store,w_scan_store,w_cov_store,w_subj_store);
//   List all_lambdas_samp = List::create(lambda_pop_store,lambda_scan_store,lambda_cov_store,lambda_subj_store);
//   List all_alphas_samp = List::create(alpha_pop_store,alpha_scan_store,alpha_cov_store,alpha_subj_store);
//   List all_taus_samp = List::create(tau_pop_store,tau_scan_store,tau_cov_store,tau_subj_store);
// 
//   List allparams_samp = List::create(all_gammas_samp, all_ws_samp, all_lambdas_samp, all_alphas_samp, all_taus_samp,
//                                         Zeta_store, rho_store, ssq_store);
//   return allparams_samp;
// }
