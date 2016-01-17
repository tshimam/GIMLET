#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h> 

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double ldcorCpp(const arma::vec dx, const arma::vec dy, const arma::vec w){
  int n = w.size();
  double f = 0;
  double s1 = 0;
  double s2 = 0;
  double s3 = 0;
  double s2a = 0;
  double s2b = 0;
  double s1x = 0;
  double s1y = 0;
  double s2x = 0;
  double s2y = 0;
  double s3x = 0;
  double s3y = 0;
  arma::vec edx(n);
  arma::vec edy(n);
  edx.zeros();
  edy.zeros();
  int k = 0;
  for (int i = 0; i < (n - 1); i++) {
    for (int j = i + 1; j < n; j++) {
      f = w(i) * w(j);
      s1 = s1 + dx(k) * dy(k) * f;
      s1x = s1x + dx(k) * dx(k) * f;
      s1y = s1y + dy(k) * dy(k) * f;
      edx(i) = edx(i) + dx(k) * w(j);
      edx(j) = edx(j) + dx(k) * w(i);
      edy(i) = edy(i) + dy(k) * w(j);
      edy(j) = edy(j) + dy(k) * w(i);
      k++;
    }
  }
  for (int i = 0; i < n; i++) {
    s3 = s3 + edx(i) * edy(i) * w(i);
    s2a = s2a + edx(i) * w(i);
    s2b = s2b + edy(i) * w(i);
    s3x = s3x + edx(i) * edx(i) * w(i);
    s3y = s3y + edy(i) * edy(i) * w(i);
  }
  s1 = 2 * s1;
  s1x = 2 * s1x;
  s1y = 2 * s1y;
  s2 = s2a * s2b;
  s2x = s2a * s2a;
  s2y = s2b * s2b;
  double r = sqrt((s1 + s2 - 2*s3) / sqrt((s1x + s2x - 2 * s3x) * (s1y + s2y - 2 * s3y)));
  return r;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec distCpp(const arma::mat X, const int p){
  int n = X.n_rows;
  int m = n * (n-1) / 2;
  arma::vec d(m);
  d.zeros();
  int k = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if(i < j) {
        d(k) = norm(X.row(i) - X.row(j),p);
        k++;
      }
    }
  }
  return d;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec gaukerCpp(const arma::mat X, const arma::rowvec mu, double lambda){
  int n = X.n_rows;
  double tker;
  arma::vec ker(n);
  for(int i = 0; i < n; i++){
    tker = norm(X.row(i) - mu, 1);
    ker(i) = exp(-(2.5 * tker * tker) / (2.0 * lambda));
  }
  return ker;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec neikerCpp(const arma::mat X, const arma::rowvec mu, int r){
  int n = X.n_rows;
  arma::vec tker(n);
  arma::vec ker(n);
  arma::uvec id(n);
  for(int i = 0; i < n; i++){
    tker(i) = norm(X.row(i) - mu, 1);
    ker(i) = 0;
  }
  id = sort_index(tker);
  for(int i = 0; i < r; i++){
    ker(id(i)) = 1;
  }
  return ker;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List gimletCpp(const arma::mat X, const arma::mat Y, const arma::mat Z,
		     const int r, const double delta, const int bootn, const int type){
  int n = X.n_rows;
  int ln = round(n*delta);
  int p = Y.n_cols;
  int q = Z.n_cols;
  double Zmin;
  double Zmax;
  double lambda;
  double S;
  arma::vec normz(n);
  arma::vec dx = distCpp(X,1);
  arma::vec dy;
  arma::vec tmpC(2*r);
  arma::vec dd(n);
  arma::vec ww(n);
  arma::vec gp(p);
  arma::uvec id(n);
  arma::mat nZ(n,q);
  arma::mat M(2*r,q);
  arma::mat W(n,2*r);
  arma::mat C(2*r,p);
  arma::mat ep(2*r,p);
  arma::mat nullW(n,2*r);
  arma::mat nullS(p,bootn,fill::zeros);
  arma::cube nullC(2*r,p,bootn,fill::zeros);
  arma::vec nullCC(2*r);
  for(int i = 0; i < q; i++){
    Zmin = min(Z.col(i));
    Zmax = max(Z.col(i));
    nZ.col(i) = (Z.col(i) - Zmin) / (Zmax - Zmin);
  }
  for(int i = 0; i < n; i++){
    normz(i) = norm(nZ.row(i),1);
  }
  id = sort_index(normz);
  int count = 0;
  for(int i = 0; i < r; i++){
    M.row(count) = nZ.row(id(i));
    count++;
  }
  for(int i = n - r; i < n; i++){
    M.row(count) = nZ.row(id(i));
    count++;
  }
  for(int i = 0; i < 2*r; i++){
    for(int j = 0; j < n; j++){
      dd(j) = norm(nZ.row(j) - M.row(i), 1);
    }
    dd = sort(dd);
    lambda = dd(ln - 1);
    if(lambda == 0){
      lambda = 1e-3;
    }
    if(type == 1){
      ww = gaukerCpp(nZ, M.row(i), lambda);
    } else {
      ww = neikerCpp(nZ, M.row(i), ln);
    }
    W.col(i) = ww / sum(ww);
  }
  ep.zeros();
  gp.zeros();
  for(int i = 0; i < p; i++){
    dy = distCpp(Y.col(i),1);
    tmpC.zeros();
    for(int j = 0; j < 2*r; j++){
      tmpC(j) = ldcorCpp(dx,dy,W.col(j));
    }
    S = mean(tmpC.subvec(r,2*r-1)) - mean(tmpC.subvec(0,r-1));
    C.col(i) = tmpC;
    for(int j = 0; j < bootn; j++){
      nullW = shuffle(W,0);
      for(int k = 0; k < 2*r; k++){
        nullC(k,i,j) = ldcorCpp(dx,dy,nullW.col(k));
        nullCC(k) = nullC(k,i,j);
        if(C(k,i) < nullC(k,i,j)){
          ep(k,i)++;
        }
      }
      nullS(i,j) = mean(nullCC.subvec(r,2*r-1)) - mean(nullCC.subvec(0,r-1));
      if(S < nullS(i,j)){
        gp(i)++;
      }
    }
  }
  ep = (ep + 1) / (bootn + 1);
  gp = (gp + 1) / (bootn + 1);
  return List::create(Named("modulator") = M,
		      Named("ldcor") = C,
		      Named("global.score") = S,
		      Named("each.p.value") = ep,
                      Named("global.p.value") = gp,
                      Named("null.global.score") = nullS,
                      Named("null.each.ldcor") = nullC);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List ldcorAllCpp(const arma::mat X, const arma::mat Y, const arma::mat Z,
		       const arma::mat M, const double delta, const int bootn, const int type){
  int n = X.n_rows;
  int ln = round(n*delta);
  int p = Y.n_cols;
  int q = Z.n_cols;
  int r = M.n_rows;
  double Zmin;
  double Zmax;
  double lambda;
  arma::vec normz(n);
  arma::vec dx = distCpp(X,1);
  arma::vec dy;
  arma::vec tmpC(r);
  arma::vec dd(n);
  arma::vec ww(n);
  arma::vec gp(p);
  arma::uvec id(n);
  arma::mat nZ(n,q);
  arma::mat W(n,r);
  arma::mat C(r,p);
  arma::mat ep(r,p);
  arma::mat nullW(n,r);
  arma::cube nullC(r,p,bootn,fill::zeros);
  arma::vec nullCC(r);
  for(int i = 0; i < q; i++){
    Zmin = min(Z.col(i));
    Zmax = max(Z.col(i));
    nZ.col(i) = (Z.col(i) - Zmin) / (Zmax - Zmin);
  }
  for(int i = 0; i < r; i++){
    for(int j = 0; j < n; j++){
      dd(j) = norm(nZ.row(j) - M.row(i), 1);
    }
    dd = sort(dd);
    lambda = dd(ln - 1);
    if(lambda == 0){
      lambda = 1e-3;
    }
    if(type == 1){
      ww = gaukerCpp(nZ, M.row(i), lambda);
    } else {
      ww = neikerCpp(nZ, M.row(i), ln);
    }
    W.col(i) = ww / sum(ww);
  }
  ep.zeros();
  for(int i = 0; i < p; i++){
    dy = distCpp(Y.col(i),1);
    tmpC.zeros();
    for(int j = 0; j < r; j++){
      tmpC(j) = ldcorCpp(dx,dy,W.col(j));
    }
    C.col(i) = tmpC;
    for(int j = 0; j < bootn; j++){
      nullW = shuffle(W,0);
      for(int k = 0; k < r; k++){
        nullC(k,i,j) = ldcorCpp(dx,dy,nullW.col(k));
        nullCC(k) = nullC(k,i,j);
        if(C(k,i) < nullC(k,i,j)){
          ep(k,i)++;
        }
      }
    }
  }
  ep = (ep + 1) / (bootn + 1);
  return List::create(Named("ldcor") = C,
		      Named("each.p.value") = ep,
                      Named("null.each.ldcor") = nullC);
}
