#include <RcppArmadillo.h>
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
const int r, const double delta, const int nperm){
  int n = X.n_rows;
  int ln = round(n*delta);
  int q = Z.n_cols;
  double Zmin;
  double Zmax;
  double lambda;
  arma::vec S(1);
  arma::vec nS(1);
  arma::vec normZ(n);
  arma::vec dx = distCpp(X,1);
  arma::vec dy = distCpp(Y,1);
  arma::vec dd(n);
  arma::vec ww(n);
  double gp;
  arma::uvec id(n);
  arma::mat nZ(n,q);
  arma::mat M(n,q);
  arma::mat W(n,n);
  arma::vec C(2*r);
  arma::vec ep(2*r);
  arma::mat nullW(n,2*r);
  arma::vec nullS(nperm,fill::zeros);
  arma::mat nullC(2*r,nperm,fill::zeros);
  arma::vec nullCC(2*r);
  for(int i = 0; i < q; i++){
    Zmin = min(Z.col(i));
    Zmax = max(Z.col(i));
    nZ.col(i) = (Z.col(i) - Zmin) / (Zmax - Zmin);
  }
  for(int i = 0; i < n; i++){
    normZ(i) = norm(nZ.row(i),1);
  }
  id = sort_index(normZ);
  for(int i = 0; i < n; i++){
    M.row(i) = nZ.row(id(i));
  }
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      dd(j) = norm(nZ.row(j) - M.row(i), 1);
    }
    dd = sort(dd);
    lambda = dd(ln - 1);
    if(lambda == 0){
      lambda = 1e-3;
    }
    ww = gaukerCpp(nZ, M.row(i), lambda);
    W.col(i) = ww / sum(ww);
  }
  ep.zeros();
  gp = 0;
  C.zeros();
  int count = 0;
  for(int j = 0; j < r; j++){
    C(count) = ldcorCpp(dx,dy,W.col(j));
    count ++;
  }
  for(int j = n - r; j < n; j++){
    C(count) = ldcorCpp(dx,dy,W.col(j));
    count ++;
  }
  S(0) = mean(C.subvec(r,2*r-1)) - mean(C.subvec(0,r-1));
  S = abs(S);
  for(int j = 0; j < nperm; j++){
    nullW = shuffle(W,0);
    count = 0;
    for(int k = 0; k < r; k++){
      nullC(count,j) = ldcorCpp(dx,dy,nullW.col(k));
      nullCC(count) = nullC(count,j);
      if(C(count) < nullC(count,j)){
        ep(count)++;
      }
      count ++;
    }
    for(int k = n - r; k < n; k++){
      nullC(count,j) = ldcorCpp(dx,dy,nullW.col(k));
      nullCC(count) = nullC(count,j);
      if(C(count) < nullC(count,j)){
        ep(count)++;
      }
      count ++;
    }
    nS(0) = mean(nullCC.subvec(r,2*r-1)) - mean(nullCC.subvec(0,r-1));
    nS = abs(nS);
    nullS(j) = nS(0);
    if(S(0) < nS(0)){
      gp++;
    }
  }
  ep = (ep + 1) / (nperm + 1);
  gp = (gp + 1) / (nperm + 1);
  return List::create(Named("modulator") = M,
                      Named("weights") = W,
                      Named("ldcor") = C,
                      Named("global.score") = S,
                      Named("each.p.value") = ep,
                      Named("global.p.value") = gp,
                      Named("null.global.score") = nullS,
                      Named("null.each.ldcor") = nullC);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List gimletCpp2(const arma::mat X, const arma::mat Y, const arma::rowvec Z, const int nperm){
  int n = X.n_rows;
  double sumZ1;
  double sumZ0;
  arma::vec S(1);
  arma::vec nS(1);
  arma::vec dx = distCpp(X,1);
  arma::vec dy = distCpp(Y,1);
  arma::vec tmpC(2);
  arma::vec dd(n);
  arma::vec ww(n);
  arma::vec gp(1);
  arma::mat W(n,2);
  arma::vec C(2);
  arma::vec ep(2);
  arma::mat nullW(n,2);
  arma::vec nullS(nperm,fill::zeros);
  arma::mat nullC(2,nperm,fill::zeros);
  arma::vec nullCC(2);
  ep.zeros();
  gp.zeros();
  
  sumZ1 = sum(Z==1);
  sumZ0 = sum(Z==0);
  for(int i = 0; i < n; i++){
    W(i,0) = 0;
    W(i,1) = 0;
    if(Z(i) == 1){
      W(i,1) = 1/sumZ1;
    } else if(Z(i) == 0){
      W(i,0) = 1/sumZ0;
    }
  }
  tmpC(0) = ldcorCpp(dx,dy,W.col(0));
  tmpC(1) = ldcorCpp(dx,dy,W.col(1));
  S(0) = tmpC(1) - tmpC(0);
  S = abs(S);
  C = tmpC;
  for(int j = 0; j < nperm; j++){
    nullW = shuffle(W,0);
    nullC(0,j) = ldcorCpp(dx,dy,nullW.col(0));
    nullCC(0) = nullC(0,j);
    if(C(0) < nullC(0,j)){
      ep(0)++;
    }
    nullC(1,j) = ldcorCpp(dx,dy,nullW.col(1));
    nullCC(1) = nullC(1,j);
    if(C(1) < nullC(1,j)){
      ep(1)++;
    }
    nS(0) = nullCC(1) - nullCC(0);
    nS = abs(nS);
    nullS(j) = nS(0);
    if(S(0) < nS(0)){
      gp(0)++;
    }
  }
  ep = (ep + 1) / (nperm + 1);
  gp = (gp + 1) / (nperm + 1);
  return List::create(Named("weights") = W,
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
const arma::mat M, const double delta, const int nperm){
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
  arma::cube nullC(r,p,nperm,fill::zeros);
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
    ww = gaukerCpp(nZ, M.row(i), lambda);
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
    for(int j = 0; j < nperm; j++){
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
  ep = (ep + 1) / (nperm + 1);
  return List::create(Named("ldcor") = C,
                      Named("each.p.value") = ep,
                      Named("null.each.ldcor") = nullC);
}
