
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector qD_MLE(NumericVector q,NumericVector ai){
  const int length = q.size();
  const int S = ai.size();
  NumericVector Q(length);
  NumericVector temp(S);
  for(int j = 0; j<length;j++){
    for(int i = 0 ; i<S;i++){
      temp[i] = pow(ai[i],q[j]);
    }
    Q[j] = (1-sum(temp)) / (q[j]-1);
  }
  return Q;
}

// [[Rcpp::export]]
NumericVector qDSub2(double q, int f1, double A, const int n){
  NumericVector tmp(n);
  for(int r=0;r<=(n-1);r++){
    tmp[r] = Rf_choose(q-1, r)*pow(A-1, r);
  }
  return tmp;
}

// [[Rcpp::export]]
NumericVector qDSub(double q, NumericVector Xi,const int n){
  const int len = Xi.size();
  NumericVector Out(len);
  for(int k=0;k<=(len-1);k++){
    int z = Xi[k];
    NumericVector tmp(n-z+1);
    for(int l=0;l<=(n-z);l++){
      tmp[l] = Rf_choose(l-q,l)*exp(Rf_lchoose(n-l-1, z-1)-Rf_lchoose(n, z));
    }
    Out[k] = sum(tmp);
  }
  return Out;
}




