// [[Rcpp::depends(RcppArmadillo)]]
#include <fstream>
#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

using namespace Rcpp;
using namespace std;
using namespace arma;


class Params
{
public:
   arma::colvec beta    ;
   arma::colvec std_err    ;
   arma::colvec resi       ;
   arma::mat henssen;
};

Params fastLm(const arma::mat& X, const arma::colvec& y) { 
    int n = X.n_rows, k = X.n_cols;    
    arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
    //arma::colvec res  = y - X*coef; // residuals 
    // std.errors of coefficients
    //double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k); 
    arma::mat henssen = arma::pinv(arma::trans(X)*X);    
    Params results;
    results.henssen = henssen;// * var(res);
    results.beta = coef;
    return results;
}
// [[Rcpp::export]]
arma::mat calynew(arma::mat ynew,arma::mat vl) {
    int nobs = ynew.n_rows;
    int nt = ynew.n_cols;
    arma::mat Y(nobs*nt,1);
    arma::mat deltv(nt,nt);
    for (int i = 0; i < nobs; ++i)
    {      
        deltv = vl.rows(nt*i,nt*(i+1)-1);   
        Y.rows(nt*i,nt*(i+1)-1) = deltv * trans(ynew.row(i));
    }
    return Y;
}
 // [[Rcpp::export]]
arma::mat fast_mvlmmthreads(int thread,arma::vec Y,arma::mat gnew,arma::mat vl,arma::vec sth) {
    int gnewst;
    int st,end,range;  
    gnewst = sth(thread-1);
    end = sth(thread);
    range = sth(thread) -sth(thread-1);
    int nobs =gnew.n_rows;
    int nmar = range;
    int nt = vl.n_cols;
   
    //arma::mat X(nobs*nt,nt*2);
    arma::mat X(nobs*nt,nt);
    arma::mat result(nmar,(nt*2)+2);
    
    arma::mat deltv(nt,nt);
    //string outfile = Rcpp::as<string>(outfile_in);
    size_t ncount;
    char *cp;
    
    //ofstream writefile (outfile.c_str(), ofstream::out);
    //X.cols(0,nt-1) = Gnewb0;
  
   int count=0;
    for(size_t i=gnewst; i<end; ++i) {   
       
        for (int j = 0; j < nobs; ++j){
          //X(span(nt*j,nt*(j+1)-1),span(0,nt-1)) = vl.rows(nt*j,nt*(j+1)-1) * gnew(j,i);
          X.rows(nt*j,nt*(j+1)-1) = vl.rows(nt*j,nt*(j+1)-1) * gnew(j,i);
        }
        Params rgres;
        rgres = fastLm(X, Y);
        arma::mat vbeta= rgres.henssen;//(span(nt,2*nt-1),span(nt,2*nt-1));
        arma::mat coef(1,nt);
        coef = rgres.beta.tail_rows(nt);
        arma::mat wald = trans(coef)*inv(vbeta)*coef;
        result(count,(2*nt))=wald(0,0);
        result(count,(2*nt)+1)=count+gnewst;
        for(int jj = 0;jj < nt;jj++){
            result(count,jj)=coef(jj);
            result(count,jj+nt)=vbeta(jj,jj);
        }
        count++;
    }
    return result;
 }

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP MatMult(arma::mat A, arma::mat B){
    arma::mat C = A * B;
    return Rcpp::wrap(C);
}
// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;
    
    return Rcpp::wrap(C);
}
// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;
    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
Eigen::MatrixXd CLegendre(Eigen::VectorXd x,int s){
  int size = x.size();
  Eigen::VectorXd at(size);
  Eigen::VectorXd Le0(size);
  Eigen::VectorXd Le1(size);
  Eigen::MatrixXd Le(size,s);
  Le0 = Le0.setOnes(size);
  at = 2*(x.array() - x.minCoeff())/(x.maxCoeff()-x.minCoeff())-1;
  Le1 = at;
  Le.col(0)=at;
  if(s==1){
      return Le;
  }
  float con;
  float count=1.0;
  for(int i=1;i<s;i++) {
     con = (2*count+1)/(count+1);
     Le.col(i) =con*at;
     Le.col(i).array() *= Le1.array();
     con=count/(count+1);
     Le.col(i).array() -= con*Le0.array();
     Le0 = Le1;
     Le1 = Le.col(i);
     count += 1.0;
  }
    return Le;
}
