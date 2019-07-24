#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//' Multiply a number by two
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
int timesTwo(int x) {
  return x * 2;
}

//' Sigmoid function
//'
//' @param a A Numeric Vector. original distance
//' @param n n for Sigmoid function.
//' @param k k for Sigmoid function.
//' @export
// [[Rcpp::export]]
NumericVector Sigmoid(NumericVector a, double n, double k) {
  int num = a.size();
  NumericVector b(num);
  for(int i = 0; i < num; i++)
  {
    b[i] = a[i]/(1 + 1/pow(2.72,(n*(a[i] - k))));
  }
  return b;
}

//' EucNorm function
//'
//' @param a A Numeric Vector.
//' @export
// [[Rcpp::export]]
double EucNorm(NumericVector a) {
  int n = a.size();
  double out = 0;

  for(int i = 0; i < n; i++)
  {
    out += pow(a[i], 2.0);
  }
  out = sqrt(out);
  return out;
}

//' Scale distance using a sigmoid function. The main purpose is to reduce the effects of small distances (those distances are most likely from random noise of ADT signals)
//'
//' @param data data matrix of ADT.
//' @param n n for Sigmoid function.
//' @param k k for Sigmoid function.
//' @export
// [[Rcpp::export]]
NumericMatrix scaleDistCpp(NumericMatrix data, double n, double k) {
  int dim_cell = data.ncol();
  NumericMatrix temp_dist_matrix(dim_cell, dim_cell);

  for(int i = 0; i < dim_cell - 1; i++)
  {
    for(int j = i + 1; j < dim_cell; j++)
    {
      NumericVector d1 = data( _ , i);
      NumericVector d2 = data( _ , j);
      NumericVector diff = abs(d1 - d2);
      diff = Sigmoid(diff,n,k);
      double value = EucNorm(diff);
      temp_dist_matrix(i,j) = value;
      temp_dist_matrix(j,i) = value;
    }
  }
  return temp_dist_matrix;
}


