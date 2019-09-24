#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

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


//' Scale distance using a sigmoid function. The main purpose is to reduce the effects of small distances (those distances are most likely from random noise of ADT signals)
//' This is an updated function. This function will distingrish distances between negative and positive from those within negative/positive groups. Only distances within negative/positive groups will be scaled.
//' This function may take longer time than normal one...
//'
//' @param data data matrix of ADT.
//' @param n n for Sigmoid function.
//' @param k k for Sigmoid function.
//' @param c c is a vector that contains constant value that seaprate negative and positive for each ADT feature. By default, c = 1,1,...,1
//' @export
// [[Rcpp::export]]
NumericMatrix scaleDistUpdateCpp(NumericMatrix data, double n, double k, NumericVector c) {
  int dim_cell = data.ncol();
  NumericMatrix temp_dist_matrix(dim_cell, dim_cell);

  for(int i = 0; i < dim_cell - 1; i++)
  {
    for(int j = i + 1; j < dim_cell; j++)
    {
      NumericVector d1 = data( _ , i);
      NumericVector d2 = data( _ , j);
      NumericVector diff = abs(d1 - d2);
      NumericVector diff_scale = Sigmoid(diff,n,k);

      NumericVector min_value = pmin(d1,d2);
      NumericVector max_value = pmax(d1,d2);
      for(int m = 0; m < c.size(); m++)
      {
        if((min_value[m] < c[m]) && (max_value[m] > c[m]))
        {
          diff_scale[m] = diff[m];
        }
      }

      double value = EucNorm(diff_scale);
      temp_dist_matrix(i,j) = value;
      temp_dist_matrix(j,i) = value;
    }
  }
  return temp_dist_matrix;
}



//' objective function for gradient descnet method
//'
//' @param alpha current value of parameter alpha.
//' @param X distance vector.
//' @param Y distance vector.
//' @export
// [[Rcpp::export]]
double objectiveFunctionCpp(double alpha, NumericVector X, NumericVector Y) {
  double n = X.size();
  NumericVector Z = X * alpha - Y;
  double diff = sum(pow(Z, 2))/(2*n);
  return diff;
}

//' gradient function for gradient descnet method
//'
//' @param alpha current value of parameter alpha.
//' @param X distance vector.
//' @param Y distance vector.
//' @export
// [[Rcpp::export]]
double gradientFunctionCpp(double alpha, NumericVector X, NumericVector Y) {
  double n = X.size();
  double diff = sum(X * alpha - Y)/n;
  return diff;
}

//' gradient descnet method
//'
//' @param X distance vector.
//' @param Y distance vector.
//' @param alpha initial value of parameter alpha.
//' @param learning_rate learning rate of GD method
//' @param low_threshold the low threshold of GD
//' @param max_iter maximum iterations
//' @export
// [[Rcpp::export]]
double gradientDescentCpp(NumericVector X, NumericVector Y, double alpha, double learning_rate, double low_threshold, int max_iter) {
  double gradient;
  double alpha_last;
  double n = X.size();
  for(int iter = 1; iter < max_iter; iter++)
  {
    // calculate gradient
    gradient = sum(X * alpha - Y)/n;

    alpha_last = alpha;
    alpha = alpha - learning_rate * gradient;
    //Rcout << "iter =" << iter << "\talpha = " << alpha << "\talpha last = " << alpha_last << "\tgradient = " << gradient << "\tlearning rate = " << learning_rate << "\n";
    if(abs(gradient) < low_threshold)
    {
      break;
    }
    if((alpha < 0)||(alpha > 1))
    {
      if(abs(alpha_last) < abs(alpha)) {
        learning_rate = learning_rate/2;
        alpha_last = 0;
        alpha = 0;
      }
    }
  }
  return alpha;
}


