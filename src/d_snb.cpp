#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>


#define printArma(x) print(NumericVector(x.begin(), x.end()))
#define printArmaI(x, i) print(NumericVector(x.begin(), x.begin()+i))
#define printArmaRow(x, i) print(NumericVector(x.begin_row(i), x.end_row(i)))
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

NumericVector repeatVector(NumericVector& x, const int dest_size) {
  NumericVector res(dest_size);
  int n = x.length();
  for(int i = 0; i < dest_size; i++) {
    res[i] = x[i%n];
  }
  return res;
}

//' d_snb
//' @name d_snb
//' @param x numeric value or vector
//' @param size_param,prob_param,mu_param,size,prob,mu
//'     numeric value or vector,
//'     atleast two of these parameters must be specified
//' @param infinite integer value specifying the maximum number
//'     of terms to compute in the series, default: 100000
//' @param log_v,log whether the results should be log values,
//'     default: false
//' @importFrom Rcpp evalCpp
//' @useDynLib stochprofML
//' @export
// [[Rcpp::export]]
NumericVector d_snb(NumericVector& x,
                    Nullable<NumericVector> size_param = R_NilValue,
                    Nullable<NumericVector> prob_param = R_NilValue,
                    Nullable<NumericVector> mu_param = R_NilValue,
                    const int& infinite = 100000, const bool& log_v = false) {

    int n = x.length(), param_check = 0;
    NumericVector res(n);
    NumericVector mu, l, prob;

    // check for null parameters
    if(size_param.isNotNull()) {
      l = NumericVector(size_param);
      param_check++;
    }
    if(prob_param.isNotNull()) {
      prob = NumericVector(prob_param);
      param_check++;
    }
    if(mu_param.isNotNull()) {
      mu = NumericVector(mu_param);
      param_check++;
    }

    if (param_check < 2) {
      warning("Two values among mu, size and prob must be provided");
      return NumericVector(0);
    }
    int m = std::max({ l.length(), prob.length(), mu.length() });

    // inflate non-null parameters
    if(mu_param.isNotNull()) {
      mu = repeatVector(mu, m);
    }
    if(size_param.isNotNull()) {
      l = repeatVector(l, m);
    }
    if(prob_param.isNotNull()) {
       prob = repeatVector(prob, m);
    }


    // initialise null parameters
    if(mu_param.isNull()) {
      mu = l/prob - l;
    }
    if(size_param.isNull()) {
      l = (prob * mu) / (1 - prob);
    }
    if(prob_param.isNull()) {
      prob = l/(l + mu);
    }

    if( max(prob) - min(prob) == 0 ) {
      for(int i = 0; i < n; i++) {
        res[i] = R::dnbinom(x[i], sum(l), mean(prob), log_v);
      }
      return res;
    } else {
      arma::rowvec alpha = l;
      arma::rowvec p = prob, q = 1 - prob;
      double p1 = max(p);
      arma::rowvec Q = ((1 - p1) * p)/(q * p1);

      double l_R2 = arma::sum(-alpha % arma::log(1/Q));

      arma::colvec l_delta_6002(infinite);
      l_delta_6002.fill(NA_REAL);
      l_delta_6002[0] = log(1.0e-305);

      arma::mat A = arma::repmat(alpha, infinite, 1);
      arma::mat E = arma::repmat(1-Q, infinite, 1);

      // column-wise operations would be a faster approach
      for(int i = 0; i < infinite; i++) {
        // lambda expressions, awesome :)
        E.row(i).for_each( [&i](double& val) { val = std::pow(val, ((double)i+1)); } );
      }

      arma::mat ae(A.n_rows, A.n_cols);
      ae = A % E;

      // compute row sums
      arma::colvec L_XI_300 = arma::log(arma::sum(ae, 1) * 1.0e305);

      int i = 1, count = 2;
      while(i < count) {
        arma::vec ks2 = arma::linspace(1, i, i);
        arma::colvec X = L_XI_300.head(i) + arma::flipud(l_delta_6002.head(i));

        double k = floor(max(X)/700);
        int kp = k > 0;
        if(arma::max(X) > 700) {
          l_delta_6002(i) = log(  sum(  exp(k * log(1e-305) + X)  * 1e-305 )  )
            + (k - 1) * kp * log(1e305) + log(1e305) - log(i);
        } else {
          l_delta_6002(i) = log(1e-305) + log(sum(exp(X))) - log(i);
        }
        count = std::min(count + 3 * (l_delta_6002[i] > l_delta_6002[i-1]), infinite);
        i++;
      }


      l_delta_6002 = l_delta_6002.head(i);

      arma::colvec x_vec = x;
      double alphaSum = arma::sum(alpha);
      x_vec.for_each( [ &l_delta_6002, &alphaSum, &p1, &l_R2 ](double& val) {
        arma::colvec seqAlong = alphaSum + arma::regspace<arma::colvec>(0, 1, l_delta_6002.n_elem-1);
        seqAlong.for_each( [&val, &p1] (double& seqVal) { seqVal = R::dnbinom(val, seqVal, p1, 1); } );
        arma::colvec h1 = l_delta_6002 + seqAlong;
        // double k2 = (arma::max(h1) > -700) ? -1 :  std::ceil(arma::max(h1) / 700);
        double k2 = std::ceil(arma::max(h1) / 700);
        val = log(arma::sum(arma::exp(k2 * log(1e-305) + h1))) + k2 * log(1e305) + l_R2 + log(1e305);
      } );

      if(log_v)
        return NumericVector(x_vec.begin(), x_vec.end());
      else {
        arma::colvec t = arma::exp(x_vec);
        return NumericVector(t.begin(), t.end());
      }
    }
}
