#ifndef __HDPCPP_CPPESTEP_H__
#define __HDPCPP_CPPESTEP_H__

void hpbcpp_cppestep(Eigen::VectorXd& eta,
		     arma::mat& betas,
		     Rcpp::NumericVector& doc_ctv,
		     Rcpp::NumericVector& muv,
		     Rcpp::NumericMatrix& siginvm,
		     Rcpp::NumericVector& sigmaentropym,
		     Rcpp::NumericVector& bound_out,
		     Rcpp::NumericMatrix& lambda_out,
		     Rcpp::NumericMatrix& sigma_ss_out,
		     arma::cube& beta_ss_out,
		     int i,
		     int aspect,
		     arma::uvec words);

#endif
