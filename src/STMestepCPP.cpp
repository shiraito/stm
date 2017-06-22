// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <RcppNumerical.h>
#include "hpbcpp_cppestep.h"

using namespace Numer;

class STMfun: public MFuncGrad
{
private:
  const Rcpp::NumericVector mu_i;
  const arma::mat betas;
  const Rcpp::NumericVector doc_ct;
  const Rcpp::NumericMatrix siginvm;
public:
  STMfun(const Rcpp::NumericVector _mu_i, const arma::mat _beta_i
	 ,const Rcpp::NumericVector _doc_ct
	 ,const Rcpp::NumericMatrix _siginv) : mu_i(_mu_i), betas(_beta_i), doc_ct(_doc_ct), siginvm(_siginv) {}
  
  double f_grad(Constvec& eta, Refvec grad)
  {
    arma::vec etas(eta.size());
    arma::vec doc_cts = Rcpp::as<arma::vec>(doc_ct);
    arma::vec mus = Rcpp::as<arma::vec>(mu_i);
    arma::mat siginvs = Rcpp::as<arma::mat>(siginvm);

    arma::rowvec expeta(etas.size()+1); 
    expeta.fill(1);
    int neta = etas.n_elem;
    for(int j=0; j <neta;  j++){
      expeta(j) = exp(etas(j));
    }
    double ndoc = sum(doc_cts);
    double part1 = arma::as_scalar(log(expeta*betas)*doc_cts - ndoc*log(sum(expeta)));
    arma::vec diff = etas - mus;
    double part2 = .5*arma::as_scalar(diff.t()*siginvs*diff);
    double out = part2 - part1;

    arma::mat betas_changed = betas;
    betas_changed.each_col() % trans(expeta);
    arma::vec part1_grad = betas_changed*(doc_cts/arma::trans(sum(betas_changed,0))) - (sum(doc_cts)/sum(expeta))*trans(expeta);
    arma::vec part2_grad = siginvs*(etas - mus);
    part1_grad.shed_row(neta);

    for (int j = 0; j < part1_grad.size(); ++j) {
      grad[j] = part2_grad(j) - part1_grad(j);
    }
    
    return out;
  }
};

// [[Rcpp::export]]
Rcpp::List estep_loop(Rcpp::List documents
		      ,Rcpp::IntegerVector beta_index
		      ,Rcpp::NumericMatrix lambda
		      ,Rcpp::NumericMatrix mu
		      ,bool update_mu
		      ,Rcpp::List beta
		      ,Rcpp::NumericMatrix siginv
		      ,Rcpp::NumericVector sigmaentropy
		      ,int N
		      ,int V
		      ,int K
		      ,int A
		      ){

  using namespace Rcpp;
  using namespace arma;

  Rcpp::NumericVector bound(N);
  Rcpp::NumericMatrix sigma_ss(K-1, K-1);

  cube beta_ss(K, V, A);
  beta_ss.fill(0);

  // if parallelization is desired, add pragma statements here
  for (int i = 0; i < N; ++i) {
    SEXP l_doc = documents[i];
    NumericMatrix doc(l_doc);
    uvec words = as<uvec>(wrap(doc( 0, _)));
    words = words - ones<uvec>(words.n_elem);
    int aspect = beta_index[i] - 1;
    NumericVector init = lambda( i, _);
    Eigen::Map<Eigen::VectorXd> eta_map = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(init);
    Eigen::VectorXd eta(eta_map);
    NumericVector mu_i;
    if (!update_mu) {
      mu_i = as<NumericVector>(mu);
    } else {
      mu_i = mu( _, i);
    }
    SEXP l_beta = beta[aspect];
    mat beta_mat = as<mat>(l_beta);
    mat beta_i = beta_mat.cols(words);

    NumericVector doc_ct = doc( 1, _);
    int Ndoc = sum(doc_ct);
    
    STMfun nll(mu_i, beta_i, doc_ct, siginv);

    double fopt;

    int status = optim_lbfgs(nll, eta, fopt);

    hpbcpp_cppestep(eta, beta_i, doc_ct, mu_i,
    		    siginv, sigmaentropy,
    		    bound, lambda, sigma_ss,
    		    beta_ss, i, aspect, words);

  }

  return Rcpp::List::create(
			    Rcpp::Named("sigma_ss") = sigma_ss,
			    Rcpp::Named("beta_ss") = beta_ss,
			    Rcpp::Named("bound") = bound,
			    Rcpp::Named("lambda") = lambda
			    );
  
}
