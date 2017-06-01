// [[Rcpp::depends(RcppArmadillo,RcppEigen)]]

#include <time.h>
#include "RcppArmadillo.h"
#include "RcppEigen.h"

// [[Rcpp::depends(RcppArmadillo,RcppEigen)]]

// // [[Rcpp::export]]
// SEXP hpbcpp(SEXP eta,
//             SEXP beta,
//             SEXP doc_ct,
//             SEXP mu,
//             SEXP siginv,
//             SEXP sigmaentropy){
 
//    Rcpp::NumericVector etav(eta); 
//    arma::vec etas(etav.begin(), etav.size(), false);
//    Rcpp::NumericMatrix betam(beta);
//    arma::mat betas(betam.begin(), betam.nrow(), betam.ncol());
//    Rcpp::NumericVector doc_ctv(doc_ct);
//    arma::vec doc_cts(doc_ctv.begin(), doc_ctv.size(), false);
//    Rcpp::NumericVector muv(mu);
//    arma::vec mus(muv.begin(), muv.size(), false);
//    Rcpp::NumericMatrix siginvm(siginv);
//    arma::mat siginvs(siginvm.begin(), siginvm.nrow(), siginvm.ncol(), false);
//    Rcpp::NumericVector sigmaentropym(sigmaentropy);
//    arma::vec entropy(sigmaentropym);

//    //Performance Nots from 3/6/2015
//    //  I tried a few different variants and benchmarked this one as roughly twice as
//    //  fast as the R code for a K=100 problem.  Key to performance was not creating
//    //  too many objects and being selective in how things were flagged as triangular.
//    //  Some additional notes in the code below.
//    //
//    //  Some things this doesn't have or I haven't tried
//    //  - I didn't tweak the arguments much.  sigmaentropy is a double, and I'm still
//    //    passing beta in the same way.  I tried doing a ", false" for beta but it didn't
//    //    change much so I left it the same as in gradient.  
//    //  - I tried treating the factors for doc_cts and colSums(EB) as a diagonal matrix- much slower.
   
//    //  Haven't Tried/Done
//    //  - each_row() might be much slower (not sure but arma is column order).  Maybe transpose in place?
//    //  - depending on costs there are some really minor calculations that could be precomputed: 
//    //     - sum(doc_ct)
//    //     - sqrt(doc_ct)
   
//    //  More on passing by reference here:
//    //  - Hypothetically we could alter beta (because hessian is last thing we do) however down
//    //    the road we may want to explore treating nonPD hessians by optimization at which point
//    //    we would need it again.

   
//    arma::colvec expeta(etas.size()+1); 
//    expeta.fill(1);
//    int neta = etas.size(); 
//    for(int j=0; j <neta;  j++){
//      expeta(j) = exp(etas(j));
//    }
//    arma::vec theta = expeta/sum(expeta);


//    // profiling
//    clock_t start = clock();
//    clock_t end = clock();
   
//    //create a new version of the matrix so we can mess with it
//    arma::mat EB(betam.begin(), betam.nrow(), betam.ncol());
//    //multiply each column by expeta
//    EB.each_col() %= expeta; //this should be fastest as its column-major ordering
  
//    //divide out by the column sums
//    EB.each_row() %= arma::trans(sqrt(doc_cts))/sum(EB,0);
    
//    //Combine the pieces of the Hessian which are matrices
//    start = clock();
//    double sum_doc_cts = sum(doc_cts);
//    end = clock();
//    Rcpp::Rcout << "double sum_doc_cts = sum(doc_cts): " << (double)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;
//    start = clock();
//    arma::mat hess = EB*EB.t();
//    end = clock();
//    Rcpp::Rcout << "arma::mat hess = EB*EB.t(): " << (double)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;
//    start = clock();
//    arma::mat crossprod_theta = theta*theta.t();
//    end = clock();
//    Rcpp::Rcout << "arma::mat crossprod_theta = theta*theta.t(): " << (double)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;
//    start = clock();
//    crossprod_theta *= sum_doc_cts;
//    end = clock();
//    Rcpp::Rcout << "crossprod_theta *= sum_doc_cts: " << (double)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;
//    start = clock();
//    hess -= crossprod_theta;
//    end = clock();
//    Rcpp::Rcout << "hess -= crossprod_theta: " << (double)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl << std::endl;
  
//    //we don't need EB any more so we turn it into phi
//    EB.each_row() %= arma::trans(sqrt(doc_cts));

//    //Now alter just the diagonal of the Hessian
//    hess.diag() -= sum(EB,1) - sum(doc_cts)*theta;
//    //Drop the last row and column
//    hess.shed_row(neta);
//    hess.shed_col(neta);
//    //Now we can add in siginv
//    hess = hess + siginvs;
//    //At this point the Hessian is complete.

   
//    //This next bit of code is from http://arma.sourceforge.net/docs.html#logging
//    //It basically keeps arma from printing errors from chol to the console.
//    std::ostream nullstream(0);
//    arma::set_stream_err2(nullstream);
   
//    ////
//    //Invert via cholesky decomposition
//    ////
//    //Start by initializing an object
//    arma::mat nu = arma::mat(hess.n_rows, hess.n_rows);
//    //This version of chol generates a boolean which tells us if it failed.
//    bool worked = arma::chol(nu,hess);
//    if(!worked) {
//      //It failed!  Oh Nos.
//      // So the matrix wasn't positive definite.  In practice this means that it hasn't
//      // converged probably along some minor aspect of the dimension.

//      //Here we make it positive definite through diagonal dominance
//      arma::vec dvec = hess.diag();
//      //find the magnitude of the diagonal 
//      arma::vec magnitudes = sum(abs(hess), 1) - abs(dvec);
//      //iterate over each row and set the minimum value of the diagonal to be the magnitude of the other terms
//      int Km1 = dvec.size();
//      for(int j=0; j < Km1;  j++){
//        if(arma::as_scalar(dvec(j)) < arma::as_scalar(magnitudes(j))) dvec(j) = magnitudes(j); //enforce diagonal dominance 
//      }
//      //overwrite the diagonal of the hessian with our new object
//      hess.diag() = dvec;
//      //that was sufficient to ensure positive definiteness so we now do cholesky
//      nu = arma::chol(hess);

//    }
//    //compute 1/2 the determinant from the cholesky decomposition
//    double detTerm = -sum(log(nu.diag()));

//    //Now finish constructing nu
//    nu = arma::inv(arma::trimatu(nu));
//    nu = nu * nu.t(); //trimatu doesn't do anything for multiplication so it would just be timesink to signal here.
   
//    //Precompute the difference since we use it twice
//    arma::vec diff = etas - mus;
//    //Now generate the bound and make it a scalar
//    double bound = arma::as_scalar(log(arma::trans(theta)*betas)*doc_cts + detTerm - .5*diff.t()*siginvs*diff - entropy);

//    // Generate a return list that mimics the R output
//    return Rcpp::List::create(
//         Rcpp::Named("phis") = EB,
//         Rcpp::Named("eta") = Rcpp::List::create(Rcpp::Named("lambda")=etas, Rcpp::Named("nu")=nu),
//         Rcpp::Named("bound") = bound
//         );
// }


// [[Rcpp::export]]
Rcpp::List estep_loop(Rcpp::List documents
		      ,Rcpp::IntegerVector beta_index
		      ,Rcpp::NumericMatrix lambda_old
		      ,Rcpp::NumericMatrix mu
		      ,bool update_mu
		      ,Rcpp::List beta
		      ){

  using namespace Rcpp;
  using namespace arma;

  int N = documents.size();

  // if parallelization is desired, add pragma statements here
  for (int i = 0; i < N; ++i) {
    SEXP l_doc = documents[i];
    NumericMatrix doc(l_doc);
    uvec words = as<uvec>(doc( i, _));
    int aspect = beta_index[i];
    NumericVector init = lambda_old( i, _);
    NumericVector mu_i = mu( _, i);
    SEXP l_beta = beta[aspect];
    mat beta_mat = as<mat>(l_beta);
    mat beta_i = beta_mat.cols(words);
    
  }
  
}
