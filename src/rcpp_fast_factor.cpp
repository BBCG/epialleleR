#include <Rcpp.h>
using namespace Rcpp;

// Fast factor from https://gallery.rcpp.org/articles/fast-factor-generation/
// Written by Kevin Ushey
// I won't need it when rcpp_merge_ends will be vectorised

template <int RTYPE>
IntegerVector rcpp_fast_factor_template( const Vector<RTYPE>& x ) {
  Vector<RTYPE> levs = sort_unique(x);
  IntegerVector out = match(x, levs);
  out.attr("levels") = as<CharacterVector>(levs);
  out.attr("class") = "factor";
  return out;
}

// [[Rcpp::export]]
SEXP rcpp_fast_factor( SEXP x ) {
  switch( TYPEOF(x) ) {
  case INTSXP: return rcpp_fast_factor_template<INTSXP>(x);
  case REALSXP: return rcpp_fast_factor_template<REALSXP>(x);
  case STRSXP: return rcpp_fast_factor_template<STRSXP>(x);
  }
  return R_NilValue;
}

// test code in R
//

/*** R
rcpp_fast_factor(c("s","p","e","e","e","e","d"))
*/

// Sourcing:
// Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_fast_factor.cpp")