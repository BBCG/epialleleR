#include <Rcpp.h>
// using namespace Rcpp;

// Parses XM tags and outputs average beta value according to context.
//

// fast, vectorised
// [[Rcpp::export("rcpp_get_xm_beta")]]
std::vector<double> rcpp_get_xm_beta(Rcpp::DataFrame &df,                       // BAM data
                                     std::string ctx_meth,                      // methylated context string, e.g. "XZ". NON-EMPTY
                                     std::string ctx_unmeth)                    // unmethylated context string, e.g. "xz". NON-EMPTY
{
  // Rcpp::CharacterVector xm = bam["XM"];                                         // template XM
  Rcpp::XPtr<std::vector<std::string>> xm((SEXP)df.attr("xm_xptr"));            // merged refspaced template XMs, as a pointer to std::vector<std::string>
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  std::vector<double> res (xm->size(), 0);
  for (unsigned int x=0; x<xm->size(); x++) {
    // checking for the interrupt
    if ((x & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();
    
    unsigned int ascii_map [128] = {0};
    std::for_each(xm->at(templid[x]).begin(), xm->at(templid[x]).end(), [&ascii_map] (unsigned int const &c) {
      ascii_map[c]++;
    });
    unsigned int n_ctx_meth = 0;
    std::for_each(ctx_meth.begin(), ctx_meth.end(), [&n_ctx_meth, &ascii_map] (unsigned int const &c) {
      n_ctx_meth += ascii_map[c];
    });
    unsigned int n_ctx_unmeth = 0;
    std::for_each(ctx_unmeth.begin(), ctx_unmeth.end(), [&n_ctx_unmeth, &ascii_map] (unsigned int const &c) {
      n_ctx_unmeth += ascii_map[c];
    });
    unsigned int n_ctx_all = n_ctx_meth + n_ctx_unmeth;
    if (n_ctx_all==0) n_ctx_all=1;
    res[x] = (double)n_ctx_meth / n_ctx_all;
  }
  
  return res;
}


// test code in R
//

/*** R
# ctx.meth   <- "ZX"
# ctx.unmeth <- "zx"
# microbenchmark::microbenchmark(rcpp_get_xm_beta(bam, ctx.meth, ctx.unmeth), times=10)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_get_xm_beta.cpp")

// #############################################################################
