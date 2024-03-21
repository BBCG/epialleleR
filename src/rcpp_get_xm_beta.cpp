#include <Rcpp.h>
#include "epialleleR.h"
// using namespace Rcpp;

// Parses XM tags and outputs average beta value according to context.
//

// fast, vectorised
// [[Rcpp::export("rcpp_get_xm_beta")]]
std::vector<double> rcpp_get_xm_beta(Rcpp::DataFrame &df,                       // BAM data
                                     const std::string ctx_meth,                // methylated context string, e.g. "XZ". NON-EMPTY
                                     const std::string ctx_unmeth)              // unmethylated context string, e.g. "xz". NON-EMPTY
{
  Rcpp::XPtr<std::vector<std::string>> seqxm((SEXP)df.attr("seqxm_xptr"));      // merged refspaced packed template SEQXMs, as a pointer to std::vector<std::string>
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  std::vector<double> res (seqxm->size(), 0);
  for (unsigned int x=0; x<seqxm->size(); x++) {
    // checking for the interrupt
    if ((x & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();
    
    unsigned int ctx_map[16] = {0};
    const char* seqxm_x = seqxm->at(templid[x]).c_str();                        // seqxm->at(templid[x]) is a reference to a corresponding SEQXM string
    const unsigned int size_x = seqxm->at(templid[x]).size();                   // length of the current read
    for (unsigned int i=0; i<size_x; i++) {                                     // char by char - it's faster this way than using std::string in the cycle
      ctx_map[seqxm_x[i] & 15]++;                                               // extract lower 4 bits (XM) and count them;
    }
    
    unsigned int n_ctx_meth = 0;
    std::for_each(ctx_meth.begin(), ctx_meth.end(), [&n_ctx_meth, &ctx_map] (unsigned int const &c) {
      n_ctx_meth += ctx_map[ctx_to_idx(c)];
    });
    unsigned int n_ctx_unmeth = 0;
    std::for_each(ctx_unmeth.begin(), ctx_unmeth.end(), [&n_ctx_unmeth, &ctx_map] (unsigned int const &c) {
      n_ctx_unmeth += ctx_map[ctx_to_idx(c)];
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
