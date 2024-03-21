#include <Rcpp.h>
#include "epialleleR.h"
// using namespace Rcpp;

// Read thresholding
// Output: bool vector with "true" for reads passing/above thresholding criteria
//
// This one would def benefit from:
// [ ] SIMD
// [ ] fewer branches
// [x] FALSE as a default

// thresholding, vectorised, ascii-based
// [[Rcpp::export("rcpp_threshold_reads")]]
std::vector<bool> rcpp_threshold_reads(Rcpp::DataFrame &df,                     // BAM data
                                       const std::string ctx_meth,              // methylated context string, e.g. "XZ". NON-EMPTY
                                       const std::string ctx_unmeth,            // unmethylated context string, e.g. "xz". NON-EMPTY
                                       const std::string ooctx_meth,            // methylated out-of-context string, e.g. "HU". Can be empty
                                       const std::string ooctx_unmeth,          // unmethylated out-of-context string, e.g. "hu". Can be empty
                                       const unsigned int min_n_ctx,            // minimum number of context bases in xm field
                                       const double min_ctx_meth_frac,          // minimum fraction of methylated to total context bases (min context beta value)
                                       const double max_ooctx_meth_frac)        // maximum fraction of methylated to total out-of-context bases (max out-of-context beta value)
{
  Rcpp::XPtr<std::vector<std::string>> seqxm((SEXP)df.attr("seqxm_xptr"));      // merged refspaced packed template SEQXMs, as a pointer to std::vector<std::string>
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  std::vector<bool> res (seqxm->size(), false);
  for (unsigned int x=0; x<seqxm->size(); x++) {
    // checking for the interrupt
    if ((x & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();
    
    unsigned int ctx_map[16] = {0};
    const char* seqxm_x = seqxm->at(templid[x]).c_str();                        // seqxm->at(templid[x]) is a reference to a corresponding SEQXM string
    const unsigned int size_x = seqxm->at(templid[x]).size();                   // length of the current read
    for (unsigned int i=0; i<size_x; i++) {                                     // char by char - it's faster this way than using std::string in the cycle
      ctx_map[unpack_ctx_idx(seqxm_x[i])]++;                                    // extract lower 4 bits (XM) and count them;
    }
    
    unsigned int n_ctx_meth = 0;
    std::for_each(ctx_meth.begin(), ctx_meth.end(), [&n_ctx_meth, &ctx_map] (unsigned int const &c) {
      n_ctx_meth += ctx_map[ctx_to_idx(c)];
    });
    if (n_ctx_meth==0) continue;                                                // next read if no methylated context bases
    
    unsigned int n_ctx_unmeth = 0;
    std::for_each(ctx_unmeth.begin(), ctx_unmeth.end(), [&n_ctx_unmeth, &ctx_map] (unsigned int const &c) {
      n_ctx_unmeth += ctx_map[ctx_to_idx(c)];
    });
    unsigned int n_ctx_all = n_ctx_meth + n_ctx_unmeth;
    if (n_ctx_all<min_n_ctx) continue;                                          // next read if total number of context bases is less than min_n_ctx
    
    double ctx_meth_frac = (double)n_ctx_meth / n_ctx_all;
    if (ctx_meth_frac<min_ctx_meth_frac) continue;                              // next read if average context beta is less than min_ctx_meth_frac
    
    unsigned int n_ooctx_meth = 0;
    std::for_each(ooctx_meth.begin(), ooctx_meth.end(), [&n_ooctx_meth, &ctx_map] (unsigned int const &c) {
      n_ooctx_meth += ctx_map[ctx_to_idx(c)];
    });
    if (n_ooctx_meth>0) {
      unsigned int n_ooctx_unmeth = 0;
      std::for_each(ooctx_unmeth.begin(), ooctx_unmeth.end(), [&n_ooctx_unmeth, &ctx_map] (unsigned int const &c) {
        n_ooctx_unmeth += ctx_map[ctx_to_idx(c)];
      });
      
      unsigned int n_ooctx_all = n_ooctx_meth + n_ooctx_unmeth;
      double ooctx_meth_frac = (double)n_ooctx_meth / n_ooctx_all;
      if (ooctx_meth_frac>max_ooctx_meth_frac) continue;                        // next read if average out-of-context beta is higher than max_ooctx_meth_frac
    }
    
    res[x] = true;                                                              // read has passed all thresholds
  }
  
  return res;
}


// test code in R
//

/*** R
bam          <- bam
ctx.meth     <- "Z"
ctx.unmeth   <- "z"
ooctx.meth   <- "XH"
ooctx.unmeth <- "xh"
min.n.ctx    <- 2
min.ctx.meth.frac   <- 0.5
max.ooctx.meth.frac <- 0.1
microbenchmark::microbenchmark(rcpp_threshold_reads(bam, ctx.meth, ctx.unmeth, ooctx.meth, ooctx.unmeth, min.n.ctx, min.ctx.meth.frac, max.ooctx.meth.frac), times=10)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_threshold_reads.cpp")

// #############################################################################
