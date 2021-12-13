#include <Rcpp.h>
// using namespace Rcpp;

// Read thresholding
// Output: bool vector with "true" for reads passing/above thresholding criteria
//
// This one would def benefit from:
// [ ] SIMD
// [ ] fewer branches
// [ ] FALSE as a default

// thresholding, vectorised, ascii-based
// [[Rcpp::export("rcpp_threshold_reads")]]
std::vector<bool> rcpp_threshold_reads(Rcpp::DataFrame &df,                     // BAM data
                                       std::string ctx_meth,                    // methylated context string, e.g. "XZ". NON-EMPTY
                                       std::string ctx_unmeth,                  // unmethylated context string, e.g. "xz". NON-EMPTY
                                       std::string ooctx_meth,                  // methylated out-of-context string, e.g. "HU". Can be empty
                                       std::string ooctx_unmeth,                // unmethylated out-of-context string, e.g. "hu". Can be empty
                                       unsigned int min_n_ctx,                  // minimum number of context bases in xm field
                                       double min_ctx_meth_frac,                // minimum fraction of methylated to total context bases (min context beta value)
                                       double max_ooctx_meth_frac)              // maximum fraction of methylated to total out-of-context bases (max out-of-context beta value)
{
  // Rcpp::CharacterVector xm = df["XM"];                                          // merged refspaced template XMs
  Rcpp::XPtr<std::vector<std::string>> xm((SEXP)df.attr("xm_xptr"));            // merged refspaced template XMs, as a pointer to std::vector<std::string>
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  std::vector<bool> res (xm->size(), true);
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
    if (n_ctx_meth==0) {
      res[x] = false;
      continue;
    }
    
    unsigned int n_ctx_unmeth = 0;
    std::for_each(ctx_unmeth.begin(), ctx_unmeth.end(), [&n_ctx_unmeth, &ascii_map] (unsigned int const &c) {
      n_ctx_unmeth += ascii_map[c];
    });
    unsigned int n_ctx_all = n_ctx_meth + n_ctx_unmeth;
    if (n_ctx_all<min_n_ctx) {
      res[x] = false;
      continue;
    }
    double ctx_meth_frac = (double)n_ctx_meth / n_ctx_all;
    if (ctx_meth_frac<min_ctx_meth_frac) {
      res[x] = false;
      continue;
    }
    
    unsigned int n_ooctx_meth = 0;
    std::for_each(ooctx_meth.begin(), ooctx_meth.end(), [&n_ooctx_meth, &ascii_map] (unsigned int const &c) {
      n_ooctx_meth += ascii_map[c];
    });
    if (n_ooctx_meth>0) {
      unsigned int n_ooctx_unmeth = 0;
      std::for_each(ooctx_unmeth.begin(), ooctx_unmeth.end(), [&n_ooctx_unmeth, &ascii_map] (unsigned int const &c) {
        n_ooctx_unmeth += ascii_map[c];
      });
      
      unsigned int n_ooctx_all = n_ooctx_meth + n_ooctx_unmeth;
      double ooctx_meth_frac = (double)n_ooctx_meth / n_ooctx_all;
      if (ooctx_meth_frac>max_ooctx_meth_frac) {
        res[x] = false;
      }
    }
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
