#include <Rcpp.h>
// using namespace Rcpp;

// Read thresholding
// Output: bool array with "true" for reads passing/above thresholding criteria
//

// thresholding, vectorised, ascii-based
// [[Rcpp::export("rcpp_threshold_reads")]]
std::vector<bool> rcpp_threshold_reads(std::vector<std::string> xm,  // merged normalised BAM XM fields
                                       std::string ctx_meth,         // methylated context string, e.g. "XZ". NON-EMPTY
                                       std::string ctx_unmeth,       // unmethylated context string, e.g. "xz". NON-EMPTY
                                       std::string ooctx_meth,       // methylated out-of-context string, e.g. "HU". Can be empty
                                       std::string ooctx_unmeth,     // unmethylated out-of-context string, e.g. "hu". Can be empty
                                       int min_n_ctx,                // minimum number of context bases in xm field
                                       double min_ctx_meth_frac,     // minimum fraction of methylated to total context bases (min context beta value)
                                       double max_ooctx_meth_frac)   // maximum fraction of methylated to total out-of-context bases (max out-of-context beta value)
{
  std::vector<bool> res (xm.size(), true);
  for (unsigned int x=0; x<xm.size(); x++) {
    // checking for the interrupt
    if (x & 1048575 == 0) Rcpp::checkUserInterrupt();
    
    unsigned int ascii_map [128] = {0};
    std::for_each(xm[x].begin(), xm[x].end(), [&ascii_map] (char const &c) {
      ascii_map[c]++;
    });
    
    unsigned int n_ctx_meth = 0;
    std::for_each(ctx_meth.begin(), ctx_meth.end(), [&n_ctx_meth, &ascii_map] (char const &c) {
      n_ctx_meth += ascii_map[c];
    });
    if (n_ctx_meth==0) {
      res[x] = false;
      continue;
    }
    
    unsigned int n_ctx_unmeth = 0;
    std::for_each(ctx_unmeth.begin(), ctx_unmeth.end(), [&n_ctx_unmeth, &ascii_map] (char const &c) {
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
    std::for_each(ooctx_meth.begin(), ooctx_meth.end(), [&n_ooctx_meth, &ascii_map] (char const &c) {
      n_ooctx_meth += ascii_map[c];
    });
    if (n_ooctx_meth>0) {
      unsigned int n_ooctx_unmeth = 0;
      std::for_each(ooctx_unmeth.begin(), ooctx_unmeth.end(), [&n_ooctx_unmeth, &ascii_map] (char const &c) {
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
rcpp_threshold_reads(c("..z..Zh..Zh..z..","..z..Zh..zh..z..","..z..Zh..ZH..z..","..x..Zh..xh..x.."),
                            "Z", "z", "XHU", "xhu", 2, 0.5, 0.1)

xm           <- bam$XM
ctx.meth     <- "Z"
ctx.unmeth   <- "z"
ooctx.meth   <- "XH"
ooctx.unmeth <- "xh"
min.n.ctx    <- 2
min.ctx.meth.frac   <- 0.5
max.ooctx.meth.frac <- 0.1
microbenchmark::microbenchmark(rcpp_threshold_reads(xm, ctx.meth, ctx.unmeth, ooctx.meth, ooctx.unmeth, min.n.ctx, min.ctx.meth.frac, max.ooctx.meth.frac), times=10)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_threshold_reads.cpp")

// #############################################################################
