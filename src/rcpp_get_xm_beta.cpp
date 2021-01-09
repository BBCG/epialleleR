#include <Rcpp.h>
using namespace Rcpp;

// Parses XM tags and outputs average beta value according to context.
//

// fast, vectorised
// [[Rcpp::export("rcpp_get_xm_beta")]]
std::vector<double> rcpp_get_xm_beta(std::vector<std::string> xm,  // merged normalised BAM XM fields
                                     std::string ctx_meth,         // methylated context string, e.g. "XZ". NON-EMPTY
                                     std::string ctx_unmeth)       // unmethylated context string, e.g. "xz". NON-EMPTY
{
  std::vector<double> res (xm.size(), 0);
  for (unsigned int x=0; x<xm.size(); x++) {
    unsigned int ascii_map [128] = {0};
    std::for_each(xm[x].begin(), xm[x].end(), [&ascii_map] (char const &c) {
      ascii_map[c]++;
    });
    unsigned int n_ctx_meth = 0;
    std::for_each(ctx_meth.begin(), ctx_meth.end(), [&n_ctx_meth, &ascii_map] (char const &c) {
      n_ctx_meth += ascii_map[c];
    });
    unsigned int n_ctx_unmeth = 0;
    std::for_each(ctx_unmeth.begin(), ctx_unmeth.end(), [&n_ctx_unmeth, &ascii_map] (char const &c) {
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
rcpp_get_xm_beta(c("..z..Zh..Zh..z..","..z..Zh..zh..z..","..z..Zh..XH..z..","..x..Zh..xh..x..","..h..Hh..hh..H.."), "ZX", "zx")
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_get_xm_beta.cpp")

// #############################################################################
