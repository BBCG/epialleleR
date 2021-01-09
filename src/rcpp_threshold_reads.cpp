#include <Rcpp.h>
using namespace Rcpp;

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
    unsigned int n_ooctx_meth = 0;
    std::for_each(ooctx_meth.begin(), ooctx_meth.end(), [&n_ooctx_meth, &ascii_map] (char const &c) {
      n_ooctx_meth += ascii_map[c];
    });
    unsigned int n_ooctx_unmeth = 0;
    std::for_each(ooctx_unmeth.begin(), ooctx_unmeth.end(), [&n_ooctx_unmeth, &ascii_map] (char const &c) {
      n_ooctx_unmeth += ascii_map[c];
    });
    
    unsigned int n_ctx_all = n_ctx_meth + n_ctx_unmeth;
    if (n_ctx_all==0) n_ctx_all=1;
    double ctx_meth_frac = (double)n_ctx_meth / n_ctx_all;
    
    unsigned int n_ooctx_all = n_ooctx_meth + n_ooctx_unmeth;
    if (n_ooctx_all==0) n_ooctx_all=1;
    double ooctx_meth_frac = (double)n_ooctx_meth / n_ooctx_all;
    
    if (n_ctx_all<min_n_ctx)
      res[x] = false;
    if (ctx_meth_frac<min_ctx_meth_frac)
      res[x] = false;
    if (ooctx_meth_frac>max_ooctx_meth_frac)
      res[x] = false;
  }
  return res;
}


// test code in R
//

/*** R
rcpp_threshold_reads(c("..z..Zh..Zh..z..","..z..Zh..zh..z..","..z..Zh..ZH..z..","..x..Zh..xh..x.."),
                            "Z", "z", "XHU", "xhu", 2, 0.5, 0.1)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_threshold_reads.cpp")

// #############################################################################


// // rcpp_count_any_of helper function. Returns number of elements that are equal
// // to one of the characters in the given character sequence
// // [[Rcpp::export]]
// int rcpp_count_any_of(std::string &str, std::string &seq) {
//   int res=0;
//   int pos = str.find_first_of(seq, 0);
//   while (pos != std::string::npos) {
//     res++;
//     pos = str.find_first_of(seq, pos+1);
//   }
//   return res;
// }

// test code in R
//
// rcpp_count_any_of("..z..Zh..Zh..z..", "Zz")
// rcpp_count_any_of("..z..Zh..Zh..z..", "Hh")
  

// // filtering, vectorised - slower than R!
// // [[Rcpp::export]]
// std::vector<bool> rcpp_filter_reads(std::vector<std::string> xm,  // merged normalised BAM XM fields
//                                            std::string ctx_meth,         // methylated context string, e.g. "XZ". NON-EMPTY
//                                            std::string ctx_unmeth,       // unmethylated context string, e.g. "xz". NON-EMPTY
//                                            std::string ooctx_meth,       // methylated out-of-context string, e.g. "HU". Can be empty
//                                            std::string ooctx_unmeth,     // unmethylated out-of-context string, e.g. "hu". Can be empty
//                                            int min_n_ctx,                // minimum number of context bases in xm field
//                                            double min_ctx_meth_frac,     // minimum fraction of methylated to total context bases
//                                            double max_ooctx_meth_frac    // maximum fraction of methylated to total out-of-context bases
//                                            ) {
//   std::vector<bool> res (xm.size(), true);
//   for (int x=0; x<xm.size(); x++) {
//     int n_ctx_meth = rcpp_count_any_of(xm[x], ctx_meth);                 // count meth context bases
//     int n_ctx_unmeth = rcpp_count_any_of(xm[x], ctx_unmeth);             // count unmeth context bases
//     int n_ctx_all = n_ctx_meth+n_ctx_unmeth;
//     if (n_ctx_all<min_n_ctx) {                                           // if too few context bases -> FALSE
//       res[x] = false;
//     } else if (min_ctx_meth_frac>0 &&                                    // if we filter by VEF
//                n_ctx_meth/(double)n_ctx_all<min_ctx_meth_frac) {         // AND VEF is too low -> FALSE
//       res[x] = false;
//     } else if (max_ooctx_meth_frac<1) {                                  // if we filter by out-of-context methylation
//       int n_ooctx_meth = rcpp_count_any_of(xm[x], ooctx_meth);           // count meth out-of-context bases
//       int n_ooctx_unmeth = rcpp_count_any_of(xm[x], ooctx_unmeth);       // count unmeth out-of-context bases
//       int n_ooctx_all = n_ooctx_meth+n_ooctx_unmeth;
//       if (n_ooctx_all>0 &&
//           n_ooctx_meth/(double)n_ooctx_all>max_ooctx_meth_frac) {        // if out-of-context methylation is too high -> FALSE
//         res[x] = false;
//       }
//     }
//   }
//   return res;
// }
