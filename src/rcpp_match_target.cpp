#include <Rcpp.h>
// using namespace Rcpp;

// Matches reads to targets by start *or* end plus/minus tolerance (amplicons)
// OR overlap (capture).
// Return value: 1-based target index or NA for non-matched.
// Only first match is taken.
//
// Currently this code doesn't take advantage from the fact that both reads and
// targets can be sorted. Thus it is possible to speed it up quite a bit, taking
// the VCF code as an example

// MATCH AMPLICON BY POSITION
// fast, vectorised
// [[Rcpp::export("rcpp_match_amplicon")]]
std::vector<int> rcpp_match_amplicon(std::vector<std::string> read_chr,  // chr of reads
                                     std::vector<int> read_start,        // start pos of reads
                                     std::vector<int> read_end,          // end pos of reads
                                     std::vector<std::string> ampl_chr,  // chr of amplicons
                                     std::vector<int> ampl_start,        // start pos of amplicons
                                     std::vector<int> ampl_end,          // end pos of amplicons
                                     int tolerance)                      // coordinate tolerance
{
  std::vector<int> res (read_start.size(), NA_INTEGER);
  for (unsigned int x=0; x<read_start.size(); x++) {
    // checking for the interrupt
    if ((x & 1048575) == 0) Rcpp::checkUserInterrupt();
    
    for (unsigned int i=0; i<ampl_start.size(); i++) {
      if ((std::abs(read_start[x]-ampl_start[i]) <= tolerance ||
           std::abs(read_end[x]-ampl_end[i]) <= tolerance) &&
           read_chr[x] == ampl_chr[i]) {
        res[x] = i+1;
        break;
      }
    }
  }
  
  return res;
}


// MATCH CAPTURE BY OVERLAP
// fast, vectorised
// [[Rcpp::export("rcpp_match_capture")]]
std::vector<int> rcpp_match_capture(std::vector<std::string> read_chr,  // chr of reads
                                    std::vector<int> read_start,        // start pos of reads
                                    std::vector<int> read_end,          // end pos of reads
                                    std::vector<std::string> capt_chr,  // chr of capture targets
                                    std::vector<int> capt_start,        // start pos of capture targets
                                    std::vector<int> capt_end,          // end pos of capture targets
                                    signed int min_overlap)             // min overlap of reads and capture targets
{
  std::vector<int> res (read_start.size(), NA_INTEGER);
  for (unsigned int x=0; x<read_start.size(); x++) {
    // checking for the interrupt
    if (x & 1048575 == 0) Rcpp::checkUserInterrupt();
    
    for (unsigned int i=0; i<capt_start.size(); i++) {
      signed int overlap = std::min(read_end[x], capt_end[i]) - std::max(read_start[x], capt_start[i]) + 1;
      if (overlap >= min_overlap && read_chr[x] == capt_chr[i]) {
        res[x] = i+1;
        break;
      }
    }
  }
  
  return res;
}



// test code in R
//

/*** R
rcpp_match_amplicon(c("chr17", "chr17", "chr17", "chr17", "chr17", "chr17", "chr16", "chr17", "chr17", "chr17" ),
                    c(43125171,43125624,43125172,43125270,43124861,43125270,43125171,43125172,43124861,43125624),
                    c(43125550,43126026,43125551,43125640,43125249,43125640,43125550,43125550,43125249,43126026),
                    c("chr17", "chr17", "chr17", "chr17", "chr17" ),
                    c(43125624,43125270,43125171,43124861,43125270),
                    c(43126026,43125640,43125550,43125249,43125550), 0)
rcpp_match_capture(c("chr1",  "chr2",   "chr3",    "chr1",  "chr1",  "chr1",  "chr1",  "chr1",  "chr1",  "chr2" ),
                   c(3067647, 47401863, 100707761, 3067600, 3067600, 3067600, 3067500, 3067703, 3067705, 3067647),
                   c(3067703, 47402069, 100708296, 3067650, 3067647, 3067645, 3067803, 3067803, 3067803, 3067703),
                   c("chr1",  "chr2",   "chr3"),
                   c(3067647, 47401863, 100707761),
                   c(3067703, 47402069, 100708296), 1)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_match_target.cpp")

// #############################################################################
