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
std::vector<int> rcpp_match_amplicon(Rcpp::DataFrame &df,                       // BAM data
                                     Rcpp::DataFrame &bed,                      // BED data
                                     int tolerance)                             // coordinate tolerance
{
  Rcpp::IntegerVector read_chr = df["rname"];                                   // template rname
  Rcpp::IntegerVector read_start = df["start"];                                 // template start
  // Rcpp::CharacterVector xm = bam["XM"];                                         // template XM
  Rcpp::XPtr<std::vector<std::string>> xm((SEXP)df.attr("xm_xptr"));            // merged refspaced template XMs, as a pointer to std::vector<std::string>
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  Rcpp::IntegerVector ampl_chr = bed["seqnames"];                               // BED rname
  Rcpp::IntegerVector ampl_start = bed["start"];                                // BED start
  Rcpp::IntegerVector ampl_end = bed["end"];                                    // BED start
  
  std::vector<int> res (read_start.size(), NA_INTEGER);
  for (unsigned int x=0; x<read_start.size(); x++) {
    // checking for the interrupt
    if ((x & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();
    
    int read_end = read_start[x] + xm->at(templid[x]).size() - 1;
    for (unsigned int i=0; i<ampl_start.size(); i++) {
      if ((read_chr[x] == ampl_chr[i]) &&
          ((std::abs(read_start[x] - ampl_start[i]) <= tolerance) ||
           (std::abs(read_end - ampl_end[i]) <= tolerance))) {
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
std::vector<int> rcpp_match_capture(Rcpp::DataFrame &df,                        // BAM data
                                    Rcpp::DataFrame &bed,                       // BED data
                                    signed int min_overlap)                     // min overlap of reads and capture targets
{
  Rcpp::IntegerVector read_chr = df["rname"];                                   // template rname
  Rcpp::IntegerVector read_start = df["start"];                                 // template start
  // Rcpp::CharacterVector xm = bam["XM"];                                         // template XM
  Rcpp::XPtr<std::vector<std::string>> xm((SEXP)df.attr("xm_xptr"));            // merged refspaced template XMs, as a pointer to std::vector<std::string>
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  Rcpp::IntegerVector capt_chr = bed["seqnames"];                               // BED rname
  Rcpp::IntegerVector capt_start = bed["start"];                                // BED start
  Rcpp::IntegerVector capt_end = bed["end"];                                    // BED start
  
  std::vector<int> res (read_start.size(), NA_INTEGER);
  for (unsigned int x=0; x<read_start.size(); x++) {
    // checking for the interrupt
    if ((x & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();
    
    int read_end = read_start[x] + xm->at(templid[x]).size() - 1;
    for (unsigned int i=0; i<capt_start.size(); i++) {
      signed int overlap = std::min(read_end, capt_end[i]) - std::max(read_start[x], capt_start[i]) + 1;
      if ((read_chr[x] == capt_chr[i]) && (overlap >= min_overlap)) {
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
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_match_target.cpp")

// #############################################################################
