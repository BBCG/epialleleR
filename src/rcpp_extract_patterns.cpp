#include <Rcpp.h>
// using namespace Rcpp;

// Scans trough reads and extracts methylation patterns from the reads that
// overlap target area. Clips overhangs if necessary.
// Return value: a data frame with the following columns:
// 1) pattern id (FNV hash)
// 1) 

// fast, vectorised
// [[Rcpp::export("rcpp_extract_patterns")]]
Rcpp::DataFrame rcpp_extract_patterns(Rcpp::DataFrame &df,                      // data frame with BAM data
                                      int target_rname,                         // 
                                      int target_start,                         // 
                                      int target_end,                           // 
                                      signed int min_overlap,                   // min overlap of reads and capture targets
                                      std::string &ctx,                         // context string for bases to include
                                      bool clip,                                // clip the matched reads to target area
                                      unsigned int reverse_offset) {            // decrease reverse strand coordinates by this value: 0 for CHH, 1 for CpG, 2 for CHG
  Rcpp::IntegerVector rname   = df["rname"];                                    // template rname
  Rcpp::IntegerVector strand  = df["strand"];                                   // template strand
  Rcpp::IntegerVector start   = df["start"];                                    // template start
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  Rcpp::XPtr<std::vector<std::string>> xm((SEXP)df.attr("xm_xptr"));            // merged refspaced template XMs, as a pointer to std::vector<std::string>
  
  std::map<int, std::vector<int>> pat;                                          // per-pattern counts
  std::vector<uint64_t> patid;                                                  // FNV hashes of patterns
  
  // // rnames are sorted, should've been taking an advantage of it...
  // Rcpp::IntegerVector::const_iterator lower = std::lower_bound(rname.begin(), rname.end(), target_rname);
  // Rcpp::IntegerVector::const_iterator upper = std::upper_bound(lower, rname.end(), target_rname);
  // for (unsigned int x=lower-rname.begin(); x<upper-rname.begin(); x++) {
  
  for (unsigned int x=0; x<rname.size(); x++) {   
    // checking for the interrupt
    if ((x & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();
    
    if (rname[x]==target_rname) {
      const unsigned int size_x = xm->at(templid[x]).size();                    // length of the current read
      const unsigned int start_x = start[x];                                    // start position of the current read
      const unsigned int end_x = start_x + size_x - 1;                          // end position of the current read
      const signed int overlap = std::min(end_x, target_end) - std::max(start_x, target_start) + 1; // overlap with target
      if (overlap>=min_overlap) {                                               // if overlaps the target
        const char* xm_x = xm->at(templid[x]).c_str();                          // xm->at(templid[x]) is a reference to a corresponding XM string
        const unsigned int offset_x = strand[x]==2 ? reverse_offset : 0;        // offset coordinates of reverse strand for symmetric methylation
        const unsigned int begin_i = clip ? (start_x - std::max(start_x, target_start)) : 0;
        const unsigned int end_i = clip ? overlap : size_x;                     // ^ clip the XM?
        for (unsigned int i=begin_i; i<end_i; i++) {                            // char by char - it's faster this way than using std::string in the cycle
          if (ctx.find(xm_x[i]) != std::string::npos) {
            const unsigned int pos = start_x + i - offset_x;
          }
        }
      }
    }

  }
  
  
  std::vector<int> x1(4,6) ;
  std::vector<int> x2(4,2) ;
  std::vector<int> x3(4,3) ;
  
  std::map< int, std::vector<int> > m ;
  m.emplace(11, x1);
  m.emplace(44, x2);
  m.emplace(33, x3);
  
  return Rcpp::wrap( m ) ;
}










// test code in R
//

/*** R
rcpp_extract_patterns(bam, 47, 43124861, 43125261, 1, "zZ", TRUE, 1)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_cx_report.cpp")

// #############################################################################
