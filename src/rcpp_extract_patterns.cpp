#include <Rcpp.h>
// using namespace Rcpp;

// Scans trough reads and extracts methylation patterns from the reads that
// overlap target area. Clips overhangs if necessary.
// Return value: a data frame with the following columns:
// 1) pattern id (FNV hash)
// 1) 

// [[Rcpp::plugins(cpp17)]]

// fast, vectorised
// [[Rcpp::export("rcpp_extract_patterns")]]
Rcpp::DataFrame rcpp_extract_patterns(Rcpp::DataFrame &df,                      // data frame with BAM data
                                      unsigned int target_rname,                // target chromosome
                                      unsigned int target_start,                // target start
                                      unsigned int target_end,                  // target end
                                      signed int min_overlap,                   // min overlap of reads and capture targets
                                      std::string &ctx,                         // context string for bases to include
                                      bool clip,                                // clip the matched reads to target area
                                      unsigned int reverse_offset) {            // decrease reverse strand coordinates by this value: 0 for CHH, 1 for CpG, 2 for CHG
  Rcpp::IntegerVector rname   = df["rname"];                                    // template rname
  Rcpp::IntegerVector strand  = df["strand"];                                   // template strand
  Rcpp::IntegerVector start   = df["start"];                                    // template start
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  Rcpp::XPtr<std::vector<std::string>> xm((SEXP)df.attr("xm_xptr"));            // merged refspaced template XMs, as a pointer to std::vector<std::string>
  
  // http://www.isthe.com/chongo/tech/comp/fnv/
#define fnv_add(hash, pointer, size) {     /* hash starts with offset_basis */ \
  for (unsigned int idx=0; idx<size; idx++) {        /* cycle through bytes */ \
    hash ^= *(pointer+idx);                       /* hash xor octet_of_data */ \
    hash *= 1099511628211u;                             /* hash * FNV_prime */ \
  }                                                                            \
}
  
  unsigned int npat = 0;                                                        // pattern counter
  typedef int T_key;                                                            // map key
  typedef std::vector<int> T_val;                                               // map value
  std::map<T_key, T_val> pat;                                                   // per-pattern methylation counts
  std::map<T_key, T_val>::iterator hint = pat.begin();                          // map iterator
  std::vector<int> pat_strand, pat_start, pat_end, pat_nbase;                   // pattern strands, starts, ends, number of bases within context
  std::vector<double> pat_beta;                                                 // pattern betas
  std::vector<std::string> pat_fnv;                                             // FNV-1a hashes of patterns
  
  // // rnames are sorted, should've been taking an advantage of it...
  // Rcpp::IntegerVector::const_iterator lower = std::lower_bound(rname.begin(), rname.end(), target_rname);
  // Rcpp::IntegerVector::const_iterator upper = std::upper_bound(lower, rname.end(), target_rname);
  // for (unsigned int x=lower-rname.begin(); x<upper-rname.begin(); x++) {
  
  for (unsigned int x=0; x<rname.size(); x++) {   
    // checking for the interrupt
    if ((x & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();                         // check for interrupt
    
    if (rname[x]==target_rname) {
      const unsigned int size_x = xm->at(templid[x]).size();                    // length of the current read
      const unsigned int start_x = start[x];                                    // start position of the current read
      const unsigned int end_x = start_x + size_x - 1;                          // end position of the current read
      const unsigned int over_start_x = std::max(start_x, target_start);        // start of overlapped area
      const unsigned int over_end_x = std::min(end_x, target_end);              // end of overlapped area
      const signed int overlap = over_end_x - over_start_x + 1;                 // overlap with target
      if (overlap>=min_overlap) {                                               // if overlaps the target
        const char* xm_x = xm->at(templid[x]).c_str();                          // xm->at(templid[x]) is a reference to a corresponding XM string
        const unsigned int offset_x = strand[x]==2 ? reverse_offset : 0;        // offset coordinates of reverse strand for symmetric methylation
        const unsigned int begin_i = clip ? (over_start_x - start_x) : 0;       // clip the XM?
        const unsigned int end_i = clip ? overlap : size_x;                     // clip the XM?
        unsigned int meth=0, total=0;                                           // counters for methylated and total within context
        uint64_t fnv=14695981039346656037u;                                    // FNV-1a hash of current pattern
        for (unsigned int i=begin_i; i<end_i; i++) {                            // char by char - it's faster this way than using std::string in the cycle
          if (ctx.find(xm_x[i]) != std::string::npos) {                         // if base is within context
            const unsigned int pos = start_x + i - offset_x;                    // position of the base
            hint = pat.try_emplace(hint, pos, T_val (npat, NA_INTEGER));        // check if this position is already included, emplace if not
            const unsigned int base = !(xm_x[i] & 0x20);                        // 0 for lowercase, 1 for uppercase
            hint->second.push_back(base);                                       // save base by position
            meth += base;                                                       // methylated?++
            total++;                                                            // total++
            fnv_add(fnv, reinterpret_cast<const char*>(&pos), sizeof(unsigned int));
            fnv_add(fnv, xm_x+i, sizeof(char));
          }
        }
        npat++;                                                                 // patterns++
        for (auto it=pat.begin(); it!=pat.end(); it++) {
          it->second.resize(npat, NA_INTEGER);                                  // top up all vectors with NA values
        }
        pat_strand.push_back(strand[x]);
        pat_start.push_back(start_x+begin_i);
        pat_end.push_back(start_x+end_i-1);
        pat_nbase.push_back(total);
        pat_beta.push_back((double)meth/total);
        char fnv_str[17] = {0};
        snprintf(fnv_str, 17, "%.16lX", fnv);
        pat_fnv.emplace_back(fnv_str, 16);
      }
    }
  }
  
  Rcpp::IntegerVector pat_rname_w(npat, (int)target_rname);                     // rname factor (size, value)
  pat_rname_w.attr("class") = rname.attr("class");
  pat_rname_w.attr("levels") = rname.attr("levels");
  
  Rcpp::IntegerVector pat_strand_w = Rcpp::wrap(pat_strand);                    // strand factor
  pat_strand_w.attr("class") = strand.attr("class");
  pat_strand_w.attr("levels") = strand.attr("levels");
  
  Rcpp::DataFrame res = Rcpp::wrap(pat);
  res.push_front(pat_fnv, "pattern");
  res.push_front(pat_beta, "beta");
  res.push_front(pat_nbase, "nbase");
  res.push_front(pat_end, "end");
  res.push_front(pat_start, "start");
  res.push_front(pat_strand_w, "strand");
  res.push_front(pat_rname_w, "seqnames");
  
  return(res) ;
}



// test code in R
//

/*** R
z <- data.table::data.table(rcpp_extract_patterns(bam, 47, 43124861, 43125249, 1, "zZ", TRUE, 0))
z
z[, c(lapply(.SD, unique), .N), by=pattern, .SDcols=grep("^X", colnames(z), value=TRUE)][order(-N)]
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_cx_report.cpp")

// #############################################################################
