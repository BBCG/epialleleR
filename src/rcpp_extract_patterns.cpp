#include <Rcpp.h>
#include <boost/container/flat_map.hpp>
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
                                      double min_ctx_freq,                      // minimum frequency of observed context at position
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
#define ctx_to_idx(c) ((c+2)>>2) & 15
  
  // consts, vars, typedefs
  const uint64_t offset_basis=14695981039346656037u;                            // FNV-1a offset basis
  unsigned int npat = 0;                                                        // pattern counter
  typedef int T_key;                                                            // map key
  typedef std::vector<int> T_val;                                               // map value
  boost::container::flat_map<T_key, T_key> pos_map;                             // all positions to figure out valid ones
  boost::container::flat_map<T_key, T_key>::iterator pos_hint = pos_map.begin();// position map iterator
  std::map<T_key, T_val> pat;                                                   // per-pattern methylation counts
  std::vector<int> pat_strand, pat_start, pat_end, pat_nbase;                   // pattern strands, starts, ends, number of bases within context
  std::vector<double> pat_beta;                                                 // pattern betas
  std::vector<std::string> pat_fnv;                                             // FNV-1a hashes of patterns
  
  pos_map.reserve(std::max((int)(target_end-target_start), 0xFFFF));
  
  // filling the context map
  unsigned int ctx_map [128] = {0};
  std::for_each(ctx.begin(), ctx.end(), [&ctx_map] (unsigned int const &c) {
    ctx_map[c] = 1;
  });
  
  // // rnames are sorted, should've been taking an advantage of it...
  // Rcpp::IntegerVector::const_iterator lower = std::lower_bound(rname.begin(), rname.end(), target_rname);
  // Rcpp::IntegerVector::const_iterator upper = std::upper_bound(lower, rname.end(), target_rname);
  // for (unsigned int x=lower-rname.begin(); x<upper-rname.begin(); x++) {
  
  // first - find valid positions
  for (unsigned int x=0; x<rname.size(); x++) {   
    // checking for the interrupt
    if ((x & 0xFFFF) == 0) Rcpp::checkUserInterrupt();                          // check for interrupt
    
    if (rname[x]==(int)target_rname) {
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
        for (unsigned int i=begin_i; i<end_i; i++) {                            // char by char - it's faster this way than using std::string in the cycle
          if (ctx_map[(int)xm_x[i]]) {                                          // if base is within context
            const unsigned int pos = start_x + i - offset_x;                    // position of the base
            pos_hint = pos_map.try_emplace(pos_hint, pos, 0);                   // check if this position is already included, emplace if not
            pos_hint->second++;                                                 // position++
          }
        }
        npat++;                                                                 // patterns++, to know how many
      }
    }
  }
  
  for (auto it=pos_map.begin(); it!=pos_map.end(); it++) {
    if ((double)it->second/npat >= min_ctx_freq)                                // if position frequency in patterns is higher than the min
      pat.emplace(it->first, T_val (npat, NA_INTEGER));                         // add position to the pattern map
  }
  
  npat = 0;
  for (unsigned int x=0; x<rname.size(); x++) {   
    // checking for the interrupt
    if ((x & 0xFFFF) == 0) Rcpp::checkUserInterrupt();                          // check for interrupt
    
    if (rname[x]==(int)target_rname) {
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
        unsigned int meth = 0, total = 0;                                       // counters for methylated and total within context
        uint64_t fnv = offset_basis;                                            // FNV-1a hash of current pattern
        for (unsigned int i=begin_i; i<end_i; i++) {                            // char by char - it's faster this way than using std::string in the cycle
          if (ctx_map[(int)xm_x[i]]) {                                          // if base is within context
            const unsigned int pos = start_x + i - offset_x;                    // position of the base
            auto hint = pat.find(pos);                                          // find and check if this position is already included
            if (hint != pat.end()) {
              const unsigned int base = ctx_to_idx(xm_x[i]);                    // rcpp_cx_report for details
              hint->second[npat] = base;                                        // save base by position
              meth += !(base & 8);                                              // methylated + (0 for lowercase, 1 for uppercase)
              total++;                                                          // total++
              fnv_add(fnv, reinterpret_cast<const char*>(&pos), sizeof(pos));   // FNV-1a: add int position
              fnv_add(fnv, xm_x+i, sizeof(char));                               // FNV-1a: add char base
            }
          }
        }
        
        if (fnv != offset_basis) {                                              // only if nonempty pattern
          npat++;                                                               // patterns++
          pat_strand.push_back(strand[x]);                                      // push strand
          pat_start.push_back(start_x+begin_i);                                 // push start
          pat_end.push_back(start_x+end_i-1);                                   // push end
          pat_nbase.push_back(total);                                           // push total
          pat_beta.push_back((double)meth/total);                               // push beta
          char fnv_str[17] = {0};
          snprintf(fnv_str, 17, "%.16lX", fnv);
          pat_fnv.emplace_back(fnv_str, 16);                                    // push FNV-1a hash
        }
      }
    }
  }
  
  if (!npat) return Rcpp::DataFrame::create();                                  // return empty DataFrame if no patterns were found
  
  Rcpp::IntegerVector pat_rname_w(npat, (int)target_rname);                     // rname factor (size, value)
  pat_rname_w.attr("class") = rname.attr("class");
  pat_rname_w.attr("levels") = rname.attr("levels");
  
  Rcpp::IntegerVector pat_strand_w = Rcpp::wrap(pat_strand);                    // strand factor
  pat_strand_w.attr("class") = strand.attr("class");
  pat_strand_w.attr("levels") = strand.attr("levels");
  
  for (auto it=pat.begin(); it!=pat.end(); it++) {                              // we reserved more, cutting off NA values now
    it->second.resize(npat);
  }
  Rcpp::DataFrame res = Rcpp::wrap(pat);                                        // wrap pattern map into DataFrame
  
  Rcpp::CharacterVector contexts = Rcpp::CharacterVector::create(               // base contexts
    NA_STRING,"H",NA_STRING,NA_STRING,NA_STRING,"X","Z",NA_STRING,
    NA_STRING,"h",NA_STRING,NA_STRING,NA_STRING,"x","z",NA_STRING
  );
  for (int i=0; i<res.length(); i++) {                                          // context factor for every column
    Rcpp::IntegerVector column = res[i];
    column.attr("class") = "factor";
    column.attr("levels") = contexts;
  }
  
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
bam <- preprocessBam(bam.file=system.file("extdata", "amplicon010meth.bam", package="epialleleR"))
z <- data.table::data.table(rcpp_extract_patterns(bam, 47, 43124861, 43125249, 1, "zZ", 0.1, TRUE, 0))
z[, c(lapply(.SD, unique), .N), by=pattern, .SDcols=grep("^X", colnames(z), value=TRUE)][order(-N)]
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_cx_report.cpp")

// #############################################################################
