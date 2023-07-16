#include <Rcpp.h>
#include <array>
#include <boost/container/flat_map.hpp>
#include "epialleleR.h"

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(BH)]]

// WTF report
// PRE-SORTED DATASET IS A REQUIREMENT.
// 
// Parses XM tags and calculates WTF. Outputs only if most frequent context
// is observed in more than 50% of reads (not including +-) and is within ctx
// string parameter.
// Output report is a data.frame with six columns and rows for every cytosine:
// rname (factor), strand (factor), pos, ctx, cov, wtf
// 
// 1) all XM positions counted in int[16]: index is equal to char+2>>2&00001111
// 2) when gap in reads or another chr - spit map to res, clear map
// 3) spit if within context and same context in more than 50% of the reads
// 
// ctx_to_idx conversion is described in epialleleR.h file
// 

// WTF numerator and denominator lookup tables are precomputed using Sk(n, k)
//
// Triangular sequence, nth element
uint64_t T(uint64_t n) {
  return n*(n+1)/2;
}
// WTF numerator or denominator, sum S of all possible WTF combinations
// of length i from 1 to n, times i
uint64_t S(uint64_t n)
{
  if (n<2) return n;
  return S(n-1) + T(n);
}
// WTF numerator or denominator, sum S of all possible WTF combinations
// of length i from 1 to k (where k<n), times i
// [[Rcpp::export]]
uint64_t Sk(uint64_t n, uint64_t k) {
  if (n<=k) return S(n);
  return S(k) + (n-k)*T(k);
}


// [[Rcpp::export]]
Rcpp::DataFrame rcpp_wtf_report(Rcpp::DataFrame &df,                            // data frame with BAM data
                                std::string ctx,                                // context string for bases to report,
                                bool discont,                                   // if WTF calculation should be discontinuous
                                int k)                                          // limit for l in WTF formula
{
  // walking trough bunch of reads <- filling the map
  // pos<<2|strand -> {0: rname,  1: pos,       2: 'H',  3: numer,  4: denom, 5: 'U',  6: 'X',  7: 'Z',
  //                   8: strand, 9: coverage, 10: 'h', 11: '+-',  12: '.',  13: 'u', 14: 'x', 15: 'z'}
  // boost::container::flat_map<uint64_t, std::array<uint64_t,16>>
  
  Rcpp::IntegerVector rname   = df["rname"];                                    // template rname
  Rcpp::IntegerVector strand  = df["strand"];                                   // template strand
  Rcpp::IntegerVector start   = df["start"];                                    // template start
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  Rcpp::XPtr<std::vector<std::string>> xm((SEXP)df.attr("xm_xptr"));            // merged refspaced template XMs, as a pointer to std::vector<std::string>
  
  // main typedefs
  typedef uint64_t T_key;                                                       // {62bit:pos, 2bit:strand}
  typedef std::array<uint64_t, 16> T_val;                                       // {0:rname, 1:pos, 8:strand, 9:coverage, 3:numerator, 4:denominator, and 10 more for 11 valid chars}
  typedef boost::container::flat_map<T_key, T_val> T_wtf_map;                   // attaboy
  
// macros
#define spit_results {                            /* save aggregated counts */ \
  for (T_wtf_map::iterator it=wtf_map.begin(); it!=wtf_map.end(); it++) {      \
    it->second[9] /= 2;                               /* halve the coverage */ \
    if (it->second[12] > it->second[9]) continue;     /* skip if most are . */ \
    else if ((it->second[2] + it->second[10]) > it->second[9])                 \
      max_freq_idx=2;                                                  /* H */ \
    else if ((it->second[6] + it->second[14]) > it->second[9])                 \
      max_freq_idx=6;                                                  /* X */ \
    else if ((it->second[7] + it->second[15]) > it->second[9])                 \
      max_freq_idx=7;                                                  /* Z */ \
    else continue;                                 /* skip if none is > 50% */ \
    if (ctx_map[max_freq_idx]) {                           /* if within ctx */ \
      res_strand.push_back(it->second[8]);                        /* strand */ \
      res_pos.push_back(it->second[1]);                              /* pos */ \
      res_ctx.push_back(max_freq_idx);                           /* context */ \
      res_cov.push_back(it->second[max_freq_idx] +                  /* meth */ \
                        it->second[max_freq_idx | 8]);            /* unmeth */ \
      res_wtf.push_back((double)it->second[3]/it->second[4]);        /* WTF */ \
    }                                                                          \
  }                                                                            \
  res_rname.resize(res_strand.size(), map_val[0]);           /* same rname! */ \
  max_pos=0;                                                                   \
  wtf_map.clear();                                                             \
  hint = wtf_map.end();                                                        \
};

  // array of contexts to print
  unsigned int ctx_map [16] = {0};
  std::for_each(ctx.begin(), ctx.end(), [&ctx_map] (unsigned int const &c) {
    ctx_map[ctx_to_idx(c)]=1;
  });
  
  // precomputed WTF numerator lookup table
  const size_t wtf_lookup_len = 4096;
  uint64_t wtf_lookup[wtf_lookup_len] = {0};
  if (k<=0) k=wtf_lookup_len;
  for (size_t n=0; n<wtf_lookup_len; n++) {
    wtf_lookup[n] = Sk(n, k);                                                   // filling the WTF values for faster computations
  }
  
  // WTF numerator buffer for current XM
  size_t num_buf_len = 8192;                                                    // maximum, though expandable length of numerator buffer
  uint64_t *num_buf  = (uint64_t*) malloc(num_buf_len * sizeof(uint64_t));      // numerator buffer

  // result
  std::vector<int> res_rname, res_strand, res_pos, res_ctx, res_cov;
  std::vector<double> res_wtf;
  size_t nitems = std::min(rname.size()*pow(ctx.size()<<2,2), 3e+9);
  res_rname.reserve(nitems); res_strand.reserve(nitems);
  res_pos.reserve(nitems); res_ctx.reserve(nitems);
  res_cov.reserve(nitems); res_wtf.reserve(nitems);
  
  // iterating over XM vector, saving the results when necessary
  T_wtf_map wtf_map;
  T_wtf_map::iterator hint;
  T_val map_val = {0};
  int max_pos = 0;
  unsigned int max_freq_idx;
  
  wtf_map.reserve(100000);                                                      // reserving helps?
  for (unsigned int x=0; x<rname.size(); x++) {
    // checking for the interrupt
    if ((x & 0xFFFF) == 0) Rcpp::checkUserInterrupt();                          // every ~65k reads
    
    const int start_x = start[x];                                               // start of the current read
    if ((start_x>max_pos) || ((uint64_t)rname[x]!=map_val[0])) {                // if current position is further downstream or another reference
      spit_results;
      map_val[0] = rname[x];
    }
    map_val[8] = strand[x];
    const char* xm_x = xm->at(templid[x]).c_str();                              // xm->at(templid[x]) is a reference to a corresponding XM string
    const unsigned int size_x = xm->at(templid[x]).size();                      // length of the current read
    
    // first, prefill WTF numerator buffer in first pass of XM
    if (num_buf_len < size_x) {
      num_buf_len = size_x;                                                     // new size
      num_buf  = (uint64_t*) realloc(num_buf, num_buf_len * sizeof(uint64_t));  // expand numerator buffer
      if (num_buf==NULL) Rcpp::stop("Unable to allocate memory");               // check memory allocation
    }
    std::memset(num_buf,  0, size_x * sizeof(uint64_t));                        // clean the buffer
    size_t mh_start = 0, mh_end = 0, mh_size = 0, h_size = 0;                   // start, end and size of the current methylated stretch (number of ctx bases); total size of haplotype
    uint64_t mh_sum = 0;                                                        // sum of WTF numerators for all methylated stretches, used to calculate continuous WTF that spans entire read pair
    for (unsigned int i=0; i<size_x; i++) {                                     // first pass to compute local WTF values, char by char
      const unsigned int base_idx = ctx_to_idx(xm_x[i]);                        // index of current base context; see the table in epialleleR.h
      if (ctx_map[base_idx]) {                                                  // if within context
        h_size++;                                                               // haplotype size++
        if (base_idx<8) {                                                       // if uppercase (methylated stretch started/continues)
          if (!mh_size) mh_start = i;                                           // store start position of methylated stretch
          mh_end = i;                                                           // store end position of methylated stretch
          mh_size++;                                                            // methylated stretch size++
        } else if (mh_size) {                                                   // if lowercase and after non-0-length methylated stretch
          std::fill(num_buf+mh_start, num_buf+mh_end+1, wtf_lookup[mh_size]);   // set values to Sk(mh_size, k) within methylated stretch
          mh_sum += wtf_lookup[mh_size];                                        // sum of numerators
          mh_size = 0;                                                          // reset the size
        }
      }
    }
    if (mh_size) {                                                              // save last non-0-length methylated stretch
      std::fill(num_buf+mh_start, num_buf+mh_end+1, wtf_lookup[mh_size]);
      mh_sum += wtf_lookup[mh_size];                                            // sum of numerators
    }
    if (!discont) std::fill_n(num_buf, size_x, mh_sum);                         // fill entire buffer if continuous WTF (same WTF value for entire read pair)
    
    // second, walk through XM once again, filling the map
    for (unsigned int i=0; i<size_x; i++) {                                     // char by char - it's faster this way than using std::string in the cycle
      const unsigned int idx_to_increase = ctx_to_idx(xm_x[i]);                 // index of context; see the table in epialleleR.h
      if (idx_to_increase==11) continue;                                        // skip +-
      map_val[1] = start_x+i;                                                   // current position
      const T_key map_key = ((T_key)map_val[1] << 2) | map_val[8];
      hint = wtf_map.try_emplace(hint, map_key, map_val);
      hint->second[idx_to_increase]++;
      hint->second[9]++;                                                        // total coverage
      hint->second[3] += num_buf[i];                                            // WTF numerator
      hint->second[4] += wtf_lookup[h_size];                                    // WTF denominator
    }
    if ((uint64_t)max_pos<map_val[1]) max_pos=map_val[1];                       // last position of C in wtf_map
  }
  spit_results;
  
  Rcpp::DataFrame res = Rcpp::DataFrame::create(                                // final CX report
    Rcpp::Named("rname") = res_rname,                                           // numeric ids (factor) for reference names
    Rcpp::Named("strand") = res_strand,                                         // numeric ids (factor) for reference strands
    Rcpp::Named("pos") = res_pos,                                               // position of cytosine
    Rcpp::Named("context") = res_ctx,                                           // cytosine context
    Rcpp::Named("coverage") = res_cov,                                          // cytosine context
    Rcpp::Named("wtf") = res_wtf                                                // number of unmethylated
  );
  
  Rcpp::IntegerVector col_rname = res["rname"];                                 // making rname a factor
  col_rname.attr("class") = "factor";
  col_rname.attr("levels") = rname.attr("levels");
  
  Rcpp::IntegerVector col_strand = res["strand"];;                              // making strand a factor
  col_strand.attr("class") = "factor";
  col_strand.attr("levels") = strand.attr("levels");
  
  Rcpp::CharacterVector contexts = Rcpp::CharacterVector::create(               // base contexts
    "NA1","CHH","NA3","NA4","NA5","CHG","CG"
  );
  Rcpp::IntegerVector col_context = res["context"];;                            // making context a factor
  col_context.attr("class") = "factor";
  col_context.attr("levels") = contexts;
  
  free(num_buf);                                                                // free manually allocated memory
  
  return res;
}


// test code in R
//

/*** R
### WTF calculations for a stretch of n mCpGs over a window of k
# for mCpG stretches of length n, sum S of all possible WTF combinations
# (of length i from 1 to n, times i) equals to:
# n   S   T
# 1   1   1
# 2   4   3
# 3  10   6
# 4  20  10
# 5  35  15
# 6  56  21
# 7  84  28
# 8 120  36
#
# which is in effect, Fibonacci-like sequence on top of triangular sequence
# T{n} = n(n+1)/2; S{n} = T{n} + S{n-1}
### 
T <- function (n) {n*(n+1)/2}
S <- function (n) {
  if (n<2) n
  else S(n-1) + T(n)
}
sapply(1:10, S)
###
# for mCpG stretches of length n, sum S of all possible WTF combinations
# (of length i from 1 to k [where k<n], times i) equals to:
#   k=2 k=3 k=4 k=5 k=6 k=7 k=8
# n   S   S   S   S   S   S   S
# 1   1   1   1   1   1   1   1
# 2   4   4   4   4   4   4   4
# 3   7  10  10  10  10  10  10
# 4  10  16  20  20  20  20  20
# 5  13  22  30  35  35  35  35
# 6  16  28  40  50  56  56  56
# 7  19  34  50  65  77  84  84
# 8  22  40  60  80  98 112 120
###
Sk <- function (n, k) {
  if (n<=k) S(n)
  else S(k) + (n-k)*T(k)
}
matrix(sapply(1:10, function (k) lapply(1:10, Sk, k=k)),
       nrow=10, dimnames=list(n=1:10, k=1:10))
###
#
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_wtf_report.cpp")
