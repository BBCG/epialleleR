#include <Rcpp.h>
#include <array>
#include <boost/container/flat_map.hpp>
#include "epialleleR.h"

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(BH)]]

// Linearized MHL report
// PRE-SORTED DATASET IS A REQUIREMENT.
// 
// Parses XM tags and calculates linearized MHL (lMHL). Outputs only if the most
// frequent context is observed in more than 50% of reads (not including +-)
// and is within ctx string parameter.
// Output report is a data.frame with seven columns and rows for every cytosine:
// rname (factor), strand (factor), pos, ctx (factor), cov, hlen, mhl
// 
// 1) all XM positions counted in int[16]: index is equal to char+2>>2&00001111
// 2) when gap in reads or another chr - spit map to res, clear map
// 3) spit if within context and same context in more than 50% of the reads
// 
// ctx_to_idx conversion is described in epialleleR.h file
// 

// lMHL numerator and denominator lookup tables are precomputed using nrS(n)
//
// Triangular sequence, nth element
// uint64_t T(uint64_t n) {
//   return (n*(n+1))/2;
// }
// lMHL numerator or denominator, sum S of all possible lMHL combinations
// of length i from 1 to n, times i
// uint64_t S(uint64_t n)
// {
//   if (n<2) return n;
//   return S(n-1) + T(n);
// }
// S can be simplified to non-recursive versions
uint64_t nrS(uint64_t n)
{
  if (n<2) return n;
  return (n*(n+1)*(n+2))/6;
}

// [[Rcpp::export]]
Rcpp::DataFrame rcpp_mhl_report(Rcpp::DataFrame &df,                            // data frame with BAM data
                                const std::string ctx,                          // context string for bases to report,
                                int hmax,                                       // maximum length of a computation window (limit for l in lMHL formula)
                                const int hmin,                                 // ignore haplotypes smaller than hmin
                                const double max_ooctx_meth_frac)               // maximum fraction of methylated to total out-of-context bases (max out-of-context beta value)
{
  // walking trough bunch of reads <- filling the map
  // pos -> { 0: rname,   1: pos,       2: 'H',  3: numer,  4: denom,  5: 'U',  6: 'X',  7: 'Z',  # + strand
  //          8: h_size,  9: coverage, 10: 'h', 11: '+-',  12: '.',   13: 'u', 14: 'x', 15: 'z',  # + strand
  //         16: ''    , 17: '',       18: 'H', 19: numer, 20: denom, 21: 'U', 22: 'X', 23: 'Z',  # - strand
  //         24: h_size, 25: coverage, 26: 'h', 27: '+-',  28: '.',   29: 'u', 30: 'x', 31: 'z'}  # - strand
  // boost::container::flat_map<uint64_t, std::array<uint64_t,32>>
  
  Rcpp::IntegerVector rname   = df["rname"];                                    // template rname
  Rcpp::IntegerVector strand  = df["strand"];                                   // template strand
  Rcpp::IntegerVector start   = df["start"];                                    // template start
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  Rcpp::XPtr<std::vector<std::string>> seqxm((SEXP)df.attr("seqxm_xptr"));      // merged refspaced packed template SEQXMs, as a pointer to std::vector<std::string>
  
  // main typedefs
  typedef uint64_t T_key;                                                       // {64bit:pos}
  typedef std::array<uint64_t, 32> T_val;                                       // {0:rname, 1:pos, 8,25:h_size, 9,25:coverage, 3,19:numerator, 4,20:denominator, and 10 more for 11 valid chars * two strands}
  typedef boost::container::flat_map<T_key, T_val> T_mhl_map;                   // attaboy
  
// macros
#define spit_results {                                                                           /* save aggregated counts */ \
  for (T_mhl_map::iterator it=mhl_map.begin(); it!=mhl_map.end(); it++) {                                                     \
    for (int s=0; s<2; s++) {                                                                      /* iterate over strands */ \
      str_shft = s<<4;                                                               /* strand shift: 0 for F and 16 for R */ \
      if (it->second[9+str_shft]==0) continue;                                                      /* skip if not covered */ \
      it->second[9+str_shft] /= 2;                                                                   /* halve the coverage */ \
      if (it->second[12+str_shft] > it->second[9+str_shft]) continue;                                /* skip if most are . */ \
      else if ((it->second[2+str_shft] + it->second[10+str_shft]) > it->second[9+str_shft])                                   \
        max_freq_idx=2;                                                                                               /* H */ \
      else if ((it->second[6+str_shft] + it->second[14+str_shft]) > it->second[9+str_shft])                                   \
        max_freq_idx=6;                                                                                               /* X */ \
      else if ((it->second[7+str_shft] + it->second[15+str_shft]) > it->second[9+str_shft])                                   \
        max_freq_idx=7;                                                                                               /* Z */ \
      else continue;                                                                              /* skip if none is > 50% */ \
      if (ctx_map[max_freq_idx]) {                                                                        /* if within ctx */ \
        res_strand.push_back(s+1);                                                                               /* strand */ \
        res_pos.push_back(it->first);                                                                               /* pos */ \
        res_ctx.push_back(max_freq_idx);                                                                        /* context */ \
        const int cov = it->second[max_freq_idx+str_shft] + it->second[(max_freq_idx+str_shft) | 8];                          \
        res_cov.push_back(cov);                                                                           /* meth + unmeth */ \
        res_hlen.push_back((double)it->second[8+str_shft]/cov);                                        /* average hap size */ \
        res_mhl.push_back((double)it->second[3+str_shft]/it->second[4+str_shft]);                                  /* lMHL */ \
      }                                                                                                                       \
    }                                                                                                                         \
  }                                                                                                                           \
  res_rname.resize(res_strand.size(), map_val[0]);                                                          /* same rname! */ \
  max_pos=0;                                                                                                                  \
  mhl_map.clear();                                                                                                            \
  hint = mhl_map.end();                                                                                                       \
};

  // array of contexts to print
  unsigned int ctx_map [16] = {0};
  std::for_each(ctx.begin(), ctx.end(), [&ctx_map] (unsigned int const &c) {
    ctx_map[ctx_to_idx(c)]=1;
  });
  
  // precomputed lMHL numerator lookup table
  const size_t mhl_lookup_len = 65536;
  uint64_t mhl_lookup[mhl_lookup_len] = {0};
  hmax = (hmax>0) ? std::min((size_t)hmax, mhl_lookup_len) : mhl_lookup_len ;   // number of context bases is always in range [1; mhl_lookup_len]
  for (size_t n=0; n<(size_t)hmax; n++) {
    mhl_lookup[n] = nrS(n);                                                     // filling the lMHL values for faster computations
  }
  std::fill_n(mhl_lookup+hmax, mhl_lookup_len-hmax, nrS(hmax));                 // if hmax < mhl_lookup_len - fill the rest of the lookup table with it
  
  // lMHL numerator buffer for current XM
  size_t num_buf_len = 8192;                                                    // maximum, though expandable length of numerator buffer
  uint64_t *num_buf  = (uint64_t*) malloc(num_buf_len * sizeof(uint64_t));      // numerator buffer

  // result
  std::vector<int> res_rname, res_strand, res_pos, res_ctx, res_cov;
  std::vector<double> res_hlen, res_mhl;
  size_t nitems = std::min(rname.size()*pow(ctx.size()<<2,2), 3e+9);
  res_rname.reserve(nitems); res_strand.reserve(nitems);
  res_pos.reserve(nitems); res_ctx.reserve(nitems);
  res_cov.reserve(nitems); res_hlen.reserve(nitems); res_mhl.reserve(nitems);
  
  // iterating over XM vector, saving the results when necessary
  T_mhl_map mhl_map;
  T_mhl_map::iterator hint;
  T_val map_val = {0};
  int max_pos = 0;
  unsigned int max_freq_idx, str_shft;
  
  mhl_map.reserve(100000);                                                      // reserving helps?
  for (unsigned int x=0; x<rname.size(); x++) {
    // checking for the interrupt
    if ((x & 0xFFFF) == 0) Rcpp::checkUserInterrupt();                          // every ~65k reads
    
    const int start_x = start[x];                                               // start of the current read
    if ((start_x>max_pos) || ((uint64_t)rname[x]!=map_val[0])) {                // if current position is further downstream or another reference
      spit_results;
      map_val[0] = rname[x];
    }
    str_shft = (strand[x]-1)<<4;                                                // strand shift: 0 for F and 16 for R
    const char* seqxm_x = seqxm->at(templid[x]).c_str();                        // seqxm->at(templid[x]) is a reference to a corresponding SEQXM string
    const unsigned int size_x = seqxm->at(templid[x]).size();                   // length of the current read
    
    // first, prefill lMHL numerator buffer in first pass of XM
    if (num_buf_len < size_x) {
      num_buf_len = size_x;                                                     // new size
      num_buf  = (uint64_t*) realloc(num_buf, num_buf_len * sizeof(uint64_t));  // expand numerator buffer
      if (!num_buf) Rcpp::stop("Unable to allocate memory");                    // check memory allocation
    }
    std::memset(num_buf,  0, size_x * sizeof(uint64_t));                        // clean the buffer
    size_t mh_start = 0, mh_end = 0, mh_size = 0, h_size = 0;                   // start, end and size of the current methylated stretch (number of ctx bases); total size of haplotype
    size_t ooctx_map [16] = {0};                                                // counter array for methylated and unmethylated out-of-context bases
    for (unsigned int i=0; i<size_x; i++) {                                     // first pass to compute local lMHL values, char by char
      const unsigned int base_idx = unpack_ctx_idx(seqxm_x[i]);                 // index of current base context; see the table in epialleleR.h
      if (ctx_map[base_idx]) {                                                  // if within context
        h_size++;                                                               // haplotype size++
        if (base_idx<8) {                                                       // if uppercase (methylated stretch started/continues)
          if (!mh_size) mh_start = i;                                           // store start position of methylated stretch
          mh_end = i;                                                           // store end position of methylated stretch
          mh_size++;                                                            // methylated stretch size++
        } else if (mh_size) {                                                   // if lowercase and after non-0-length methylated stretch
          std::fill(num_buf+mh_start, num_buf+mh_end+1, mhl_lookup[mh_size]);   // set values to nrS(mh_size) within methylated stretch
          mh_size = 0;                                                          // reset the size
        }
      } else {                                                                  // only if out-of-context
        ooctx_map[base_idx]++;                                                  // increment the counter of that symbol
      }
    }
    size_t ooctx_meth = ooctx_map[2]+ooctx_map[5]+ooctx_map[6]+ooctx_map[7];    // sum of o-o-ctx methylated
    size_t ooctx_unmeth = ooctx_map[10]+ooctx_map[13]+ooctx_map[14]+ooctx_map[15];  // sum of o-o-ctx unmethylated
    double ooctx_meth_frac = (double)ooctx_meth / (ooctx_meth+ooctx_unmeth);    // fraction of o-o-ctx methylated
    if ((int)h_size<hmin || ooctx_meth_frac>max_ooctx_meth_frac) continue;      // skip read if haplotype is smaller than hmin or too many o-o-ctx meth bases
    if (mh_size) {                                                              // save last non-0-length methylated stretch
      std::fill(num_buf+mh_start, num_buf+mh_end+1, mhl_lookup[mh_size]);
    }
    
    // second, walk through XM once again, filling the map
    for (unsigned int i=0; i<size_x; i++) {                                     // char by char - it's faster this way than using std::string in the cycle
      const unsigned int idx_to_increase = unpack_ctx_idx(seqxm_x[i]);          // index of context; see the table in epialleleR.h
      if (idx_to_increase==11) continue;                                        // skip +-
      map_val[1] = start_x+i;                                                   // current position
      hint = mhl_map.try_emplace(hint, (T_key)(map_val[1]), map_val);
      hint->second[idx_to_increase+str_shft]++;
      hint->second[9+str_shft]++;                                               // total coverage
      hint->second[8+str_shft] += h_size;                                       // sum haplotype sizes
      hint->second[3+str_shft] += num_buf[i];                                   // lMHL numerator
      hint->second[4+str_shft] += mhl_lookup[h_size];                           // lMHL denominator
    }
    if ((uint64_t)max_pos<map_val[1]) max_pos=map_val[1];                       // last position of C in mhl_map
  }
  spit_results;
  
  Rcpp::DataFrame res = Rcpp::DataFrame::create(                                // final CX report
    Rcpp::Named("rname") = res_rname,                                           // numeric ids (factor) for reference names
    Rcpp::Named("strand") = res_strand,                                         // numeric ids (factor) for reference strands
    Rcpp::Named("pos") = res_pos,                                               // position of cytosine
    Rcpp::Named("context") = res_ctx,                                           // cytosine context
    Rcpp::Named("coverage") = res_cov,                                          // cytosine coverage
    Rcpp::Named("length") = res_hlen,                                           // average haplotype length
    Rcpp::Named("lmhl") = res_mhl                                               // lMHL value
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
### lMHL calculations for a stretch of n mCpGs
# for mCpG stretches of length n, sum S of all possible lMHL combinations
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
nrS <- function (n) {
  if (n<2) n
  else n*(n+1)*(n+2)/6
}
sapply(1:30, S)
sapply(1:30, nrS)
microbenchmark::microbenchmark(sapply(1:200, S), sapply(1:200, nrS))

### lMHL calculations for a stretch of n mCpGs over a window of k.
### The following is not used, as it doesn't provide proper granularity for
### long-range sequencing
#
# for mCpG stretches of length n, sum S of all possible lMHL combinations
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
nrSk <- function (n, k) {
  if (n<=k) nrS(n)
  else nrS(k) + (n-k)*T(k)
}
matrix(sapply(1:10, function (k) lapply(1:10, Sk, k=k)),
       nrow=10, dimnames=list(n=1:10, k=1:10))
matrix(sapply(1:10, function (k) lapply(1:10, nrSk, k=k)),
       nrow=10, dimnames=list(n=1:10, k=1:10))
###
#
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_mhl_report.cpp")
