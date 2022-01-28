#include <Rcpp.h>
#include <array>
#include <boost/container/flat_map.hpp>

// Was using C++17 for std::map::try_emplace
// Very much required for boost::container::flat_map - 10x speed gain
// [[Rcpp::plugins(cpp17)]]
// 
// Using boost::container::flat_map as it is >3x faster than any of std::*
// [[Rcpp::depends(BH)]]
// sparse_map - https://www.codeproject.com/Articles/866996/Fast-Implementations-of-Maps-with-Integer-Keys-in
// can actually be the fastest map for the task, however I didn't try it yet
// because it requires me to include non-standard header and change the way I
// emplace new items.


// CX report, vectorised, summarising, context-aware, linearly scalable
// PRE-SORTED DATASET IS A REQUIREMENT.
// 
// Parses XM tags and outputs summarised CX report only if most frequent context
// is observed in more than 50% of reads (not including +-) and is within ctx
// string parameter.
// Output report is a data.frame with six columns and rows for every cytosine:
// rname (factor), strand (factor), pos, ctx (char), meth, unmeth
// 
// 1) all XM positions counted in int[16]: index is equal to char+2>>2&00001111
// 2) when gap in reads or another chr - spit map to res, clear map
// 3) spit if within context and same context in more than 50% of the reads
// 
// Here's the ctx_to_idx convertion:
// ctx  bin       +2        >>2&15  idx
// +    00101011  00101101  1011    11
// -    00101101  00101111  1011    11
// .    00101110  00110000  1100    12
// H    01001000  01001010  0010    2
// U    01010101  01010111  0101    5
// X    01011000  01011010  0110    6
// Z    01011010  01011100  0111    7
// h    01101000  01101010  1010    10
// u    01110101  01110111  1101    13
// x    01111000  01111010  1110    14
// z    01111010  01111100  1111    15
// 
// [[Rcpp::export("rcpp_cx_report")]]
Rcpp::DataFrame rcpp_cx_report(Rcpp::DataFrame &df,                             // data frame with BAM data
                               Rcpp::LogicalVector &pass,                       // does it pass the threshold
                               std::string ctx)                                 // context string for bases to report
{
  // walking trough bunch of reads <- filling the map
  // pos<<2|strand -> {0: rname,  1: pos,       2: 'H',  3: '',    4: '',   5: 'U',  6: 'X',  7: 'Z',
  //                   8: strand, 9: coverage, 10: 'h', 11: '+-', 12: '.', 13: 'u', 14: 'x', 15: 'z'}
  // boost::container::flat_map<uint64_t, std::array<int,16>>
  
  Rcpp::IntegerVector rname   = df["rname"];                                    // template rname
  Rcpp::IntegerVector strand  = df["strand"];                                   // template strand
  Rcpp::IntegerVector start   = df["start"];                                    // template start
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  Rcpp::XPtr<std::vector<std::string>> xm((SEXP)df.attr("xm_xptr"));            // merged refspaced template XMs, as a pointer to std::vector<std::string>
  
  // main typedefs
  typedef uint64_t T_key;                                                       // {62bit:pos, 2bit:strand}
  typedef std::array<int,16> T_val;                                             // {0:rname, 1:pos, 8:strand, 9:coverage, and 10 more for 11 valid chars}
  typedef boost::container::flat_map<T_key, T_val> T_cx_fmap;                   // attaboy
  
  // macros
  #define ctx_to_idx(c) ((c+2)>>2) & 15
  #define spit_results {                                                       \
    for (it=cx_map.begin(); it!=cx_map.end(); it++) {                          \
      it->second[9] /= 2;                             /* halve the coverage */ \
      if (it->second[12] > it->second[9]) continue;   /* skip if most are . */ \
      else if ((it->second[2] + it->second[10]) > it->second[9])               \
        max_freq_idx=2;                                                /* H */ \
      else if ((it->second[6] + it->second[14]) > it->second[9])               \
        max_freq_idx=6;                                                /* X */ \
      else if ((it->second[7] + it->second[15]) > it->second[9])               \
        max_freq_idx=7;                                                /* Z */ \
      else continue;                               /* skip if none is > 50% */ \
      if (ctx_map[max_freq_idx]) {                         /* if within ctx */ \
        res_rname.push_back(it->second[0]);                        /* rname */ \
        res_strand.push_back(it->second[8]);                      /* strand */ \
        res_pos.push_back(it->second[1]);                            /* pos */ \
        res_ctx.push_back(max_freq_idx);                         /* context */ \
        res_meth.push_back(it->second[max_freq_idx]);               /* meth */ \
        res_unmeth.push_back(it->second[max_freq_idx | 8]);       /* unmeth */ \
      }                                                                        \
    }                                                                          \
    max_pos=0;                                                                 \
    cx_map.clear();                                                            \
    hint = cx_map.end();                                                       \
  }

  unsigned int ctx_map [16] = {0};                                              // array of contexts to print
  std::for_each(ctx.begin(), ctx.end(), [&ctx_map] (unsigned int const &c) {
    ctx_map[ctx_to_idx(c)]=1;
  });

  // result
  std::vector<int> res_rname, res_strand, res_pos, res_ctx, res_meth, res_unmeth;
  size_t nitems = std::min(rname.size()*pow(ctx.size()<<2,2), 1e+9);
  res_rname.reserve(nitems); res_strand.reserve(nitems);
  res_pos.reserve(nitems); res_ctx.reserve(nitems);
  res_meth.reserve(nitems); res_unmeth.reserve(nitems);
  
  // iterating over XM vector, saving the results when necessary
  T_cx_fmap cx_map;
  T_cx_fmap::iterator it, hint;
  T_key map_key;
  T_val map_val = {0};
  int max_pos;
  unsigned int idx_to_increase, max_freq_idx, pass_x;
  
  cx_map.reserve(100000);                                                       // reserving helps?
  for (unsigned int x=0; x<rname.size(); x++) {
    // checking for the interrupt
    if ((x & 0xFFFF) == 0) Rcpp::checkUserInterrupt();                          // every ~65k reads
    
    if ((start[x]>max_pos) || (rname[x]!=map_val[0])) {                         // if current position is further downstream or another reference
      spit_results;
      map_val[0] = rname[x];
    }
    map_val[8] = strand[x];
    pass_x = (!pass[x])<<3;                                                     // should we lowercase this XM (TRUE==0, FALSE==8)
    const char* xm_x = xm->at(templid[x]).c_str();                              // xm->at(templid[x]) is a reference to a corresponding XM string
    const unsigned int size_x = xm->at(templid[x]).size();
    for (unsigned int i=0; i<size_x; i++) {                                     // char by char - it's faster this way than using std::string in the cycle
      idx_to_increase = ctx_to_idx(xm_x[i]);                                    // see the table above
      if (idx_to_increase==11) continue;                                        // skip +-
      idx_to_increase |= pass_x;                                                // if not pass - lowercase
      map_val[1] = start[x]+i;
      map_key = ((T_key)map_val[1] << 2) | map_val[8];
      hint = cx_map.try_emplace(hint, map_key, map_val);
      hint->second[idx_to_increase]++;
      hint->second[9]++;                                                        // total coverage
    }
    max_pos=max_pos<map_val[1]?map_val[1]:max_pos;                              // last position of C in cx_map
  }
  spit_results;
  
  Rcpp::DataFrame res = Rcpp::DataFrame::create(                                // final CX report
    Rcpp::Named("rname") = res_rname,                                           // numeric ids (factor) for reference names
    Rcpp::Named("strand") = res_strand,                                         // numeric ids (factor) for reference strands
    Rcpp::Named("pos") = res_pos,                                               // position of cytosine
    Rcpp::Named("context") = res_ctx,                                           // cytosine context
    Rcpp::Named("meth") = res_meth,                                             // number of methylated
    Rcpp::Named("unmeth") = res_unmeth                                          // number of unmethylated
  );
  
  return res;
}


// test code in R
//

/*** R
# microbenchmark::microbenchmark(rcpp_cx_report(rname, strand, start, xm, pass, ctx), times=10)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_cx_report.cpp")

// #############################################################################
// ## branched spit
//  #define spit_results {                                                       \
//    for (it=cx_map.begin(); it!=cx_map.end(); it++) {                          \
//      it->second[9] /= 2;                              /* half the coverage */ \
//      if (it->second[12] > it->second[9]) continue;   /* skip if most are . */ \
//      max_freq_idx = ((it->second[2]+it->second[10]) > it->second[9])*2+ /*H*/ \
//                     ((it->second[6]+it->second[14]) > it->second[9])*6+ /*X*/ \
//                     ((it->second[7]+it->second[15]) > it->second[9])*7; /*Z*/ \
//      if (ctx_map[max_freq_idx]) {                         /* if within ctx */ \
//        res_rname.push_back(it->second[0]);                        /* rname */ \
//        res_strand.push_back(it->second[8]);                      /* strand */ \
//        res_pos.push_back(it->second[1]);                            /* pos */ \
//        res_ctx.push_back(max_freq_idx);                         /* context */ \
//        res_meth.push_back(it->second[max_freq_idx]);               /* meth */ \
//        res_unmeth.push_back(it->second[max_freq_idx | 8]);       /* unmeth */ \
//      }                                                                        \
//    }                                                                          \
//    max_pos=0;                                                                 \
//    cx_map.clear();                                                            \
//    hint = cx_map.end();                                                       \
//  }
