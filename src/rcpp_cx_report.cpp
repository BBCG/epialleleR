#include <Rcpp.h>
#include <array>
#include <boost/container/flat_map.hpp>
// using namespace Rcpp;

// Was using C++17 for std::map::try_emplace
// Most likely not required for boost::container::flat_map
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
// Output report is an int array with six fields for every cytosine:
// rname (factor), pos, strand (factor), ctx (char), meth, unmeth
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
std::vector<int> rcpp_cx_report(std::vector<int> rname,            // int value for factorised BAM rname field
                                std::vector<int> strand,           // int value for factorised BAM strand field
                                std::vector<int> start,            // read start = min(pos,mpos) of BAM fields
                                std::vector<std::string> xm,       // merged normalised BAM XM fields
                                std::vector<bool> pass,            // boolean: if read has passed filtering
                                std::string ctx)                   // context string for bases to report
{
  // walking trough bunch of reads <- filling the map
  // pos<<2|strand -> {0: rname,  1: pos,       2: 'H',  3: '',    4: '',   5: 'U',  6: 'X',  7: 'Z',
  //                   8: strand, 9: coverage, 10: 'h', 11: '+-', 12: '.', 13: 'u', 14: 'x', 15: 'z'}
  // boost::container::flat_map<uint64_t, std::array<int,16>>
  
  // main typedefs
  typedef uint64_t T_key;                                          // {62bit:pos, 2bit:strand}
  typedef std::array<int,16> T_val;                                // {0:rname, 1:pos, 8:strand, 9:coverage, and 10 more for 11 valid chars}
  typedef boost::container::flat_map<T_key, T_val> T_cx_fmap;      // attaboy
  
  // macros
  #define ctx_to_idx(c) ((c+2)>>2) & 15
  #define spit_results {                                                       \
    for (it=cx_map.begin(); it!=cx_map.end(); it++) {                          \
      it->second[9] -= it->second[11];                 /* discarding the +- */ \
      if (it->second[12]*2 > it->second[9])           /* skip if most are . */ \
        continue;                                                              \
      else if ((it->second[2] + it->second[10])*2 > it->second[9])             \
        max_freq_ctx='H', max_freq_idx=2;                                      \
      else if ((it->second[6] + it->second[14])*2 > it->second[9])             \
        max_freq_ctx='X', max_freq_idx=6;                                      \
      else if ((it->second[7] + it->second[15])*2 > it->second[9])             \
        max_freq_ctx='Z', max_freq_idx=7;                                      \
      else continue;                               /* skip if none is > 50% */ \
      if (ctx.find(max_freq_ctx)!=std::string::npos) {     /* if within ctx */ \
        res.push_back(it->second[0]);                              /* rname */ \
        res.push_back(it->second[1]);                                /* pos */ \
        res.push_back(it->second[8]);                             /* strand */ \
        res.push_back(max_freq_ctx);                             /* context */ \
        res.push_back(it->second[max_freq_idx]);                    /* meth */ \
        max_freq_idx |= 8;                                     /* lowercase */ \
        res.push_back(it->second[max_freq_idx]);                  /* unmeth */ \
      }                                                                        \
    }                                                                          \
    cx_map.clear();                                                            \
    hint = cx_map.end();                                                       \
  }

  // result
  // { (rname, pos, strand, ctx, meth, unmeth) * n }
  std::vector<int> res;
  // reserving space makes it only slow:
  // res.reserve(rname.size()*pow(ctx.size(),2));

  // iterating over XM vector, saving the results when necessary
  T_cx_fmap cx_map;
  T_cx_fmap::iterator it, hint;
  T_key map_key;
  T_val map_val = {0};
  unsigned int idx_to_increase, max_freq_ctx, max_freq_idx, pass_x;
  
  cx_map.reserve(100000);                                          // reserving helps
  for (unsigned int x=0; x<rname.size(); x++) {
    // checking for the interrupt
    if ((x & 1048575) == 0) Rcpp::checkUserInterrupt();
    
    if (start[x]>map_val[1] || rname[x]!=map_val[0]) {
      spit_results;
    }
    map_val[0] = rname[x];
    map_val[8] = strand[x];
    pass_x = pass[x]?0:8;                                          // should we lowercase this XM
    for (unsigned int i=0; i<xm[x].size(); i++) {
      map_val[1] = start[x]+i;
      map_key = ((T_key)map_val[1] << 2) | map_val[8];
      idx_to_increase = ctx_to_idx(xm[x][i]);                      // see the table above
      idx_to_increase |= pass_x;                                   // if not pass - lowercase
      hint = cx_map.try_emplace(hint, map_key, map_val);
      hint->second[idx_to_increase]++;
      hint->second[9]++;                                           // total coverage
    }
  }
  spit_results;
  
  return res;
}


// test code in R
//

/*** R
matrix(rcpp_cx_report(c(2,2), c(1,1), c(3689466,3689466),
                      c("z......xh....Z.......Z......z...Zx...x..xh.hhhh..xhh.Zxh......xhh.......xh..x..x..x..z....z......z......Z",
                        "z......xh....Z.......Z......z...Zx...x..xh.hhhh..xhh.Zxh......xhh.......xh..x..x..x..z....z......z......Z"),
                      c(TRUE,FALSE), "zZ"),
       ncol=6, byrow=TRUE)

rname  <- as.integer(bam$rname)
strand <- as.integer(bam$strand)
start  <- bam$start
xm     <- bam$XM
pass   <- rep(TRUE, length(rname))
ctx    <- "ZXHzxh" # "Zz" # 
# microbenchmark::microbenchmark(rcpp_cx_report(rname, strand, start, xm, pass, ctx), times=10)
#
system.time(res <- data.table::data.table(matrix(rcpp_cx_report(rname, strand, start, xm, pass, ctx), ncol=6, byrow=TRUE,
                                                 dimnames=list(NULL, c("rname","pos","strand","context","meth","unmeth")))) )
dim(res)

*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_cx_report.cpp")

// #############################################################################
// 
// Optimised previous flat_map version, which 
//   DOES NOT CHECK CORRECTNESS OF THE CONTEXT and
//   IS POORLY SCALABLE BEYOND 2e+6 READS FROM WHOLE-GENOME APPROACHES
// is pasted below:
// 
// // [[Rcpp::export("rcpp_cx_report")]]
// std::vector<int> rcpp_cx_report(std::vector<int> rname,            // int value for factorised BAM rname field
//                                 std::vector<int> strand,           // int value for factorised BAM strand field
//                                 std::vector<int> start,            // read start = min(pos,mpos) of BAM fields
//                                 std::vector<std::string> xm,       // merged normalised BAM XM fields
//                                 std::vector<bool> pass,            // boolean: if read has passed filtering
//                                 std::string ctx)                   // context string for bases to report
// {
//   // walking trough all reads <- filling the map
//   // uint64_t -> {rname, pos, strand, context, meth, unmeth}
//   // boost::container::flat_map<uint64_t, std::array<int,6>>
//   
//   // main typedefs
//   typedef uint64_t T_key;                                          // {16bit:rname, 32bit:pos, 8bit:strand, 8bit:context
//   typedef std::array<int,6> T_val;                                 // {0:rname, 1:pos, 2:strand, 3:context, 4:meth, 5:unmeth}
//   typedef boost::container::flat_map<T_key, T_val> T_cx_fmap;      // attaboy
//   // typedef std::map<T_key, T_val> T_cx_map;                      // 2x slower than boost::container::flat_map
//   // typedef std::unordered_map<T_key, T_val> T_cx_umap;           // not used since slower than ordered std::map even w/o sorting
//   
//   // iterating over XM vector
//   // for 1'618'360 reads, all pass, "zxhZXH": cycling without try_emplace takes 1.12s
//   // compared to 3.25s with try_emplace but w/o it->second[idx_to_increase]++
//   // it->second[idx_to_increase]++ takes no time, most time is taken by try_emplace
//   // hint is correct most of the time, I believe
//   T_cx_fmap cx_map;
//   T_cx_fmap::iterator it = cx_map.end();
//   T_key map_key = 0;
//   T_val map_val = {0, 0, 0, 0, 0, 0};
//   unsigned int idx_to_increase = 0;
//   cx_map.reserve(rname.size()*pow(ctx.size(),2));                  // exclusive to flat_map, up to 20% faster
//   for (unsigned int x=0; x<rname.size(); x++) {
//     map_val[0] = rname[x];
//     map_val[2] = strand[x];
//     int found = xm[x].find_first_of(ctx, 0);
//     while (found != std::string::npos) {
//       map_val[1] = start[x] + found;
//       idx_to_increase = (pass[x] && xm[x][found]>='A' && xm[x][found]<='Z') ? 4 : 5;
//       map_val[3] = xm[x][found] | 32;                              // lowercasing reported context
//       map_key = ((T_key)map_val[0] << 48) |
//                 ((T_key)map_val[1] << 16) |
//                 ((T_key)map_val[2] << 8 ) |
//                 ((T_key)map_val[3]);
//       it = cx_map.try_emplace(it, map_key, map_val);
//       it->second[idx_to_increase]++;
//       found = xm[x].find_first_of(ctx, found+1);
//     }
//   }
//   
//   // walking through the map, returning the int vector of
//   // { (rname, pos, strand, ctx, meth, unmeth) * n }
//   // there must be a better and more efficient way - iter+insert takes ~5-10% time
//   std::vector<int> res;
//   res.reserve(cx_map.size()*6);
//   for (it=cx_map.begin(); it!=cx_map.end(); it++) {
//     res.insert(res.end(), it->second.begin(), it->second.end());
//   }
//   return res;
//   
// }
//