#include <Rcpp.h>
#include <boost/container/flat_map.hpp>
using namespace Rcpp;

// Using C++17 for std::map::try_emplace
// Not sure if required for boost::container::flat_map
// [[Rcpp::plugins(cpp17)]]

// Using boost::container::flat_map here as it is >3x faster than any of std::*
// [[Rcpp::depends(BH)]]
// sparse_map - https://www.codeproject.com/Articles/866996/Fast-Implementations-of-Maps-with-Integer-Keys-in
// can actually be the fastest map for the task, however I didn't try it yet
// because it hequires me to include non-standard header and change the way I
// emplace new items. Besides, looks like genome-assisted counting is 
// inevitable, and it might be faster for very large (whole genome) data.
// For the smaller (any MiSeq file in any context) flat_map is acceptable.

// Parses XM tags and outputs SUMMARISED CX report according to filter/context.
// Output report is an int array with six fields for every cytosine:
// rname (factor), pos, strand (factor), ctx (char), meth, unmeth


// CX report, vectorised, SUMMARISING by rname+pos+strand+context
// [[Rcpp::export("rcpp_cx_report")]]
std::vector<int> rcpp_cx_report(std::vector<int> rname,            // int value for factorised BAM rname field
                                std::vector<int> strand,           // int value for factorised BAM strand field
                                std::vector<int> start,            // read start = min(pos,mpos) of BAM fields
                                std::vector<std::string> xm,       // merged normalised BAM XM fields
                                std::vector<bool> pass,            // boolean: if read has passed filtering
                                std::string ctx)                   // context string for bases to report
{
  // walking trough all reads <- filling the map
  // uint64_t -> {rname, pos, strand, context, meth, unmeth}
  // boost::container::flat_map<uint64_t, std::array<int,6>>
  
  // main typedefs
  typedef uint64_t T_key;                                          // {16bit:rname, 32bit:pos, 8bit:strand, 8bit:context
  typedef std::array<int,6> T_val;                                 // {0:rname, 1:pos, 2:strand, 3:context, 4:meth, 5:unmeth}
  typedef boost::container::flat_map<T_key, T_val> T_cx_fmap;      // attaboy
  // typedef std::map<T_key, T_val> T_cx_map;                      // 2x slower than boost::container::flat_map
  // typedef std::unordered_map<T_key, T_val> T_cx_umap;           // not used since slower than ordered std::map even w/o sorting
  
  // iterating over XM vector
  // for 1'618'360 reads, all pass, "zxhZXH": cycling without try_emplace takes 1.12s
  // compared to 3.25s with try_emplace but w/o it->second[idx_to_increase]++
  // it->second[idx_to_increase]++ takes no time, most time is taken by try_emplace
  // hint is correct most of the time, I believe
  T_cx_fmap cx_map;
  T_cx_fmap::iterator it = cx_map.end();
  T_key map_key = 0;
  T_val map_val = {0, 0, 0, 0, 0, 0};
  unsigned int idx_to_increase = 0;
  cx_map.reserve(rname.size()*pow(ctx.size(),2));                  // exclusive to flat_map, up to 20% faster
  for (unsigned int x=0; x<rname.size(); x++) {
    map_val[0] = rname[x];
    map_val[2] = strand[x];
    int found = xm[x].find_first_of(ctx, 0);
    while (found != std::string::npos) {
      map_val[1] = start[x] + found;
      idx_to_increase = (pass[x] && xm[x][found]>='A' && xm[x][found]<='Z') ? 4 : 5;
      map_val[3] = xm[x][found] | 32;                              // lowercasing reported context
      map_key = ((T_key)map_val[0] << 48) |
                ((T_key)map_val[1] << 16) |
                ((T_key)map_val[2] << 8 ) |
                ((T_key)map_val[3]);
      it = cx_map.try_emplace(it, map_key, map_val);
      it->second[idx_to_increase]++;
      found = xm[x].find_first_of(ctx, found+1);
    }
  }
  
  // walking through the map, returning the int vector of
  // { (rname, pos, strand, ctx, meth, unmeth) * n }
  // there must be a better and more efficient way - iter+insert takes ~5-10% time
  std::vector<int> res;
  res.reserve(cx_map.size()*6);
  for (it=cx_map.begin(); it!=cx_map.end(); it++) {
    res.insert(res.end(), it->second.begin(), it->second.end());
  }
  return res;
  
}


// TODO: there will be another function
// 1) all XM positions counted in int[16]: index is equal to char+2>>2&00001111
// 2) when gap in reads or another chr - spit map to res, clear map
// 3) find max - spit if within context
// Here's the array
// char  bin       +2        >>2&15  dec0b
// +     00101011  00101101  1011    11
// -     00101101  00101111  1011    11
// .     00101110  00110000  1100    12
// H     01001000  01001010  0010    2
// U     01010101  01010111  0101    5
// X     01011000  01011010  0110    6
// Z     01011010  01011100  0111    7
// h     01101000  01101010  1010    10
// u     01110101  01110111  1101    13
// x     01111000  01111010  1110    14
// z     01111010  01111100  1111    15
// 


// test code in R
//

/*** R
matrix(rcpp_cx_report(c(2,2), c(1,1), c(3689466,3689466),
                      c("z......xh....Z.......Z......z...Zx...x..xh.hhhh..xhh.Zxh......xhh.......xh..x..x..x..z....z......z......Z",
                        "z......xh....Z.......Z......z...Zx...x..xh.hhhh..xhh.Zxh......xhh.......xh..x..x..x..z....z......z......Z"),
                      c(TRUE,FALSE), "zZ"),
       ncol=6, byrow=TRUE)

.bam   <- bam.processed[order(bam.processed$rname, bam.processed$start),]
rname  <- as.integer(.bam$rname)
strand <- as.integer(.bam$strand)
start  <- .bam$start
xm     <- .bam$XM
pass   <- rep(TRUE, length(rname))
ctx    <- "zxhZXH"
# microbenchmark::microbenchmark(rcpp_cx_report(rname, strand, start, xm, pass, ctx), times=10)
# 
res <- matrix(rcpp_cx_report(rname, strand, start, xm, pass, ctx), ncol=6, byrow=TRUE)
dim(res)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_cx_report.cpp")

// #############################################################################

// For some reason on 1618360 reads, 2011500 positions to report in zZ context
// unordered_map was slower than std::map even without sorting
// 
// UNORDERED_MAP code snippets:
// 
// struct key_hasher {
//   std::size_t operator()(const T_key & k) const {
//     std::size_t h = (k[1]<<8) | k[0];
//     return(h);
//   }
// };
// typedef std::unordered_map<T_key, T_val, key_hasher> T_cx_umap;
// 
// 
// // sorting the unsorted_map
// std::vector<T_key> keys;
// keys.reserve(cx_map.size());
// for (auto& it : cx_map) {
//   keys.push_back(it.first);
// }
// std::sort(keys.begin(), keys.end());
// 
// 
// #############################################################################
// 
// // old rcpp_parse_xm
// std::vector<int> res;
// for (int x=0; x<rname.size(); x++) {
//   std::vector<int> row = {rname[x], strand[x], 0, 0, 0, 0};       // rname, strand, pos, ctx, meth, unmeth
//   for (int i=0; i<ctx.size(); i++) {
//     row[3] = ctx[i] | 32;                                         // lowercasing reported context
//     bool is_meth = (pass[x] && ctx[i]>='A' && ctx[i]<='Z') ? true : false;
//     int found = xm[x].find(ctx[i], 0);
//     while (found != std::string::npos) {
//       row[2] = start[x] + found;
//       row[4] = is_meth ? 1 : 0;
//       row[5] = is_meth ? 0 : 1;
//       res.insert(res.end(), row.begin(), row.end());
//       found = xm[x].find(ctx[i], found+1);
//     }
//   }
// }
// return res;
// 
// #############################################################################
// 
// Perfectly optimised std version, but still 2x slower than flat_map:
// // [[Rcpp::export("rcpp_cx_report")]]
// std::vector<int> rcpp_cx_report(std::vector<int> rname,            // int value for factorised BAM rname field
//                                 std::vector<int> strand,           // int value for factorised BAM strand field
//                                 std::vector<int> start,            // read start = min(pos,mpos) of BAM fields
//                                 std::vector<std::string> xm,       // merged normalised BAM XM fields
//                                 std::vector<bool> pass,            // boolean: if read has passed filtering
//                                 std::string ctx)                   // context string for bases to report
// {
//   // walking trough all reads <- filling the std::map
//   // uint64_t -> {rname, pos, strand, context, meth, unmeth}
//   // std::map<uint64_t, std::array<int,6>>
//   
//   // main typedefs
//   typedef uint64_t T_key;                                          // {16bit:rname, 32bit:pos, 8bit:strand, 8bit:context
//   typedef std::array<int,6> T_val;                                 // {0:rname, 1:pos, 2:strand, 3:context, 4:meth, 5:unmeth}
//   typedef std::map<T_key, T_val> T_cx_map;
//   // typedef std::unordered_map<T_key, T_val> T_cx_umap;           // not used since slower even w/o sorting
//   
//   // iterating over XM vector
//   // for 1'618'360 reads, all pass, "zxhZXH": cycling without try_emplace takes 1.10s
//   // compared to 7.80s with try_emplace but w/o it->second[idx_to_increase]++
//   // it->second[idx_to_increase]++ takes no time, most time is taken by try_emplace
//   // even if hint is correct most of the time (I believe)
//   T_cx_map cx_map;
//   T_cx_map::iterator it = cx_map.end();
//   T_key map_key = 0;
//   T_val map_val = {0, 0, 0, 0, 0, 0};
//   unsigned int idx_to_increase = 0;
//   for (unsigned int x=0; x<rname.size(); x++) {
//     map_val[0] = rname[x];
//     map_val[2] = strand[x];
//     for (unsigned int i=0; i<ctx.size(); i++) {
//       map_val[3] = ctx[i] | 32;                                     // lowercasing reported context
//       idx_to_increase = (pass[x] && ctx[i]>='A' && ctx[i]<='Z') ? 4 : 5;
//       int found = xm[x].find(ctx[i], 0);
//       while (found != std::string::npos) {
//         map_val[1] = start[x] + found;
//         map_key = ((T_key)map_val[0] << 48) |
//           ((T_key)map_val[1] << 16) |
//           ((T_key)map_val[2] << 8 ) |
//           ((T_key)map_val[3]);
//         it = cx_map.try_emplace(it, map_key, map_val);
//         it->second[idx_to_increase]++;
//         it++;
//         found = xm[x].find(ctx[i], found+1);
//       }
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