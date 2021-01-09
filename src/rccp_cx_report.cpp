#include <Rcpp.h>
using namespace Rcpp;

// Using C++17 for std::map::try_emplace
// [[Rcpp::plugins(cpp17)]]

// Parses XM tags and outputs SUMMARISED CX report according to filter/context.
// Output report is an int array with six fields for every cytosine:
// rname (factor), pos, strand (factor), ctx (char), meth (0|1), unmeth (0|1)

// CX report, vectorised, SUMMARISING by rname+pos+strand+context
// [[Rcpp::export("rcpp_cx_report")]]
std::vector<int> rcpp_cx_report(std::vector<int> rname,       // int value for factorised BAM rname field
                                std::vector<int> strand,      // int value for factorised BAM strand field
                                std::vector<int> start,       // read start = min(pos,mpos) of BAM fields
                                std::vector<std::string> xm,  // merged normalised BAM XM fields
                                std::vector<bool> pass,       // boolean: if read has passed filtering
                                std::string ctx)              // context string for bases to report
{
  // walking trough all reads <- filling the std::map
  // {rname, pos, strand, context} -> {meth, unmeth}
  // std::map<std::array<int,4>, std::array<int,2>> OR
  // std::unordered_map<std::array<int,4>, std::array<int,2>>
  // with custom hashing: pos<<8 | chr
  
  // main typedefs
  typedef std::array<unsigned int,4> T_key;                   // {0:rname, 1:pos, 2:strand, 3:context
  typedef std::array<unsigned int,2> T_val;                   // {0:meth, 1:unmeth}
  typedef std::map<T_key, T_val> T_cx_map;
  struct key_hasher {
    std::size_t operator()(const T_key & k) const {
      std::size_t h = (k[1]<<8) | k[0];
      return(h);
    }
  };
  typedef std::unordered_map<T_key, T_val, key_hasher> T_cx_umap;

  // iterating over XM vector
  T_cx_map cx_map;
  T_key map_key = {0, 0, 0, 0};
  T_val map_val = {0, 0};
  for (unsigned int x=0; x<rname.size(); x++) {
    map_key[0] = rname[x];
    map_key[2] = strand[x];
    for (unsigned int i=0; i<ctx.size(); i++) {
      map_key[3] = ctx[i] | 32;                               // lowercasing reported context
      unsigned char idx_to_increase = (pass[x] && ctx[i]>='A' && ctx[i]<='Z') ? 0 : 1;
      int found = xm[x].find(ctx[i], 0);
      while (found != std::string::npos) {
        map_key[1] = start[x] + found;
        auto [it, exists] = cx_map.try_emplace(map_key, map_val);
        it->second[idx_to_increase]++;
        found = xm[x].find(ctx[i], found+1);
      }
    }
  }

  // walking through the map, returning the int vector of
  // { (rname, pos, strand, ctx, meth, unmeth) * n }
  // there must be a better and more efficient way - iter+insert takes ~10% time
  T_cx_map::iterator it;
  std::vector<int> res;
  res.reserve(cx_map.size()*6);
  for (it=cx_map.begin(); it!=cx_map.end(); it++) {
    res.insert(res.end(), it->first.begin(), it->first.end());
    res.insert(res.end(), it->second.begin(), it->second.end());
  }
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

n <- 1000
rname  <- as.integer(.bam.proc$rname) # rep(c(1:20), 5*n) # 20*5
strand <- as.integer(.bam.proc$strand) # rep(c(1,2), 50*n) # 2*50
start  <- .bam.proc$start # rep(c(1:(n*20)), 5) + 1000000
xm     <- .bam.proc$XM # rep("z......xh....Z.......Z......z...Zx...x..xh.hhhh..xhh.Zxh......xhh.......xh..x..x..x..z....z......z......Z", 100*n)
pass   <- .bam.proc$pass # rep(c(TRUE,FALSE), 50*n)
ctx    <- "zxhZXH"
# microbenchmark::microbenchmark(rcpp_cx_report(rname, strand, start, xm, pass, ctx), times=10)
# 
res <- matrix(rcpp_cx_report(rname, strand, start, xm, pass, ctx), ncol=6, byrow=TRUE)
dim(res)
# 
# simple rcpp_parse_xm mean time: 32-35 ms for n=1000
# ordered std::map takes 132 ms for n=1000 (130 ms if memory was reserved for res )
# std::unordered_map takes 71 ms for n=1000 (with memory reserved)
# 
# real data - 1006 capture (1'618'360 reads, 2'011'500 positions to report in zZ context):
# 10 times, mo mapping (parse_xm), mean time: 0.83s
# 10 times, std::unordered_map, mean time: 2.80s
# 10 times, ordered std::map, mean time: 2.41s
# report by "generateCytosineReport" with rcpp_parse_xm+summarise takes 6-8s
# report by "generateCytosineReport" with rcpp_cx_report+summarise takes 10-12s
# report by "generateCytosineReport" with rcpp_cx_report takes 2.9-3.0s
# 
# real data - 1006 capture (1'618'360 reads, 18'765'434 positions to report in zxhZXH context):
# 10 times, std::unordered_map, mean time: 20.1s
# 10 times, ordered std::map, mean time: 12.6s
# report by "generateCytosineReport" with rcpp_parse_xm+summarise takes 66-72s"
# report by "generateCytosineReport" with rcpp_parse_xm+summarise takes 90-97s"
# report by "generateCytosineReport" with rcpp_cx_report takes 23-25s with unordered_map and 15-17s with std::map
# 
# real data - WHIP050-D501-D701 amplicon (167'241 reads, 27'903 positions to report in zZ context):
# 10 times, std::unordered_map, mean time: 137ms
# 10 times, ordered std::map, mean time: 200ms
# 
# real data - WHIP050-D501-D701 amplicon (167'241 reads, 409'589 positions to report in zxhZXH context):
# 10 times, std::unordered_map, mean time: 1.12s
# 10 times, ordered std::map, mean time: 1.43s
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
