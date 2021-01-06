#include <Rcpp.h>
using namespace Rcpp;

// Parses XM tags and outputs CX report according to filter and context.
// Output report is an int array with six fields for every cytosine:
// rname (factor), strand (factor), pos, ctx (char), meth (0|1), unmeth (0|1)

// CX report, vectorised, the main one now
// [[Rcpp::export]]
std::vector<int> rcpp_parse_xm_vector(std::vector<int> rname,       // int value for factorised BAM rname field
                                      std::vector<int> strand,      // int value for factorised BAM strand field
                                      std::vector<int> start,       // read start = min(pos,mpos) of BAM fields
                                      std::vector<std::string> xm,  // merged normalised BAM XM fields
                                      std::vector<bool> pass,       // boolean: if read has passed filtering
                                      std::string ctx) {            // context string for bases to report
  std::vector<int> res;                                             // { (rname, strand, pos, ctx, meth, unmeth) * n }
  for (int x=0; x<rname.size(); x++) {
    std::vector<int> row = {rname[x], strand[x], 0, 0, 0, 0};       // rname, strand, pos, ctx, meth, unmeth
    for (int i=0; i<ctx.size(); i++) {
      row[3] = ctx[i];
      bool is_meth = (pass[x] && ctx[i]>='A' && ctx[i]<='Z') ? true : false;
      int found = xm[x].find(ctx[i], 0);
      while (found != std::string::npos) {
        row[2] = start[x] + found;
        row[4] = is_meth ? 1 : 0;
        row[5] = is_meth ? 0 : 1;
        res.insert(res.end(), row.begin(), row.end());
        found = xm[x].find(ctx[i], found+1);
      } 
    }
  }
  return res;
}


// test code in R
//

/*** R
matrix(rcpp_parse_xm_vector(c(2,2), c(1,1), c(3689466,3689466),
                            c("z......xh....Z.......Z......z...Zx...x..xh.hhhh..xhh.Zxh......xhh.......xh..x..x..x..............z......Z",
                              "z......xh....Z.......Z......z...Zx...x..xh.hhhh..xhh.Zxh......xhh.......xh..x..x..x..............z......Z"),
                            c(TRUE,FALSE), "zZ"),
       ncol=6, byrow=TRUE)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_parse_xm.cpp")

// #############################################################################


// // [[Rcpp::export]]
// char rcpp_base_to_lower(char base) {
//   return (base>='A' && base<='Z') ? base+32 : base;
// }
// 
// // [[Rcpp::export]]
// bool rcpp_base_is_meth(char base) {
//   return (base>='A' && base<='Z') ? true : false;
// }
// 
// sapply(c("u","U","h","H","x","X","z","Z",".","-","+"), rcpp_base_to_lower)
// sapply(c("u","U","h","H","x","X","z","Z",".","-","+"), rcpp_base_is_meth)


// // simple CX report
// // [[Rcpp::export]]
// std::vector<int> rcpp_parse_xm(int rname, int strand, int start, std::string xm, bool pass, std::string ctx) {
//   std::vector<int> res; // { (rname, strand, pos, ctx, meth, unmeth) * n }
//   std::vector<int> row = {rname, strand, 0, 0, 0, 0}; // rname, strand, pos, ctx, meth, unmeth
//   for (int i=0; i<ctx.size(); i++) {
//     row[3] = ctx[i];
//     bool is_meth = (pass && ctx[i]>='A' && ctx[i]<='Z') ? true : false;
//     int found = xm.find(ctx[i], 0);
//     while (found != std::string::npos) {
//       row[2] = start + found;
//       row[4] = is_meth ? 1 : 0;
//       row[5] = is_meth ? 0 : 1;
//       res.insert(res.end(), row.begin(), row.end());
//       found = xm.find(ctx[i], found+1);
//     } 
//   }
//   return res;
// }
// 
// test code in R
//
// matrix(rcpp_parse_xm(2, 1, 3689466, "z......xh....Z.......Z......z...Zx...x..xh.hhhh..xhh.Zxh......xhh.......xh..x..x..x..............z......Z", TRUE, "zZ"), ncol=6, byrow=TRUE)
// matrix(rcpp_parse_xm(2, 1, 3689466, "z......xh....Z.......Z......z...Zx...x..xh.hhhh..xhh.Zxh......xhh.......xh..x..x..x..............z......Z", FALSE, "zZ"), ncol=6, byrow=TRUE)
  
  
