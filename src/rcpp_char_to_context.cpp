#include <Rcpp.h>
// using namespace Rcpp;

// converts XM char code to context string
//

// fast, vectorised
// [[Rcpp::export("rcpp_char_to_context")]]
std::vector<std::string> rcpp_char_to_context(std::vector<unsigned char> ctx)   // char
{
  std::map <unsigned char, std::string> ctx_map = { {'z',"CG"},  {'Z',"CG"},  {'x',"CHG"}, {'X',"CHG"}, {'h',"CHH"}, {'H',"CHH"}, {'u',"CN"},  {'U',"CN"} };
  std::vector<std::string> res (ctx.size());
  for (unsigned int x=0; x<ctx.size(); x++) {
    // checking for the interrupt
    if ((x & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();
    
    res[x] = ctx_map[ctx[x]];
  }
  return res;
}

// test code in R
//

/*** R
# z <- rep(utf8ToInt("ZHUX.zhux-+"),100000)
# microbenchmark::microbenchmark(rcpp_char_to_context(z))
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_char_to_context.cpp")

// #############################################################################
