#include <Rcpp.h>
using namespace Rcpp;

// converts XM char code to context string
//

// fast, vectorised
// [[Rcpp::export("rcpp_char_to_context_vector")]]
std::vector<std::string> rcpp_char_to_context_vector(std::vector<int> ctx) {
  std::map <int, std::string> ctx_map = {
    {'z',"CG"},  {'Z',"CG"},  {'x',"CHG"}, {'X',"CHG"},
    {'h',"CHH"}, {'H',"CHH"}, {'u',"UNK"}, {'U',"UNK"}
  };
  std::vector<std::string> res (ctx.size());
  for (int x=0; x<ctx.size(); x++) {
    res[x] = ctx_map[ctx[x]];
  }
  return res;
}


// test code in R
//

/*** R
rcpp_char_to_context_vector(utf8ToInt("ZHUX.zhux-+"))
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_char_to_context.cpp")

// #############################################################################