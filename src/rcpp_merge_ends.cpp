#include <Rcpp.h>
using namespace Rcpp;

// Fast merge of paired end sequences by means of std::basic_string and
// std::copy. Accepts as parameters five vectors of equal length (bam entries).
// Output: full insert seq
//


// 'std::copy' version, vectorised, the main one now
// [[Rcpp::export("rcpp_merge_ends_vector")]]
std::vector<std::string> rcpp_merge_ends_vector(std::vector<int> read1_pos,         // BAM pos field for read 1
                                                std::vector<std::string> read1_seq, // normalised BAM XM field for read 1
                                                std::vector<int> read2_pos,         // BAM pos field for read 2
                                                std::vector<std::string> read2_seq, // normalised BAM XM field for read 2
                                                std::vector<int> isize) {           // BAM isize field
  std::vector<std::string> res;
  for (int x=0; x<read1_pos.size(); x++) {
    int start=read1_pos[x]<read2_pos[x]?read1_pos[x]:read2_pos[x];
    std::string insert(abs(isize[x]),'.');
    std::copy(read2_seq[x].begin(), read2_seq[x].end(), insert.begin()+read2_pos[x]-start);
    std::copy(read1_seq[x].begin(), read1_seq[x].end(), insert.begin()+read1_pos[x]-start);
    res.push_back(insert);
  }
  return res;
}


// test code in R
//

/*** R
rcpp_merge_ends_vector(c(31094072,156777322,9960854,70350904),
                       c("A........x.....hhh.........hh...x.hhh.hh.h.h..xh......x.z.......zxhhh.hhhh..z.h..Z.........x.....h.B",
                         "Ezxh.h....h.z.h.z..h.hh......h..h...z.h...z...z..x........h....h..xhh...z...hh....x...h..z.h........F",
                         "I...h.............h...hh......hhx...hx.......x...h..h......hh.h.....h.hh..x...h..............h.x....J",
                         "M........h.h.h..h.......Z....hx......h.h.................h........h.x..xZ...Z.Z.....h.h.H.hx......hhN"),
                       c(31094061,156777252,9960907,70351011),
                       c("C...h......x........x.....hhh.........hh...x.hhh.hh.h.h..xh......x.z.......zxhhh.hhhh..z.h..Z.......D",
                         "G...z.z.hh..x.z..zxhh.zxh.zx.zx.zx...zxhh..z.h....x...z.z....h.z.h.h.z.zxh.h....h.z.h.z..h.hh......hH",
                         "K.....hh.h.....h.hh..x...h..............h.x....hh.h...h.........h..x.....h..x..............hh......L",
                         "O.Z.....xZ....Z...x....h.Z.....x.......Z.....h.hhx...h...h.........hhxZ..h.h............h..........P"),
                       c(-111,-171,153,207)) 
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_merge_ends.cpp")

// #############################################################################


// // 'std::copy' version, non-vectorised
// // [[Rcpp::export]]
// std::string rcpp_merge_ends(std::vector<int> flag, std::vector<int> pos, std::vector<int> isize, std::vector<std::string> seq) {
//   int first_read_idx=(flag[0]&128)>>7;
//   int second_read_idx=first_read_idx^1;
//   int start=pos[0]<pos[1]?pos[0]:pos[1];
//   std::string insert(abs(isize[1]),'.');
//   std::copy(seq[second_read_idx].begin(), seq[second_read_idx].end(), insert.begin()+pos[second_read_idx]-start);
//   std::copy(seq[first_read_idx].begin(), seq[first_read_idx].end(), insert.begin()+pos[first_read_idx]-start);
//   return insert;
// }

// test code in R
//
// rcpp_merge_ends(c(83,163), c(31094072,31094061), c(-111,111), c("A........x.....hhh.........hh...x.hhh.hh.h.h..xh......x.z.......zxhhh.hhhh..z.h..Z.........x.....h.B",
//                   "C...h......x........x.....hhh.........hh...x.hhh.hh.h.h..xh......x.z.......zxhhh.hhhh..z.h..Z.......D"))
//   rcpp_merge_ends(c(83,163), c(156777322,156777252), c(-171,171), c("Ezxh.h....h.z.h.z..h.hh......h..h...z.h...z...z..x........h....h..xhh...z...hh....x...h..z.h........F",
//                     "G...z.z.hh..x.z..zxhh.zxh.zx.zx.zx...zxhh..z.h....x...z.z....h.z.h.h.z.zxh.h....h.z.h.z..h.hh......hH"))
//   rcpp_merge_ends(c(99,147), c(9960854,9960907), c(153,-153), c("I...h.............h...hh......hhx...hx.......x...h..h......hh.h.....h.hh..x...h..............h.x....J",
//                     "K.....hh.h.....h.hh..x...h..............h.x....hh.h...h.........h..x.....h..x..............hh......L"))
//   rcpp_merge_ends(c(99,147), c(70350904,70351011), c(207,-207), c("M........h.h.h..h.......Z....hx......h.h.................h........h.x..xZ...Z.Z.....h.h.H.hx......hhN",
//                     "O.Z.....xZ....Z...x....h.Z.....x.......Z.....h.hhx...h...h.........hhxZ..h.h............h..........P"))
  




// // 'std::basic_string::replace' version, no use
// // [[Rcpp::export]]
// std::string rcpp_merge_ends(std::vector<int> flag, std::vector<int> pos, std::vector<int> isize, std::vector<std::string> seq) {
//   int first_read_idx=(flag[0]&128)>>7;
//   int second_read_idx=first_read_idx^1;
//   int start=pos[0]<pos[1]?pos[0]:pos[1];
//   std::string insert(abs(isize[1]),'.');
//   insert.replace(abs(pos[second_read_idx]-start), seq[second_read_idx].size(), seq[second_read_idx]);
//   insert.replace(abs(pos[first_read_idx]-start), seq[first_read_idx].size(), seq[first_read_idx]);
//   return(insert);
// }


