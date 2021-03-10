#include <Rcpp.h>
// using namespace Rcpp;

// Fast merge of paired end sequences by means of std::basic_string and
// std::copy. Accepts as parameters five vectors of equal length (bam entries).
// Output: full insert seq
//


// 'std::copy' version, vectorised, the main one now
// [[Rcpp::export("rcpp_merge_ends")]]
std::vector<std::string> rcpp_merge_ends(std::vector<int> read1_pos,         // BAM pos field for read 1
                                         std::vector<std::string> read1_seq, // normalised BAM XM field for read 1
                                         std::vector<int> read2_pos,         // BAM pos field for read 2
                                         std::vector<std::string> read2_seq, // normalised BAM XM field for read 2
                                         std::vector<int> isize,             // BAM isize field
                                         char gap)                           // char to fill gaps with
{
  std::vector<std::string> res (read1_seq.size());
  
  for (unsigned int x=0; x<read1_pos.size(); x++) {
    // checking for the interrupt
    if ((x & 1048575) == 0) Rcpp::checkUserInterrupt();
    
    if (read1_seq[x].size()==abs(isize[x])) {
      res[x] = read1_seq[x];
      continue;
    }
    unsigned int start=read1_pos[x]<read2_pos[x] ? read1_pos[x] : read2_pos[x];
    res[x] = std::string (abs(isize[x]), gap);
    std::copy(read2_seq[x].begin(), read2_seq[x].end(), res[x].begin()+read2_pos[x]-start);
    std::copy(read1_seq[x].begin(), read1_seq[x].end(), res[x].begin()+read1_pos[x]-start);
  }
  return res;
}


// test code in R
//

/*** R
rcpp_merge_ends(c(31094072,156777322,9960854,70350904),
                c("A........x.....hhh.........hh...x.hhh.hh.h.h..xh......x.z.......zxhhh.hhhh..z.h..Z.........x.....h.B",
                  "Ezxh.h....h.z.h.z..h.hh......h..h...z.h...z...z..x........h....h..xhh...z...hh....x...h..z.h........F",
                  "I...h.............h...hh......hhx...hx.......x...h..h......hh.h.....h.hh..x...h..............h.x....J",
                  "M........h.h.h..h.......Z....hx......h.h.................h........h.x..xZ...Z.Z.....h.h.H.hx......hhN"),
                c(31094061,156777252,9960907,70351011),
                c("C...h......x........x.....hhh.........hh...x.hhh.hh.h.h..xh......x.z.......zxhhh.hhhh..z.h..Z.......D",
                  "G...z.z.hh..x.z..zxhh.zxh.zx.zx.zx...zxhh..z.h....x...z.z....h.z.h.h.z.zxh.h....h.z.h.z..h.hh......hH",
                  "K.....hh.h.....h.hh..x...h..............h.x....hh.h...h.........h..x.....h..x..............hh......L",
                  "O.Z.....xZ....Z...x....h.Z.....x.......Z.....h.hhx...h...h.........hhxZ..h.h............h..........P"),
                c(-111,-171,153,207), '.')

n <- 100000
read1_pos <- rep(c(31094072,156777322,9960854,70350904), n)
read1_seq <- rep(c("A........x.....hhh.........hh...x.hhh.hh.h.h..xh......x.z.......zxhhh.hhhh..z.h..Z.........x.....h.B",
                   "Ezxh.h....h.z.h.z..h.hh......h..h...z.h...z...z..x........h....h..xhh...z...hh....x...h..z.h........F",
                   "I...h.............h...hh......hhx...hx.......x...h..h......hh.h.....h.hh..x...h..............h.x....J",
                   "M........h.h.h..h.......Z....hx......h.h.................h........h.x..xZ...Z.Z.....h.h.H.hx......hhN"), n)
read2_pos <- rep(c(31094061,156777252,9960907,70351011), n)
read2_seq <- rep(c("C...h......x........x.....hhh.........hh...x.hhh.hh.h.h..xh......x.z.......zxhhh.hhhh..z.h..Z.......D",
                   "G...z.z.hh..x.z..zxhh.zxh.zx.zx.zx...zxhh..z.h....x...z.z....h.z.h.h.z.zxh.h....h.z.h.z..h.hh......hH",
                   "K.....hh.h.....h.hh..x...h..............h.x....hh.h...h.........h..x.....h..x..............hh......L",
                   "O.Z.....xZ....Z...x....h.Z.....x.......Z.....h.hhx...h...h.........hhxZ..h.h............h..........P"), n)
isize     <- rep(c(-111,-171,153,207), n)
gap       <- '.'


read1_pos <- bam.dt[isfirst==TRUE]$pos
read1_seq <- bam.dt[isfirst==TRUE]$XM.norm
read2_pos <- bam.dt[isfirst!=TRUE]$pos
read2_seq <- bam.dt[isfirst!=TRUE]$XM.norm
isize     <- bam.dt[isfirst==TRUE]$isize
gap       <- '-'
z <- rcpp_merge_ends(read1_pos, read1_seq, read2_pos, read2_seq, isize, gap)
# microbenchmark::microbenchmark(rcpp_merge_ends(read1_pos, read1_seq, read2_pos, read2_seq, isize, gap), times=10)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_merge_ends.cpp")

// #############################################################################
