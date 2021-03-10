#include <Rcpp.h>
// using namespace Rcpp;

// Applies CIGAR to the sequence given (such as XM tag)
// Outpus query strings in a reference space

// Specifications are described in https://samtools.github.io/hts-specs/SAMv1.pdf


// fast vectorised alternative to GenomicAlignments::sequenceLayer
// [[Rcpp::export("rcpp_apply_cigar")]]
std::vector<std::string> rcpp_apply_cigar(std::vector<std::string> cigar,       // CIGAR strings
                                          std::vector<std::string> query)       // query (XM) strings
{
  std::vector<std::string> res(query.size());

  // iterating over input
  for (unsigned int x=0; x<cigar.size(); x++) {
    // checking for the interrupt
    if ((x & 1048575) == 0) Rcpp::checkUserInterrupt();
      
    // move to first nondigit
    int found = 0;
    while (found<cigar[x].size() & cigar[x][found]>='0' & cigar[x][found]<='9')
      found++;
    
    // quick copy if simple match
    if (found==cigar[x].size()-1 & cigar[x][found]=='M') {
      res[x] = query[x];
      continue;
    }

    int prev = 0;                                                               // previous cigar_op position
    int query_pos = 0;                                                          // current position within the query string

    while (found<cigar[x].size()) {
      int cigar_oplen = std::stoi(cigar[x].substr(prev, found-prev));

      // applying cigar_op
      switch(cigar[x][found]) {
        case 'M' :
        case '=' :
        case 'X' :
          res[x].append(query[x], query_pos, cigar_oplen);
          query_pos += cigar_oplen;
          break;
        case 'I' :
          query_pos += cigar_oplen;
          break;
        case 'D' :
        case 'N' :
          res[x] += std::string(cigar_oplen,'-');
          break;
        case 'S' :
          query_pos += cigar_oplen;
          break;
        case 'H' :
        case 'P' :
          break;
        default :
          Rcpp::stop("Unknown CIGAR operation at BAM line ", x);
      }

      // next op
      prev = ++found;
      while (found<cigar[x].size() & cigar[x][found]>='0' & cigar[x][found]<='9')
        found++;
    }

    // copy the rest as GenomicAlignments::sequenceLayer does
    if (query_pos < query[x].size()-1)
      res[x].append(query[x], query_pos, std::string::npos); // query[x].size()-query_pos
  }

  return res;
}


// test code in R
//

/*** R
rcpp_apply_cigar(c("10M2D4I12M","1M2I223M","226M"),
                 c("abcdefghijklmnopqrstuvwxyz",
                   "....h......x......Z.........x...........xh..x................x.......h.......h.....x.....Z...x.Z......h.....Z..h.....x.......h.......x....................................x........xh.h..hhh.Z...x..x.....Z..x..h.......xhhh.h....",
                   "..h.....x...........x..x...z.h.h.Z.....x............Zx.hh...x..h.hh..xhhh.....Z..x......x...Z.z.hhh.hhhh..Zx.h.h.hh...hhh........h..........x.Z.hhh.hhhh..Zx.h.h.h....hhhh.....Z..x......x...Z.Z.hhh.hhhh..Zx.h.h.hh...hhh........"))

cigar <- bam.dt$cigar
query <- bam.dt$XM
system.time(xm.my  <- rcpp_apply_cigar(cigar,query))
system.time(xm.ref <- as.character(GenomicAlignments::sequenceLayer(Biostrings::BStringSet(query, use.names=FALSE), cigar), use.names=FALSE))
identical(xm.my, xm.ref)
# microbenchmark::microbenchmark(rcpp_apply_cigar(cigar,query), times=10)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_cx_report.cpp")

// #############################################################################
