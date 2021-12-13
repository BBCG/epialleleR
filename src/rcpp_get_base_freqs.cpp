#include <Rcpp.h>
// using namespace Rcpp;

// Matches reads with given 1-base positions (VCF) and returns base frequencies.
// FUNCTION ASSUMES THAT BOTH READS AND VCF ENTRIES ARE SORTED!
// Return value: std::array<int,32> for each position in VCF
// Array indices are calculated as (char&7 | pass?15:0 | strand?31:0),
// thus they appear as following: U+*8, U-*8, M+*8, M-*8
//
// char  &7   dec0b
// Aa    001  1
// Cc    011  3
// Gg    111  7
// Nn    110  6
// Tt    100  4

// MATCH VCF ENTRIES, RETURN BASE FREQS
// fast, vectorised
// [[Rcpp::export("rcpp_get_base_freqs")]]
Rcpp::NumericMatrix rcpp_get_base_freqs(Rcpp::DataFrame &df,                    // BAM data
                                        std::vector<bool> pass,                 // read passes the threshold?
                                        Rcpp::DataFrame &vcf)                   // VCF data
{
  Rcpp::IntegerVector read_rname = df["rname"];                                 // template rname
  Rcpp::IntegerVector read_strand = df["strand"];                               // template strand
  Rcpp::IntegerVector read_start = df["start"];                                 // template start
  // Rcpp::CharacterVector read_seq = bam["seq"];                                  // template SEQ
  Rcpp::XPtr<std::vector<std::string>> read_seq((SEXP)df.attr("seq_xptr"));     // merged refspaced template SEQ, as a pointer to std::vector<std::string>
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  Rcpp::IntegerVector vcf_chr = vcf["seqnames"];                                // VCF rname
  Rcpp::IntegerVector vcf_pos = vcf["start"];                                   // VCF start
  
  Rcpp::NumericMatrix res(vcf_pos.size(),32);
  
  int cur_vcf=0;
  for (unsigned int x=0; x<read_start.size(); x++) {
    // checking for the interrupt
    if ((x & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();
    
    int read_end = read_start[x] + read_seq->at(templid[x]).size() - 1;
    for (unsigned int i=cur_vcf; i<vcf_pos.size(); i++) {
      if (vcf_chr[i]<read_rname[x] || 
          (vcf_chr[i]==read_rname[x] && vcf_pos[i]<read_start[x])) {            // skip VCF if before read
        cur_vcf=i;
        continue;
      }
      if (vcf_chr[i]>read_rname[x] ||
          (vcf_chr[i]==read_rname[x] && vcf_pos[i]>read_end)) {                 // start over from cur_vcf for new read
        break;
      }
      if (vcf_chr[i]==read_rname[x] &&
          vcf_pos[i]>=read_start[x] && vcf_pos[i]<=read_end) {                  // match found
        int idx = read_seq->at(templid[x])[vcf_pos[i]-read_start[x]] & 7;
        if (read_strand[x]==2) idx |= 8;
        if (pass[x]) idx |= 16;
        res(i,idx)++;
      }
    }
  }

  return res;
}


// test code in R
//

/*** R
# microbenchmark::microbenchmark(rcpp_get_base_freqs(bam, pass, vcf.dt), times=10)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_get_base_freqs.cpp")

// #############################################################################
