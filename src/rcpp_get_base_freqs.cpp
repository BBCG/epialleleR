#include <Rcpp.h>
#include <htslib/hts.h>
#include "epialleleR.h"

// Matches reads with given 1-base positions (VCF) and returns base frequencies.
// FUNCTION ASSUMES THAT BOTH READS AND VCF ENTRIES ARE SORTED!
// Return value: std::array<int,20> for each position in VCF
// Array indices are 0123 for ACGT, and 4 for N and extended IUPAC,
// thus they appear as following: U+*5, U-*5, M+*5, M-*5


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
  Rcpp::XPtr<std::vector<std::string>> seqxm((SEXP)df.attr("seqxm_xptr"));      // merged refspaced packed template SEQXMs, as a pointer to std::vector<std::string>
  Rcpp::IntegerVector templid = df["templid"];                                  // template id, effectively holds indexes of corresponding std::string in std::vector
  
  Rcpp::IntegerVector vcf_chr = vcf["seqnames"];                                // VCF rname
  Rcpp::IntegerVector vcf_pos = vcf["start"];                                   // VCF start
  
  Rcpp::NumericMatrix res(vcf_pos.size(),20);
  
  int cur_vcf=0;
  for (unsigned int x=0; x<read_start.size(); x++) {
    // checking for the interrupt
    if ((x & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();
    
    const int read_rname_x = read_rname[x];
    const int read_start_x = read_start[x];
    const int read_end_x = read_start_x + seqxm->at(templid[x]).size() - 1;
    for (unsigned int i=cur_vcf; i<vcf_pos.size(); i++) {
      const int vcf_chr_i = vcf_chr[i];
      const int vcf_pos_i = vcf_pos[i];
      if (vcf_chr_i<read_rname_x || 
          (vcf_chr_i==read_rname_x && vcf_pos_i<read_start_x)) {                // skip VCF if before read
        cur_vcf=i;
        continue;
      }
      if (vcf_chr_i>read_rname_x ||
          (vcf_chr_i==read_rname_x && vcf_pos_i>read_end_x)) {                  // start over from cur_vcf for new read
        break;
      }
      if (vcf_chr_i==read_rname_x &&
          vcf_pos_i>=read_start_x && vcf_pos_i<=read_end_x) {                   // match found
        int idx = seq_nt16_int[unpack_seq_idx(seqxm->at(templid[x])[vcf_pos_i-read_start_x])]; // index of a base, [0;4]
        idx += (read_strand[x]-1) * 5;                                          // shift by 5 if '-' strand (==2)
        idx += ((bool)(pass[x])) * 10;                                          // shift by 10 if pass==TRUE (==1)
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
