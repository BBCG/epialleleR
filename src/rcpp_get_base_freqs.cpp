#include <Rcpp.h>
using namespace Rcpp;

// Matches reads with given 1-base positions (VCF) and returns base frequences.
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
NumericMatrix rcpp_get_base_freqs(std::vector<int> read_rname,          // factor for read chr
                                  std::vector<int> read_strand,         // factor for read strand
                                  std::vector<int> read_start,          // start pos of reads
                                  std::vector<int> read_end,            // end pos of reads
                                  std::vector<std::string> read_seq,    // seq of reads
                                  std::vector<bool> pass,               // read passes the threshold?
                                  std::vector<int> vcf_chr,             // factor for VCF entries chr
                                  std::vector<int> vcf_pos)             // pos of VCF entries
{
  NumericMatrix res(vcf_pos.size(),32);
  
  int cur_read=0;
  for (unsigned int i=0; i<vcf_pos.size(); i++) {
    for (unsigned int x=cur_read; x<read_start.size(); x++) {
      if (read_rname[x]<vcf_chr[i] || read_end[x]<vcf_pos[i]) {         // skip reads if VCF is ahead
        cur_read=x;
        continue;
      }
      if (read_rname[x]>vcf_chr[i] || read_start[x]>vcf_pos[i])         // start over from cur_read for new VCF
        break;
      if (read_rname[x]==vcf_chr[i] && read_start[x]<=vcf_pos[i] &&
          read_end[x]>=vcf_pos[i]) {                                    // match found
        int idx = read_seq[x][vcf_pos[i]-read_start[x]] & 7;
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
rcpp_get_base_freqs(c(2,2,2,2,2),
                    c(1,2,1,2,2),
                    c(3067650,3067650,3067651,3067652,3067655),
                    c(3067788,3067858,3067921,3067928,3067812),
                    c("GTTGTGTTTTATTTTTTTTAATGTGGAGGTGTGGAGTTGTATGTGTGGGTTTTAGTGGAGTCGTTATAGGTTTTATTATATAATAAAGGGTAGGGAGGGTAAGGTTAGGAGTTTTGTGGGGTGGTGTAGAGATTTTGTA",
                      "GTTGTGTTTTATTTTTTTTAATGTGGAGGTGTGGAGTTGTATGTGTGGGTTTTAGTGGAGTTGTTATAGGTTTTATTATATAATAAAGGGTAGGGAGGGTNNNNNNNNGAGTTTTGTGGGGTGGTGTAGAGATTTTGTATTGTGGGTGTGATTTTTGTTTTGGGTTTAGGTTGATTTTGTGGTAGTTTTAAGTTTTGTTGTTGGGTAGG",
                       "TTGTGTTTTATTTTTTTTAATGTGGAGGTGTGGAGTTGTATGTGTGGGTTTTAGTGGAGTTGTTATAGGTTTTATTATATAATAAAGGGTAGGGAGGGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGATTTTGTGGTAGTTTTAAGTTTTGTTGTTGGGTAGGATTTTGGTGTTTTGGAAGGGTGTGGGGTAGTAATTTGGTTAGGGTAGTGTTTTGGTTTTTTTT",
                        "TACATTTCACCTCCTTTAACACAAAAACACAAAATTACACATATAAATCTCAATAAAACCACCACAAATCTTATTACACAACAAAAAACAAAAAAAACAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACAACAACTCCAAATCCCATTACCAAACAAAACTTCAATACTCCAAAAAAATACAAAACAACAACCTAACCAAAACAACACCTCAATTTCCTTCCAATTCA",
                           "GTTTTATTTTTTTTAATGTGGAGGTGCGGAGTTGTATGTGTGGGTTTTAGTGGAGTCGTTATAGGTTTTATTATATAATAAAGGGTAGGGAGGGTAAGGTTAGGAGTTTTGTGGGGAGGTGTAGAGATTTTGTATTGTGGGTGTGATTTTTGTTTTGG"),
                    # ** *                                                       *              * *                *        *                                        *
                    c(TRUE,TRUE,TRUE,FALSE,FALSE),
                    c(2,2,2,2,2,2,2,2,2,10),
                    c(3067649,3067650,3067652,3067708,3067723,3067725,3067742,3067751,3067792,3067751))

read_rname  <- as.integer(bam$rname)
read_strand <- as.integer(bam$strand)
read_start  <- bam$start
read_end    <- bam$start+bam$width-1
read_seq    <- bam$seq
pass        <- bam$pass
vcf_chr     <- as.integer(seqnames(vcf.ranges))
vcf_pos     <- start(vcf.ranges)

z <- rcpp_get_base_freqs(read_rname, read_strand, read_start, read_end, read_seq, pass, vcf_chr, vcf_pos)
# microbenchmark::microbenchmark(rcpp_get_base_freqs(read_rname, read_strand, read_start, read_end, read_seq, pass, vcf_chr, vcf_pos), times=10)

*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_get_base_freqs.cpp")

// #############################################################################
