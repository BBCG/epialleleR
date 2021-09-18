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
Rcpp::NumericMatrix rcpp_get_base_freqs(std::vector<int> read_rname,          // factor for read chr
                                        std::vector<int> read_strand,         // factor for read strand
                                        std::vector<int> read_start,          // start pos of reads
                                        std::vector<int> read_end,            // end pos of reads
                                        std::vector<std::string> read_seq,    // seq of reads
                                        std::vector<bool> pass,               // read passes the threshold?
                                        std::vector<int> vcf_chr,             // factor for VCF entries chr
                                        std::vector<int> vcf_pos)             // pos of VCF entries
{
  Rcpp::NumericMatrix res(vcf_pos.size(),32);
  
  int cur_vcf=0;
  for (unsigned int x=0; x<read_start.size(); x++) {
    // checking for the interrupt
    if ((x & 1048575) == 0) Rcpp::checkUserInterrupt();
    
    for (unsigned int i=cur_vcf; i<vcf_pos.size(); i++) {
      if (vcf_chr[i]<read_rname[x] || 
          (vcf_chr[i]==read_rname[x] && vcf_pos[i]<read_start[x])) {    // skip VCF if before read
        cur_vcf=i;
        continue;
      }
      if (vcf_chr[i]>read_rname[x] ||
          (vcf_chr[i]==read_rname[x] && vcf_pos[i]>read_end[x])) {      // start over from cur_vcf for new read
        break;
      }
      if (vcf_chr[i]==read_rname[x] &&
          vcf_pos[i]>=read_start[x] && vcf_pos[i]<=read_end[x]) {       // match found
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
freqs <-
rcpp_get_base_freqs(c(2,2,2,2,2),
                    c(1,2,1,2,2),
                    c(3067650,3067650,3067651,3067652,3067655),
                    c(3067788,3067858,3067740,3067928,3067812),
                    c("GTTGTGTTTTATTTTTTTTAATGTGGAGGTGTGGAGTTGTATGTGTGGGTTTTAGTGGAGTCGTTATAGGTTTTATTATATAATAAAGGGTAGGGAGGGTAAGGTTAGGAGTTTTGTGGGGTGGTGTAGAGATTTTGTA",
                      "GTTGTGTTTTATTTTTTTTAATGTGGAGGTGTGGAGTTGTATGTGTGGGTTTTAGTGGAGTTGTTATAGGTTTTATTATATAATAAAGGGTAGGGAGGGTNNNNNNNNGAGTTTTGTGGGGTGGTGTAGAGATTTTGTATTGTGGGTGTGATTTTTGTTTTGGGTTTAGGTTGATTTTGTGGTAGTTTTAAGTTTTGTTGTTGGGTAGG",
                       "TTGTGTTTTATTTTTTTTAATGTGGAGGTGTGGAGTTGTATGTGTGGGTTTTAGTGGAGTTGTTATAGGTTTTATTATATAATAAAGGGT", #AGGGAGGGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGATTTTGTGGTAGTTTTAAGTTTTGTTGTTGGGTAGGATTTTGGTGTTTTGGAAGGGTGTGGGGTAGTAATTTGGTTAGGGTAGTGTTTTGGTTTTTTTT
                        "TACATTTCACCTCCTTTAACACAAAAACACAAAATTACACATATAAATCTCAATAAAACCACCACAAATCTTATTACACAACAAAAAACAAAAAAAACAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACAACAACTCCAAATCCCATTACCAAACAAAACTTCAATACTCCAAAAAAATACAAAACAACAACCTAACCAAAACAACACCTCAATTTCCTTCCAATTCA",
                           "GTTTTATTTTTTTTAATGTGGAGGTGCGGAGTTGTATGTGTGGGTTTTAGTGGAGTCGTTATAGGTTTTATTATATAATAAAGGGTAGGGAGGGTAAGGTTAGGAGTTTTGTGGGGAGGTGTAGAGATTTTGTATTGTGGGTGTGATTTTTGTTTTGG"),
                    # ** *                                                       *              * *                *        *                                        *
                    c(TRUE,TRUE,TRUE,FALSE,FALSE),
                    c(2,2,2,2,2,2,2,2,2,10),
                    c(3067649,3067650,3067652,3067708,3067723,3067725,3067742,3067751,3067792,3067751))
matrix(freqs, ncol=32, dimnames=list(c(),c("","U+A","","U+C","U+T","","U+N","U+G",
                                           "","U-A","","U-C","U-T","","U-N","U-G",
                                           "","M+A","","M+C","M+T","","M+N","M+G",
                                           "","M-A","","M-C","M-T","","M-N","M-G")))

bam.slice   <- bam.processed #[rname=="chr17" & start<43125475 & start+width>43125475]
read_rname  <- as.integer(bam.slice$rname)
read_strand <- as.integer(bam.slice$strand)
read_start  <- bam.slice$start
read_end    <- bam.slice$start+bam.slice$width-1
read_seq    <- bam.slice$seq
pass        <- if(is.null(bam.slice$pass)) rep(TRUE,nrow(bam.slice)) else bam.slice$pass
vcf_chr     <- as.integer(seqnames(vgranges))
vcf_pos     <- start(vgranges)

z <- rcpp_get_base_freqs(read_rname, read_strand, read_start, read_end, read_seq, pass, vcf_chr, vcf_pos)
microbenchmark::microbenchmark(rcpp_get_base_freqs(read_rname, read_strand, read_start, read_end, read_seq, pass, vcf_chr, vcf_pos), times=10)

*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_get_base_freqs.cpp")

// #############################################################################
