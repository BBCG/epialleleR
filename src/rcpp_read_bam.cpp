#include <stdlib.h>
#include <Rcpp.h>
#include <htslib/sam.h>

// [[Rcpp::depends(Rhtslib)]]

// an attempt to read and preprocess BAM in place
// [[Rcpp::export("rcpp_read_bam")]]
int rcpp_read_bam (std::string fn                          // file name
)
{
  samFile *bam_fp = hts_open(fn.c_str(), "r");                                  // open file
  bam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);                                    // read file header
  bam1_t *bam = bam_init1();                                                    // create BAM alignment structure
  
  while(sam_read1(bam_fp, bam_hdr, bam) > 0) {   // 
    
  }
  
  bam_destroy1(bam);
  sam_close(bam_fp);
  
  return(0);
}






// test code in R
//

/*** R
rcpp_read_bam("~/work/packages/epialleleR/inst/extdata/amplicon000meth.bam")
# microbenchmark::microbenchmark(rcpp_apply_cigar(cigar,query,gap), times=10)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_read_bam.cpp")

// #############################################################################