// #include <stdlib.h>
#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
// #include <htslib/bgzf.h>
// #include <zlib.h>

// [[Rcpp::depends(Rhtslib)]]
// [[Rcpp::depends(zlibbioc)]]

// an attempt to read and preprocess BAM in place
// [[Rcpp::export("rcpp_read_bam")]]
int rcpp_read_bam (std::string fn                          // file name
)
{
  htsFile *bam_fp = hts_open(fn.c_str(), "r");                                  // open file
  bam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);                                    // read file header
  bam1_t *bam = bam_init1();                                                    // create BAM alignment structure

  int i=0;
  
  while( sam_read1(bam_fp, bam_hdr, bam) > 0 ) {   //
    i++;
  }

  bam_destroy1(bam);
  hts_close(bam_fp);
  
  return(i);
}






// test code in R
// Sourcing doesn't work on OS X

/*** R
rcpp_read_bam("/Users/oleksii.nikolaienko/work/packages/epialleleR/inst/extdata/amplicon000meth.bam")
# microbenchmark::microbenchmark(rcpp_apply_cigar(cigar,query,gap), times=10)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_read_bam.cpp")

// #############################################################################