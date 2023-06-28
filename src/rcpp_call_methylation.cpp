#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rhtslib)]]

// Makes methylation calls using either genomic sequence or MM/ML tags
// and writes them in XM tag
//
// Returns simple statistics

// this one makes calls based on genomic sequence in the absence of MM/ML tags
// [[Rcpp::export]]
Rcpp::List rcpp_call_methylation_genome (std::string in_fn,                     // input BAM file name
                                         std::string out_fn,                    // output BAM file name
                                         Rcpp::List genome)                     // genome object
{
  
  // https://github.com/samtools/htslib/blob/e86126bc45f6cf6b1c3d67e13de96aeaa5e58805/test/test_view.c#L389
  
  
  
  // wrap and return the results
  Rcpp::List res = Rcpp::List::create(                                          // final List
    Rcpp::Named("n") = 0                                                        // numeric for number of reads with MM tags
  );
  
  return(res);
}


// #############################################################################
// test code and sourcing don't work on OS X
/*** R
setwd("~/work/packages/epialleleR/")
devtools::document()
devtools::load_all()
system.time(genome <- rcpp_read_genome("/scratch/ref/DRAGEN/hg38_plus_lambda_ChrY_PAR_masked.fa.gz"))
*/
// #############################################################################