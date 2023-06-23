#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

// [[Rcpp::depends(Rhtslib)]]

// Makes methylation calls using either genomic sequence or MM/ML tags
// and writes them in XM tag
//
// Returns simple statistics

// [[Rcpp::export]]
Rcpp::List rcpp_call_methylation_genome (std::string in_fn,                     // input BAM file name
                                         std::string out_fn)                    // output BAM file name
{
  
  
  
  
  
  // wrap and return the results
  Rcpp::List res = Rcpp::List::create(                                          // final List
    Rcpp::Named("n") = 0                                                        // numeric for number of reads with MM tags
  );
                                           
                                           return(res);
}


// #############################################################################
// test code and sourcing don't work on OS X
/*** R
*/
// #############################################################################