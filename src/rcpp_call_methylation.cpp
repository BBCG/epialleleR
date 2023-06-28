#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

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
                                         Rcpp::List genome,                     // genome object
                                         int nthreads)                          // HTSlib threads, >0 for multiple
{
  // file IO
  htsFile *in_fp = hts_open(in_fn.c_str(), "r");                                // try open input file
  if (in_fp==NULL) Rcpp::stop("Unable to open input BAM file for reading");     // fall back if error
  htsFile *out_fp = hts_open(out_fn.c_str(), "wb");                             // try open input file
  if (out_fp==NULL) Rcpp::stop("Unable to open output BAM file for writing");   // fall back if error
  // shared thread pool
  htsThreadPool thread_pool = {NULL, 0};                                        // thread pool cuts time by 30%
  if (nthreads>0) {
    thread_pool.pool = hts_tpool_init(nthreads);                                // when initiated for >0 threads
    hts_set_opt(in_fp,  HTS_OPT_THREAD_POOL, &thread_pool);                     // and bound to the input
    hts_set_opt(out_fp, HTS_OPT_THREAD_POOL, &thread_pool);                     // and output file pointers
  }
  // https://github.com/samtools/htslib/blob/e86126bc45f6cf6b1c3d67e13de96aeaa5e58805/test/test_view.c#L389
  
  bam_hdr_t *in_hdr = sam_hdr_read(in_fp);                                      // try read input file header
  if (in_hdr==NULL) Rcpp::stop("Unable to read input BAM header");              // fall back if error
  if (sam_hdr_write(out_fp, in_hdr) < 0) Rcpp::stop("Unable to write header");  // try write output file header
  
  bam1_t *in_rec = bam_init1();                                                 // create BAM alignment structure
  
  // vars
  int nrecs = 0;                                                                // counters: BAM records
  int max_query_width = 1024;                                                   // max query width, expanded if necessary
  int query_width = max_query_width;                                            // query width
  uint8_t *xm = (uint8_t*) malloc(max_query_width * sizeof(uint8_t));           // XM array
  
  while(sam_read1(in_fp, in_hdr, in_rec) > 0) {                                 // read rec by rec
    nrecs++;                                                                    // BAM alignment records ++
    if ((nrecs & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();                     // every ~1M reads check for the interrupt
    
    query_width = abs(in_rec->core.l_qseq);                                     // query width
    if (query_width > max_query_width) {                                        // if sequence is longer than XM holder
      max_query_width = query_width;                                            // new max
      xm = (uint8_t *) realloc(xm, max_query_width);                            // expand xm holder
      if (xm==NULL) Rcpp::stop("No memory for BAM record #%i", nrecs);          // check memory allocation
    }
    
    // meth calling
    
    
    
    if (sam_write1(out_fp, in_hdr, in_rec) < 0) Rcpp::stop("Unable to write BAM"); // write rec
  }

  // cleaning
  bam_destroy1(in_rec);                                                         // clean BAM alignment structure 
  hts_close(in_fp);                                                             // close input BAM file
  hts_close(out_fp);                                                            // close output BAM file
  if (thread_pool.pool) hts_tpool_destroy(thread_pool.pool);                    // free thread pool
  free(xm);                                                                     // and free manually allocated memory
  
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
system.time(genome <- rcpp_read_genome("/Users/oleksii.nikolaienko/work/data/hg38/hg38_plus_lambda_ChrY_PAR_masked.fa.gz"))
inout <- c(system.file("extdata", "test", "paired-name-xm.bam", package="epialleleR"), "/tmp/out.bam")
rcpp_call_methylation_genome(inout[1], inout[2], genome, 0)
sapply(inout, file.size)
*/
// #############################################################################