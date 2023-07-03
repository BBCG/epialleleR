#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/bgzf.h>
#include <htslib/thread_pool.h>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(Rhtslib)]]

// Reads genomic sequences.
// Takes as an input either .fa or bgzipped .fa.gz (by means of HTSlib).
//
// Returns a list with:
// 1) field "rid"      - numeric ids of reference sequences
// 2) field "rname"    - names of reference sequences
// 3) field "rlen"     - lengths of reference sequences
// 4) attribute "rseq" - XPtr to std::vector of std::string with genomic sequences


// main sub that performs the reading
// [[Rcpp::export]]
Rcpp::List rcpp_read_genome (std::string fn,                                    // input: a name of (optionally bgzipped and/or indexed) FASTA file
                             int nthreads)                                      // HTSlib threads, >0 for multiple
{
  // containers
  std::vector<uint64_t> rid;                                                    // numeric ids of reference sequences
  std::vector<std::string> rname;                                               // names of reference sequences
  std::vector<uint64_t> rlen;                                                   // lengths of reference sequences
  std::vector<std::string>* rseq = new std::vector<std::string>;                // reference sequences themselves
  
  // file IO
  faidx_t *faidx = fai_load(fn.c_str());                                        // FASTA index
  if (faidx==NULL) Rcpp::stop("Unable to open FASTA index for reading");        // fall back if error
  
  hts_tpool *tpool = NULL;                                                      // thread pool works using a dirty hack
  if (nthreads>0) {
    tpool = hts_tpool_init(nthreads);
    bgzf_thread_pool(*(BGZF **)faidx, tpool, 0);                                // conversion of faidx_t to BGZF because it's first in the struct
  }

  // vars
  int flen = 0;                                                                 // fetched sequence length

  // fetch sequences
  for (size_t i=0; i<(unsigned int)faidx_nseq(faidx); i++) {
    rid.push_back(i);
    const char *name = faidx_iseq(faidx, i);
    int64_t length = faidx_seq_len(faidx, name);
    
    rname.push_back(name);
    rlen.push_back(length);
    
    char *sequence = faidx_fetch_seq(faidx, name, 0, length-1, &flen);          // try fetch sequence
    if (length!=flen) Rcpp::stop("Corrupted FASTA index. Delete and try again");// if fetched bytes differ from expected
    
    // encodeContextLookup(sequence, length);                                      // encode context
    
    rseq->emplace_back((const char*) sequence, length);                         // emplace encoded
    free(sequence);                                                             // free allocated
  }
  
  fai_destroy(faidx);                                                           // free allocated
  if (tpool) hts_tpool_destroy(tpool);                                          // free thread pool which isn't working yet
  
  // wrap and return the results
  Rcpp::List res = Rcpp::List::create(                                          // final List
    Rcpp::Named("rid") = rid,                                                   // numeric ids of reference sequences
    Rcpp::Named("rname") = rname,                                               // names of reference sequences
    Rcpp::Named("rlen") = rlen                                                  // lengths of reference sequences
  );
  
  Rcpp::XPtr<std::vector<std::string>> rseq_xptr(rseq, true);
  res.attr("rseq_xptr") = rseq_xptr;                                            // external pointer to sequences
  
  return(res);
}


// #############################################################################
// test code and sourcing don't work on OS X
/*** R
system.time(genome <- rcpp_read_genome("/scratch/ref/DRAGEN/hg38_plus_lambda_ChrY_PAR_masked.fa.gz", 1))
*/
// #############################################################################