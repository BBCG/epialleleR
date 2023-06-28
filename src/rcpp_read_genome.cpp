#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/bgzf.h>
#include <htslib/thread_pool.h>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(Rhtslib)]]

// Reads genomic sequences, converts bases to IUPAC codes (HTSlib seq_nt16_str),
// adds cytosine context information to the output.
// Takes as an input either .fa or .fa.gz (by means of Boost).
//
// Returns a list with:
// 1) field "rid"      - numeric ids of reference sequences
// 2) field "rname"    - names of reference sequences
// 3) field "rlen"     - lengths of reference sequences
// 4) attribute "rseq" - vector of std::string in which every byte is:
//    a) first (leftmost) 4 bits   - IUPAC code of nucleotide
//    b) next (middle) 2 bits      - cytosine context in + strand
//    c) last (rightmost) 2 bits   - cytosine context in - strand
//
// Nucleotide base is encoded using HTSlib's seq_nt16_table:
//   https://github.com/samtools/htslib/blob/develop/hts.c
// Cytosine context is encoded as follows:
//   .=0, h=1, x=2, z=3 


// genomic sequence-to-context lookup tables
const unsigned char triad_forward_context[512] = {
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  16,  16,   0,  16,  16,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  16,  16,   0,  16,  16,
  0,  16,   0,  16,  16,   0,  16,  16,   0,   0,   0,   0,   0,   0,   0,   0,
  0,  16,   0,  16,  16,   0,  16,  16,   0,  16,   0,  16,  16,   0,  16,  16,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  36,   0,  36,  36,   0,  36,  40,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  36,   0,  36,  36,   0,  36,  40,
  0,  36,   0,  36,  36,   0,  36,  40,   0,   0,   0,   0,   0,   0,   0,   0,
  0,  36,   0,  36,  36,   0,  36,  40,   0,  44,   0,  44,  44,   0,  44,  44,
  0,   0,   0,   0,   0,   0,   0,   0,   0, 128,   0, 128, 128,   0, 128, 128,
  0,   0,   0,   0,   0,   0,   0,   0,   0, 128,   0, 128, 128,   0, 128, 128,
  0, 128,   0, 128, 128,   0, 128, 128,   0,   0,   0,   0,   0,   0,   0,   0,
  0, 128,   0, 128, 128,   0, 128, 128,   0, 128,   0, 128, 128,   0, 128, 128,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0, 240,   0, 240, 240,   0, 240, 240,
  0,   0,   0,   0,   0,   0,   0,   0,   0, 240,   0, 240, 240,   0, 240, 240,
  0, 240,   0, 240, 240,   0, 240, 240,   0,   0,   0,   0,   0,   0,   0,   0,
  0, 240,   0, 240, 240,   0, 240, 240,   0, 240,   0, 240, 240,   0, 240, 240,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  64,   0,  64,  64,   0,  64,  64,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  64,   0,  64,  64,   0,  64,  64,
  0,  64,   0,  64,  64,   0,  64,  64,   0,   0,   0,   0,   0,   0,   0,   0,
  0,  64,   0,  64,  64,   0,  64,  64,   0,  64,   0,  64,  64,   0,  64,  64
};
const unsigned char triad_reverse_context[512] = {
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  65,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  67,
  0,  16,   0,  32, 128,   0, 240,  65,   0,   0,   0,   0,   0,   0,   0,   0,
  0,  16,   0,  32, 128,   0, 240,  65,   0,  16,   0,  32, 128,   0, 240,  65,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  66,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  67,
  0,  16,   0,  32, 128,   0, 240,  66,   0,   0,   0,   0,   0,   0,   0,   0,
  0,  16,   0,  32, 128,   0, 240,  66,   0,  16,   0,  32, 128,   0, 240,  66,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  65,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  67,
  0,  16,   0,  32, 128,   0, 240,  65,   0,   0,   0,   0,   0,   0,   0,   0,
  0,  16,   0,  32, 128,   0, 240,  65,   0,  16,   0,  32, 128,   0, 240,  65,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  65,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  67,
  0,  16,   0,  32, 128,   0, 240,  65,   0,   0,   0,   0,   0,   0,   0,   0,
  0,  16,   0,  32, 128,   0, 240,  65,   0,  16,   0,  32, 128,   0, 240,  65,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  65,
  0,   0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  67,
  0,  16,   0,  32, 128,   0, 240,  65,   0,   0,   0,   0,   0,   0,   0,   0,
  0,  16,   0,  32, 128,   0, 240,  65,   0,  16,   0,  32, 128,   0, 240,  65
};

// Encodes sequence and cytosine context using lookup tables
inline int encodeContextLookup (char *source, size_t size)
{
  const size_t buf_size = 64;                                                   // size of chunks
  const size_t buf_overlap = 4;                                                 // overlap = 2nt+2nt
  char seq[buf_size];                                                           // source seq
  char end[buf_size];                                                           // tail seq
  char out[buf_size];                                                           // encoded seq + context
  unsigned int idx;                                                             // index for lookup
  
// macros
#define main_loop {                           /* single-pass encode context */ \
  for (size_t i=0; i<buf_size-2; i++) {                    /* loop over buf */ \
    idx = ((seq[i] &7) << 6) | ((seq[i+1] &7) << 3) | (seq[i+2] &7); /* idx */ \
    out[i+2] = triad_reverse_context[idx];        /* seq + context of third */ \
    out[i] = out[i] | triad_forward_context[idx]; /* seq + context of first */ \
  }                                                                            \
};

  memset(end, 'N', buf_size);                                                   // prepare tail buffer
  std::memcpy(end, source+size-buf_size+2, buf_size-2);                         // read last chunk: "...NN"
  
  memset(seq, 'N', buf_size);                                                   // prepare seq buffer
  std::memcpy(seq+2, source, buf_size-2);                                       // read first chunk: "NN..."
  for (size_t b=1; b <= size/(buf_size-buf_overlap); b++) {                     // chunk by chunk
    memset(out, 0, sizeof out);                                                 // clean destination
    main_loop;                                                                  // one-pass encode
    std::memcpy(seq, source+b*(buf_size-buf_overlap)-2, buf_size);              // read another chunk
    std::memcpy(source+(b-1)*(buf_size-buf_overlap), out+2, buf_size-buf_overlap); // copy to source
  }
  
  // now, the tail
  memset(out, 0, sizeof out);                                                   // clean destination
  std::memcpy(seq, end, buf_size);                                              // copy tail to seq buffer
  main_loop;                                                                    // one-pass encode
  std::memcpy(source+size-buf_size+buf_overlap, out+2, buf_size-buf_overlap);   // copy to source
  
  return (size);                                                                // number of bases processed
}



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
  
  // hts_tpool *tpool = NULL;                                                      // thread pool not working yet
  // if (nthreads>0) {
  //   tpool = hts_tpool_init(nthreads);
  //   bgzf_thread_pool(faidx->bgzf, tpool, 0);
  // }
  // 
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
    
    encodeContextLookup(sequence, length);                                      // encode context
    
    rseq->emplace_back((const char*) sequence, length);                         // emplace encoded
    free(sequence);                                                             // free allocated
  }
  
  fai_destroy(faidx);                                                           // free allocated
  // if (tpool) hts_tpool_destroy(tpool);                                          // free thread pool which isn't working yet
  
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
system.time(genome <- rcpp_read_genome("/scratch/ref/DRAGEN/hg38_plus_lambda_ChrY_PAR_masked.fa.gz"))
*/
// #############################################################################