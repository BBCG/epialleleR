#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>

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

// [[Rcpp::export]]
Rcpp::List rcpp_read_genome (std::string fn)                                    // input (optionally [b]gzipped and/or indexed) FASTA file name
{
  // containers
  std::vector<uint64_t> rid;                                                    // numeric ids of reference sequences
  std::vector<std::string> rname;                                               // names of reference sequences
  std::vector<uint64_t> rlen;                                                   // lengths of reference sequences
  std::vector<std::string>* rseq = new std::vector<std::string>;                // reference sequences themselves
  
  // file IO
  faidx_t *faidx = fai_load(fn.c_str());                                        // FASTA index
  if (faidx==NULL) Rcpp::stop("Unable to open FASTA index for reading");        // fall back if error
  
  // vars
  int flen = 0;                                                                 // fetched sequence length
  unsigned int C0 = 0, C1 = 0, G1 = 0, G2 = 0;                                  // C in position 0 and 1, G in position 1 and 2
  unsigned int ctx_plus = 0, ctx_minus = 0;                                     // cytosine contexts of 0+ and 2- bases
  
  // fetch sequences
  for (size_t i=0; i<faidx_nseq(faidx); i++) {
    rid.push_back(i);
    const char *name = faidx_iseq(faidx, i);
    uint64_t length = faidx_seq_len(faidx, name);
    
    rname.push_back(name);
    rlen.push_back(length);
    
    char *sequence = faidx_fetch_seq(faidx, name, 0, length-1, &flen);          // try fetch sequence
    if (length!=flen) Rcpp::stop("Corrupted FASTA index. Delete and try again");// if fetched bytes differ from expected
    // Rcpp::Rcout << "Fetched " << flen << std::endl;
    
    // convert to IUPAC
    for (size_t j=0; j<length; j++) {
      sequence[j] = seq_nt16_table[ (unsigned char) sequence[j] ] << 4 ;
    }
    
    // get context
    for (size_t j=0; j<length-2; j++) {
      C0 = (sequence[j]   == 0b00100000);                                       // IUPAC 'C' << 4
      C1 = (sequence[j+1] == 0b00100000);
      G1 = (sequence[j+1] == 0b01000000);                                       // IUPAC 'G' << 4
      G2 = (sequence[j+2] == 0b01000000);
      ctx_plus  = ( C0 << (G1 | G2) ) | (C0 & G1) ;
      ctx_minus = ( G2 << (C0 | C1) ) | (G2 & C1) ;
      sequence[j] = sequence[j]   | ( (unsigned char) ctx_plus << 2 );
      sequence[j] = sequence[j+2] | ( (unsigned char) ctx_minus );
    }
    
    rseq->emplace_back((const char*) sequence, length);
    free(sequence);
  }
  
  
  fai_destroy(faidx);
  
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
*/
// #############################################################################