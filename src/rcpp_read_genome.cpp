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

// genomic sequence-to-context lookup index
#define triad_to_idx(t) ( (t[0]&7)<<6 + (t[1]&7)<<3 + (t[2]&7) )

// genomic sequence-to-context lookup tables
const unsigned char triad_forward_context[512] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  16,  16,   0,  16,  16,   0,
    0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  16,  16,   0,  16,  16,   0,
   16,   0,  16,  16,   0,  16,  16,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
   16,   0,  16,  16,   0,  16,  16,   0,  16,   0,  16,  16,   0,  16,  16,   0, 
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   
    0,   0,   0,   0,   0,   0,   0,   0,  36,   0,  36,  36,   0,  36,  40,   0,  
    0,   0,   0,   0,   0,   0,   0,   0,  36,   0,  36,  36,   0,  36,  40,   0,  
   36,   0,  36,  36,   0,  36,  40,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
   36,   0,  36,  36,   0,  36,  40,   0,  44,   0,  44,  44,   0,  44,  44,   0,   
    0,   0,   0,   0,   0,   0,   0,   0, 128,   0, 128, 128,   0, 128, 128,   0,   
    0,   0,   0,   0,   0,   0,   0,   0, 128,   0, 128, 128,   0, 128, 128,   0, 
  128,   0, 128, 128,   0, 128, 128,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  128,   0, 128, 128,   0, 128, 128,   0, 128,   0, 128, 128,   0, 128, 128,   0,  
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
    0,   0,   0,   0,   0,   0,   0,   0, 240,   0, 240, 240,   0, 240, 240,   0, 
    0,   0,   0,   0,   0,   0,   0,   0, 240,   0, 240, 240,   0, 240, 240,   0, 
  240,   0, 240, 240,   0, 240, 240,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  240,   0, 240, 240,   0, 240, 240,   0, 240,   0, 240, 240,   0, 240, 240,   0, 
    0,   0,   0,   0,   0,   0,   0,   0,  64,   0,  64,  64,   0,  64,  64,   0,  
    0,   0,   0,   0,   0,   0,   0,   0,  64,   0,  64,  64,   0,  64,  64,   0,  
   64,   0,  64,  64,   0,  64,  64,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
   64,   0,  64,  64,   0,  64,  64,   0,  64,   0,  64,  64,   0,  64,  64,   0
};
const unsigned char triad_reverse_context[512] = {
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
    0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  65,   0, 
    0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  65,   0,
   16,   0,  32, 128,   0, 240,  65,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
   16,   0,  32, 128,   0, 240,  65,   0,  16,   0,  32, 128,   0, 240,  67,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  66,   0, 
    0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  66,   0,
   16,   0,  32, 128,   0, 240,  66,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
   16,   0,  32, 128,   0, 240,  66,   0,  16,   0,  32, 128,   0, 240,  67,   0,
    0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  65,   0,
    0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  65,   0,
   16,   0,  32, 128,   0, 240,  65,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   16,   0,  32, 128,   0, 240,  65,   0,  16,   0,  32, 128,   0, 240,  67,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  65,   0,
    0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  65,   0, 
   16,   0,  32, 128,   0, 240,  65,   0,   0,   0,   0,   0,   0,   0,   0,   0, 
   16,   0,  32, 128,   0, 240,  65,   0,  16,   0,  32, 128,   0, 240,  67,   0, 
    0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  65,   0,
    0,   0,   0,   0,   0,   0,   0,   0,  16,   0,  32, 128,   0, 240,  65,   0,
   16,   0,  32, 128,   0, 240,  65,   0,   0,   0,   0,   0,   0,   0,   0,   0,
   16,   0,  32, 128,   0, 240,  65,   0,  16,   0,  32, 128,   0, 240,  67,   0
};

// Encodes sequence and cytosine context using a simple loop
// Try rewrite this using lookup table or SIMD intrinsics
inline int encodeContextLoop (char *source, size_t size)
{
  const size_t buf_size = 16;
  const size_t buf_overlap = 4;
  char seq[buf_size];
  char isC[buf_size];
  char isG[buf_size];
  char Fctx[buf_size];
  char Rctx[buf_size];
  
  memset(seq, 'N', 2);
  // read first batch
  std::memcpy(seq+2, source, buf_size-2);
  // batch by batch
  for (size_t b=1; b <= size/(buf_size-buf_overlap); b++) {
    // IUPAC encode
    for (size_t i=0; i<buf_size; i++) {
      seq[i] = seq_nt16_table[ (unsigned char) seq[i] ];
    }
    // check if Cs or Gs
    for (size_t i=0; i<buf_size; i++) {
      isC[i] = seq[i]==seq_nt16_table['C'];                                     // is C?
      isG[i] = seq[i]==seq_nt16_table['G'];                                     // is G?
    }
    // calculate context
    for (size_t i=0; i<buf_size-2; i++) {
      Fctx[i] = ( isC[i] << (isG[i+1] | isG[i+2]) ) | (isC[i] & isG[i+1]);
      Rctx[i+2] = ( isG[i+2] << (isC[i] | isC[i+1]) ) | (isC[i+1] & isG[i+2]);
    }
    // pack sequence and context
    for (size_t i=buf_overlap/2; i<buf_size-buf_overlap/2; i++) {
      seq[i] = (seq[i]<<4) | (Fctx[i]<<2) | (Rctx[i]);
    }
    // copy to destination
    std::memcpy(source+(b-1)*(buf_size-buf_overlap), seq+2, buf_size-buf_overlap);
    // read another chunk
    std::memcpy(seq, source+b*(buf_size-buf_overlap)-2, buf_size);
  }
  // number of bases processed
  return (size/(buf_size-buf_overlap))*(buf_size-buf_overlap);
}


// Encodes sequence and cytosine context using lookup tables
inline int encodeContextLookup (char *source, size_t size)
{
  const size_t buf_size = 16;
  const size_t buf_overlap = 4;
  char seq[buf_size];
  char out[buf_size];
  
  memset(seq, 'N', 2);
  // read first batch
  std::memcpy(seq+2, source, buf_size-2);
  // batch by batch
  for (size_t b=1; b <= size/(buf_size-buf_overlap); b++) {
    // clean out
    memset(out, 0, sizeof out);
    // one-pass encode
    for (size_t i=0; i<buf_size-2; i++) {
      out[i] = out[i] | triad_forward_context[triad_to_idx(seq+i)];
      out[i+2] = out[i+2] | triad_reverse_context[triad_to_idx(seq+i)];
    }
    // copy to destination
    std::memcpy(source+(b-1)*(buf_size-buf_overlap), out+2, buf_size-buf_overlap);
    // read another chunk
    std::memcpy(seq, source+b*(buf_size-buf_overlap)-2, buf_size);
  }
  // number of bases processed
  return (size/(buf_size-buf_overlap))*(buf_size-buf_overlap);
}



// temp sub to see the seq and context
int decodeContext (const char *p, size_t s) {
  char ctx_map[] = ".hxz";
  const size_t buf_size = 1024;
  char seq[buf_size];
  char Fctx[buf_size];
  char Rctx[buf_size];
  for (size_t i=0; i<std::min(s, buf_size); i++) {
    seq[i]  = seq_nt16_str[(unsigned char) (p[i] & 0b11110000) >> 4];
    Fctx[i] = ctx_map[(unsigned char) (p[i] & 0b00001100) >> 2];
    Rctx[i] = ctx_map[(unsigned char) (p[i] & 0b00000011) ];
  }
  
  Rcpp::Rcout << " Seq:" << std::string(seq,  std::min(s, buf_size)) << std::endl;
  Rcpp::Rcout << "Fctx:" << std::string(Fctx, std::min(s, buf_size)) << std::endl;
  Rcpp::Rcout << "Rctx:" << std::string(Rctx, std::min(s, buf_size)) << std::endl;
  
  return 0;
}


// [[Rcpp::export]]
Rcpp::List rcpp_read_genome (std::string fn)                                    // input: a name of (optionally bgzipped and/or indexed) FASTA file 
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

  // fetch sequences
  for (size_t i=0; i<(unsigned int)faidx_nseq(faidx); i++) {
    rid.push_back(i);
    const char *name = faidx_iseq(faidx, i);
    int64_t length = faidx_seq_len(faidx, name);
    
    rname.push_back(name);
    rlen.push_back(length);
    
    char *sequence = faidx_fetch_seq(faidx, name, 0, length-1, &flen);          // try fetch sequence
    if (length!=flen) Rcpp::stop("Corrupted FASTA index. Delete and try again");// if fetched bytes differ from expected
    
    // save the tail
    char tail[24];
    memset(tail+20, 'N', 4);
    memcpy(tail, sequence+length-20, 20);
    
    // // check before packing
    // Rcpp::Rcout << "Fetched " << flen << std::endl;
    // Rcpp::Rcout << " BGN:" << std::string(sequence, 64) << std::endl;
    // Rcpp::Rcout << " END:" << std::string(sequence+length-64, 64) << std::endl;

    // encode context
    encodeContextLoop(sequence, length);                                        // sequence first
    encodeContextLoop(tail, sizeof tail);                                       // then tail
    memcpy(sequence+length-16, tail+4, 16);                                     // put 16 last at their place
    
    // // check after packing
    // decodeContext(sequence, 64);
    // decodeContext(sequence+length-64, 64);
    
    // emplace
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



// temp sub to validate context conversion
// [[Rcpp::export]]
int genometest ()
{
  char full_seq[] = "GTGACAGCCGCCCTTGGGAGACGACGGCGTCTGCAACCAGCAGCCTCCAAAGGGTGCAGCCAGGAGGCTCAGCTTGTCCGCCTCCGGGGCTCGGGGCTAAG"; //ACTGACTGCCGGCCGGCCCGGGGCTCGGGGAG";
  char corrFctx[] = "....x..xz.hhh........z..z..z..x..h..hx..x..hh.hh........x..hx......h.x..h....xz.hh.xz....h.z....h....";
  char corrRctx[] = "h.h...x..z.....hhh.h..z..zx.z...x......x..x........hhh.h..x...xh.hh....x...x...z.....zxhh...zxhh....h";
  
  size_t sizeof_full_seq = sizeof full_seq-1;
  
  Rcpp::Rcout << "  In:" << std::string(full_seq, sizeof full_seq) << std::endl;
  Rcpp::Rcout << "Fctx:" << std::string(corrFctx, sizeof corrFctx) << std::endl;
  Rcpp::Rcout << "Rctx:" << std::string(corrRctx, sizeof corrRctx) << std::endl;
  
  char tail[24];
  memset(tail+20, 'N', 4);
  memcpy(tail, full_seq + sizeof_full_seq - 20, 20);
  Rcpp::Rcout << "Tail:" << std::string(tail, sizeof tail) << std::endl;
  
  int n = encodeContextLoop(full_seq, sizeof_full_seq);
  Rcpp::Rcout << "n=" << n << std::endl;
  
  decodeContext(full_seq, sizeof_full_seq);
  
  encodeContextLoop(tail, sizeof tail);
  memcpy(full_seq + sizeof_full_seq - 16, tail+4, 16);
  
  decodeContext(full_seq, sizeof_full_seq);
  
  return 0;
}



// temp sub to testread genome with Boost
// [[Rcpp::export]]
int boosttest ()
{
  return 0;
}





// #############################################################################
// test code and sourcing don't work on OS X
/*** R
setwd("~/work/packages/epialleleR/")
devtools::document()
devtools::load_all()
genometest()
system.time(genome <- rcpp_read_genome("/scratch/ref/DRAGEN/hg38_plus_lambda_ChrY_PAR_masked.fa.gz"))
*/
// #############################################################################