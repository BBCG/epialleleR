#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rhtslib)]]

// Makes methylation calls using either genomic sequence [done] or MM/ML tags
// [pending] and writes them in XM tag.
// Relies on the presence of XG (Illumina/Bismark) or YC (bwa-meth) genome
// conversion tags.
//
// This calling is read-based: the cytosine context of the reference sequence 
// is calculated for every read in case there's a substitution in
// the query upstream/downstream of the context and if it changes the context.
// E.g., if ACGT in the reference becomes ACNT, then my context becomes .h..
// This seems to be identical to Illumina DRAGEN. The only difference is that
// I don't use Unknown context for CNN or NNG, I call them as h.. and ..h here.
//
// Returns simple statistics on records parsed and calls made.


// genomic sequence-to-context lookup tables
const unsigned char triad_forward_context[512] = {
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', 'h', '.', 'h', 'h', '.', 'h', 'x',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', 'h', '.', 'h', 'h', '.', 'h', 'x',
  '.', 'h', '.', 'h', 'h', '.', 'h', 'x', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', 'h', '.', 'h', 'h', '.', 'h', 'x', '.', 'z', '.', 'z', 'z', '.', 'z', 'z',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'
};
const unsigned char triad_reverse_context[512] = {
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'h',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'z',
  '.', '.', '.', '.', '.', '.', '.', 'h', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', 'h', '.', '.', '.', '.', '.', '.', '.', 'h',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'x',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'z',
  '.', '.', '.', '.', '.', '.', '.', 'x', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', 'x', '.', '.', '.', '.', '.', '.', '.', 'x',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'h',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'z',
  '.', '.', '.', '.', '.', '.', '.', 'h', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', 'h', '.', '.', '.', '.', '.', '.', '.', 'h',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'h',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'z',
  '.', '.', '.', '.', '.', '.', '.', 'h', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', 'h', '.', '.', '.', '.', '.', '.', '.', 'h',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'h',
  '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'z',
  '.', '.', '.', '.', '.', '.', '.', 'h', '.', '.', '.', '.', '.', '.', '.', '.',
  '.', '.', '.', '.', '.', '.', '.', 'h', '.', '.', '.', '.', '.', '.', '.', 'h'
};


// this one makes calls based on genomic sequence in the absence of MM/ML tags
// [[Rcpp::export]]
Rcpp::List rcpp_call_methylation_genome (std::string in_fn,                     // input BAM file name
                                         std::string out_fn,                    // output BAM file name
                                         Rcpp::List &genome,                    // genome object (list+XPtr)
                                         std::string tag,                       // what tag to read genome strand from (XG/YC)
                                         int nthreads)                          // HTSlib threads, >0 for multiple
{
  // genome data
  std::vector<std::string> rname = genome["rname"];                             // reference sequence names
  std::vector<uint64_t> rlen = genome["rlen"];                                  // reference sequence lengths
  Rcpp::XPtr<std::vector<std::string>> rseq((SEXP)genome.attr("rseq_xptr"));    // packed reference sequence + context, as a pointer to std::vector<std::string>
  
  // file IO
  htsFile *in_fp = hts_open(in_fn.c_str(), "r");                                // try open input file
  if (in_fp==NULL) Rcpp::stop("Unable to open input BAM file for reading");     // fall back if error
  htsFile *out_fp = hts_open(out_fn.c_str(), "wb");                             // try open output file
  if (out_fp==NULL) Rcpp::stop("Unable to open output BAM file for writing");   // fall back if error
  // shared thread pool
  htsThreadPool thread_pool = {NULL, 0};                                        // thread pool cuts time by 30%
  if (nthreads>0) {
    thread_pool.pool = hts_tpool_init(nthreads);                                // when initiated for >0 threads
    hts_set_opt(in_fp,  HTS_OPT_THREAD_POOL, &thread_pool);                     // and bound to the input
    hts_set_opt(out_fp, HTS_OPT_THREAD_POOL, &thread_pool);                     // and output file pointers
  }
  
  bam_hdr_t *in_hdr = sam_hdr_read(in_fp);                                      // try read input file header
  if (in_hdr==NULL) Rcpp::stop("Unable to read input BAM header");              // fall back if error
  if (sam_hdr_write(out_fp, in_hdr) < 0) Rcpp::stop("Unable to write header");  // try write output file header
  
  // vars
  int nrecs = 0, ncalled = 0;                                                   // counters: BAM records, records with methylation called
  int max_query_width = 1024;                                                   // max query width, expanded if necessary
  int query_width = max_query_width;                                            // query width
  char *rs = (char*) malloc((max_query_width+4) * sizeof(char));                // sequence of the reference plus 2x2nt on sides
  char *xm = (char*) malloc(max_query_width * sizeof(char));                    // XM array
  
  // compare (+remap?) reference sequences in the genome and in the BAM header
  // I don't do remap atm, let's see if order of reference sequences in genome
  // and BAM is ever different...
  for (int i=0; i<in_hdr->n_targets; i++) {                                     // for all BAM header refseqs
    if ((in_hdr->target_len[i] != rlen[i]) ||                                   // if any of length
        (rname[i].compare(in_hdr->target_name[i])!=0))                          // or name are different
      Rcpp::stop("BAM reference sequence doesn't match the provided genome sequence"); // freak out
  }
  
  bam1_t *in_rec = bam_init1();                                                 // create BAM alignment structure
  while(sam_read1(in_fp, in_hdr, in_rec) > 0) {                                 // read rec by rec
    nrecs++;                                                                    // BAM alignment records ++
    if ((nrecs & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();                     // every ~1M reads check for the interrupt
    
    char *record_strand = (char*) bam_aux_get(in_rec, tag.c_str());             // genome strand ("ZCT" or "ZGA")
    char *record_xm = (char*) bam_aux_get(in_rec, "XM");                        // methylation string (XM)
    if ((in_rec->core.flag & BAM_FUNMAP) ||                                     // if unmapped
        (record_strand==NULL) ||                                                // or genome strand is unknown
        (record_xm!=NULL))                                                      // or XM is already present
      goto writeout;                                                            // don't do anything, just write out
    
    {
      query_width = abs(in_rec->core.l_qseq);                                   // query width
      if (query_width > max_query_width) {                                      // if sequence is longer than XM holder
        max_query_width = query_width;                                          // new max
        rs = (char *) realloc(rs, (max_query_width+4) * sizeof(char));          // expand rs holder
        xm = (char *) realloc(xm, max_query_width * sizeof(char));              // expand xm holder
        if ((rs==NULL) || (xm==NULL)) Rcpp::stop("No memory for BAM record #%i", nrecs); // check memory allocation
      }

      // apply CIGAR to reference seq (convert from reference to query space)
      const char *refseq = rseq->at(in_rec->core.tid).c_str() + in_rec->core.pos; // reference sequence
      uint32_t n_cigar = in_rec->core.n_cigar;                                  // number of CIGAR operations
      uint32_t *record_cigar = bam_get_cigar(in_rec);                           // CIGAR array
      uint32_t ref_pos = 0;                                                     // starting position in reference array
      uint32_t dest_pos = 2;                                                    // starting position in destination (query space) array
      for (size_t i=0; i<n_cigar; i++) {                                        // op by op
        uint32_t cigar_op = bam_cigar_op(record_cigar[i]);                      // CIGAR operation
        uint32_t cigar_oplen = bam_cigar_oplen(record_cigar[i]);                // CIGAR operation length
        switch(cigar_op) {
        case BAM_CMATCH :                                                       // 'M', 0, consumes both query and reference
        case BAM_CEQUAL :                                                       // '=', 7, consumes both query and reference
          memcpy(rs+dest_pos, refseq+ref_pos, cigar_oplen);                     // will be turned into a context later
          ref_pos += cigar_oplen;
          dest_pos += cigar_oplen;
          break;
        case BAM_CDIFF :                                                        // 'X', 8, consumes both query and reference
          memset(rs+dest_pos, 'N', cigar_oplen);                                // we don't actually know what was inserted
          ref_pos += cigar_oplen;
          dest_pos += cigar_oplen;
          break;
        case BAM_CINS :                                                         // 'I', 1, consumes query only
        case BAM_CSOFT_CLIP :                                                   // 'S', 4, consumes query only
          memset(rs+dest_pos, 'N', cigar_oplen);                                // we don't actually know what was inserted
          dest_pos += cigar_oplen;
          break;
        case BAM_CDEL :                                                         // 'D', 2, consumes reference only
        case BAM_CREF_SKIP :                                                    // 'N', 3, consumes reference only
          ref_pos += cigar_oplen;
          break;
        case BAM_CHARD_CLIP :                                                   // 'H', 5
        case BAM_CPAD :                                                         // 'P', 6
        case BAM_CBACK :
          break;
        default :
          Rcpp::stop("Unknown CIGAR operation for BAM entry %s",                // unknown CIGAR operation
                     bam_get_qname(in_rec));
        }
      }
      rs[0] = in_rec->core.pos>=2 ? refseq[-2] : 'N';                           // -2 base of reference sequence in front of query
      rs[1] = in_rec->core.pos>=1 ? refseq[-1] : 'N';                           // -1 base of reference sequence in front of query
      int bases_left = rseq->at(in_rec->core.tid).size() - in_rec->core.pos - ref_pos; // bases from the last reference position till its end
      rs[query_width+2] = bases_left >= 1 ? refseq[ref_pos+0] : 'N';            // 1st next base of reference sequence beyond query
      rs[query_width+3] = bases_left >= 2 ? refseq[ref_pos+1] : 'N';            // 2st next base of reference sequence beyond query
      
      uint8_t *record_pseq = bam_get_seq(in_rec);                               // packed sequence string (4 bit per base)
      int context_shift = (record_strand[1] == 'C') ? 2 : 0;                    // start making genomic context from 2nd base for fwd strand, and 0th base for reverse
      const unsigned char* context_map = context_shift ? triad_forward_context : triad_reverse_context; // lookup table to use
      for (int i=0; i<query_width; i++) {
        unsigned int idx = ((rs[i+context_shift] &7) << 6) |                    // last 3 bits of the first base << 6
                           ((rs[i+1+context_shift] &7) << 3) |                  // last 3 bits of the second base << 3
                           (rs[i+2+context_shift] &7);                          // last 3 bits of the third base
        xm[i] = context_map[idx];                                               // look up context
        
        // *** the actual methylation calling starts here ***
        if (xm[i]!='.') {                                                       // if it's a hxz
          if (seq_nt16_str[bam_seqi(record_pseq,i)]==record_strand[1]) {        // if query base is the first char of genome conversion: C,G
            xm[i] &= 0b11011111;                                                // uppercase the context char
          } else if (seq_nt16_str[bam_seqi(record_pseq,i)]!=record_strand[2]) { // if query base is NOT the second char of genome conversion: T,A
            xm[i]='.';                                                          // clear the context
          }                                                                     // leave lowercase otherwise
        }
        // *** the actual methylation calling ends here ***
        
      }
      
      // if ((record_xm!=NULL) && (memcmp(xm, record_xm+1, query_width)!=0)) {     // check differences with available XM
      //   Rcpp::Rcout << nrecs << ": S|" << std::string(rs, query_width+4) << std::endl;
      //   Rcpp::Rcout << nrecs << ": R|  " << std::string(record_xm+1, query_width) << std::endl;
      //   Rcpp::Rcout << nrecs << ": C|  " << std::string(xm, query_width) << std::endl;
      // }
      
      // if (record_xm==NULL)
      bam_aux_update_str(in_rec, "XM", query_width, xm);                        // since XM tag is absent, add it
      ncalled++;                                                                // successfully called
    }
    
    writeout:
    if ((tag!="XG") && (record_strand!=NULL))                                   // if XG tag is absent but genome strand is known
      bam_aux_append(in_rec, "XG", 'Z', 2, (const uint8_t*)record_strand+1);    // add it
    if (sam_write1(out_fp, in_hdr, in_rec) < 0) Rcpp::stop("Unable to write BAM"); // write record
  }

  // cleaning
  bam_destroy1(in_rec);                                                         // clean BAM alignment structure 
  hts_close(in_fp);                                                             // close input BAM file
  hts_close(out_fp);                                                            // close output BAM file
  if (thread_pool.pool) hts_tpool_destroy(thread_pool.pool);                    // free thread pool
  free(rs);                                                                     // and free manually allocated memory
  free(xm);
  
  // wrap and return the results
  Rcpp::List res = Rcpp::List::create(                                          // final List
    Rcpp::Named("nrecs") = nrecs,                                               // number of BAM records
    Rcpp::Named("ncalled") = ncalled                                            // number of XM tags successfully called
  );
  
  return(res);
}


// #############################################################################
// test code and sourcing don't work on OS X
/*** R
setwd("~/work/packages/epialleleR/")
devtools::document()
devtools::load_all()
system.time(genome <- rcpp_read_genome("/scratch/ref/DRAGEN/hg38_plus_lambda_ChrY_PAR_masked.fa.gz", 1))
system.time(genome <- rcpp_read_genome("/Users/oleksii.nikolaienko/work/data/hg38/hg38_plus_lambda_ChrY_PAR_masked.fa.gz", 1))
inout <- c(system.file("extdata", "test", "paired-name-xm.bam", package="epialleleR"), "/tmp/out.bam")
rcpp_call_methylation_genome(inout[1], inout[2], genome, "XG", 8)
sapply(inout, file.size)
*/
// #############################################################################