#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include "epialleleR.h"

// [[Rcpp::depends(Rhtslib)]]

// An optimised attempt to read and preprocess BAM in place. To do:
// [ ] OpenMP SIMD?
// [+] HTSlib threads
// [+] rec_seq_rs and rec_xm_rs as char*
// [?] reverse QNAME
// [ ] free resources on interrupt

// SHORT-READ PAIRED-END BAM

// [[Rcpp::export]]
Rcpp::DataFrame rcpp_read_bam_paired (std::string fn,                           // file name
                                      int min_mapq,                             // min read mapping quality
                                      int min_baseq,                            // min base quality
                                      bool skip_duplicates,                     // skip marked duplicates
                                      int trim5,                                // trim bases from 5'
                                      int trim3,                                // trim bases from 3'
                                      int nthreads)                             // HTSlib threads, >0 for multiple
{
  // constants
  int max_qname_width = 1024;                                                   // max QNAME length, not expanded yet, ever error-prone?
  int max_templ_width = 8192;                                                   // max insert size, expanded if necessary
  min_baseq = min_baseq - (min_baseq>0);                                        // decrease base quality by one to include bases with QUAL==min_baseq
  
  // file IO
  htsFile *bam_fp = hts_open(fn.c_str(), "r");                                  // try open file
  if (!bam_fp) Rcpp::stop("Unable to open BAM file for reading");               // fall back if error
  htsThreadPool thread_pool = {NULL, 0};                                        // thread pool cuts time by 30%
  if (nthreads>0) {
    thread_pool.pool = hts_tpool_init(nthreads);                                // when initiated for >0 threads
    hts_set_opt(bam_fp, HTS_OPT_THREAD_POOL, &thread_pool);                     // and bound to the file pointer
  }
  bam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);                                    // try read file header
  if (!bam_hdr) Rcpp::stop("Unable to read BAM header");                        // fall back if error  
  bam1_t *bam_rec = bam_init1();                                                // create BAM alignment structure
  
  // main containers
  std::vector<std::string>* seq = new std::vector<std::string>;                 // SEQ
  std::vector<std::string>* xm = new std::vector<std::string>;                  // XM
  std::vector<int> rname, strand, start;                                        // id for RNAME, id for CT==1/GA==2, POS
  int nrecs = 0, ntempls = 0;                                                   // counters: BAM records, templates (consecutive proper read pairs)
  
  // reserve some memory
  rname.reserve(0xFFFFF); strand.reserve(0xFFFFF); start.reserve(0xFFFFF); 
  seq->reserve(0xFFFFF); xm->reserve(0xFFFFF);
  
  // template holders
  char *templ_qname = (char*) malloc(max_qname_width * sizeof(char));           // template QNAME
  uint8_t *templ_qual_rs = (uint8_t*) malloc(max_templ_width * sizeof(uint8_t));// template QUAL array
  uint8_t *templ_seq_rs  = (uint8_t*) malloc(max_templ_width * sizeof(uint8_t));// template SEQ array
  uint8_t *templ_xm_rs   = (uint8_t*) malloc(max_templ_width * sizeof(uint8_t));// template XM array
  int templ_rname = 0, templ_start = 0, templ_strand = 0, templ_width = 0;      // template RNAME, POS, STRAND, ISIZE
  
  #define push_template {                                       /* pushing template data to vectors */ \
    rname.push_back(templ_rname + 1);                                                    /* RNAME+1 */ \
    strand.push_back(templ_strand);                                                       /* STRAND */ \
    start.push_back(templ_start + trim5 + 1);                                              /* POS+1 */ \
    seq->emplace_back((const char*) templ_seq_rs + trim5, templ_width - (trim5+trim3));      /* SEQ */ \
    xm->emplace_back( (const char*) templ_xm_rs  + trim5, templ_width - (trim5+trim3));       /* XM */ \
    ntempls++;                                                                                /* +1 */ \
  }
  
  // process alignments
  while( sam_read1(bam_fp, bam_hdr, bam_rec) > 0 ) {                            // rec by rec
    nrecs++;                                                                    // BAM alignment records ++
    if ((nrecs & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();                     // every ~1M reads check for the interrupt
    
    if ((bam_rec->core.qual < min_mapq) ||                                      // skip if mapping quality < min.mapq
        (!(bam_rec->core.flag & BAM_FPROPER_PAIR)) ||                           // or if not a proper pair
        (skip_duplicates && (bam_rec->core.flag & BAM_FDUP))) continue;         // or if record is an optical/PCR duplicate
    
    char *rec_strand = (char*) bam_aux_get(bam_rec, "XG");                      // genome strand
    char *rec_xm = (char*) bam_aux_get(bam_rec, "XM");                          // methylation string
    if (!rec_strand || !rec_xm) continue;                                       // skip if no XM/XG tags (no methylation info available)
    
    // check if not the same template (QNAME)
    if ((strcmp(templ_qname, bam_get_qname(bam_rec)) != 0)) {                
      // store previous template if it's a valid record
      if (templ_strand!=0) push_template;                                       // templ_strand is 0 for empty records (very start of BAM and/or when invalid records are in front of BAM)
      
      // initialize new template
      strcpy(templ_qname, bam_get_qname(bam_rec));                              // store template QNAME
      templ_rname = bam_rec->core.tid;                                          // store template RNAME
      templ_start = bam_rec->core.pos < bam_rec->core.mpos ?                    // smallest of POS,MPOS is a start
        bam_rec->core.pos : bam_rec->core.mpos;
      templ_width = abs(bam_rec->core.isize);                                   // template ISIZE
      templ_strand = ( rec_strand[1] == 'C' ) ? 1 : 2 ;                         // STRAND is 1 if "ZCT"/"+", 2 if "ZGA"/"-"
      
      // resize containers if necessary
      if (templ_width > max_templ_width) {
        max_templ_width = templ_width;                                          // expand template holders
        templ_qual_rs = (uint8_t *) realloc(templ_qual_rs, max_templ_width);
        templ_seq_rs  = (uint8_t *) realloc(templ_seq_rs,  max_templ_width);
        templ_xm_rs   = (uint8_t *) realloc(templ_xm_rs,   max_templ_width);
        if (!templ_qual_rs || !templ_seq_rs || !templ_xm_rs) Rcpp::stop("Unable to allocate memory for BAM record #%i", nrecs); // check memory allocation
      }
      std::memset(templ_qual_rs, (uint8_t) min_baseq, templ_width);             // clean template holders
      std::memset(templ_seq_rs, 'N', templ_width);
      std::memset(templ_xm_rs,  '-', templ_width);
      // std::fill_n(templ_qual_rs, templ_width, (uint8_t) min_baseq);             // clean template holders
      // std::fill_n(templ_seq_rs,  templ_width, 'N');
      // std::fill_n(templ_xm_rs,   templ_width, '-');
     }
    
    // add another read to the template
    // source containers
    uint8_t *rec_qual = bam_get_qual(bam_rec);                                  // quality string (Phred scale with no +33 offset)
    rec_xm++;                                                                   // remove leading 'Z' from XM string
    uint8_t *rec_pseq = bam_get_seq(bam_rec);                                   // packed sequence string (4 bit per base)
    
    // apply CIGAR
    uint32_t n_cigar = bam_rec->core.n_cigar;                                   // number of CIGAR operations
    uint32_t *rec_cigar = bam_get_cigar(bam_rec);                               // CIGAR array
    uint32_t query_pos = 0;                                                     // starting position in query array
    uint32_t dest_pos = bam_rec->core.pos - templ_start;                        // starting position in destination array
    for (size_t i=0; i<n_cigar; i++) {                                          // op by op
      uint32_t cigar_op = bam_cigar_op(rec_cigar[i]);                           // CIGAR operation
      uint32_t cigar_oplen = bam_cigar_oplen(rec_cigar[i]);                     // CIGAR operation length
      switch(cigar_op) {
        case BAM_CMATCH :                                                       // 'M', 0
        case BAM_CEQUAL :                                                       // '=', 7
        case BAM_CDIFF :                                                        // 'X', 8
          for (size_t j=0; j<cigar_oplen; j++) {
            if (rec_qual[query_pos+j] > templ_qual_rs[dest_pos+j]) {
              templ_qual_rs[dest_pos+j] = rec_qual[query_pos+j];
              templ_seq_rs[dest_pos+j] = seq_nt16_str[bam_seqi(rec_pseq,query_pos+j)];
              templ_xm_rs[dest_pos+j] = rec_xm[query_pos+j];
            }
          }
          query_pos += cigar_oplen;
          dest_pos += cigar_oplen;
          break;
        case BAM_CINS :                                                         // 'I', 1
        case BAM_CSOFT_CLIP :                                                   // 'S', 4
          query_pos += cigar_oplen;
          break;
        case BAM_CDEL :                                                         // 'D', 2
        case BAM_CREF_SKIP :                                                    // 'N', 3
          dest_pos += cigar_oplen;
          break;
        case BAM_CHARD_CLIP :                                                   // 'H', 5
        case BAM_CPAD :                                                         // 'P', 6
        case BAM_CBACK :
          break;
        default :
          Rcpp::stop("Unknown CIGAR operation for BAM entry %s",                // unknown CIGAR operation
                     bam_get_qname(bam_rec));
      }
    }
  }
  
  // push last, yet unsaved template (no empty files enter this function)
  push_template;
  
  // cleaning
  bam_destroy1(bam_rec);                                                        // clean BAM alignment structure 
  hts_close(bam_fp);                                                            // close BAM file
  if (thread_pool.pool) hts_tpool_destroy(thread_pool.pool);                    // free thread pool
  free(templ_qname);                                                            // and free manually allocated memory
  free(templ_qual_rs);
  free(templ_seq_rs);
  free(templ_xm_rs);
  
  // wrap and return the results
  Rcpp::DataFrame res = Rcpp::DataFrame::create(                                // final DF
    Rcpp::Named("rname") = rname,                                               // numeric ids (factor) for reference names
    Rcpp::Named("strand") = strand,                                             // numeric ids (factor) for reference strands
    Rcpp::Named("start") = start                                                // start positions of reads
  );
  
  // factor levels
  std::vector<std::string> chromosomes (                                        // vector of reference names
      bam_hdr->target_name, bam_hdr->target_name + bam_hdr->n_targets);
  std::vector<std::string> strands = {"+", "-"};
  
  Rcpp::IntegerVector col_rname = res["rname"];                                 // make rname a factor
  col_rname.attr("class") = "factor";
  col_rname.attr("levels") = chromosomes;
  
  Rcpp::IntegerVector col_strand = res["strand"];                               // make strand a factor
  col_strand.attr("class") = "factor";
  col_strand.attr("levels") = strands;
  
  Rcpp::XPtr<std::vector<std::string>> seq_xptr(seq, true);
  res.attr("seq_xptr") = seq_xptr;                                              // external pointer to sequences
  Rcpp::XPtr<std::vector<std::string>> xm_xptr(xm, true);
  res.attr("xm_xptr") = xm_xptr;                                                // external pointer to methylation strings
  
  res.attr("nrecs") = nrecs;                                                    // number of records in BAM file
  res.attr("npushed") = ntempls;                                                // number of templates pushed to data.frame
  
  return(res);
}

// #############################################################################

// SHORT-READ SINGLE-END BAM

// [[Rcpp::export]]
Rcpp::DataFrame rcpp_read_bam_single (std::string fn,                           // file name
                                      int min_mapq,                             // min read mapping quality
                                      int min_baseq,                            // min base quality
                                      bool skip_duplicates,                     // skip marked duplicates
                                      int trim5,                                // trim bases from 5'
                                      int trim3,                                // trim bases from 3'
                                      int nthreads)                             // HTSlib threads, >0 for multiple
{
  // constants
  int max_record_width  = 1024;                                                 // max record width, expanded if necessary
  
  // file IO
  htsFile *bam_fp = hts_open(fn.c_str(), "r");                                  // try open file
  if (!bam_fp) Rcpp::stop("Unable to open BAM file for reading");               // fall back if error
  htsThreadPool thread_pool = {NULL, 0};                                        // thread pool cuts time by 30%
  if (nthreads>0) {
    thread_pool.pool = hts_tpool_init(nthreads);                                // when initiated for >0 threads
    hts_set_opt(bam_fp, HTS_OPT_THREAD_POOL, &thread_pool);                     // and bound to the file pointer
  }
  bam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);                                    // try read file header
  if (!bam_hdr) Rcpp::stop("Unable to read BAM header");                        // fall back if error  
  bam1_t *bam_rec = bam_init1();                                                // create BAM alignment structure
  
  // main containers
  std::vector<std::string>* seq = new std::vector<std::string>;                 // SEQ
  std::vector<std::string>* xm = new std::vector<std::string>;                  // XM
  std::vector<int> rname, strand, start;                                        // id for RNAME, id for CT==1/GA==2, POS
  int nrecs = 0, npushed = 0;                                                   // counters: BAM records read, BAM records pushed to data.table
  
  // reserve some memory
  rname.reserve(0xFFFFF); strand.reserve(0xFFFFF); start.reserve(0xFFFFF); 
  seq->reserve(0xFFFFF); xm->reserve(0xFFFFF);
  
  // read holders
  int record_width = max_record_width;                                          // record ISIZE/TLEN
  uint8_t *record_seq_rs  = (uint8_t*) malloc(record_width * sizeof(uint8_t));  // record SEQ array
  uint8_t *record_xm_rs   = (uint8_t*) malloc(record_width * sizeof(uint8_t));  // record XM array
  
  // process alignments
  while( sam_read1(bam_fp, bam_hdr, bam_rec) > 0 ) {                            // rec by rec
    nrecs++;                                                                    // BAM alignment records ++
    if ((nrecs & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();                     // every ~1M reads check for the interrupt

    if ((bam_rec->core.qual < min_mapq) ||                                      // skip if mapping quality < min.mapq
        (skip_duplicates && (bam_rec->core.flag & BAM_FDUP))) continue;         // or if record is an optical/PCR duplicate
    
    char *record_strand = (char*) bam_aux_get(bam_rec, "XG");                   // genome strand
    char *record_xm = (char*) bam_aux_get(bam_rec, "XM");                       // methylation string
    if (!record_strand || !record_xm) continue;                                 // skip if no XM/XG tags (no methylation info available)
    
    // get record sequence, XM, quality string
    uint8_t *record_qual = bam_get_qual(bam_rec);                               // quality string (Phred scale with no +33 offset)
    record_xm++;                                                                // remove leading 'Z' from XM string
    uint8_t *record_pseq = bam_get_seq(bam_rec);                                // packed sequence string (4 bit per base)
    
    // get CIGAR
    uint32_t n_cigar = bam_rec->core.n_cigar;                                   // number of CIGAR operations
    uint32_t *record_cigar = bam_get_cigar(bam_rec);                            // CIGAR array
    record_width = bam_cigar2rlen(n_cigar, record_cigar);                       // reference length for the current query
    
    // resize containers if necessary
    if (record_width > max_record_width) {
      max_record_width = record_width;                                          // expand template holders
      record_seq_rs  = (uint8_t *) realloc(record_seq_rs, record_width * sizeof(uint8_t));
      record_xm_rs   = (uint8_t *) realloc(record_xm_rs,  record_width * sizeof(uint8_t));
      if (!record_seq_rs || !record_xm_rs) Rcpp::stop("Unable to allocate memory for BAM record #%i", nrecs); // check memory allocation
    }
    
    // prepare for the new record
    std::memset(record_seq_rs, 'N', record_width);
    std::memset(record_xm_rs,  '-', record_width);
    
    // apply CIGAR
    uint32_t query_pos = 0;                                                     // starting position in query array
    uint32_t dest_pos = 0;                                                      // starting position in destination array
    for (size_t i=0; i<n_cigar; i++) {                                          // op by op
      uint32_t cigar_op = bam_cigar_op(record_cigar[i]);                        // CIGAR operation
      uint32_t cigar_oplen = bam_cigar_oplen(record_cigar[i]);                  // CIGAR operation length
      switch(cigar_op) {
      case BAM_CMATCH :                                                         // 'M', 0
      case BAM_CEQUAL :                                                         // '=', 7
      case BAM_CDIFF :                                                          // 'X', 8
        for (size_t j=0; j<cigar_oplen; j++) {
          if (record_qual[query_pos+j] >= min_baseq) {
            record_seq_rs[dest_pos+j] = seq_nt16_str[bam_seqi(record_pseq,query_pos+j)];
            record_xm_rs[dest_pos+j] = record_xm[query_pos+j];
          }
        }
        query_pos += cigar_oplen;
        dest_pos += cigar_oplen;
        break;
      case BAM_CINS :                                                           // 'I', 1
      case BAM_CSOFT_CLIP :                                                     // 'S', 4
        query_pos += cigar_oplen;
        break;
      case BAM_CDEL :                                                           // 'D', 2
      case BAM_CREF_SKIP :                                                      // 'N', 3
        dest_pos += cigar_oplen;
        break;
      case BAM_CHARD_CLIP :                                                     // 'H', 5
      case BAM_CPAD :                                                           // 'P', 6
      case BAM_CBACK :
        break;
      default :
        Rcpp::stop("Unknown CIGAR operation for BAM entry %s",                  // unknown CIGAR operation
                   bam_get_qname(bam_rec));
      }
    }
    
    // pushing record data to vectors
    rname.push_back(bam_rec->core.tid + 1);                                     // RNAME+1 
    strand.push_back(( record_strand[1] == 'C' ) ? 1 : 2);                      // STRAND is 1 if "ZCT"/"+", 2 if "ZGA"/"-"
    start.push_back(bam_rec->core.pos + trim5 +1);                              // POS+1 
    seq->emplace_back((const char*) record_seq_rs + trim5, dest_pos - (trim5+trim3)); // SEQ 
    xm->emplace_back( (const char*) record_xm_rs + trim5,  dest_pos - (trim5+trim3)); // XM 
    npushed++;                                                                  // +1 
  }
  
  // cleaning
  bam_destroy1(bam_rec);                                                        // clean BAM alignment structure 
  hts_close(bam_fp);                                                            // close BAM file
  if (thread_pool.pool) hts_tpool_destroy(thread_pool.pool);                    // free thread pool
  free(record_seq_rs);                                                          // and free manually allocated memory
  free(record_xm_rs);
  
  // wrap and return the results
  Rcpp::DataFrame res = Rcpp::DataFrame::create(                                // final DF
    Rcpp::Named("rname") = rname,                                               // numeric ids (factor) for reference names
    Rcpp::Named("strand") = strand,                                             // numeric ids (factor) for reference strands
    Rcpp::Named("start") = start                                                // start positions of reads
  );
  
  // factor levels
  std::vector<std::string> chromosomes (                                        // vector of reference names
      bam_hdr->target_name, bam_hdr->target_name + bam_hdr->n_targets);
  std::vector<std::string> strands = {"+", "-"};
  
  Rcpp::IntegerVector col_rname = res["rname"];                                 // make rname a factor
  col_rname.attr("class") = "factor";
  col_rname.attr("levels") = chromosomes;
  
  Rcpp::IntegerVector col_strand = res["strand"];                               // make strand a factor
  col_strand.attr("class") = "factor";
  col_strand.attr("levels") = strands;
  
  Rcpp::XPtr<std::vector<std::string>> seq_xptr(seq, true);
  res.attr("seq_xptr") = seq_xptr;                                              // external pointer to sequences
  Rcpp::XPtr<std::vector<std::string>> xm_xptr(xm, true);
  res.attr("xm_xptr") = xm_xptr;                                                // external pointer to methylation strings
  
  res.attr("nrecs") = nrecs;                                                    // number of records in BAM file
  res.attr("npushed") = npushed;                                                // number of records pushed to data.frame
  
  return(res);
}


// #############################################################################

// LONG-READ SINGLE-END BAM
// https://samtools.github.io/hts-specs/SAMv1.pdf
// https://samtools.github.io/hts-specs/SAMtags.pdf

// EXPLAIN THE LOGIC IN DOCS - NAMELY, DIFF BETWEEN SHORT-READ (BISULFITE) AND
// LONG-READ (NATIVE) SEQUENCING METHYLATION CALLING (GENOME IS ABSOLUTELY
// NECESSARY FOR THE FIRST BUT NOT THE SECOND)

// [[Rcpp::export]]
Rcpp::DataFrame rcpp_read_bam_mm (std::string fn,                               // file name
                                  int min_mapq,                                 // min read mapping quality
                                  int min_baseq,                                // min base quality
                                  int min_prob,                                 // min probability of 5mC modification
                                  bool highest_prob,                            // consider only if 5mC probability is the highest of all mods at particular pos
                                  bool skip_duplicates,                         // skip marked duplicates
                                  int trim5,                                    // trim bases from 5'
                                  int trim3,                                    // trim bases from 3'
                                  int nthreads)                                 // HTSlib threads, >0 for multiple
{
  // constants
  int max_query_width   = 1024;                                                 // max NON-refspaced query width, expanded if necessary
  int max_record_width  = 1024;                                                 // max refspaced record width, expanded if necessary
  const int max_nmods   = 16;                                                   // allow MAX 16 modifications per base

  // file IO
  htsFile *bam_fp = hts_open(fn.c_str(), "r");                                  // try open file
  if (!bam_fp) Rcpp::stop("Unable to open BAM file for reading");               // fall back if error
  htsThreadPool thread_pool = {NULL, 0};                                        // thread pool cuts time by 30%
  if (nthreads>0) {
    thread_pool.pool = hts_tpool_init(nthreads);                                // when initiated for >0 threads
    hts_set_opt(bam_fp, HTS_OPT_THREAD_POOL, &thread_pool);                     // and bound to the file pointer
  }
  bam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);                                    // try read file header
  if (!bam_hdr) Rcpp::stop("Unable to read BAM header");                        // fall back if error
  bam1_t *bam_rec = bam_init1();                                                // create BAM alignment structure
  
  // base modifications
  hts_base_mod_state *mod_state = hts_base_mod_state_alloc();                   // allocate space for base modification states
  hts_base_mod base_mods[max_nmods];                                            // allocate an array of MAX 16 possible modifications per base
  int mod_pos = 0, nmods = 0;                                                   // position of modified base in the query, number of modifications at that base

  // main containers
  std::vector<std::string>* seq = new std::vector<std::string>;                 // SEQ
  std::vector<std::string>* xm = new std::vector<std::string>;                  // XM
  std::vector<int> rname, strand, start;                                        // id for RNAME, id for CT==1/GA==2, POS
  int nrecs = 0, npushed = 0;                                                   // counters: BAM records read, BAM records pushed to data.table

  // reserve some memory
  rname.reserve(0xFFFFF); strand.reserve(0xFFFFF); start.reserve(0xFFFFF);
  seq->reserve(0xFFFFF); xm->reserve(0xFFFFF);

  // read holders
  int query_width = max_query_width;                                            // NON-refspaced query length
  uint8_t *record_seq = (uint8_t*) malloc((query_width+4) * sizeof(uint8_t));   // NON-refspaced query SEQ array, plus NN at the end
  uint8_t *record_xm  = (uint8_t*) malloc( query_width    * sizeof(uint8_t));   // NON-refspaced query XM array
  int record_width = max_record_width;                                          // refspaced record ISIZE/TLEN
  uint8_t *record_seq_rs  = (uint8_t*) malloc(record_width * sizeof(uint8_t));  // refspaced record SEQ array
  uint8_t *record_xm_rs   = (uint8_t*) malloc(record_width * sizeof(uint8_t));  // refspaced record XM array

  // process alignments
  while( sam_read1(bam_fp, bam_hdr, bam_rec) > 0 ) {                            // rec by rec
    nrecs++;                                                                    // BAM alignment records ++
    if ((nrecs & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();                     // every ~1M reads check for the interrupt

    if ((bam_rec->core.qual < min_mapq) ||                                      // skip if mapping quality < min.mapq
        (skip_duplicates && (bam_rec->core.flag & BAM_FDUP))) continue;         // or if record is an optical/PCR duplicate

    int record_strand = ! ( bam_rec->core.flag & BAM_FREVERSE );                // genome strand, 0 if forward, 1 if reverse

    // get record sequence, quality string
    uint8_t *record_qual = bam_get_qual(bam_rec);                               // quality string (Phred scale with no +33 offset)
    uint8_t *record_pseq = bam_get_seq(bam_rec);                                // packed sequence string (4 bit per base)

    // get CIGAR
    uint32_t n_cigar = bam_rec->core.n_cigar;                                   // number of CIGAR operations
    uint32_t *record_cigar = bam_get_cigar(bam_rec);                            // CIGAR array
    query_width = abs(bam_rec->core.l_qseq);                                    // NON-refspaced query width
    record_width = bam_cigar2rlen(n_cigar, record_cigar);                       // reference length for the current query (refspaced)

    // resize containers if necessary
    if (query_width > max_query_width) {
      max_query_width = query_width;                                            // expand template holders
      record_seq = (uint8_t *) realloc(record_seq, (query_width+4) * sizeof(uint8_t));
      record_xm  = (uint8_t *) realloc(record_xm,   query_width    * sizeof(uint8_t));
      if (!record_seq || !record_xm) Rcpp::stop("Unable to allocate memory for BAM record #%i", nrecs); // check memory allocation

    }
    if (record_width > max_record_width) {
      max_record_width = record_width;                                          // expand template holders
      record_seq_rs = (uint8_t *) realloc(record_seq_rs, record_width * sizeof(uint8_t));
      record_xm_rs  = (uint8_t *) realloc(record_xm_rs,  record_width * sizeof(uint8_t));
      if (!record_seq_rs || !record_xm_rs) Rcpp::stop("Unable to allocate memory for BAM record #%i", nrecs); // check memory allocation
    }

    // prepare for the new record
    std::memset(record_seq_rs, 'N', record_width);
    std::memset(record_xm_rs,  '-', record_width);
    
    // unpack the sequence string
    for (int i=0; i<query_width; i++) {
      record_seq[i+2] = seq_nt16_str[bam_seqi(record_pseq,i)];
    }
    std::memset(record_seq, 'N', 2);
    std::memset(record_seq+query_width, 'N', 2);

    // create a non-refspaced context string
    int context_shift = record_strand ? 0 : 2;                                  // start making genomic context from 2nd base for fwd strand, and 0th base for reverse
    const unsigned char* context_map = context_shift ? triad_forward_context : triad_reverse_context; // lookup table to use
    for (int i=0; i<query_width; i++) {
      record_xm[i] = triad_to_ctx((record_seq+context_shift+i), context_map);   // look up context
    }
    
    // parse base modifications: any location not reported is implicitly
    // assumed to contain no modification
    bam_parse_basemod(bam_rec, mod_state);                                      // fill the mod_state structure
    while ((nmods = bam_next_basemod(bam_rec, mod_state, base_mods, max_nmods, &mod_pos)) > 0) { // cycle through modified bases
      if (record_strand) mod_pos = query_width - mod_pos - 1;                   // if mapped to reverse strand, then position is on revcomplemented query
      Rcpp::Rcout << record_xm[mod_pos];
      if (record_xm[mod_pos]<'A') continue;                                     // skip if this position is not a 'hxz'
      unsigned int meth_prob = 0, max_other_prob = 0;                           // scaled probabilities: methylation mod, maximum of any other mods of a current base
      for (int i=0; i<nmods; i++) {                                             // cycle through all mods of a current base
        if (base_mods[i].modified_base=='m' || base_mods[i].modified_base==-27551) { // if it's a 5mC
          meth_prob = (unsigned int) base_mods[i].qual;                         // -1 (unknown probability) becomes MAX_INT (highest probability)
        } else if (max_other_prob < (unsigned int) base_mods[i].qual) {         // if NOT a 5mC and probability is higher than max_other_prob
            max_other_prob = (unsigned int) base_mods[i].qual;                  // record the highest probability of other modifications
        }
      }
      if (meth_prob>0 &&                                                        // if there is a 5mC modification
          meth_prob>=min_prob &&                                                // and its probability is not less than min_prob
          (!highest_prob || meth_prob>max_other_prob)) {                        // and its probability is either highest or highest_prob==FALSE
        record_xm[mod_pos] &= 0b11011111;                                       // uppercase the context char
      }
    }
    Rcpp::Rcout << "\n";

    // apply CIGAR
    uint32_t query_pos = 0;                                                     // starting position in query array
    uint32_t dest_pos = 0;                                                      // starting position in destination array
    for (size_t i=0; i<n_cigar; i++) {                                          // op by op
      uint32_t cigar_op = bam_cigar_op(record_cigar[i]);                        // CIGAR operation
      uint32_t cigar_oplen = bam_cigar_oplen(record_cigar[i]);                  // CIGAR operation length
      switch(cigar_op) {
      case BAM_CMATCH :                                                         // 'M', 0
      case BAM_CEQUAL :                                                         // '=', 7
      case BAM_CDIFF :                                                          // 'X', 8
        for (size_t j=0; j<cigar_oplen; j++) {
          if (record_qual[query_pos+j] >= min_baseq) {
            record_seq_rs[dest_pos+j] = record_seq[query_pos+2+j];
            record_xm_rs[dest_pos+j] = record_xm[query_pos+j];
          }
        }
        query_pos += cigar_oplen;
        dest_pos += cigar_oplen;
        break;
      case BAM_CINS :                                                           // 'I', 1
      case BAM_CSOFT_CLIP :                                                     // 'S', 4
        query_pos += cigar_oplen;
        break;
      case BAM_CDEL :                                                           // 'D', 2
      case BAM_CREF_SKIP :                                                      // 'N', 3
        dest_pos += cigar_oplen;
        break;
      case BAM_CHARD_CLIP :                                                     // 'H', 5
      case BAM_CPAD :                                                           // 'P', 6
      case BAM_CBACK :
        break;
      default :
        Rcpp::stop("Unknown CIGAR operation for BAM entry %s",                  // unknown CIGAR operation
                   bam_get_qname(bam_rec));
      }
    }

    // pushing record data to vectors
    rname.push_back(bam_rec->core.tid + 1);                                     // RNAME+1
    strand.push_back(record_strand + 1);                                        // STRAND is 1 if "CT"/"+", 2 if "GA"/"-"
    start.push_back(bam_rec->core.pos + trim5 + 1);                             // POS+1
    seq->emplace_back((const char*) record_seq_rs + trim5, dest_pos - (trim5+trim3)); // SEQ
    xm->emplace_back( (const char*) record_xm_rs + trim5,  dest_pos - (trim5+trim3)); // XM
    npushed++;                                                                  // +1

    Rcpp::Rcout << "dest_pos:" << dest_pos << ", record_width:" << record_width << ", query_width:" << query_width << ", strand:" << strand.back() << "\n";
    Rcpp::Rcout << seq->back()/*.substr(0, std::min(record_width, 100))*/ << "\n";
    Rcpp::Rcout << xm->back()/*.substr(0, std::min(record_width, 100))*/ << "\n\n";
  }

  // cleaning
  hts_base_mod_state_free(mod_state);                                           // free base modification state structure
  bam_destroy1(bam_rec);                                                        // free BAM alignment structure
  hts_close(bam_fp);                                                            // close BAM file
  if (thread_pool.pool) hts_tpool_destroy(thread_pool.pool);                    // free thread pool
  free(record_seq_rs);                                                          // and free manually allocated memory
  free(record_xm_rs);
  free(record_seq);
  free(record_xm);

  // wrap and return the results
  Rcpp::DataFrame res = Rcpp::DataFrame::create(                                // final DF
    Rcpp::Named("rname") = rname,                                               // numeric ids (factor) for reference names
    Rcpp::Named("strand") = strand,                                             // numeric ids (factor) for reference strands
    Rcpp::Named("start") = start                                                // start positions of reads
  );

  // factor levels
  std::vector<std::string> chromosomes (                                        // vector of reference names
      bam_hdr->target_name, bam_hdr->target_name + bam_hdr->n_targets);
  std::vector<std::string> strands = {"+", "-"};

  Rcpp::IntegerVector col_rname = res["rname"];                                 // make rname a factor
  col_rname.attr("class") = "factor";
  col_rname.attr("levels") = chromosomes;

  Rcpp::IntegerVector col_strand = res["strand"];                               // make strand a factor
  col_strand.attr("class") = "factor";
  col_strand.attr("levels") = strands;

  Rcpp::XPtr<std::vector<std::string>> seq_xptr(seq, true);
  res.attr("seq_xptr") = seq_xptr;                                              // external pointer to sequences
  Rcpp::XPtr<std::vector<std::string>> xm_xptr(xm, true);
  res.attr("xm_xptr") = xm_xptr;                                                // external pointer to methylation strings

  res.attr("nrecs") = nrecs;                                                    // number of records in BAM file
  res.attr("npushed") = npushed;                                                // number of records pushed to data.frame

  return(res);
}




// #############################################################################
// test code and sourcing don't work on OS X
/*** R
*/
// #############################################################################