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
  std::vector<std::string>* seqxm = new std::vector<std::string>;               // SEQXM, leftmost 4 bits are SEQ and rightmost 4 are XM
  std::vector<int> rname, strand, start;                                        // id for RNAME, id for CT==1/GA==2, POS
  int nrecs = 0, ntempls = 0;                                                   // counters: BAM records, templates (consecutive proper read pairs)
  
  // reserve some memory
  rname.reserve(0xFFFFF); strand.reserve(0xFFFFF); start.reserve(0xFFFFF); 
  seqxm->reserve(0xFFFFF);
  
  // template holders
  char *templ_qname = (char*) malloc(max_qname_width * sizeof(char));           // template QNAME
  uint8_t *templ_qual_rs  = (uint8_t*)malloc(max_templ_width * sizeof(uint8_t));// template QUAL array
  uint8_t *templ_seqxm_rs = (uint8_t*)malloc(max_templ_width * sizeof(uint8_t));// template SEQXM array
  int templ_rname = 0, templ_start = 0, templ_strand = 0, templ_width = 0;      // template RNAME, POS, STRAND, ISIZE
  
  #define push_template {                                        /* pushing template data to vectors */ \
    rname.push_back(templ_rname + 1);                                                     /* RNAME+1 */ \
    strand.push_back(templ_strand);                                                        /* STRAND */ \
    start.push_back(templ_start + trim5 + 1);                                               /* POS+1 */ \
    seqxm->emplace_back((const char*) templ_seqxm_rs + trim5, templ_width - (trim5+trim3)); /* SEQXM */ \
    ntempls++;                                                                                 /* +1 */ \
  }
  
  // process alignments
  while( sam_read1(bam_fp, bam_hdr, bam_rec) > 0 ) {                            // rec by rec
    nrecs++;                                                                    // BAM alignment records ++
    if ((nrecs & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();                     // every ~1M reads check for the interrupt
    
    if ((bam_rec->core.flag & BAM_FUNMAP) ||                                    // skip if unmapped
        (!(bam_rec->core.flag & BAM_FPROPER_PAIR)) ||                           // or if not a proper pair
        (bam_rec->core.qual < min_mapq) ||                                      // or if mapping quality < min.mapq
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
      templ_strand = 2 - (rec_strand[1] == 'C');                                // STRAND is 1 if "ZCT"/"+", 2 if "ZGA"/"-"
      
      // resize containers if necessary
      if (templ_width > max_templ_width) {
        max_templ_width = templ_width;                                          // expand template holders
        templ_qual_rs  = (uint8_t *) realloc(templ_qual_rs,  max_templ_width);
        templ_seqxm_rs = (uint8_t *) realloc(templ_seqxm_rs, max_templ_width);
        if (!templ_qual_rs || !templ_seqxm_rs) Rcpp::stop("Unable to allocate memory for BAM record #%i", nrecs); // check memory allocation
      }
      std::memset(templ_qual_rs, (uint8_t) min_baseq, templ_width);             // clean template holders
      std::memset(templ_seqxm_rs, 0b11111011, templ_width);                     // fill SEQXM with 'N-', i.e., '15,11'
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
              templ_seqxm_rs[dest_pos+j] = bam_seqi_shifted(rec_pseq,query_pos+j) | ctx_to_idx(rec_xm[query_pos+j]);
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
  free(templ_seqxm_rs);
  
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
  
  Rcpp::XPtr<std::vector<std::string>> seqxm_xptr(seqxm, true);
  res.attr("seqxm_xptr") = seqxm_xptr;                                          // external pointer to packed sequences + methylation strings
  
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
  std::vector<std::string>* seqxm = new std::vector<std::string>;               // SEQXM, leftmost 4 bits are SEQ and rightmost 4 are XM
  std::vector<int> rname, strand, start;                                        // id for RNAME, id for CT==1/GA==2, POS
  int nrecs = 0, npushed = 0;                                                   // counters: BAM records read, BAM records pushed to data.table
  
  // reserve some memory
  rname.reserve(0xFFFFF); strand.reserve(0xFFFFF); start.reserve(0xFFFFF); 
  seqxm->reserve(0xFFFFF);
  
  // read holders
  int record_width = max_record_width;                                          // record ISIZE/TLEN
  uint8_t *record_seqxm_rs  = (uint8_t*) malloc(record_width * sizeof(uint8_t));// record SEQXM array
  
  // process alignments
  while( sam_read1(bam_fp, bam_hdr, bam_rec) > 0 ) {                            // rec by rec
    nrecs++;                                                                    // BAM alignment records ++
    if ((nrecs & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();                     // every ~1M reads check for the interrupt

    if ((bam_rec->core.flag & BAM_FUNMAP) ||                                    // skip if unmapped
        (bam_rec->core.qual < min_mapq) ||                                      // or if mapping quality < min.mapq
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
      record_seqxm_rs = (uint8_t *) realloc(record_seqxm_rs, record_width * sizeof(uint8_t));
      if (!record_seqxm_rs) Rcpp::stop("Unable to allocate memory for BAM record #%i", nrecs); // check memory allocation
    }
    
    // prepare for the new record
    std::memset(record_seqxm_rs, 0b11111011, record_width);                     // fill SEQXM with 'N-', i.e., '15,11'
    
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
            record_seqxm_rs[dest_pos+j] = bam_seqi_shifted(record_pseq,query_pos+j) | ctx_to_idx(record_xm[query_pos+j]);
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
    seqxm->emplace_back((const char*) record_seqxm_rs + trim5, dest_pos - (trim5+trim3)); // SEQXM
    npushed++;                                                                  // +1 
  }
  
  // cleaning
  bam_destroy1(bam_rec);                                                        // clean BAM alignment structure 
  hts_close(bam_fp);                                                            // close BAM file
  if (thread_pool.pool) hts_tpool_destroy(thread_pool.pool);                    // free thread pool
  free(record_seqxm_rs);                                                        // and free manually allocated memory
  
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
  
  Rcpp::XPtr<std::vector<std::string>> seqxm_xptr(seqxm, true);
  res.attr("seqxm_xptr") = seqxm_xptr;                                          // external pointer to packed sequences + methylation strings
  
  res.attr("nrecs") = nrecs;                                                    // number of records in BAM file
  res.attr("npushed") = npushed;                                                // number of records pushed to data.frame
  
  return(res);
}


// #############################################################################

// LONG-READ SINGLE-END BAM
// https://samtools.github.io/hts-specs/SAMv1.pdf
// https://samtools.github.io/hts-specs/SAMtags.pdf
// https://github.com/samtools/hts-specs/tree/master/test/SAMtags
// https://github.com/samtools/htslib/tree/develop/test/base_mods

// EXPLAIN THE LOGIC IN DOCS - NAMELY, DIFF BETWEEN SHORT-READ (BISULFITE) AND
// LONG-READ (NATIVE) SEQUENCING METHYLATION CALLING (GENOME IS ABSOLUTELY
// NECESSARY FOR THE FIRST BUT NOT THE SECOND)

// SEQXM packing is not efficient here. Would be great to rewrite methylation
// calling to allow all IUPAC bases, which will use 4096-byte context lookup
// tables, HTSlib codes for bases and, therefore, will save some ops by
// avoiding unnecessary conversions

// [[Rcpp::export]]
Rcpp::DataFrame rcpp_read_bam_mm_single (std::string fn,                        // file name
                                         int min_mapq,                          // min read mapping quality
                                         int min_baseq,                         // min base quality
                                         int min_prob,                          // min probability of 5mC modification
                                         bool highest_prob,                     // consider only if 5mC probability is the highest of all mods at particular pos
                                         bool skip_duplicates,                  // skip marked duplicates
                                         int trim5,                             // trim bases from 5'
                                         int trim3,                             // trim bases from 3'
                                         int nthreads)                          // HTSlib threads, >0 for multiple
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
  std::vector<std::string>* seqxm = new std::vector<std::string>;               // SEQXM, leftmost 4 bits are SEQ and rightmost 4 are XM
  std::vector<int> rname, strand, start;                                        // id for RNAME, id for CT==1/GA==2, POS
  int nrecs = 0, npushed = 0;                                                   // counters: BAM records read, BAM records pushed to data.table

  // reserve some memory
  rname.reserve(0xFFFFF); strand.reserve(0xFFFFF); start.reserve(0xFFFFF);
  seqxm->reserve(0xFFFFF);

  // read holders
  int query_width = max_query_width;                                            // NON-refspaced query length
  uint8_t *record_seq = (uint8_t*) malloc(sizeof(uint8_t) *(query_width+4));    // NON-refspaced query SEQ array, plus NN at the end
  uint8_t *record_xm[2];                                                        // NON-refspaced query 2D XM array for both strands
  for (int s=0; s<2; s++) record_xm[s] = (uint8_t*) malloc(query_width * sizeof(uint8_t)); // allocate memory for query 2D XM array
  int record_width = max_record_width;                                          // refspaced record ISIZE/TLEN
  uint8_t *record_seqxm_rs[2];                                                  // refspaced record 2D SEQXM array for both strands
  for (int s=0; s<2; s++) record_seqxm_rs[s] = (uint8_t*) malloc(record_width * sizeof(uint8_t)); // allocate memory for record 2D SEQXM array
  
  // process alignments
  while( sam_read1(bam_fp, bam_hdr, bam_rec) > 0 ) {                            // rec by rec
    nrecs++;                                                                    // BAM alignment records ++
    if ((nrecs & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();                     // every ~1M reads check for the interrupt

    if ((bam_rec->core.flag & BAM_FUNMAP) ||                                    // skip if unmapped
        (bam_rec->core.qual < min_mapq) ||                                      // or if mapping quality < min.mapq
        (skip_duplicates && (bam_rec->core.flag & BAM_FDUP))) continue;         // or if record is an optical/PCR duplicate

    int record_strand = (bool) (bam_rec->core.flag & BAM_FREVERSE);             // genome strand, 0 if forward, 1 if reverse

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
      record_seq = (uint8_t *) realloc(record_seq, (query_width + 4) * sizeof(uint8_t));
      for(int s=0; s<2; s++) record_xm[s] = (uint8_t*) realloc(record_xm[s], query_width * sizeof(uint8_t));
      if (!record_seq || !record_xm[0] || !record_xm[1]) Rcpp::stop("Unable to allocate memory for BAM record #%i", nrecs); // check memory allocation
    }
    if (record_width > max_record_width) {
      max_record_width = record_width;                                          // expand template holders
      for(int s=0; s<2; s++) record_seqxm_rs[s] = (uint8_t*) realloc(record_seqxm_rs[s], record_width * sizeof(uint8_t));
      if (!record_seqxm_rs[0] || !record_seqxm_rs[1]) Rcpp::stop("Unable to allocate memory for BAM record #%i", nrecs); // check memory allocation
    }

    // prepare for the new record
    std::memset(record_seqxm_rs[0], 0b11111011, record_width);                  // fill refspaced SEQXM holder for forward strand with 'N-'
    std::memset(record_seqxm_rs[1], 0b11111011, record_width);                  // fill refspaced SEQXM holder for reverse strand with 'N-'
    
    // unpack the sequence string, restore flanking NN's
    for (int i=0; i<query_width; i++) {
      record_seq[i+2] = seq_nt16_str[bam_seqi(record_pseq,i)];
    }
    std::memset(record_seq, 'N', 2);
    std::memset(record_seq+query_width+2, 'N', 2);

    // create a non-refspaced context strings for both strands
    for (int i=0; i<query_width; i++) {
      record_xm[0][i] = triad_to_ctx((record_seq+i+2), triad_forward_context);   // look up fwd context
      record_xm[1][i] = triad_to_ctx((record_seq+i),   triad_reverse_context);   // look up rev context
    }
    
    // BOTH STRANDS CAN HAVE OVERLAPPING MODS ('C+m' and 'G-m')!
    // parse base modifications: any location not reported is implicitly
    // assumed to contain no modification
    bool strand_has_mods[2] = {0, 0};                                           // bool array with TRUE if XM for fwd==[0] or rev==[1] strand has base mods
    bam_parse_basemod(bam_rec, mod_state);                                      // fill the mod_state structure
    while ((nmods = bam_next_basemod(bam_rec, mod_state, base_mods, max_nmods, &mod_pos)) > 0) { // cycle through modified bases
      int ismeth[2] = {0, 0},                                                   // bool array for pos having meth mod at fwd==[0] or rev==[1] strand;
        meth_prob[2] = {-2, -2},                                                // int array for probability of meth mod at fwd==[0] or rev==[1] strand;
        max_other_prob[2] = {-2, -2};                                           // int array for probability of any other mod at fwd==[0] or rev==[1] strand;
      for (int i=0; i<nmods; i++) {                                             // cycle through all mods of a current base
        if (base_mods[i].modified_base=='m' || base_mods[i].modified_base==-27551) { // if it's a 5mC (and any of 'C+m' or 'G-m')
          ismeth[base_mods[i].strand] = 1;                                      // base has meth mod
          meth_prob[base_mods[i].strand] = base_mods[i].qual;                   // record meth prob
        } else if (max_other_prob[base_mods[i].strand] < base_mods[i].qual) {   // if NOT a 5mC and probability is higher than max_other_prob
            max_other_prob[base_mods[i].strand] = base_mods[i].qual;            // record the highest probability of other modifications
        }
      }
      for (int s=0; s<2; s++) {                                                 // as the same pos can have mods on both strands, cycle through strands and apply modification to relevant context string
        int ctx_strand = abs(record_strand - s);                                // have to flip the context strand for revcomplemented query (because mods are always on NON-revcomp)
        if (ismeth[s] &&                                                        // if there is a C+m or G-m modification
            meth_prob[s]>=min_prob &&                                           // and its probability is not less than min_prob
            (!highest_prob || meth_prob[s]>max_other_prob[s]) &&                // and its probability is either highest or highest_prob==FALSE
            record_xm[ctx_strand][mod_pos]>'A') {                               // and its not a '.-'
          record_xm[ctx_strand][mod_pos] &= 0b11011111;                         // uppercase the context char
          strand_has_mods[ctx_strand] = 1;                                      // record that the context string for this strand has mods
        }
      }
    }
    
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
          if (record_qual[query_pos+j] >= min_baseq) {                          // apply CIGAR op and pack SEQ + XM simultaneously
            const uint8_t seq_idx = (seq_nt16_table[record_seq[query_pos+2+j]]) << 4;
            record_seqxm_rs[0][dest_pos+j] = seq_idx | ctx_to_idx(record_xm[0][query_pos+j]);
            record_seqxm_rs[1][dest_pos+j] = seq_idx | ctx_to_idx(record_xm[1][query_pos+j]);
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

    // pushing record data to vectors, once for the record strand (even if there are no 'C+m') and once again if the other strand has mods too (has 'G-m')
    strand_has_mods[record_strand] = 1;                                         // always push at least one context string
    for (int s=0; s<2; s++) {
      if (strand_has_mods[s]) {
        rname.push_back(bam_rec->core.tid + 1);                                 // RNAME+1
        strand.push_back(s + 1);                                                // STRAND is 1 if "CT"/"+", 2 if "GA"/"-"
        start.push_back(bam_rec->core.pos + trim5 + 1);                         // POS+1
        seqxm->emplace_back( (const char*) record_seqxm_rs[s] + trim5, dest_pos - (trim5+trim3)); // SEQXM
        npushed++;                                                              // +1
      }
    }
  }

  // cleaning
  hts_base_mod_state_free(mod_state);                                           // free base modification state structure
  bam_destroy1(bam_rec);                                                        // free BAM alignment structure
  hts_close(bam_fp);                                                            // close BAM file
  if (thread_pool.pool) hts_tpool_destroy(thread_pool.pool);                    // free thread pool
  for(int i=0; i<2; i++) free(record_seqxm_rs[i]);
  free(record_seq);
  for(int i=0; i<2; i++) free(record_xm[i]);

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

  Rcpp::XPtr<std::vector<std::string>> seqxm_xptr(seqxm, true);
  res.attr("seqxm_xptr") = seqxm_xptr;                                          // external pointer to sequences

  res.attr("nrecs") = nrecs;                                                    // number of records in BAM file
  res.attr("npushed") = npushed;                                                // number of records pushed to data.frame

  return(res);
}




// #############################################################################
// test code and sourcing don't work on OS X
/*** R
*/
// #############################################################################