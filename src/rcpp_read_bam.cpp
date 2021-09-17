#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

// [[Rcpp::depends(Rhtslib)]]

// An optimised attempt to read and preprocess BAM in place. To do:
// [ ] OpenMP SIMD?
// [+] HTSlib threads
// [+] rec_seq_rs and rec_xm_rs as char*
// [?] reverse QNAME
// [ ] free resources on interrupt

// [[Rcpp::export]]
Rcpp::DataFrame rcpp_read_bam_paired (std::string fn,                           // file name
                                      int min_mapq,                             // min read mapping quality
                                      int min_baseq,                            // min base quality
                                      bool skip_duplicates,                     // skip marked duplicates
                                      int nthreads)                             // HTSlib threads, >1 for multiple
{
  // constants
  size_t max_qname_width = 256;                                                 // max QNAME length, not expanded yet, ever error-prone?
  size_t max_templ_width = 1024;                                                // max insert size, expanded if necessary
  
  // file IO
  htsFile *bam_fp = hts_open(fn.c_str(), "r");                                  // try open file
  if (bam_fp==NULL) Rcpp::stop("Unable to open BAM file for reading");          // fall back if error
  htsThreadPool thread_pool = {NULL, 0};                                        // thread pool cuts time by 30%
  if (nthreads>1) {
    thread_pool.pool = hts_tpool_init(nthreads);                                // when initiated for 2 threads
    hts_set_opt(bam_fp, HTS_OPT_THREAD_POOL, &thread_pool);                     // and bound to the file pointer
  }
  bam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);                                    // try read file header
  if (bam_hdr==NULL) Rcpp::stop("Unable to read BAM header");                   // fall back if error  
  bam1_t *bam_rec = bam_init1();                                                // create BAM alignment structure
  
  // main containers
  std::vector<std::string> seq, xm;                                             // SEQ, XM
  std::vector<int> rname, strand, start;                                        // id for RNAME, id for CT==1/GA==2, POS
  int nrecs = 0, ntempls = 0;                                                   // counters: BAM records, templates (read pairs)
  
  // reserve some memory
  rname.reserve(10000); strand.reserve(10000); start.reserve(10000); 
  seq.reserve(10000); xm.reserve(10000);
  
  // template holders
  char *templ_qname = (char*) malloc(max_qname_width * sizeof(char));           // template QNAME
  char *rec_qname   = (char*) malloc(max_qname_width * sizeof(char));           // current read QNAME
  uint8_t *templ_qual_rs = (uint8_t*) malloc(max_templ_width * sizeof(uint8_t));// template QUAL array
  uint8_t *templ_seq_rs  = (uint8_t*) malloc(max_templ_width * sizeof(uint8_t));// template SEQ array
  uint8_t *templ_xm_rs   = (uint8_t*) malloc(max_templ_width * sizeof(uint8_t));// template XM array
  int templ_rname = 0, templ_start = 0, templ_strand = 0, templ_width = 0;      // template RNAME, POS, STRAND, ISIZE
  
  #define unsorted (ntempls==0) || (nrecs/(ntempls>>1) < 3)                     /* TRUE if no templates OR fraction of two-read templates is <67% */
  #define push_template {               /* pushing template data to vectors */ \
    rname.push_back(templ_rname + 1);                            /* RNAME+1 */ \
    strand.push_back(templ_strand);                               /* STRAND */ \
    start.push_back(templ_start + 1);                              /* POS+1 */ \
    seq.push_back(std::string((char *) templ_seq_rs, templ_width));  /* SEQ */ \
    xm.push_back( std::string((char *) templ_xm_rs,  templ_width));   /* XM */ \
    ntempls++;                                                        /* +1 */ \
  }
  
  // process alignments
  while( sam_read1(bam_fp, bam_hdr, bam_rec) > 0 ) {                            // rec by rec
    nrecs++;                                                                    // BAM alignment records ++
    if ((nrecs & 0xFFFFF) == 0) {                                               // every ~1M reads
      Rcpp::checkUserInterrupt();                                               // checking for the interrupt
      if (unsorted) break;                                                      // break out if seemingly unsorted
    }
    if ((bam_rec->core.qual < min_mapq) ||                                      // skip if mapping quality < min.mapq
        (!(bam_rec->core.flag & BAM_FPROPER_PAIR)) ||                           // or if not a proper pair
        (skip_duplicates && (bam_rec->core.flag & BAM_FDUP))) continue;         // or if record is an optical/PCR duplicate
    
    // save qname
    strcpy(rec_qname, bam_get_qname(bam_rec));
    // std::reverse(rec_qname, rec_qname + bam_rec->core.l_qname                // reversing it to make strcmp faster:
    //                - bam_rec->core.l_extranul - 1);                          // doesn't make any sense now, because compare only once
    
    // check if not the same template (QNAME)
    if ((strcmp(templ_qname, rec_qname) != 0)) {                
      // store previous template if not first record
      if (nrecs!=1) push_template;
      
      // initialize new template
      strcpy(templ_qname, rec_qname);                                           // store template QNAME
      templ_rname = bam_rec->core.tid;                                          // store template RNAME
      templ_start = bam_rec->core.pos < bam_rec->core.mpos ?                    // smallest of POS,MPOS is a start
        bam_rec->core.pos : bam_rec->core.mpos;
      templ_width = abs(bam_rec->core.isize);                                   // template ISIZE
      char *rec_strand = (char*) bam_aux_get(bam_rec, "XG");                    // genome strand
      if (rec_strand==NULL)                                                     // fall back if absent
        Rcpp::stop("Genome strand is missing for BAM record #%i", nrecs);
      templ_strand = ( rec_strand[1] == 'C' ) ? 1 : 2 ;                         // STRAND is 1 if "ZCT"/"+", 2 if "ZGA"/"-"
      
      // resize containers if necessary
      if (templ_width > max_templ_width) {
        max_templ_width = templ_width;                                          // expand template holders
        templ_qual_rs = (uint8_t *) realloc(templ_qual_rs, max_templ_width);
        templ_seq_rs  = (uint8_t *) realloc(templ_seq_rs,  max_templ_width);
        templ_xm_rs   = (uint8_t *) realloc(templ_xm_rs,   max_templ_width);
        if (templ_qual_rs==NULL || templ_seq_rs==NULL || templ_xm_rs==NULL)     // check memory allocation
          Rcpp::stop("Unable to allocate memory for BAM record #%i", nrecs);
      }
      std::fill_n(templ_qual_rs, templ_width, (uint8_t) min_baseq);             // clean template holders
      std::fill_n(templ_seq_rs,  templ_width, 'N');
      std::fill_n(templ_xm_rs,   templ_width, '-');
     }
    
    // add another read to the template
    // source containers
    uint8_t *rec_qual = bam_get_qual(bam_rec);                                  // quality string (Phred scale with no +33 offset)
    char *rec_xm = (char*) bam_aux_get(bam_rec, "XM");                          // methylation string
    if (rec_xm==NULL)                                                           // fall back if absent 
      Rcpp::stop("Methylation string is missing for BAM record #%i", nrecs); 
    rec_xm++;                                                                   // remove leading 'Z'
    uint8_t *rec_pseq = bam_get_seq(bam_rec);                                   // packed sequence string (4 bit per base)
    
    // apply CIGAR
    uint32_t n_cigar = bam_rec->core.n_cigar;                                   // number of CIGAR operations
    uint32_t *rec_cigar = bam_get_cigar(bam_rec);                               // CIGAR array
    uint32_t query_pos = 0;                                                     // starting position in query array
    uint32_t dest_pos = bam_rec->core.pos - templ_start;                        // starting position in destination array
    for (int i=0; i<n_cigar; i++) {                                             // op by op
      uint32_t cigar_op = bam_cigar_op(rec_cigar[i]);                           // CIGAR operation
      uint32_t cigar_oplen = bam_cigar_oplen(rec_cigar[i]);                     // CIGAR operation length
      switch(cigar_op) {
        case BAM_CMATCH :                                                       // 'M', 0
        case BAM_CEQUAL :                                                       // '=', 7
        case BAM_CDIFF :                                                        // 'X', 8
          for (int j=0; j<cigar_oplen; j++) {
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
  
  // stop if single-end or seemingly unsorted
  if (unsorted)
    Rcpp::stop("BAM seems to be predominantly single-end or not sorted by QNAME. Single-end alignments are not supported yet. If paired-end, please sort using 'samtools sort -n -o out.bam in.bam'");
  
  // push last, yet unsaved template
  push_template;
  
  // cleaning
  bam_destroy1(bam_rec);                                                        // clean BAM alignment structure 
  hts_close(bam_fp);                                                            // close BAM file
  if (thread_pool.pool) hts_tpool_destroy(thread_pool.pool);                    // free thread pool
  free(templ_qname);                                                            // and free manually allocated memory
  free(rec_qname); 
  free(templ_qual_rs);
  free(templ_seq_rs);
  free(templ_xm_rs);
  
  // wrap and return the results
  std::vector<std::string> chromosomes (                                        // vector of reference names
      bam_hdr->target_name, bam_hdr->target_name + bam_hdr->n_targets);
  std::vector<std::string> strands = {"+", "-"};
  // make rname a factor
  Rcpp::IntegerVector w_rname = Rcpp::wrap(rname);
  w_rname.attr("class") = "factor";
  w_rname.attr("levels") = chromosomes;
  // make strand a factor
  Rcpp::IntegerVector w_strand = Rcpp::wrap(strand);
  w_strand.attr("class") = "factor";
  w_strand.attr("levels") = strands;
  
  Rcpp::DataFrame res = Rcpp::DataFrame::create(                                // final DF
    Rcpp::Named("rname") = w_rname,                                             // numeric ids (factor) for reference names
    Rcpp::Named("strand") = w_strand,                                           // numeric ids (factor) for reference strands
    Rcpp::Named("start") = start,                                               // start positions of reads
    Rcpp::Named("seq") = seq,                                                   // sequences
    Rcpp::Named("XM") = xm                                                      // methylation strings
  );
  
  return(res);
}

// 
// Reads single records, trims some, requires post-trimming.
// Not efficient at all
// 
// // [[Rcpp::export("rcpp_read_bam")]]
// Rcpp::List rcpp_read_bam (std::string fn,                                       // file name
//                           int min_mapq,                                         // min read mapping quality
//                           int min_baseq,                                        // min base quality
//                           bool skip_duplicates)                                 // skip marked duplicates
// {
//   // file IO
//   htsFile *bam_fp = hts_open(fn.c_str(), "r");                                  // try open file
//   if (bam_fp==NULL) Rcpp::stop("Unable to open BAM file for reading");          // fall back if error
//   bam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);                                    // try read file header
//   if (bam_hdr==NULL) Rcpp::stop("Unable to read BAM header");                   // fall back if error  
//   bam1_t *bam_rec = bam_init1();                                                // create BAM alignment structure
// 
//   // main containers
//   std::vector<std::string> qname, seq, xm;                                      // QNAME, SEQ, XM
//   std::vector<int> flag, rname, strand, start, width;                           // FLAG, id for RNAME + 1, id for CT==1/GA==2, POS + 1, read length
//   int nrecs = 0, nreads = 0, npaired = 0, ntempls = 0;                          // counters: BAM records, quality reads, paired reads, mergeable templates
//   qname.reserve(10000); seq.reserve(10000); xm.reserve(10000);
//   flag.reserve(10000); rname.reserve(10000); strand.reserve(10000);
//   start.reserve(10000); width.reserve(10000);
// 
//   std::string rec_seq_rs, rec_xm_rs;                                            // seq and XM strings in reference space
//   rec_seq_rs.reserve(255);
//   rec_xm_rs.reserve(255);
//   
//   // process alignments
//   while( sam_read1(bam_fp, bam_hdr, bam_rec) > 0 ) {                            // rec by rec
//     nrecs++;                                                                    // BAM alignment records ++
//     if ((nrecs & 1048575) == 0) Rcpp::checkUserInterrupt();                     // checking for the interrupt
//     if ((bam_rec->core.qual < min_mapq) ||                                      // skip if mapping quality < min.mapq
//         (skip_duplicates && (bam_rec->core.flag & BAM_FDUP))) continue;         // or if record is an optical/PCR duplicate
//     nreads++;                                                                   // quality reads ++
//     
//     // reverse qname to make sort/search/compare faster
//     char *rec_qname = bam_get_qname(bam_rec);                                   // QNAME
//     std::reverse(rec_qname, rec_qname + bam_rec->core.l_qname                   // reversing it
//                    - bam_rec->core.l_extranul - 1);
//     // genome strand
//     int rec_strand = ( bam_aux_get(bam_rec, "XG")[1] == 'C' ) ? 1 : 2 ;         // 1 if "ZCT"/"+", 2 if "ZGA"/"-"
//     
//     // append to main containers
//     flag.push_back(bam_rec->core.flag);
//     qname.push_back(rec_qname);                                                 // read QNAME
//     rname.push_back(bam_rec->core.tid + 1);                                     // id of read RNAME + 1
//     strand.push_back(rec_strand);                                               // id of reference strand
//     start.push_back(bam_rec->core.pos + 1);                                     // read POS + 1
//     
//     // prepare seq and XM strings
//     uint32_t read_len = bam_rec->core.l_qseq;                                   // read length
//     uint8_t *rec_qual = bam_get_qual(bam_rec);                                  // quality string (Phred scale with no +33 offset)
//     char *rec_xm = (char*) bam_aux_get(bam_rec, "XM") + 1;                      // methylation string without leading 'Z'
//     if (rec_xm==NULL) Rcpp::stop("Methylation string is absent/incomplete");    // fall back if absent 
//     uint8_t *rec_pseq = bam_get_seq(bam_rec);                                   // packed sequence string (4 bit per base)
//     char rec_seq [read_len];                                   // unpacked sequence string
//     for (int i=0; i<read_len; i++) {                                            // iterate over 4bit nucleotides
//       if (rec_qual[i]<min_baseq) {                                              // if base quality < min.baseq
//         rec_seq[i] = 'N';                                                       // sequence char to 'N'
//         rec_xm[i] = '-';                                                        // methylation char to '-'
//       } else {                                                                  // if good quality
//         rec_seq[i] = seq_nt16_str[bam_seqi(rec_pseq,i)];                        // turn nucleotide into char (array is "=ACMGRSVTWYHKDBN")
//       }
//     }
//     
//     // apply CIGAR
//     uint32_t rec_n_cigar = bam_rec->core.n_cigar;                               // number of CIGAR operations
//     uint32_t *rec_cigar = bam_get_cigar(bam_rec);                               // CIGAR array
//     uint32_t query_pos = 0;                                                     // current position in query string
//     for (int i=0; i<rec_n_cigar; i++) {                                         // op by op
//       uint32_t cigar_op = bam_cigar_op(rec_cigar[i]);                           // CIGAR operation
//       uint32_t cigar_oplen = bam_cigar_oplen(rec_cigar[i]);                     // CIGAR operation length
//       switch(cigar_op) {
//         case BAM_CMATCH :                                                       // 'M', 0
//         case BAM_CEQUAL :                                                       // '=', 7
//         case BAM_CDIFF :                                                        // 'X', 8
//           rec_seq_rs.append(rec_seq + query_pos, cigar_oplen);
//           rec_xm_rs.append(rec_xm + query_pos, cigar_oplen);
//           query_pos += cigar_oplen;
//           break;
//         case BAM_CINS :                                                         // 'I', 1
//         case BAM_CSOFT_CLIP :                                                   // 'S', 4
//           query_pos += cigar_oplen;
//           break;
//         case BAM_CDEL :                                                         // 'D', 2
//         case BAM_CREF_SKIP :                                                    // 'N', 3
//           rec_seq_rs.append(cigar_oplen, 'N');
//           rec_xm_rs.append(cigar_oplen, '-');
//           break;
//         case BAM_CHARD_CLIP :                                                   // 'H', 5
//         case BAM_CPAD :                                                         // 'P', 6
//         case BAM_CBACK :
//           break;
//         default :
//           Rcpp::stop("Unknown CIGAR operation for BAM entry %s",                // unknown CIGAR operation
//                      bam_get_qname(bam_rec));
//       }
//       // if (query_pos < read_len-1)                                               // if not the entire string was converted to the reference space
//       //   rec_seq_rs.append(rec_seq + query_pos, read_len - query_pos);           // copy the rest as GenomicAlignments::sequenceLayer does
//     }
//     
//     
//     // trim leftmost READ2, post-trimming of rightmost READ2 is done elsewhere
//     if ((bam_rec->core.flag & BAM_FPROPER_PAIR) &&                              // if read is a proper pair
//         (bam_rec->core.flag & BAM_FREAD2)) {                                    // and it's a Read2
//       int pos_diff = bam_rec->core.mpos - bam_rec->core.pos;                    // get distance between start positions of mates
//       if (pos_diff == 0) {                                                      // if mates fully overlap
//         rec_seq_rs.clear();                                                     // then READ2 is erased
//         rec_xm_rs.clear();
//       } else if ((pos_diff > 0) && (pos_diff < rec_seq_rs.length())) {          // else if it's a leftmost read && distance between mates is less than refspaced read length
//         rec_seq_rs.erase(pos_diff);                                             // trim leftmost Read2 until rightmost Read1
//         rec_xm_rs.erase(pos_diff);
//       }
//     }
//     
//     // append to main containers
//     width.push_back(rec_seq_rs.length());                                       // read length
//     seq.push_back(rec_seq_rs);                                                  // read SEQ
//     xm.push_back(rec_xm_rs);                                                    // read XM
//     
//     // clearing temporary vars
//     rec_seq_rs.clear();
//     rec_xm_rs.clear();
//   }
// 
//   // cleaning
//   bam_destroy1(bam_rec);                                                        // clean BAM alignment structure 
//   hts_close(bam_fp);                                                            // close BAM file
//   
//   // wrap and return the results
//   std::vector<std::string> chromosomes (                                        // vector of reference names
//       bam_hdr->target_name, bam_hdr->target_name + bam_hdr->n_targets);
//   std::vector<std::string> strands = {"+", "-"};                                // vector of DNA strands
//   std::vector<int> stats = {nrecs, nreads, npaired, ntempls};                   // statistics
//   
//   Rcpp::List res = Rcpp::List::create(                                          // final list
//     Rcpp::Named("stats") = stats,                                               // BAM loading statistics
//     Rcpp::Named("chromosomes") = chromosomes,                                   // lookup table of reference names
//     Rcpp::Named("strands") = strands,                                           // lookup table of DNA strands
//     Rcpp::Named("flag") = flag,                                                 // flag
//     Rcpp::Named("qname") = qname,                                               // query names
//     Rcpp::Named("rname") = rname,                                               // numeric ids for reference names
//     Rcpp::Named("strand") = strand,                                             // numeric ids for reference strands
//     Rcpp::Named("start") = start,                                               // start positions of reads
//     Rcpp::Named("width") = width,                                               // lengths of reads
//     Rcpp::Named("seq") = seq,                                                   // sequences
//     Rcpp::Named("XM") = xm                                                      // methylation strings
//   );
//   return(res);
// }
// 
// 
// 
// // In place post-trimming of rightmost READ2
// //
// 
// // [[Rcpp::export]]
// void rcpp_posttrim_read2(Rcpp::DataFrame &df) {
//   Rcpp::IntegerVector flag = df["flag"];
//   Rcpp::IntegerVector start = df["start"];
//   Rcpp::CharacterVector qname = df["qname"];
//   Rcpp::CharacterVector seq = df["seq"];
//   Rcpp::CharacterVector xm = df["XM"];
//   int dpos, ndel;
//   
//   // row by row
//   for (size_t x=1; x<flag.size(); x++) {
//     if ((x & 1048575) == 0) Rcpp::checkUserInterrupt();                         // checking for the interrupt
//     
//     // trim rightmost READ2
//     if ((qname[x]==qname[x-1]) &&                                               // if this is a same template (df is presorted on QNAME, POS)
//         (flag[x] & BAM_FPROPER_PAIR) &&                                         // and this read is a proper pair
//         (flag[x] & BAM_FREAD2)) {                                               // and it's a READ2
//       dpos = start[x] - start[x-1];                                             // get distance between start positions of mates
//       ndel = xm[x-1].size() - dpos;                                             // number of chars to delete from the left side of READ2
//       if ((dpos > 0) && (ndel > 0)) {                                           // if it's a rightmost read && distance between mates is less than refspaced read length
//         seq[x] = Rcpp::as<std::string>(seq[x]).erase(0, ndel);                  // trim rightmost READ2 until the end of leftmost Read1
//         xm[x] = Rcpp::as<std::string>(xm[x]).erase(0, ndel);
//         start[x] = start[x-1] + xm[x-1].size();                                 // adjust new POS
//       }
//     }
//   }
// }


// #############################################################################
// test code and sourcing don't work on OS X
/*** R
*/
// #############################################################################