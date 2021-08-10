#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

// [[Rcpp::depends(Rhtslib)]]

// an optimised attempt to read and preprocess BAM in place
// [[Rcpp::export("rcpp_read_bam")]]
Rcpp::List rcpp_read_bam (std::string fn,                                       // file name
                          int min_mapq,                                         // min read mapping quality
                          int min_baseq,                                        // min base quality
                          bool skip_duplicates                                  // skip marked duplicates
)
{
  // file IO
  htsFile *bam_fp = hts_open(fn.c_str(), "r");                                  // try open file
  if (bam_fp==NULL) Rcpp::stop("Unable to open BAM file for reading");          // fall back if error
  bam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);                                    // try read file header
  if (bam_hdr==NULL) Rcpp::stop("Unable to read BAM header");                   // fall back if error  
  bam1_t *bam_rec = bam_init1();                                                // create BAM alignment structure

  // main containers
  std::vector<std::string> qname, seq, seq_rs, xm, xm_rs;
  std::vector<char> strand;
  std::vector<int> rname, start, width;
  int nreads = 0, npairs = 0;                                                   // read/pair counters
  uint16_t prev_l_qname = 0;                                                    // length of qname of the previous read
  char *prev_qname = (char*) malloc(256);                                       // qname of the previous read
  
  // process alignments
  while( sam_read1(bam_fp, bam_hdr, bam_rec) > 0 ) {                            // rec by rec
    if ((nreads & 1048575) == 0) Rcpp::checkUserInterrupt();                    // checking for the interrupt
    if ((bam_rec->core.qual < min_mapq) ||                                      // skip if mapping quality < min.mapq
        (skip_duplicates & (bam_rec->core.flag & BAM_FDUP))) continue;          // or if record is an optical/PCR duplicate
    nreads++;                                                                   // count quality reads
    
    // prepare seq and XM strings
    uint32_t read_len = bam_rec->core.l_qseq;                                   // read length
    uint8_t *rec_qual = bam_get_qual(bam_rec);                                  // quality string (Phred scale with no +33 offset)
    char *rec_xm = (char*) bam_aux_get(bam_rec, "XM") + 1;                      // methylation string without leading 'Z'
    if (rec_xm==NULL) Rcpp::stop("Methylation string is absent/incomplete");    // fall back if absent 
    uint8_t *rec_pseq = bam_get_seq(bam_rec);                                   // packed sequence string (4 bit per base)
    char *rec_seq = (char*) malloc(read_len);                                   // unpacked sequence string
    for (int i=0; i<read_len; i++) {                                            // iterate over 4bit nucleotides
      if (rec_qual[i]<min_baseq) {                                              // if base quality < min.baseq
        rec_seq[i] = 'N';                                                       // sequence char to 'N'
        rec_xm[i] = '-';                                                        // methylation char to '-'
      } else {                                                                  // if good quality
        rec_seq[i] = seq_nt16_str[bam_seqi(rec_pseq,i)];                        // turn nucleotide into char (array is "=ACMGRSVTWYHKDBN")
      }
    }
    
    // apply CIGAR
    std::string rec_seq_rs, rec_xm_rs;                                          // seq and XM strings in reference space
    uint32_t rec_n_cigar = bam_rec->core.n_cigar;                               // number of CIGAR operations
    uint32_t *rec_cigar = bam_get_cigar(bam_rec);                               // CIGAR array
    uint32_t query_pos = 0;                                                     // current position in query string
    for (int i=0; i<rec_n_cigar; i++) {                                         // op by op
      uint32_t cigar_op = bam_cigar_op(rec_cigar[i]);                           // CIGAR operation
      uint32_t cigar_oplen = bam_cigar_oplen(rec_cigar[i]);                     // CIGAR operation length
      switch(cigar_op) {
        case BAM_CMATCH :                                                       // 'M', 0
        case BAM_CEQUAL :                                                       // '=', 7
        case BAM_CDIFF :                                                        // 'X', 8
          rec_seq_rs.append(rec_seq + query_pos, cigar_oplen);
          rec_xm_rs.append(rec_xm + query_pos, cigar_oplen);
          query_pos += cigar_oplen;
          break;
        case BAM_CINS :                                                         // 'I', 1
        case BAM_CSOFT_CLIP :                                                   // 'S', 4
          query_pos += cigar_oplen;
          break;
        case BAM_CDEL :                                                         // 'D', 2
        case BAM_CREF_SKIP :                                                    // 'N', 3
          rec_seq_rs.append(cigar_oplen, 'N');
          rec_xm_rs.append(cigar_oplen, '-');
          break;
        case BAM_CHARD_CLIP :                                                   // 'H', 5
        case BAM_CPAD :                                                         // 'P', 6
        case BAM_CBACK :
          break;
        default :
          Rcpp::stop("Unknown CIGAR operation for BAM entry %s",                // unknown CIGAR operation
                     bam_get_qname(bam_rec));
      }
      // if (query_pos < read_len-1)                                               // if not the entire string was converted to the reference space
      //   rec_seq_rs.append(rec_seq + query_pos, read_len - query_pos);           // copy the rest as GenomicAlignments::sequenceLayer does
    }
    
    // trim second read
    uint16_t rec_l_qname = bam_rec->core.l_qname;
    char *rec_qname = bam_get_qname(bam_rec);
    if ((rec_l_qname==prev_l_qname) & (strcmp(rec_qname, prev_qname)==0)) {
      
    }
    
    // append to main containers
    qname.push_back(rec_qname);
    rname.push_back(bam_rec->core.tid + 1);
    seq.push_back(std::string (rec_seq, read_len));
    seq_rs.push_back(rec_seq_rs);
    xm.push_back(std::string (rec_xm, read_len));
    xm_rs.push_back(rec_xm_rs);
    
    // save qname and its length
    prev_l_qname = bam_rec->core.l_qname;
    strcpy(prev_qname, rec_qname);
  }

  // cleaning
  bam_destroy1(bam_rec);                                                        // clean BAM alignment structure 
  hts_close(bam_fp);                                                            // close BAM file
  
  // wrap and return the results
  std::vector<std::string> chromosomes (                                        // vector of reference names
      bam_hdr->target_name, bam_hdr->target_name + bam_hdr->n_targets);
  
  Rcpp::List res = Rcpp::List::create(                                          // final list
    Rcpp::Named("nreads") = nreads,                                             // number of quality reads
    Rcpp::Named("npairs") = npairs,                                             // number of read pairs
    Rcpp::Named("chromosomes") = chromosomes,                                   // reference names
    Rcpp::Named("qname") = qname,                                               // query names
    Rcpp::Named("rname") = rname,                                               // numeric ids for reference names
    Rcpp::Named("seq") = seq,                                                   // sequences
    Rcpp::Named("seq_refspace") = seq_rs,
    Rcpp::Named("XM") = xm,                                                     // methylation strings
    Rcpp::Named("XM_refspace") = xm_rs
  );
  return(res);
}


// #############################################################################
// test code and sourcing doesn't work on OS X
/*** R
*/
// #############################################################################