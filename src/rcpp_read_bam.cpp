#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

// [[Rcpp::depends(Rhtslib)]]

// An optimised attempt to read and preprocess BAM in place
//
// optimization strategies:
// [ ] OpenMP SIMD?
// [ ] HTSlib threads?
// [-] POS==MPOS && second read && proper pair - skip most of ops <- no use in reality
// [ ] rec_seq_rs and rec_xm_rs as char*
// [ ] don't do CIGAR if single M op?
// [+] reverse QNAME
// [ ] optimize cx report as well - use 8 byes per cpg and do & instead of if

// [[Rcpp::export("rcpp_read_bam")]]
Rcpp::List rcpp_read_bam (std::string fn,                                       // file name
                          int min_mapq,                                         // min read mapping quality
                          int min_baseq,                                        // min base quality
                          bool skip_duplicates)                                 // skip marked duplicates
{
  // file IO
  htsFile *bam_fp = hts_open(fn.c_str(), "r");                                  // try open file
  if (bam_fp==NULL) Rcpp::stop("Unable to open BAM file for reading");          // fall back if error
  bam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);                                    // try read file header
  if (bam_hdr==NULL) Rcpp::stop("Unable to read BAM header");                   // fall back if error  
  bam1_t *bam_rec = bam_init1();                                                // create BAM alignment structure

  // main containers
  std::vector<std::string> qname, seq, xm;                                      // QNAME, SEQ, XM
  std::vector<int> flag, rname, strand, start, width;                           // FLAG, id for RNAME + 1, id for CT==1/GA==2, POS + 1, read length
  int nrecs = 0, nreads = 0, npaired = 0, ntempls = 0;                          // counters: BAM records, quality reads, paired reads, mergeable templates
  qname.reserve(10000); seq.reserve(10000); xm.reserve(10000);
  flag.reserve(10000); rname.reserve(10000); strand.reserve(10000);
  start.reserve(10000); width.reserve(10000);

  std::string rec_seq_rs, rec_xm_rs;                                            // seq and XM strings in reference space
  rec_seq_rs.reserve(255);
  rec_xm_rs.reserve(255);
  
  // process alignments
  while( sam_read1(bam_fp, bam_hdr, bam_rec) > 0 ) {                            // rec by rec
    nrecs++;                                                                    // BAM alignment records ++
    if ((nrecs & 1048575) == 0) Rcpp::checkUserInterrupt();                     // checking for the interrupt
    if ((bam_rec->core.qual < min_mapq) ||                                      // skip if mapping quality < min.mapq
        (skip_duplicates && (bam_rec->core.flag & BAM_FDUP))) continue;         // or if record is an optical/PCR duplicate
    nreads++;                                                                   // quality reads ++
    
    // reverse qname to make sort/search/compare faster
    char *rec_qname = bam_get_qname(bam_rec);                                   // QNAME
    std::reverse(rec_qname, rec_qname + bam_rec->core.l_qname                   // reversing it
                   - bam_rec->core.l_extranul - 1);
    // genome strand
    int rec_strand = ( bam_aux_get(bam_rec, "XG")[1] == 'C' ) ? 1 : 2 ;         // 1 if "ZCT"/"+", 2 if "ZGA"/"-"
    
    // append to main containers
    flag.push_back(bam_rec->core.flag);
    qname.push_back(rec_qname);                                                 // read QNAME
    rname.push_back(bam_rec->core.tid + 1);                                     // id of read RNAME + 1
    strand.push_back(rec_strand);                                               // id of reference strand
    start.push_back(bam_rec->core.pos + 1);                                     // read POS + 1
    
    // prepare seq and XM strings
    uint32_t read_len = bam_rec->core.l_qseq;                                   // read length
    uint8_t *rec_qual = bam_get_qual(bam_rec);                                  // quality string (Phred scale with no +33 offset)
    char *rec_xm = (char*) bam_aux_get(bam_rec, "XM") + 1;                      // methylation string without leading 'Z'
    if (rec_xm==NULL) Rcpp::stop("Methylation string is absent/incomplete");    // fall back if absent 
    uint8_t *rec_pseq = bam_get_seq(bam_rec);                                   // packed sequence string (4 bit per base)
    char rec_seq [read_len];                                   // unpacked sequence string
    for (int i=0; i<read_len; i++) {                                            // iterate over 4bit nucleotides
      if (rec_qual[i]<min_baseq) {                                              // if base quality < min.baseq
        rec_seq[i] = 'N';                                                       // sequence char to 'N'
        rec_xm[i] = '-';                                                        // methylation char to '-'
      } else {                                                                  // if good quality
        rec_seq[i] = seq_nt16_str[bam_seqi(rec_pseq,i)];                        // turn nucleotide into char (array is "=ACMGRSVTWYHKDBN")
      }
    }
    
    // apply CIGAR
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
    
    
    // trim leftmost READ2, post-trimming of rightmost READ2 is done elsewhere
    if ((bam_rec->core.flag & BAM_FPROPER_PAIR) &&                              // if read is a proper pair
        (bam_rec->core.flag & BAM_FREAD2)) {                                    // and it's a Read2
      int pos_diff = bam_rec->core.mpos - bam_rec->core.pos;                    // get distance between start positions of mates
      if (pos_diff == 0) {                                                      // if mates fully overlap
        rec_seq_rs.clear();                                                     // then READ2 is erased
        rec_xm_rs.clear();
      } else if ((pos_diff > 0) && (pos_diff < rec_seq_rs.length())) {          // else if it's a leftmost read && distance between mates is less than refspaced read length
        rec_seq_rs.erase(pos_diff);                                             // trim leftmost Read2 until rightmost Read1
        rec_xm_rs.erase(pos_diff);
      }
    }
    
    // append to main containers
    width.push_back(rec_seq_rs.length());                                       // read length
    seq.push_back(rec_seq_rs);                                                  // read SEQ
    xm.push_back(rec_xm_rs);                                                    // read XM
    
    // clearing temporary vars
    rec_seq_rs.clear();
    rec_xm_rs.clear();
  }

  // cleaning
  bam_destroy1(bam_rec);                                                        // clean BAM alignment structure 
  hts_close(bam_fp);                                                            // close BAM file
  
  // wrap and return the results
  std::vector<std::string> chromosomes (                                        // vector of reference names
      bam_hdr->target_name, bam_hdr->target_name + bam_hdr->n_targets);
  std::vector<std::string> strands = {"+", "-"};                                // vector of DNA strands
  std::vector<int> stats = {nrecs, nreads, npaired, ntempls};                   // statistics
  
  Rcpp::List res = Rcpp::List::create(                                          // final list
    Rcpp::Named("stats") = stats,                                               // BAM loading statistics
    Rcpp::Named("chromosomes") = chromosomes,                                   // lookup table of reference names
    Rcpp::Named("strands") = strands,                                           // lookup table of DNA strands
    Rcpp::Named("flag") = flag,                                                 // flag
    Rcpp::Named("qname") = qname,                                               // query names
    Rcpp::Named("rname") = rname,                                               // numeric ids for reference names
    Rcpp::Named("strand") = strand,                                             // numeric ids for reference strands
    Rcpp::Named("start") = start,                                               // start positions of reads
    Rcpp::Named("width") = width,                                               // lengths of reads
    Rcpp::Named("seq") = seq,                                                   // sequences
    Rcpp::Named("XM") = xm                                                      // methylation strings
  );
  return(res);
}



// In place post-trimming of rightmost READ2
//

// [[Rcpp::export]]
void rcpp_posttrim_read2(Rcpp::DataFrame &df) {
  Rcpp::IntegerVector flag = df["flag"];
  Rcpp::IntegerVector start = df["start"];
  Rcpp::CharacterVector qname = df["qname"];
  Rcpp::CharacterVector seq = df["seq"];
  Rcpp::CharacterVector xm = df["XM"];
  
  // row by row
  for (size_t x=1; x<flag.size(); x++) {
    if ((x & 1048575) == 0) Rcpp::checkUserInterrupt();                         // checking for the interrupt
    
    // trim rightmost READ2
    if ((qname[x]==qname[x-1]) &&                                    // if this is a rightmost read (df is presorted on QNAME, POS)
        (flag[x] & BAM_FPROPER_PAIR) &&                              // and this read is a proper pair
        (flag[x] & BAM_FREAD2)) {                                    // and it's a READ2
      int dpos = start[x] - start[x-1];                    // get distance between start positions of mates
      int ndel = xm[x-1].size() - dpos;                    // number of chars to delete from the left side of READ2
      if ((dpos > 0) && (ndel > 0)) {              // ??? if it's a rightmost read && distance between mates is less than refspaced read length
        seq[x] = Rcpp::as<std::string>(seq[x]).erase(0, ndel);                                             // trim leftmost Read2 until rightmost Read1
        xm[x] = Rcpp::as<std::string>(xm[x]).erase(0, ndel);
        start[x] = start[x-1] + xm[x-1].size();
      }
    }
  }
}


// #############################################################################
// test code and sourcing don't work on OS X
/*** R
*/
// #############################################################################