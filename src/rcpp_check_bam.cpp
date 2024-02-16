#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

// [[Rcpp::depends(Rhtslib)]]

// Checks BAM before reading for:
// [+] being name-sorted
// [+] counting all AUX tags, with special interest in:
//     [+] XM: Illumina/Bismark methylation calls
//     [+] XG: Illumina/Bismark genome strand (CT or GA)
//     [+] YD: bwa-meth genome strand (CT or GA)
//     [+] ZS: bsmap genome strand (++, +-, -+, -- as <reference><read>)
//     [+] MM: BAM base modification tag
//
// Returns list with simple counts (see *main counters* and *wrap and return*)

// [[Rcpp::export]]
Rcpp::List rcpp_check_bam (std::string fn)                                      // file name
{
  // constants
  const unsigned int max_recs = 1024;                                           // max number of reads to check
  const unsigned int max_qname_width = 1024;                                    // max QNAME length, not expanded yet, ever error-prone?

  // file IO
  htsFile *bam_fp = hts_open(fn.c_str(), "r");                                  // try open file
  if (!bam_fp) Rcpp::stop("Unable to open BAM file for reading");               // fall back if error
  bam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);                                    // try read file header
  if (!bam_hdr) Rcpp::stop("Unable to read BAM header");                        // fall back if error  
  bam1_t *bam_rec = bam_init1();                                                // create BAM alignment structure
  
  // main counters
  std::map<std::string,unsigned int> aux_map = {{"nrecs",0},                    /* BAM records */                   \
                                                {"npaired",0},                  /* have BAM_FPROPER_PAIR flag */    \
                                                {"ntempls",0}};                 /* templates (consecutive pairs) */
  // template holders
  char *templ_qname = (char*) malloc(max_qname_width * sizeof(char));           // template QNAME
  
  // process alignments
  while( (sam_read1(bam_fp, bam_hdr, bam_rec) > 0) && (aux_map["nrecs"] < max_recs) ) { // rec by rec until max_recs or EOF
    aux_map["nrecs"]++;                                                         // BAM alignment records ++
    if (bam_rec->core.flag & BAM_FPROPER_PAIR) aux_map["npaired"]++;            // if is BAM_FPROPER_PAIR
    for (uint8_t *aux = bam_aux_first(bam_rec); aux; aux = bam_aux_next(bam_rec, aux)) { // cycle through all AUX fields
      aux_map[std::string(bam_aux_tag(aux), 2)]++;                              // try increment/emplace current tag
    }
    
    // check if not the same template (QNAME)
    if (strcmp(templ_qname, bam_get_qname(bam_rec)) == 0) aux_map["ntempls"]++; // if the same template (QNAME)
    strcpy(templ_qname, bam_get_qname(bam_rec));                                // store template QNAME
  }
  
  // cleaning
  bam_destroy1(bam_rec);                                                        // clean BAM alignment structure 
  sam_hdr_destroy(bam_hdr);                                                     // free header (I don't do it when rcpp_read_bam!!!)
  hts_close(bam_fp);                                                            // close BAM file
  free(templ_qname);                                                            // and free manually allocated memory
  
  // wrap and return the results
  return(Rcpp::wrap(aux_map));
}


// #############################################################################
// test code and sourcing don't work on OS X
/*** R
*/
// #############################################################################