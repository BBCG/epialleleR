#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

// [[Rcpp::depends(Rhtslib)]]

// Checks BAM before reading for:
// [+] being name-sorted
// [+] having XM: Illumina/Bismark methylation calls
// [+] having XG: Illumina/Bismark genome strand (CT or GA)
// [+] having YC: bwa-meth genome strand (CT or GA)
// [+] having MM: BAM base modification tag
//
// Returns list with simple counts (see *main counters* and *wrap and return*)

// [[Rcpp::export]]
Rcpp::List rcpp_check_bam (std::string fn)                                      // file name
{
  // constants
  int max_recs = 1024;                                                          // max number of reads to check
  int max_qname_width = 1024;                                                   // max QNAME length, not expanded yet, ever error-prone?

  // file IO
  htsFile *bam_fp = hts_open(fn.c_str(), "r");                                  // try open file
  if (bam_fp==NULL) Rcpp::stop("Unable to open BAM file for reading");          // fall back if error
  bam_hdr_t *bam_hdr = sam_hdr_read(bam_fp);                                    // try read file header
  if (bam_hdr==NULL) Rcpp::stop("Unable to read BAM header");                   // fall back if error  
  bam1_t *bam_rec = bam_init1();                                                // create BAM alignment structure
  
  // main counters
  int nrecs = 0,                                                                // BAM records, templates (consecutive pairs)
    npp = 0,                                                                    // have BAM_FPROPER_PAIR flag
    ntempls = 0,                                                                // templates (consecutive pairs)
    nxm = 0,                                                                    // have XM tags
    nxg = 0,                                                                    // have XG tags
    nyc = 0,                                                                    // have YC tags
    nmm = 0;                                                                    // have MM tags
  
  // template holders
  char *templ_qname = (char*) malloc(max_qname_width * sizeof(char));           // template QNAME
  
// process alignments
while( (sam_read1(bam_fp, bam_hdr, bam_rec) > 0) & (nrecs < max_recs) ) {      // rec by rec until max_recs or EOF
  nrecs++;                                                                      // BAM alignment records ++
  if (bam_rec->core.flag & BAM_FPROPER_PAIR) npp++;                             // if is BAM_FPROPER_PAIR
  if (bam_aux_get(bam_rec, "XM") != NULL) nxm++;                                // if has XM tag
  if (bam_aux_get(bam_rec, "XG") != NULL) nxg++;                                // if has XG tag
  if (bam_aux_get(bam_rec, "YC") != NULL) nyc++;                                // if has YC tag
  if (bam_aux_get(bam_rec, "MM") != NULL) nmm++;                                // if has MM tag
  
  // check if not the same template (QNAME)
  if (strcmp(templ_qname, bam_get_qname(bam_rec)) == 0) ntempls++;              // if the same template (QNAME)
  strcpy(templ_qname, bam_get_qname(bam_rec));                                  // store template QNAME
}

// cleaning
bam_destroy1(bam_rec);                                                          // clean BAM alignment structure 
sam_hdr_destroy(bam_hdr);                                                       // free header (I don't do it when rcpp_read_bam!!!)
hts_close(bam_fp);                                                              // close BAM file
free(templ_qname);                                                              // and free manually allocated memory

// wrap and return the results
Rcpp::List res = Rcpp::List::create(                                            // final List
  Rcpp::Named("nrecs") = nrecs,                                                 // numeric for number of records
  Rcpp::Named("npp") = npp,                                                     // numeric for number of proper paired reads
  Rcpp::Named("ntempls") = ntempls,                                             // numeric for number of consecutively paired reads
  Rcpp::Named("nxm") = nxm,                                                     // numeric for number of reads with XM tags
  Rcpp::Named("nxg") = nxg,                                                     // numeric for number of reads with XG tags
  Rcpp::Named("nyc") = nyc,                                                     // numeric for number of reads with YC tags
  Rcpp::Named("nmm") = nmm                                                      // numeric for number of reads with MM tags
);

return(res);
}


// #############################################################################
// test code and sourcing don't work on OS X
/*** R
*/
// #############################################################################