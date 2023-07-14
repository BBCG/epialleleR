#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <boost/range/algorithm.hpp>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(Rhtslib)]]

// Creates small BAM files given necessary fields and tags supplied as Lists
//
// Returns simple statistics on records written.


// Checks if c is a cigar operation by basically checking if it's not a number
// doesn't check validity of cigar_op (MIDNSHP=X)
bool is_op (char c) {return c>57;}


// [[Rcpp::export]]
int rcpp_simulate_bam (std::vector<std::string> header,                         // header records
                       Rcpp::List &fields,                                      // mandatory fields
                       Rcpp::List &i_tags,                                      // optional integer tags
                       Rcpp::List &s_tags,                                      // optional string tags
                       std::string out_fn)                                      // output BAM file name
{
  size_t nrecs = ((Rcpp::CharacterVector)fields["qname"]).length();             // number of records
  
  // file IO
  htsFile *out_fp = hts_open(out_fn.c_str(), "wb");                             // try open output file
  if (out_fp==NULL) Rcpp::stop("Unable to open output BAM file for writing");   // fall back if error
  
  // header
  bam_hdr_t *out_hdr = sam_hdr_init();                                          // try init file header
  if (out_hdr==NULL) Rcpp::stop("Unable to init BAM header");                   // fall back if error
  for (size_t i=0; i<header.size(); i++)
    sam_hdr_add_lines(out_hdr, header[i].c_str(), 0);                           // fill the header
  if (sam_hdr_write(out_fp, out_hdr) < 0) Rcpp::stop("Unable to write header"); // try write output file header
  
  // filling and writing out records
  int call_res;                                                                 // result of HTSlib call
  bam1_t *out_rec = bam_init1();                                                // create BAM alignment structure
  for (size_t i=0; i<header.size(); i++) {
    // size_t l_qname = ((Rcpp::CharacterVector)fields["qname"])[i].length();      // Length of the query name
    // size_t n_cigar = boost::range::count_if(                                    // Number of CIGAR operations
    //   ((Rcpp::CharacterVector)fields["cigar"])[i], is_op
    // );
    // size_t l_seq = ((Rcpp::CharacterVector)fields["seq"])[i].length();          // Length of the query sequence and sequence quality string
    // size_t l_aux = i_tags.size() + s_tags.size();                               // Length to be reserved for auxiliary field data
    // 
    // call_res = bam_set1(                                                        // set all mandatory fields:
    //   out_rec,                                                                  // BAM alignment structure
    //   l_qname,                                                                  // 
    //   ((Rcpp::CharacterVector)fields["qname"])[i].c_str(),                      // 
    //   ((Rcpp::IntegerVector)fields["flag"])[i],                                 // 
    //   ((Rcpp::IntegerVector)fields["tid"])[i],                                  // 
    //   ((Rcpp::IntegerVector)fields["pos"])[i],                                  // 
    //   ((Rcpp::IntegerVector)fields["mapq"])[i],                                 // 
    //   n_cigar,
    //   
    // );
    // if (call_res<0) Rcpp::stop("Unable to fill BAM record");                    // fall back on error
    
    
    
    call_res = sam_write1(out_fp, out_hdr, out_rec);                            // write record
    if (call_res<0) Rcpp::stop("Unable to write BAM record");                   // fall back on error
  }
  
  // cleaning
  bam_hdr_destroy(out_hdr);                                                     // clean BAM header
  bam_destroy1(out_rec);                                                        // clean BAM alignment structure
  hts_close(out_fp);                                                            // close output BAM file
  
  return nrecs;
}


// #############################################################################
// test code and sourcing don't work on OS X
/*** R
file <- tempfile(pattern="simulated" ,fileext=".bam")
rcpp_simulate_bam(
  c("@SQ\tSN:seq1\tLN:1575", "@SQ\tSN:seq2\tLN:1584", "@PG\tID:epialleleR\tPN:epialleleR\tVN:1.9.5\tCL:rcpp_simulate_bam()"),
  list(
    qname=paste0("r", 1:2),
    flag=c(4,4),
    tid=c(1,1),
    pos=c(1,5),
    mapq=c(60,60),
    cigar=c("M10","M10"),
    mtid=c(0,0),
    mpos=c(0,0),
    isize=c(10,10),
    seq=c("ACTGACTGAC","TGACTGACTG"),
    qual=c("AAAAAAAAAA","AAAAAAAAAA")
  ),
  list(),
  list(
    XG=c("CT","CT"),
    XM=c("ZZzzZZzzZZ","ZZzzZZzzZZ")
  ),
  file)

*/
// #############################################################################