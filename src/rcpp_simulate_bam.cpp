#include <Rcpp.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(Rhtslib)]]

// Creates small BAM files given necessary fields
// and tags, all supplied as data.frames.
//
// Returns simple statistics on records written.


// template function to convert R's Rcpp::NumericVector to any numeric type
// and add values to a BAM array tag
template<typename T>
int save_array_tag(bam1_t *b, const char tag[2], uint8_t type, const Rcpp::NumericVector &l) {
  std::vector<T> vec (l.size());                                                // new container
  for (auto i=0; i<l.size(); i++) { vec[i] = static_cast<T>(l[i]); }            // copy tag values
  int res = bam_aux_update_array(b, tag, type, vec.size(), vec.data());         // add array tag
  return res;
}

// [[Rcpp::export]]
int rcpp_simulate_bam (std::vector<std::string> header,                         // header records
                       Rcpp::DataFrame &fields,                                 // mandatory fields
                       Rcpp::DataFrame &i_tags,                                 // optional integer tags
                       Rcpp::DataFrame &f_tags,                                 // optional float tags
                       Rcpp::DataFrame &s_tags,                                 // optional string tags
                       Rcpp::DataFrame &a_tags,                                 // optional array tags
                       std::vector<std::string> a_types,                        // types ('cCsSiIf') of optional array tags
                       std::string out_fn)                                      // output BAM file name
{
  // main
  std::vector<std::string> qname = Rcpp::as<std::vector<std::string>>(fields["qname"]);   // Query name
  std::vector<uint16_t> flag = Rcpp::as<std::vector<uint16_t>>(fields["flag"]);           // Bitwise flag, a combination of the BAM_F* constants
  std::vector<int32_t> tid = Rcpp::as<std::vector<int32_t>>(fields["tid"]);               // Chromosome ID
  std::vector<hts_pos_t> pos = Rcpp::as<std::vector<hts_pos_t>>(fields["pos"]);           // 0-based leftmost coordinate
  std::vector<uint8_t> mapq = Rcpp::as<std::vector<uint8_t>>(fields["mapq"]);             // Mapping quality
  std::vector<std::string> cigar = Rcpp::as<std::vector<std::string>>(fields["cigar"]);   // CIGAR data
  std::vector<int32_t> mtid = Rcpp::as<std::vector<int32_t>>(fields["mtid"]);             // Chromosome ID of next read in template
  std::vector<hts_pos_t> mpos = Rcpp::as<std::vector<hts_pos_t>>(fields["mpos"]);         // 0-based leftmost coordinate of next read in template
  std::vector<hts_pos_t> isize = Rcpp::as<std::vector<hts_pos_t>>(fields["isize"]);       // Observed template length
  std::vector<std::string> seq = Rcpp::as<std::vector<std::string>>(fields["seq"]);       // Sequence
  std::vector<std::string> qual = Rcpp::as<std::vector<std::string>>(fields["qual"]);     // Sequence quality
  
  std::vector<std::string> i_cols = Rcpp::as<std::vector<std::string>>(i_tags.names());   // Column names for integer tags
  std::vector<std::string> f_cols = Rcpp::as<std::vector<std::string>>(f_tags.names());   // Column names for float tags
  std::vector<std::string> s_cols = Rcpp::as<std::vector<std::string>>(s_tags.names());   // Column names for string tags
  std::vector<std::string> a_cols = Rcpp::as<std::vector<std::string>>(a_tags.names());   // Column names for array tags
  
  // file IO
  htsFile *out_fp = hts_open(out_fn.c_str(), "wb");                             // try open output file
  if (!out_fp) Rcpp::stop("Unable to open output BAM file for writing");        // fall back if error
  
  // header
  bam_hdr_t *out_hdr = sam_hdr_init();                                          // try init file header
  if (!out_hdr) Rcpp::stop("Unable to init BAM header");                        // fall back if error
  for (size_t i=0; i<header.size(); i++)
    sam_hdr_add_lines(out_hdr, header[i].c_str(), 0);                           // fill the header
  if (sam_hdr_write(out_fp, out_hdr) < 0) Rcpp::stop("Unable to write header"); // try write output file header
  
  // filling and writing out records
  int call_res;                                                                 // result of HTSlib call
  bam1_t *out_rec = bam_init1();                                                // create BAM alignment structure
  uint32_t *cigar_mem = NULL;                                                   // destination uint32_t CIGAR buffer
  size_t n_cigar = 0;                                                           // allocated number of CIGAR buffer elements
  for (size_t i=0; i<qname.size(); i++) {
    call_res = sam_parse_cigar(cigar[i].c_str(), NULL, &cigar_mem, &n_cigar);   // fill CIGAR array
    if (call_res<0) Rcpp::stop("Unable to fill CIGAR array");                   // fall back on error
    
    for (char *q=(char*)qual[i].c_str(); *q!='\0'; q++) *q -= 33;               // remove offset from QUAL string
      
    call_res = bam_set1(                                                        // set all mandatory fields
      out_rec, qname[i].size(), qname[i].c_str(),
      flag[i], tid[i], pos[i], mapq[i],
      n_cigar, cigar_mem,
      mtid[i], mpos[i], isize[i],
      seq[i].size(), seq[i].c_str(), qual[i].c_str(),
      0
    );
    if (call_res<0) Rcpp::stop("Unable to fill BAM record");                    // fall back on error
    
    for (size_t c=0; c<i_cols.size(); c++) {                                    // add integer tags
      int tag = ((Rcpp::IntegerVector)(i_tags[c]))[i];                          // tag value
      call_res = bam_aux_update_int(out_rec, i_cols[c].c_str(), tag);           // add tag
      if (call_res<0) Rcpp::stop("Unable to add %s tag", i_cols[c]);            // fall back on error
    }
    
    for (size_t c=0; c<f_cols.size(); c++) {                                    // add float tags
      float tag = ((Rcpp::NumericVector)(f_tags[c]))[i];                        // tag value
      call_res = bam_aux_update_float(out_rec, f_cols[c].c_str(), tag);         // add tag
      if (call_res<0) Rcpp::stop("Unable to add %s tag", f_cols[c]);            // fall back on error
    }
    
    for (size_t c=0; c<s_cols.size(); c++) {                                    // add string tags
      std::string tag = Rcpp::as<std::string>( ((Rcpp::StringVector)(s_tags[c]))[i] );     // tag value
      call_res = bam_aux_update_str(out_rec, s_cols[c].c_str(), tag.size(), tag.c_str());  // add tag
      if (call_res<0) Rcpp::stop("Unable to add %s tag", s_cols[c]);            // fall back on error
    }
    
    for (size_t c=0; c<a_cols.size(); c++) {                                    // add array tags
      switch (a_types[c][0]) {
      case 'f':                                                                 // float array
        call_res = save_array_tag<float>(out_rec, a_cols[c].c_str(), a_types[c][0], ((Rcpp::List)(a_tags[c]))[i]);
        break;
      case 'c':                                                                 // int8_t array
        call_res = save_array_tag<int8_t>(out_rec, a_cols[c].c_str(), a_types[c][0], ((Rcpp::List)(a_tags[c]))[i]);
        break;
      case 'C':                                                                 // uint8_t array
        call_res = save_array_tag<uint8_t>(out_rec, a_cols[c].c_str(), a_types[c][0], ((Rcpp::List)(a_tags[c]))[i]);
        break;
      case 's':                                                                 // int16_t array
        call_res = save_array_tag<int16_t>(out_rec, a_cols[c].c_str(), a_types[c][0], ((Rcpp::List)(a_tags[c]))[i]);
        break;
      case 'S':                                                                 // uint16_t array
        call_res = save_array_tag<uint16_t>(out_rec, a_cols[c].c_str(), a_types[c][0], ((Rcpp::List)(a_tags[c]))[i]);
        break;
      case 'i':                                                                 // int32_t array
        call_res = save_array_tag<int32_t>(out_rec, a_cols[c].c_str(), a_types[c][0], ((Rcpp::List)(a_tags[c]))[i]);
        break;
      case 'I':                                                                 // uint32_t array
        call_res = save_array_tag<uint32_t>(out_rec, a_cols[c].c_str(), a_types[c][0], ((Rcpp::List)(a_tags[c]))[i]);
        break;
      }
      
      if (call_res<0) Rcpp::stop("Unable to add %s tag", a_cols[c]);            // fall back on error
    }
    
    call_res = sam_write1(out_fp, out_hdr, out_rec);                            // write record
    if (call_res<0) Rcpp::stop("Unable to write BAM record");                   // fall back on error
  }
  
  // cleaning
  if (cigar_mem) free(cigar_mem);                                               // CIGAR array
  bam_hdr_destroy(out_hdr);                                                     // clean BAM header
  bam_destroy1(out_rec);                                                        // clean BAM alignment structure
  hts_close(out_fp);                                                            // close output BAM file
  
  return qname.size();
}


// #############################################################################
// test code and sourcing don't work on OS X
/*** R
*/
// #############################################################################
