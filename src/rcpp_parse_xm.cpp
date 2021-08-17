#include <Rcpp.h>
#include <htslib/sam.h>

// Trims rightmost READ2 methylation strings passed by reference.
// Pairs and parses XM tags, outputs number of bases in hHxXzZuU contexts.
//

// [[Rcpp::export("rcpp_parse_xm")]]
std::vector<int> rcpp_parse_xm(std::vector<std::string> qname,                  // QNAME strings
                               std::vector<std::string> xm)                     // XM strings
{
  int ascii_map [128] = {0};                                                    // methylation call counters
  std::vector<int> res (qname.size() * 8, 0);                                   // result array
  size_t x = 0, x0 = 0;                                                         // read counter, and it's value for the first read of consecutive paired reads
  unsigned char calls [10] = {'H','U','X','Z','h','u','x','z','-','.'};         // what chars to copy to the results and/or clean after copying
  
  // macros
  #define spit_n_clean {                   /* save results, clean ascii_map */ \
    for (size_t y=x0; y<x; y++) {    /* do it for all reads in the template */ \
      size_t idx=(y<<3);                            /* 8 ints for each read */ \
      for(size_t i=0; i<8; i++) {         /* for first eight possible chars */ \
        res[idx+i] = ascii_map[calls[i]];                         /* assign */ \
      }                                                                        \
    }                                                                          \
    for(size_t i=0; i<10; i++) {              /* for all ten possible chars */ \
      ascii_map[calls[i]] = 0;                                     /* clear */ \
    }                                                                          \
  }

  // str by str
  for (x=0; x<xm.size(); x++) {
    if ((x & 1048575) == 0) Rcpp::checkUserInterrupt();                         // checking for the interrupt

    // compare QNAMEs if x>0
    if ((x>0) && (qname[x]!=qname[x-1])) {                                      // if not the first line and new QNAME
      spit_n_clean;                                                             // save previously accumulated counts
      x0 = x;                                                                   // redefine start of consecutive read group
    }
    
    // accumulate counts
    size_t l = xm[x].length();                                                  // XM length
    const char * s = xm[x].c_str();                                             // pointer to XM c-string
    for (size_t i=0; i<l; i++) {                                                // for all chars
      ascii_map[(unsigned char)s[i]]++;                                         // count them
    }
  }
  spit_n_clean;                                                                 // save counts for the final consecutive read group

  return res;
}


