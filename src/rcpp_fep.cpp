#include <Rcpp.h>
#include <htslib/kfunc.h>

// [[Rcpp::depends(Rhtslib)]]

// Computes Fisher Exact P using HTSlib's implementation

// [[Rcpp::export]]
std::vector<double> rcpp_fep (Rcpp::DataFrame &df,                              // data.table by reference with the following columns:
                              std::vector<std::string> colnames)                // four strings with column names
                              
{
  Rcpp::IntegerVector A = df[colnames[0]];                                      // A
  Rcpp::IntegerVector B = df[colnames[1]];                                      // B
  Rcpp::IntegerVector C = df[colnames[2]];                                      // C
  Rcpp::IntegerVector D = df[colnames[3]];                                      // D
  
  std::vector<double> p(A.size(), NA_REAL);
  double left, right, two;
  
  for (unsigned int x=0; x<A.size(); x++) {
    // checking for the interrupt
    if ((x & 0xFFFFF) == 0) Rcpp::checkUserInterrupt();
    
    if (!Rcpp::IntegerVector::is_na(A[x]) &&
        !Rcpp::IntegerVector::is_na(B[x]) &&
        !Rcpp::IntegerVector::is_na(C[x]) &&
        !Rcpp::IntegerVector::is_na(D[x])) {
      kt_fisher_exact(A[x], B[x], C[x], D[x], &left, &right, &two);
      p[x] = two;
    }
  }
  
  return p;
}


// test code in R
//

/*** R
d <- data.table::data.table(matrix(c(NA, 1:8095), ncol=4))
n <- c("V1", "V2", "V3", "V4");
system.time( p <- rcpp_fep(d, n) )
system.time( f <- apply(d, 1, function (x) {if (any(is.na(x))) NA else stats::fisher.test(matrix(x, nrow=2))$p.value}) )
max(abs(p-f), na.rm=TRUE)
# microbenchmark::microbenchmark(rcpp_fep(d, n), times=10)
*/

// Sourcing:
// Rcpp::sourceCpp("rcpp_fep.cpp")

// #############################################################################