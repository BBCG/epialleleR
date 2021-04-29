test_internals <- function () {
  RUnit::checkException(
    epialleleR:::.processBam(list("bam"=list()), FALSE)
  )
  
  RUnit::checkEquals(
    epialleleR:::.processBam(
      list(list(qname="a", flag=1, rname="n", strand="+", pos=1, cigar="1M",
                seq="A", tag=list(XM=".", XG="+"))), FALSE),
    data.table::data.table(qname="a", rname="n", strand="+",
                           start=1, seq="A", XM=".", width=1)
  )
  
  RUnit::checkException(
    epialleleR:::.processBam(
      list(list(qname="a", flag=1, rname="n", strand="+", pos=1, cigar="1M",
                seq="A", tag=list(XM=".", XG="+")), mpos=2, isize=1), FALSE)
  )
  
  RUnit::checkEquals(
    epialleleR:::rcpp_apply_cigar("2S4M4=2X1D4I1N6M1H1P", "abcdefghijklmnopqrstuvwxyz", "?"),
    "cdefghijkl??qrstuvwxyz"
    # as.character(
    #   GenomicAlignments::sequenceLayer(
    #     Biostrings::BStringSet("abcdefghijklmnopqrstuvwxyz"),
    #     "2S4M4=2X1D4I1N6M1H1P", D.letter="?", N.letter="?"
    #   )
    # )
  )
  
  RUnit::checkException(
    epialleleR:::rcpp_apply_cigar("2S4M4=2X1D4I1N6Z1H1P", "abcdefghijklmnopqrstuvwxyz", "?")
  )
}
