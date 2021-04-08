test_generateBedEcdf <- function () {
  amplicon.ecdfs <- generateBedEcdf(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=system.file("extdata", "amplicon.bed", package="epialleleR"),
    bed.rows=NULL, verbose=FALSE
  )
  
  RUnit::checkEquals(
    sapply(unlist(amplicon.ecdfs, use.names=FALSE), function (x) {x(0.5)}),
    c(0.916666666667, 1, 0.885245901639, 1, 0.946236559140, 1,
      0.892857142857, 1, 0.877358490566, 0.952830188679),
    tolerance=1e-08
  )
}
