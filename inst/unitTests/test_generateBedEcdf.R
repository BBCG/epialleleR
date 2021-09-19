test_generateBedEcdf <- function () {
  amplicon.ecdfs <- generateBedEcdf(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=system.file("extdata", "amplicon.bed", package="epialleleR"),
    bed.rows=c(1,2), verbose=TRUE
  )
  
  RUnit::checkEquals(
    sapply(unlist(amplicon.ecdfs, use.names=FALSE), function (x) {x(0.5)}),
    c(0.916666666667, 1, 0.885245901639, 1),
    tolerance=1e-08
  )
  
  amplicon.ecdfs <- generateBedEcdf(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=system.file("extdata", "amplicon.bed", package="epialleleR"),
    min.mapq=30, min.baseq=20,
    bed.rows=NULL, verbose=TRUE
  )
  
  RUnit::checkEquals(
    sapply(unlist(amplicon.ecdfs, use.names=FALSE), function (x) {x(0.5)}),
    c(0.916666666667, 1, 0.885245901639, 1, 0.946236559140, 1,
      0.892857142857, 1, 0.868131868132, 1),
    tolerance=1e-08
  )
}
