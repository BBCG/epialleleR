# preprocessBam

This function reads and preprocesses BAM file.

## Usage

``` r
preprocessBam(
  bam.file,
  paired = NULL,
  override.check = FALSE,
  min.mapq = 0,
  min.baseq = 0,
  min.prob = -1,
  highest.prob = TRUE,
  skip.duplicates = FALSE,
  skip.secondary = TRUE,
  skip.qcfail = TRUE,
  skip.supplementary = TRUE,
  trim = 0,
  targets = NULL,
  zero.based.targets = FALSE,
  clip.to.targets = FALSE,
  nthreads = 1,
  verbose = TRUE
)
```

## Arguments

- bam.file:

  BAM file location string.

- paired:

  boolean for expected alignment endness: \`TRUE\` for paired-end,
  \`FALSE\` for single-end, or \`NULL\` for auto detect (the default).

- override.check:

  boolean to use supplied endness (\`paired\` parameter) even if it is
  different from the autodetected one (default: FALSE).

- min.mapq:

  non-negative integer threshold for minimum read mapping quality.
  Default is 0, however a higher value (e.g., 30) is recommended.

- min.baseq:

  non-negative integer threshold for minimum nucleotide base quality.
  Default is 0, however a higher value (e.g., at least 13 as in
  \`samtools mpileup\`) is recommended.

- min.prob:

  integer threshold for minimum scaled probability of modification
  (methylation) to consider. Affects processing of long-read sequencing
  alignments only. According to SAM/BAM specification, the continuous
  base modification probability range 0.0 to 1.0 is remapped in equal
  sized portions to the discrete integers 0 to 255 inclusively. Default
  is -1, however a higher value (e.g., 178 which corresponds to a
  modification probability of ~0.7) is strongly recommended. When
  default (-1), all C+m and G-m cytosine methylation modifications
  recorded in MM/Mm tag will be included, even if ML/Ml tag with
  probabilities is absent (in such case, probability of modification
  equals -1). If a base has at least one modifications, and neither the
  probability of any of modifications nor the probability of the
  conventional base reach the \`min.prob\`, this base will be masked in
  this particular read and not counted neither as modified nor as
  conventional.

- highest.prob:

  boolean defining if methylation modification must have the highest
  probability among all modifications at a particular base to be
  considered in the analyses (default: TRUE). Affects processing
  long-read sequencing alignments only. If default (TRUE) and ML/Ml tag
  with probability scores is absent, then cytosines with more than one
  modification will be omitted (as the probability of all modifications
  will be equal).

- skip.duplicates:

  boolean defining if duplicate aligned reads should be skipped
  (default: FALSE). Option has no effect if duplicate reads were not
  marked by alignment software.

- skip.secondary:

  boolean defining if secondary alignments should be skipped (default:
  TRUE). Do not change.

- skip.qcfail:

  boolean defining if alignments failing QC should be skipped (default:
  TRUE). Do not change.

- skip.supplementary:

  boolean defining if supplementary alignments should be skipped
  (default: TRUE). Do not change.

- trim:

  non-negative integer or vector of length 2 for the number of
  nucleotide bases to be trimmed from 5' and 3' ends of a template
  (i.e., read pair for paired-end BAM or read for single-end BAM).
  Default: 0 for no trimming. Specifying \`trim=1\` will result in
  removing of a single base from both ends, while specifying
  \`trim=c(1,2)\` will result in removing of a single base from 5' end
  and 2 bases from 3' end.

- targets:

  Browser Extensible Data (BED) file location string OR object of class
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  holding genomic coordinates for regions of interest. The style of
  seqlevels of BED file/object must be the same as the style of
  seqlevels of BAM file. During targets loading, the genomic regions are
  reduced (i.e., overlapping regions are merged), and the strand
  information is discarded. These reduced regions are then used to only
  load overlapping reads (merged read pairs for paired-end sequencing).

- zero.based.targets:

  boolean defining if BED file coordinates are zero-based (default:
  FALSE).

- clip.to.targets:

  boolean defining if overlapping reads (merged read pairs in case of
  paired-end sequencing) should be clipped to retain only fragments
  overlapping with \`targets\` (default: FALSE). If a read overlaps with
  two or more targets, all overlapping fragments will be retained.

- nthreads:

  non-negative integer for the number of additional HTSlib threads to be
  used during BAM file decompression (default: 1). Two threads (and
  usually no more than four) make sense for files larger than 100 MB.

- verbose:

  boolean to report progress and timings (default: TRUE).

## Value

[`data.table`](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
object containing preprocessed BAM data.

NB: most of the BAM data is stored not in the data.table *per se*, but
as an external object linked to this data.table. Therefore,
saving/loading this data.table cannot be used to save/recover
preprocessed BAM data.

## Details

The function loads and preprocesses BAM file, saving time when multiple
analyses are to be performed on large input files. Currently, HTSlib is
used to read the data, therefore it is possible to speed up the loading
by means of HTSlib decompression threads.

This function is also called internally when BAM file location is
supplied as an input for other \`epialleleR\` methods.

\`preprocessBam\` currently allows to load both short-read (e.g.,
bisulfite) and long-read (native) sequencing alignments. Specific
requirements for these types of data are given below.

## Short-read sequencing

For preprocessing of short reads (and therefore for all reporting
methods), \`epialleleR\` requires genomic strand (XG tag) and a
methylation call string (XM tag) to be present in a BAM file - i.e.,
methylation calling must be performed after read mapping/alignment by
your software of choice. It is the case for BAM files produced by
Bismark Bisulfite Read Mapper and Methylation Caller, Illumina DRAGEN,
Illumina Cloud analysis solutions, as well as contemporary Illumina
sequencing instruments with on-board read mapping/alignment (NextSeq
1000/2000, NovaSeq X), therefore such files can be analysed without
additional steps. For alignments produced by other tools, e.g., BWA-meth
or BSMAP, methylation calling must be performed prior to BAM loading /
reporting, by means of [`callMethylation`](callMethylation.md).

## Long-read sequencing

For preprocessing of long reads, \`epialleleR\` requires presence of MM
(Mm) and ML (Ml) tags that hold information on base modifications and
related probabilities, respectively. These are standard tags described
in SAM/BAM format specification, therefore relevant tools for analysis
and alignment of long sequencing reads should be able to produce them.

## Other details

\`preprocessBam\` always tests if BAM file is paired- or single-ended
and has all the necessary tags available. It is recommended to use
\`verbose\` processing and check messages for correct identification of
alignment endness. Otherwise, if the \`paired\` parameter is set
explicitly, exception is thrown when expected endness differs from the
autodetected one. It is nevertheless possible to override the
autodetected endness and load BAM as specified in \`paired\` by setting
the \`override.check\` to TRUE.

During preprocessing of paired-end alignments, paired reads are merged
according to their base quality: nucleotide base with the highest value
in the QUAL string is taken, unless its quality is less than
\`min.baseq\`, which results in no information for that particular
position ("-"/"N"). These merged reads are then processed as a single
entity in all \`epialleleR\` methods. Due to merging, overlapping bases
in read pairs are counted only once, and the base with the highest
quality is taken.

It is currently a requirement that paired-end BAM file must be sorted by
QNAME instead of genomic location (i.e., "unsorted") to perform merging
of paired-end reads. Error message is shown if it is sorted by genomic
location, in this case please sort it by QNAME using 'samtools sort -n
-o out.bam in.bam'.

During preprocessing of single-end alignments, no read merging is
performed. Only bases with quality of at least \`min.baseq\` are
considered. Lower base quality results in no information for that
particular position ("-"/"N").

For RRBS-like protocols, it is possible to trim alignments from one or
both ends. Trimming is performed during BAM loading and will therefore
influence results of all downstream \`epialleleR\` methods. Internally,
trimming is performed at the level of a template (i.e., read pair for
paired-end BAM or individual read for single-end BAM). This ensures that
only necessary parts (real ends of sequenced fragment) are removed for
paired-end sequencing reads.

It is also possible to load only a subset of reads (read pairs) of
interest or only fragments of such reads by supplying a list of targets
(see description of \`targets\` and other related options). The
subsetting is performed without using BAM index.

## Specific considerations for long-read sequencing data

Any location not reported is implicitly assumed to contain no
modification.

According to SAM format specification, MM base modification tags are
allowed to list modifications observed not only on the original
sequenced strand (e.g., \`C+m\`) but also on the opposite strand (e.g.,
\`G-m\`). The logic of their processing is as follows (with the examples
given below):

- if an alignment record has no methylation modifications (neither
  \`C+m\`, nor \`G-m\` are present), this record is, naturally,
  considered to be a single read with no cytosines methylated

- if an alignment record has \`C+m\` modification (base modifications on
  the original sequenced strand), then this record is, naturally,
  considered to be a single read with cytosine modifications on the
  sequenced strand

- if an alignment record has \`G-m\` modification (base modifications on
  the strand opposite to sequenced), then this record is treated as two
  reads, with the original sequenced strand having no modifications,
  while the opposite strand having cytosine modifications

- if both \`C+m\` and \`G-m\` are present, then this record is treated
  as two reads, with both strands having cytosine modifications

## See also

[`preprocessGenome`](preprocessGenome.md) for preloading reference
sequences and [`callMethylation`](callMethylation.md) for methylation
calling.

[`generateCytosineReport`](generateCytosineReport.md) for methylation
statistics at the level of individual cytosines,
[`generateBedReport`](generateBedReport.md) for genomic region-based
statistics, [`generateVcfReport`](generateVcfReport.md) for evaluating
epiallele-SNV associations, [`extractPatterns`](extractPatterns.md) for
exploring methylation patterns and [`plotPatterns`](plotPatterns.md) for
pretty plotting of its output, [`generateBedEcdf`](generateBedEcdf.md)
for analysing the distribution of per-read beta values, and
\`epialleleR\` vignettes for the description of usage and sample data.

Sequence Alignment/Map [format
specifications](https://samtools.github.io/hts-specs/SAMv1.pdf),
specifications for [optional SAM
tags](https://samtools.github.io/hts-specs/SAMtags.pdf), duplicate
alignments marking by
[Samtools](http://www.htslib.org/doc/samtools-markdup.md) and [Illumina
DRAGEN Bio IT
Platform](https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/DuplicateMarking_fDG.htm).

## Examples

``` r
  # short-read sequencing
  capture.data <- preprocessBam(
    system.file("extdata", "capture.bam", package="epialleleR"),
    targets=as("chr17:43120000-43130000", "GRanges")
  )
#> Checking BAM file: 
#> short-read, paired-end, name-sorted alignment detected
#> Reading paired-end BAM file 
#> [0.013s]
  generateCytosineReport(capture.data, threshold.reads=TRUE)
#> Filtering and thresholding reads 
#> [0.000s]
#> Preparing cytosine report 
#> [0.001s]
#>      rname strand      pos context  meth unmeth
#>     <fctr> <fctr>    <int>  <fctr> <int>  <int>
#>  1:  chr17      - 43123785      CG     1      0
#>  2:  chr17      + 43124075      CG     1      0
#>  3:  chr17      + 43124077      CG     1      0
#>  4:  chr17      - 43124568      CG     0      1
#>  5:  chr17      - 43124590      CG     0      1
#>  6:  chr17      - 43124611      CG     0      1
#>  7:  chr17      - 43124625      CG     0      1
#>  8:  chr17      + 43125677      CG     0      1
#>  9:  chr17      + 43125690      CG     0      1
#> 10:  chr17      + 43125713      CG     0      1
#> 11:  chr17      + 43125745      CG     0      1
#> 12:  chr17      + 43125769      CG     0      1
#> 13:  chr17      + 43125830      CG     0      1
#> 14:  chr17      + 43125835      CG     0      1
#> 15:  chr17      + 43125851      CG     0      2
#> 16:  chr17      + 43125908      CG     0      2
#> 17:  chr17      + 43125923      CG     0      2
#> 18:  chr17      + 43125955      CG     0      2
#> 19:  chr17      + 43125957      CG     0      2
#> 20:  chr17      + 43126252      CG     0      1
#> 21:  chr17      + 43126255      CG     0      1
#> 22:  chr17      + 43126258      CG     0      1
#> 23:  chr17      + 43126264      CG     0      1
#> 24:  chr17      + 43126286      CG     0      1
#> 25:  chr17      + 43126324      CG     0      1
#> 26:  chr17      + 43126330      CG     0      1
#> 27:  chr17      + 43126332      CG     0      1
#> 28:  chr17      + 43126340      CG     0      1
#> 29:  chr17      + 43126363      CG     0      1
#> 30:  chr17      + 43126393      CG     0      1
#> 31:  chr17      + 43126408      CG     0      1
#> 32:  chr17      + 43126412      CG     0      1
#> 33:  chr17      - 43126869      CG     1      0
#> 34:  chr17      - 43126871      CG     1      0
#> 35:  chr17      - 43126876      CG     1      0
#> 36:  chr17      - 43126880      CG     1      0
#> 37:  chr17      - 43126890      CG     1      0
#> 38:  chr17      - 43126892      CG     1      0
#> 39:  chr17      - 43126907      CG     1      0
#> 40:  chr17      - 43126909      CG     1      0
#> 41:  chr17      - 43126918      CG     1      0
#> 42:  chr17      - 43126925      CG     1      0
#> 43:  chr17      - 43126933      CG     1      0
#> 44:  chr17      - 43126939      CG     1      0
#> 45:  chr17      - 43126962      CG     1      0
#> 46:  chr17      - 43126985      CG     1      0
#>      rname strand      pos context  meth unmeth
#>     <fctr> <fctr>    <int>  <fctr> <int>  <int>
  
  # long-read sequencing
  longread.data <- preprocessBam(
    system.file("extdata", "longread.bam", package="epialleleR"),
    min.mapq=30, min.baseq=20, min.prob=178
  )
#> Checking BAM file: 
#> long-read, single-end, unsorted alignment detected
#> Reading single-end BAM file 
#> [0.005s]
  generateCytosineReport(longread.data, threshold.reads=FALSE)
#> Filtering reads 
#> [0.000s]
#> Preparing cytosine report 
#> [0.027s]
#>       rname strand      pos context  meth unmeth
#>      <fctr> <fctr>    <int>  <fctr> <int>  <int>
#>   1:  chr17      - 43115270      CG     1      0
#>   2:  chr17      - 43115300      CG     1      0
#>   3:  chr17      - 43115371      CG     1      0
#>   4:  chr17      - 43115417      CG     1      0
#>   5:  chr17      - 43115427      CG     1      0
#>  ---                                            
#> 903:  chr17      + 43136994      CG     0      1
#> 904:  chr17      + 43137174      CG     1      0
#> 905:  chr17      + 43137332      CG     1      0
#> 906:  chr17      + 43137364      CG     1      0
#> 907:  chr17      + 43137391      CG     1      0
  
  # Specifics of long-read alignment processing
  out.bam <- tempfile(pattern="out-", fileext=".bam")
  
  simulateBam(
    seq=c("ACGCCATYCGCGCCA"),
    Mm=c("C+m,0,2,0;"),
    Ml=list(as.integer(c(102,128,153))),
    output.bam.file=out.bam
  )
#> Writing sample BAM 
#> [0.002s]
#> [1] 1
  generateCytosineReport(out.bam, threshold.reads=FALSE, report.context="CX")
#> Checking BAM file: 
#> long-read, single-end, unsorted alignment detected
#> Reading single-end BAM file 
#> [0.002s]
#> Filtering reads 
#> [0.000s]
#> Preparing cytosine report 
#> [0.000s]
#>     rname strand   pos context  meth unmeth
#>    <fctr> <fctr> <int>  <fctr> <int>  <int>
#> 1:   chrS      +     2      CG     1      0
#> 2:   chrS      +     4     CHH     0      1
#> 3:   chrS      +     5     CHH     0      1
#> 4:   chrS      +     9      CG     1      0
#> 5:   chrS      +    11      CG     1      0
#> 6:   chrS      +    13     CHH     0      1
#> 7:   chrS      +    14     CHH     0      1
  
  simulateBam(
    seq=c("ACGCCATYCGCGCCA"),
    Mm=c("G-m,0,0,0;"),
    Ml=list(as.integer(c(138,101,96))),
    output.bam.file=out.bam
  )
#> Writing sample BAM 
#> [0.002s]
#> [1] 1
  generateCytosineReport(out.bam, threshold.reads=FALSE, report.context="CX")
#> Checking BAM file: 
#> long-read, single-end, unsorted alignment detected
#> Reading single-end BAM file 
#> [0.001s]
#> Filtering reads 
#> [0.000s]
#> Preparing cytosine report 
#> [0.001s]
#>      rname strand   pos context  meth unmeth
#>     <fctr> <fctr> <int>  <fctr> <int>  <int>
#>  1:   chrS      +     2      CG     0      1
#>  2:   chrS      -     3      CG     1      0
#>  3:   chrS      +     4     CHH     0      1
#>  4:   chrS      +     5     CHH     0      1
#>  5:   chrS      +     9      CG     0      1
#>  6:   chrS      -    10      CG     1      0
#>  7:   chrS      +    11      CG     0      1
#>  8:   chrS      -    12      CG     1      0
#>  9:   chrS      +    13     CHH     0      1
#> 10:   chrS      +    14     CHH     0      1
  
  simulateBam(
    seq=c("ACGCCATYCGCGCCA"),
    Mm=c("C+m,0,2,0;G-m,0,0,0;"),
    Ml=list(as.integer(c(102,128,153,138,101,96))),
    output.bam.file=out.bam
  )
#> Writing sample BAM 
#> [0.003s]
#> [1] 1
  generateCytosineReport(out.bam, threshold.reads=FALSE, report.context="CX")
#> Checking BAM file: 
#> long-read, single-end, unsorted alignment detected
#> Reading single-end BAM file 
#> [0.002s]
#> Filtering reads 
#> [0.000s]
#> Preparing cytosine report 
#> [0.000s]
#>      rname strand   pos context  meth unmeth
#>     <fctr> <fctr> <int>  <fctr> <int>  <int>
#>  1:   chrS      +     2      CG     1      0
#>  2:   chrS      -     3      CG     1      0
#>  3:   chrS      +     4     CHH     0      1
#>  4:   chrS      +     5     CHH     0      1
#>  5:   chrS      +     9      CG     1      0
#>  6:   chrS      -    10      CG     1      0
#>  7:   chrS      +    11      CG     1      0
#>  8:   chrS      -    12      CG     1      0
#>  9:   chrS      +    13     CHH     0      1
#> 10:   chrS      +    14     CHH     0      1
  
```
