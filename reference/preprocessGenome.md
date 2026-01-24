# preprocessGenome

This function reads and preprocesses (optionally \`bgzip\`ped) FASTA
file with reference sequences.

## Usage

``` r
preprocessGenome(genome.file, nthreads = 1, verbose = TRUE)
```

## Arguments

- genome.file:

  reference (genomic) sequences file location string.

- nthreads:

  non-negative integer for the number of additional HTSlib threads to be
  used during file decompression (default: 1).

- verbose:

  boolean to report progress and timings (default: TRUE).

## Value

list object containing preprocessed reference sequence data.

## Details

The function loads and preprocesses reference (genomic) sequences,
saving time when methylation calling needs to be performed on multiple
BAM files. Currently, reading the data is done by means of HTSlib,
therefore it is possible to speed up the loading by means of HTSlib
decompression threads when FASTA file is compressed by \`bgzip\`.

This function is also called internally when file location is supplied
as an input for [`callMethylation`](callMethylation.md) method.

\`preprocessGenome\` checks if index file is present, and if not,
creates it automatically. It is possible and recommended to use
compressed FASTA file as an input, but the file must be compressed by
\`bgzip\` (part of samtools/HTSlib). When FASTA file is compressed,
faster loading can be achieved using (typically one) additional HTSlib
decompression thread.

During loading, both lowercase and uppercase ACGTN symbols are allowed
and correctly recognised, however all the other symbols (e.g., extended
IUPAC symbols, MRSVWYHKDB) within sequences are converted to N.

Please also note that for the purpose of methylation calling, the very
same reference genome must be used for both alignment (when BAM is
produced) and calling cytosine methylation by
[`callMethylation`](callMethylation.md) method.

## See also

[`callMethylation`](callMethylation.md) for methylation calling, and
\`epialleleR\` vignettes for the description of usage and sample data.

Block compression/decompression utility
[bgzip](http://www.htslib.org/doc/bgzip.md).

## Examples

``` r
  genome.file <- system.file("extdata", "test", "reference.fasta.gz", package="epialleleR")
  genome.data <- preprocessGenome(genome.file)
#> Reading reference genome file 
#> [0.000s]
```
