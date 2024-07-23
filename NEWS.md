# Version 1.41.2 [2024-07-22]

## Bug Fixes

 * `poolRuns()` would give an error "Error in colMeans2(oldphenodata,
   cols = numericCols, useNames = FALSE) : Argument 'x' must be a
   matrix or a vector."


# Version 1.41.1 [2024-07-20]

## Miscellaneous

 * Fix markup typo in `help("exportBins")`.

 * Add package test for `poolRuns()`.
 

# Version 1.40.0 [2024-05-01]

## Release

 * The version number was bumped for the Bioconductor release version,
     which is now Bioconductor 3.19 for R (>= 4.4.0).


# Version 1.38.0 [2023-10-25]

## Release

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.18 for R (>= 4.3.1).


# Version 1.36.0 [2023-04-26]

## Release

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.17 for R (>= 4.3.0).


# Version 1.34.0 [2022-11-01]

## Release

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.16 for R (>= 4.2.2).


# Version 1.33.1 [2022-04-27]

## Miscellaneous

 * Now the package gives an informative error message when an outdated
   version of the **future** package is used.  It requires **future**
   (>= 1.22.1).


## Bug Fixes

 * A few functions used `class(x) == "data.frame"` rather than
   `inherits(x, "data.frame")`.


# Version 1.33.0 [2022-04-26]

## Release

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.16 for R-devel.


# Version 1.32.0 [2022-04-26]

## Release

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.15 for R (>= 4.2.0).


# Version 1.31.0 [2021-10-27]

## Release

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.15 for R-devel.


# Version 1.30.0 [2021-10-27]

## Release

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.14 for R (>= 4.1.1).


# Version 1.29.6 [2021-10-23]

## Bug Fixes

 * `segmentBins()` would _report_ the sample names as `"NA"` in output
   messages if the sample name contained hyphens, or other symbols
   automatically replaced by `data.frame(..., check.names = TRUE)`.
   This was a harmless bug.


# Version 1.29.5 [2021-10-20]

## Performance

 * All internal row and column-based **matrixStats** calls now avoid
   overhead from handling row and column names.


## Miscellaneous

 * Moved **future** from Imports to Suggests.


# Version 1.29.4 [2021-10-16]

## Documentation

 * Vignette now illustrate parallelization using the `multisession`
   future strategy, instead of the deprecated `multiprocess` strategy.


## Deprecation and Defunct

 * Argument `seeds` of `segmentBins()` is defunct.  It has been
   deprecated and ignored since **QDNAseq** 1.21.3 (September 2019).


# Version 1.29.3 [2021-10-04]

## New Features

 * Now argument `logTransform` of `exportBins()` is ignored if `type =
   "calls"`.

 * Now `exportBins()` returns the pathname to the files written.


## Software Quality

 * Test code coverage was increased from 42% to 52%.

 * Add package test for `exportBins()`.


## Bug Fixes

 * `exportBins(fit, format = "seg", ...)` and `format = "vcf"` would
   merge segments with equal copy-number calls if they were
   interweaved with copy-neutral segments.

 * `exportBins(fit, format = "seg", ...)` and `format = "vcf"`
   produced an obscure error with messages `Error in dimnames(x) <-
   dn : length of dimnames [2] not equal to array extent` for samples
   with no copy-number abberations.

 * `exportBins(fit, format = "seg", file = ...)` and `format = "vcf"`
   did not respect argument `file` but instead wrote files of its own
   names to the current working directory.

 * `exportBins()` would corrupt option `scipen`.  Now it is left
   unchanged.


## Known Issues

 * `callBins()` produces warnings on `Recycling array of length 1 in
   vector- array arithmetic is deprecated. Use c() or as.vector()
   instead.` in R (>= 3.4.0).  This is a problem in the package
   **CGHcall** dependency and is something that needs to be fixed
   there.  For further details, please see
   https://github.com/tgac-vumc/CGHcall/issues/2.


# Version 1.29.2 [2021-09-22]

## Software Quality

 * Test code coverage was increased from 32% to 39%.

 * Added package tests for `binReadCounts()`.


## Bug Fixes

 * `binReadCounts()` would fail when specifying argument `chunkSize`.
   The fix was to require **future** package version 1.22.1 or newer.


# Version 1.29.1 [2021-08-26]

## Software Quality

 * Add package test for `binReadCounts()`.


# Version 1.28.0 [2021-05-19]

## Release

 * The version number was bumped for the Bioconductor release version, which is now Bioconductor 3.13 for R (>= 4.0.3).


# Version 1.26.0 [2020-10-28]

## Release

 * The version number was bumped for the Bioconductor release version, which is now Bioconductor 3.12 for R (>= 4.0.0).


# Version 1.24.0 [2019-10-30]

## Release

 * The version number was bumped for the Bioconductor release version, which is now Bioconductor 3.11 for R (>= 3.6.1).


# Version 1.22.0 [2018-10-30]

## Release

 * Bioconductor 3.10


# Version 1.21.6 [2019-09-25]

## Bug Fixes

 * Link to the 'GEM library' tool was broken.


# Version 1.21.5 [2019-09-09]

## Bug Fixes

 * `plot(..., `logTransform=TRUE`, doSegments=TRUE)` on QDNAseqSignals would position segments that were out of range incorrectly, because it forgot to take the log transform on those outliers.


# Version 1.21.4 [2019-09-06]

## Significant Changes

 * Bin annotation data files are no longer downloaded
   automatically. They are instead part of Bioconductor annotation
   packages **QDNAseq.hg19** and **QDNAseq.mm10**, which needs to be
   installed by the user.  If not installed, an error is now produced.
   The reason for this change, is that the **QDNAseq** maintainers
   will no longer host QDNAseq bin annotation files online (in the
   cloud).


# Version 1.21.3 [2019-09-04]

## New Features

 * All functions that produce verbose output gained argument
   `verbose`.  For backward compatibility, it currently defaults to
   '`verbose=TRUE'` but that may be changed to '`verbose=FALSE'` in a
   future release.


## Improvements

 * `callBins()` now respects option `QDNAseq::verbose` for controlling
   whether output from the **CGHcall** package should be relayed or
   not.

 * MEMORY: Utilize more memory-efficient **matrixStats** functions
   `colSums2()`, `colMeans2()`, etc.


## Deprecation and Defunct

 * Argument `seeds` of `segmentBins()` is deprecated because it did
   not use proper parallel random number generation (RNG).  We now
   instead rely on `future.apply::future_lapply(...,
   future.seed=TRUE)` for this.


# Version 1.21.2 [2019-09-03]

## Significant Changes

 * Package now imports the **future** and **future.apply** packages;
   previously the **future** was listed as a suggested package.

 * Package no longer depends on **BiocParallel**.

 * `binReadCounts()` now uses the future framework instead of
   **BiocParallel** for parallelization.


## Improvements

 * MEMORY: Avoiding data type coercions in more places by for instance
   making sure that vectors and matrices are initated with the values
   of the correct data type (instead of the default NA, which is a
   logical value).


## Software Quality

 * Using `future_lapply()` and `future_apply()` of the well-tested
   **future.apply** package instead of internal analogue
   implementations.

 * TESTS: Now testing numerical reproducibility also for parallel
   processing (using **future** strategies `multisession` and
   `multicore`).

 * TESTS: Now asserting numerical reproducibility of also
   `segmentBins()` and `callBins()`.


# Version 1.21.1 [2019-08-30]

## Significant Changes

 * `exportVCF()` is no longer exported. Use `exportBins(...,
   format="vcf")` instead.


# Version 1.20.0 [2019-05-02]

## Release

 * Bioconductor 3.9.


# Version 1.18.0 [2018-10-30]

## Release

 * Bioconductor 3.8.


# Version 1.16.0 [2018-04-30]

## Release

 * Bioconductor 3.7.


# Version 1.14.0 [2017-10-30]

## Release

 * Bioconductor 3.6.


# Version 1.12.0 [2017-04-23]

## Release

 * Bioconductor 3.5.


## Improvements

 * VCF and SEG file export have been implemented to allow use of
   downstream analysis tools such as Cartegenia (NGS) Bench.

 * `binReadCounts()` now supports parallel computing.

 * `calculateBlackListByRegions()` has been implemented for convient
   bin overlap calculation of any set of regions.


# Version 1.10.0 [2016-10-18]

## Release

 * Bioconductor 3.4.


# Version 1.8.0 [2016-05-04]

## Release

 * Bioconductor 3.3.


## Improvements

 * `estimateCorrection()`, `segmentBins()`, `createBins()`, and
   `calculateBlacklist()` now support parallel computing (see vignette
   for more details)

 * `callBins()` can now also use cutoffs instead of **CGHcall**.

 * `binReadCounts()` now contains parameter `pairedEnds` to specify
   when using paired-end data, so that expected variance can be
   calculated correctly

 * `segmentBins()` now allows seeds for random number generation to be
   specified for reproducibility.

 * `binReadCounts()` supports chunked processing of BAM files.

 * `estimateCorrection()` now also allows correcting for only GC
   content or mappability.


## Bug Fixes

 * `applyFilters()` and `highlightFilters()` now work properly when
   using a numerical value for parameter residual.

 * `highlightFilters()` no longer highlights entire chromosomes for
   which the residual filter is missing altogether, which matches the
   behavior of `applyFilters()`.

 * `getBinAnnotations()` now allows custom bin annotations to be
   loaded via the path parameter even when an annotation package has
   been installed.

 * phenodata files with a single variable are now handled correctly.

 * `calculateMappability()` now retains correct chromosome order even
   when `bigWigAverageOverBed` reorders them.

 * `calculateBlacklist()` now correctly handles non-integer chromosome
   names.


# Version 1.6.0 [2015-10-14]

## Release

 * Bioconductor 3.2.


## Improvements

 * `chromosomes()` no longer tries to return chromosome names as an
   integer vector, but returns a character vector instead.


## Bug Fixes

 * `plot()` and `calculateMappability()` now work also for other
   chromosome names besides numbers and `"X"`, `"Y"`, and `"MT"`.


# Version 1.4.2 [2015-08-20]

## Bug Fixes

 * `createBins()` now properly selects chromosomes when ignoring
   mitochondrial DNA. Please note that the mitochondrial DNA is only
   ignored when it is called either `"chrMT"`, `"chrM"`, `"MT"`, or
   `"M"`.


# Version 1.4.1 [2015-06-30]

## Improvements

 * `correctBins()` now filters out bins with missing loess correction
   estimate.


# Version 1.4.0 [2015-04-17]

## Release

 * Bioconductor 3.1.


# Version 1.2.4 [2015-01-21]

## Other

 * Update package maintainer.


# Version 1.2.3 [2015-01-20]

## Bug Fixes

 * `createBins()` now also works for species not supported by
   **GenomeInfoDb** (such as Canis lupus familiaris).


# Version 1.2.2 [2014-12-23]

## Improvements

 * `applyFilters()` now ignores the residual filter for chromosomes
   for which it does not exist (instead of always filtering out the
   entire chromosome).


# Version 1.2.1 [2014-11-01]

## Improvements

 * Add an asterisk in the noise measure to clarify it's not a regular
   standard deviation or variance, but first scaled with the mean (so
   that the mean is 1.0 and the relationship holds between variance
   and 1/N, where N is the average number of reads per bins).

 * Clarify the use of log2/sqrt transformation in `segmentBins()`.


# Version 1.2.0 [2014-10-14]

## Release

 * Bioconductor 3.0.


## Improvements

 * `segmentBins()` supports another transformation option besides
   log2: `sqrt(x + 3/8)`, which stabilizes the variance.

 * `plot()` can skip plotting of segments and calls by specifying
   `doSegments=FALSE` and `doCalls=FALSE`.

 * `exportBins()` also supports BED files.


## Bug Fixes

 * `exportBins()` saves base pair positions in fixed notation instead
   of scientific.


# Version 1.0.5 [2014-06-13]

## Bug Fixes

 * Gix a bug caused by package **matrixStats** changing `madDiff()`
   from an S4 to an S3 method in version 0.9.4, released on
   2014-05-23.


# Version 1.0.4 [2014-05-23]

## Bug Fixes

 * `getBinAnnotations()` fixed after being broken by a change in
   Bitbucket.


# Version 1.0.3 [2014-05-23]

## Improvements

 * Added `exportBins()` for exporting to a file.

 * Switch graphics in the vignette to PNG.


# Version 1.0.2 [2014-05-15]

## Improvements

 * `plot()` honors user-specified values for `xlab` and `xaxt`.

 * `plot()` allows omission of labels for the standard deviation and
   the number of data points.

 * Improve diagnostic messages.


# Version 1.0.1 [2014-04-17]

## Bug Fixes

 * `smoothOutlierBins()` correctly ignores bins filtered out.


# Version 1.0.0 [2014-04-14]

## Release

 * Initial release as part of Bioconductor 2.14.
