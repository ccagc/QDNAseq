#########################################################################/**
# @RdocClass QDNAseqReadCounts
#
# @alias chromosomes,QDNAseqReadCounts-method
# @alias bpstart,QDNAseqReadCounts-method
# @alias bpend,QDNAseqReadCounts-method
#
# @title "Container for QDNAseq read count data"
#
# \description{
#   @get "title"
# }
#
# \section{Assay data elements}{
#   An object of this class contains (a subset) the following elements:
#   \describe{
#     \item{\code{counts}}{(@numeric) Binned uncorrected read counts
#       as non-negative integers in \eqn{\{0,1,2,...\}}.
#       An object with only this field is created by @see "binReadCounts".}
#     \item{\code{corrected}}{(@numeric; optional) Corrected "count" signals
#       as non-negative (alway true?!?) doubles in \eqn{[0,+\infty)}.
#       This element is added after calling @see "correctBins".}
#     \item{\code{copynumber}}{(@numeric; optional) Corrected and normalized
#       "count" signals in \eqn{(-\infty,+\infty)} (mostly (always??!) non-negative).
#       This element is added after calling @see "normalizeBins".}
#   }
# }
#
# \section{Missing values}{
#   The bin data (assay data) may contain missing values.
# }
#
#
# @author
#*/#########################################################################


setClass('QDNAseqReadCounts',
  contains    = 'eSet',
  prototype   = prototype(new('VersionedBiobase',
  versions    = c(classVersion('eSet'), QDNAseq=packageVersion('QDNAseq'))))
)

# EOF
