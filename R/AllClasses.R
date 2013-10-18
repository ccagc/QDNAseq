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
# @author
#*/#########################################################################


setClass('QDNAseqReadCounts',
  contains    = 'eSet',
  prototype   = prototype(new('VersionedBiobase',
  versions    = c(classVersion('eSet'), QDNAseq=packageVersion('QDNAseq'))))
)

# EOF
