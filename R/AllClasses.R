setClass('QDNAseqSignals',
  contains    = 'eSet',
  prototype   = prototype(new('VersionedBiobase',
  versions    = c(classVersion('eSet'), QDNAseq=packageVersion('QDNAseq'))))
)

setClass('QDNAseqReadCounts',
  contains    = 'QDNAseqSignals',
  prototype   = prototype(new('VersionedBiobase',
  versions    = c(classVersion('eSet'), QDNAseq=packageVersion('QDNAseq'))))
)

setClass('QDNAseqCopyNumbers',
  contains    = 'QDNAseqSignals',
  prototype   = prototype(new('VersionedBiobase',
  versions    = c(classVersion('eSet'), QDNAseq=packageVersion('QDNAseq'))))
)

# EOF
