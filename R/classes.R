setClass('qdnaseq',
  contains    = 'eSet',
  prototype   = prototype(new('VersionedBiobase',
  versions    = c(classVersion('eSet'), qdnaseq='1.0.0')))
)

# EOF
