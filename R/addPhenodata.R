addPhenodata <- function(obj, phenofile) {
  pdata <- read.table(phenofile, header=TRUE, sep='\t', as.is=TRUE, row.names=1)
  obj[['phenodata']] <- cbind(obj[['phenodata']], pdata[rownames(obj[['phenodata']]),])
  obj
}

# EOF
