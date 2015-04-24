#########################################################################/**
# @RdocClass QDNAseqSignals
#
# @alias chromosomes,QDNAseqSignals-method
# @alias bpstart,QDNAseqSignals-method
# @alias bpend,QDNAseqSignals-method
#
# @title "A parent class for containers of QDNAseq data"
#
# \description{
#     @get "title"
# }
#
# @author "IS"
#*/#########################################################################

setMethod("binsToUse", signature=c(object="QDNAseqSignals"),
    definition=function(object) {
    
    if ("use" %in% colnames(fData(object))) {
        return(fData(object)$use)
    } else if ("filter" %in% colnames(fData(object))) {
        fData(object)$use <- fData(object)$filter
        fData(object)$filter <- NULL
        return(fData(object)$use)
    }
    rep(TRUE, times=nrow(object))
})

setMethod("binsToUse", signature=c(object="AnnotatedDataFrame"),
    definition=function(object) {
    
    if ("use" %in% colnames(object)) {
        return(object$use)
    } else if ("filter" %in% colnames(object)) {
        object$use <- object$filter
        object$filter <- NULL
        return(object$use)
    }
    rep(TRUE, times=nrow(object))
})

setMethod("chromosomes", signature=c(object="QDNAseqSignals"),
    definition=function(object) {
    
    tmp <- fData(object)$chromosome
    tmp[tmp=="X"] <- "23"
    tmp[tmp=="Y"] <- "24"
    as.integer(tmp)
})

setMethod("bpstart", signature=c(object="QDNAseqSignals"),
    definition=function(object) {
    
    fData(object)$start
})

setMethod("bpend", signature=c(object="QDNAseqSignals"),
    definition=function(object) {
    
    fData(object)$end
})

setReplaceMethod("binsToUse",
    signature=c(object="QDNAseqSignals", value="logical"),
    definition=function(object, value) {
    
    fData(object)$use <- value
    object
})

setReplaceMethod("binsToUse",
    signature=c(object="AnnotatedDataFrame", value="logical"),
    definition=function(object, value) {
    
    object$use <- value
    object
})

# EOF
