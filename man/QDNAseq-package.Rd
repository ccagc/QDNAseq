%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  999.package.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{QDNAseq-package}
\alias{QDNAseq-package}
\alias{QDNAseq}
\docType{package}

\title{Package QDNAseq}


\description{
    Quantitative DNA sequencing for chromosomal aberrations.
        The genome is divided into non-overlapping fixed-sized bins, number of
        sequence reads in each counted, adjusted with a simultaneous
        two-dimensional loess correction for sequence mappability and GC
        content, and filtered to remove spurious regions in the genome.
        Downstream steps of segmentation and calling are also implemented via
        packages DNAcopy and CGHcall, respectively.
}

\section{Details}{
    A package to detect chromosomal aberrations from whole-genome sequencing
    data. \code{\link{QDNAseqReadCounts}} and \code{\link{QDNAseqCopyNumbers}} classes are
    used as the main data structures.
}

\section{How to cite this package}{
    Whenever using this package, please cite:
    Scheinin I, Sie D, Bengtsson H, van de Wiel MA, Olshen AB, van Thuijl
HF, van Essen HF, Eijk PP, Rustenburg F, Meijer GA, Reijneveld JC,
Wesseling P, Pinkel D, Albertson DG and Ylstra B (2014). "DNA copy
number analysis of fresh and formalin-fixed specimens by shallow
whole-genome sequencing with identification and exclusion of
problematic regions in the genome assembly." _Genome Research_, *24*,
pp. 2022-2032.
}

\author{Ilari Scheinin}

\section{License}{
    This package is licensed under GPL.
}
\keyword{package}
