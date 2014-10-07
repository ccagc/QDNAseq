%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  clonality.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{clonalityTest}
\alias{clonalityTest}


\alias{clonalityTest}

\title{Clonality testing using copy number data}

\description{
      This functions will calculate the correlation score and the log
      likelihood ratio (llr) for two given data sets. The combination of these
      two measures is used to determine their clonal relationship.
}

\arguments{
      \item{obj}{A QDNAseq object containing segmented data}
      \item{sbjctLst}{Character vector containing identifiers for the
          subjects, eg. patient numbers}
      \item{labels}{Character vector containing the labels for the samples,
          eg. biopsy number}
      \item{corData}{OPTIONAL use pre calculated corData from "correlation"
           function}
      \item{llrData}{OPTIONAL use pre calculated llrData from "likelihood"
           function}
      \item{refdata}{OPTIONAL collection of clonaly unrelated data to be used
           as the negative reference set}
}

\value{
    Returns a list containing data frames for the correlation and llr scores
}

\references{
Ostrovnaya, I., Olshen, A. B., Seshan, V.E., Orlow, I., Albertson, D. G. and
Begg, C. B. (2010), A metastasis or a second independent cancer? Evaluating
the clonal origin of tumors using array copy number data.
Statistics in Medicine, 29: 1608-1621

Ostrovnaya, I. and Begg, C.  Testing Clonal Relatedness of Tumors Using Array
Comparative Genomic Hybridization: A Statistical Challenge
Clin Cancer Res March 1, 2010 16:1358-1367

Venkatraman, E. S. and Olshen, A. B. (2007). A faster circular binary
segmentation algorithm for the analysis of array CGH data.
Bioinformatics, 23:657 63.

Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. (2004). Circular
binary segmentation for the analysis of array-based DNA copy number data.
Biostatistics 5: 557-572.

}

\seealso{
    Internally, \code{\link[Clonality]{clonality.analysis}} of the \pkg{Clonality}
    package is used.
}

\author{Daoud Sie}

