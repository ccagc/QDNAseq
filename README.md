<div id="badges"><!-- pkgdown markup -->
<a href="https://bioconductor.org/packages/QDNAseq/"><img border="0" src="https://bioconductor.org/shields/years-in-bioc/QDNAseq.svg" alt="Bioconductor since badge"/></a> <a href="https://bioconductor.org/checkResults/release/bioc-LATEST/QDNAseq/"><img border="0" src="https://bioconductor.org/shields/build/release/bioc/QDNAseq.svg" alt="Bioconductor release build status"/></a> <a href="https://bioconductor.org/checkResults/devel/bioc-LATEST/QDNAseq/"><img border="0" src="https://bioconductor.org/shields/build/devel/bioc/QDNAseq.svg" alt="Bioconductor devel build status"/></a> <a href="https://github.com/ccagc/QDNAseq/actions?query=workflow%3AR-CMD-check"><img border="0" src="https://github.com/ccagc/QDNAseq/workflows/R-CMD-check/badge.svg?branch=develop" alt="Build status"/></a> <a href="https://codecov.io/gh/ccagc/QDNAseq"><img border="0" src="https://codecov.io/gh/ccagc/QDNAseq/branch/develop/graph/badge.svg" alt="Coverage Status"/></a> 
</div>

# QDNAseq: Quantitative DNA Sequencing for Chromosomal Aberrations

This repository contains source code for the R/Bioconductor package
[QDNAseq](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html), which implements the QDNAseq method (Scheinin et al., 2014).  For instructions on how to use QDNAseq, see the ['Introduction to QDNAseq'](https://bioconductor.org/packages/release/bioc/vignettes/QDNAseq/inst/doc/QDNAseq.pdf) vignette, which also installed together with  the package.  Please remember to cite Scheinin et al. (2014) whenever using QDNAseq in your research.

Please use [the official Bioconductor instructions](https://bioconductor.org/packages/release/bioc/html/QDNAseq.html) to install QDNAseq;

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("QDNAseq")
```

To analyze human data, the [QDNAseq.hg19](https://bioconductor.org/packages/release/data/experiment/html/QDNAseq.hg19.html) package must be installed;
```r
BiocManager::install("QDNAseq.hg19")
```

To analyze mouse data, the [QDNAseq.mm10](https://bioconductor.org/packages/release/data/experiment/html/QDNAseq.mm10.html) package must be installed;
```r
BiocManager::install("QDNAseq.mm10")
```

To install the devel versions of QDNAseq, QDNAseq.hg19, and QDNAseq.mm10, see the [Bioconductor devel installation instruction](https://bioconductor.org/packages/devel/bioc/html/QDNAseq.html).



## References

* Scheinin I, Sie D, Bengtsson H, van de Wiel MA, Olshen AB, van Thuijl HF, van
Essen HF, Eijk PP, Rustenburg F, Meijer GA, Reijneveld JC, Wesseling P, Pinkel
D, Albertson DG and Ylstra B. **DNA copy number analysis of fresh and
formalin-fixed specimens by shallow whole-genome sequencing with identification
and exclusion of problematic regions in the genome assembly.** *Genome
Research* **24**: 2022-2032, 2014. [PMC4248318](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4248318/)


## GitHub repository and Bioconductor repository

This GitHub repository is (manually) kept in sync with the Bioconductor [devel version](https://bioconductor.org/packages/devel/bioc/html/QDNAseq.html) of the QDNAseq package.  More precisely, the [master](https://github.com/ccagc/QDNAseq/tree/master) branch in this GitHub repository is in sync with what is on the Bioconductor git repository (e.g. `git clone https://git.bioconductor.org/packages/QDNAseq`).  The source code previous versions is available via git tags, e.g. QDNAseq [1.21.1](https://github.com/ccagc/QDNAseq/tree/1.21.1).
