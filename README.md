### QDNAseq: Quantitative DNA sequencing for chromosomal aberrations

This repository contains source code for the R/Bioconductor package QDNAseq.

#### Citing QDNAseq

To cite QDNAseq, please use:  
Scheinin I, Sie D, Bengtsson H, van de Wiel MA, Olshen AB, van Thuijl HF, van
Essen HF, Eijk PP, Rustenburg F, Meijer GA, Reijneveld JC, Wesseling P, Pinkel
D, Albertson DG and Ylstra B (2014) **DNA copy number analysis of fresh and
formalin-fixed specimens by shallow whole-genome sequencing with identification
and exclusion of problematic regions in the genome assembly.** *Genome
Research* **24**: 2022-2032

#### Repositories, Branches, and Versions

Repository [QDNAseq][github] contains two branches: *master* and *release*.
Branch *master* contains the current Bioconductor
[development version][bioc-devel], and branch *release* contains the current
[release version][bioc-release].

The *release* branch of repository [QDNAseq][github] is also mirrored in the
[QDNAseq-release][github-release] repository as branch *master*. This is
because the two repositories use the [Bioconductor Git-SVN bridge][bridge]
to mirror changes to the Bioconductor SVN repository, and the bridge can only
use branch *master* of each GitHub repository.

#### Software quality

* R CMD check status:
 - Bioconductor: <a href="http://master.bioconductor.org/checkResults/devel/bioc-LATEST/QDNAseq/">Multipleplatform build/check report</a>
 - Travis CI: <a href="https://travis-ci.org/ccagc/QDNAseq"><img src="https://travis-ci.org/ccagc/QDNAseq.svg?branch=master" alt="Build status"></a>
* Test coverage status:
 - Coveralls CI: <a href='https://coveralls.io/r/ccagc/QDNAseq?branch=devel'><img src='https://coveralls.io/repos/ccagc/QDNAseq/badge.png?branch=devel' alt='Coverage Status' /></a>

[bioc-devel]: http://bioconductor.org/packages/devel/bioc/html/QDNAseq.html
[bioc-release]: http://bioconductor.org/packages/release/bioc/html/QDNAseq.html
[bridge]: http://bioconductor.org/developers/how-to/git-svn/
[github]: https://github.com/ccagc/QDNAseq
[github-release]: https://github.com/ccagc/QDNAseq-release
