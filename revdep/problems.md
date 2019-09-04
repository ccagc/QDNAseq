# ACE

<details>

* Version: 1.2.0
* Source code: https://github.com/cran/ACE
* URL: https://github.com/tgac-vumc/ACE
* Date/Publication: 2019-05-02
* Number of recursive dependencies: 85

Run `revdep_details(,"ACE")` for more info

</details>

## In both

*   checking dependencies in R code ... NOTE
    ```
    Unexported object imported by a ':::' call: ‘QDNAseq:::sdDiffTrim’
      See the note in ?`:::` about the use of this operator.
    ```

*   checking Rd cross-references ... NOTE
    ```
    Package unavailable to check Rd xrefs: ‘corrplot’
    ```

# GeneBreak

<details>

* Version: 1.14.0
* Source code: https://github.com/cran/GeneBreak
* URL: https://github.com/stefvanlieshout/GeneBreak
* Date/Publication: 2019-05-02
* Number of recursive dependencies: 33

Run `revdep_details(,"GeneBreak")` for more info

</details>

## In both

*   checking R code for possible problems ... NOTE
    ```
    .glmbreak: no visible global function definition for ‘glm’
    .glmbreak: no visible global function definition for ‘predict’
    addGeneAnnotation,CopyNumberBreakPoints: no visible global function
      definition for ‘head’
    bpStats,CopyNumberBreakPoints: no visible global function definition
      for ‘sd’
    bpStats,CopyNumberBreakPoints: no visible global function definition
      for ‘p.adjust’
    Undefined global functions or variables:
      glm head p.adjust predict sd
    Consider adding
      importFrom("stats", "glm", "p.adjust", "predict", "sd")
      importFrom("utils", "head")
    to your NAMESPACE file.
    ```

# HiCcompare

<details>

* Version: 1.6.0
* Source code: https://github.com/cran/HiCcompare
* Date/Publication: 2019-05-02
* Number of recursive dependencies: 143

Run `revdep_details(,"HiCcompare")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Package suggested but not available for checking: 'multiHiCcompare'
    ```

*   checking installed package size ... NOTE
    ```
      installed size is  6.5Mb
      sub-directories of 1Mb or more:
        data   5.6Mb
    ```

*   checking R code for possible problems ... NOTE
    ```
    ...
    sim.other.methods: no visible binding for global variable ‘IF2’
    sim.other.methods: no visible binding for global variable ‘adj.M’
    sim.other.methods: no visible binding for global variable ‘M’
    sim_matrix: no visible binding for global variable ‘bias.slope’
    total_sum: no visible binding for global variable ‘IF2’
    total_sum: no visible binding for global variable ‘M’
    total_sum: no visible binding for global variable ‘IF1’
    total_sum: no visible binding for global variable ‘chr1’
    volcano: no visible binding for global variable ‘A’
    volcano: no visible binding for global variable ‘adj.IF1’
    volcano: no visible binding for global variable ‘adj.IF2’
    volcano: no visible binding for global variable ‘p.value’
    volcano: no visible binding for global variable ‘D’
    Undefined global functions or variables:
      A D IF IF1 IF2 M Z adj.IF1 adj.IF2 adj.M axis bias.slope bp
      centromere_locations chr1 chr2 count fold.change i j p.adj p.value
      pnorm region1 region2 start1 start2
    Consider adding
      importFrom("graphics", "axis")
      importFrom("stats", "D", "pnorm")
    to your NAMESPACE file.
    ```

# QDNAseq.hg19

<details>

* Version: 1.14.0
* Source code: https://github.com/cran/QDNAseq.hg19
* URL: https://github.com/tgac-vumc/QDNAseq.hg19
* BugReports: https://github.com/tgac-vumc/QDNAseq.hg19/issues
* Date/Publication: 2019-05-07
* Number of recursive dependencies: 33

Run `revdep_details(,"QDNAseq.hg19")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 35.4Mb
      sub-directories of 1Mb or more:
        data  35.3Mb
    ```

# QDNAseq.mm10

<details>

* Version: 1.14.0
* Source code: https://github.com/cran/QDNAseq.mm10
* URL: https://github.com/tgac-vumc/QDNAseq.mm10
* BugReports: https://github.com/tgac-vumc/QDNAseq.mm10/issues
* Date/Publication: 2019-05-07
* Number of recursive dependencies: 33

Run `revdep_details(,"QDNAseq.mm10")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is 31.3Mb
      sub-directories of 1Mb or more:
        data  31.2Mb
    ```

