# DONUTS

## Introduction
DONUTS (**D**ecomp**o**sing **n**ature and n**u**r**t**ure using GWAS **s**ummary statistics) is a novel statistic framework that can estimate direct and indirect genetic effects at the SNP level. It requires GWAS summary statistics as input, allows differential paternal and maternal effects, and accounts for GWAS sample overlap and assortative mating. DONUTS has low computational burden and can complete genome-wide analyses within seconds. We developed an R package `DONUTS` with a main function `donuts()` for the framework.

## Prerequisites
The R package was developed and tested on Linux and Mac OS environments. It should also work on Windows. The statistical computing software [R](https://www.r-project.org/) (>= 3.5.1) and the following R packages are required (older version may also work but we did not test):
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html) (>= 1.13.0)
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) (>= 0.8.3)

## Tutorial
The inputs are GWAS summary statistics. Genetic covariance intercept in [LD score regression](https://github.com/bulik/ldsc) could be used to accounts for potential sample overlap among input GWAS. Briefly, there are 2 steps:

1. run [LD score regression](https://github.com/bulik/ldsc) using any pair of the input GWAS summary statistics to get the genetic covariance intercept.
2. read the GWAS summary statistics into R as data.frame, and run `donuts()` function in our `DONUTS` package. 
    * depending on the trait and/or research purpose, 2 or 3 GWAS can be used as input and you can specify it using the `mode` argument---see details below and also in the function's help page.
    * if needed, assortative mating at each locus can be specified by the `alpha` argument within `donuts()` function.

The following GWAS may be used as input:
* GWAS-O: own phenotype ~ own genotype. This is the standard GWAS that people usually do.
* GWAS-M: offspring phenotype ~ mother's genotype
* GWAS-P: offspring phenotype ~ father's genotype
* GWAS-MP: offspring phenotype ~ parental genotype. In this GWAS, mothers and fathers from different families are pooled together. This GWAS can be used as an input only when the indirect maternal effect equals to the indirect paternal effects at each locus, or, when there are equal numbers of mothers and fathers in the GWAS.

Note, since the direct and indirect effect sizes are linear combinations of the effect sizes from the input marginal GWAS, it is critical that the phenotype in the input GWAS are on the same scale (i.e., either all standardized or all on the raw scale).

`donuts()` can take the following arguments (more details can be found in the function's help):
* `ss.own`: data.frame; GWAS-O summary statistics.
* `ss2`: data.frame; the 2nd input GWAS summary statistics. Depending on the `mode`, this could be GWAS-M, GWAS-P, or GWAS-MP.
* `ss3`: data.frame; the 3rd input GWAS summary statistics. Default is `NULL`. This is GWAS-P when `mode` = 3.
* `l12`: numeric; `ss.own` and `ss2` LDSC genetic covariance intercept.
* `l13`: numeric; `ss.own` and `ss3` LDSC genetic covariance intercept.
* `l23`: numeric; `ss2` and `ss3` LDSC genetic covariance intercept.
* `n1`: integer; sample size for `ss.own`; default is `NULL`. No need to specify if summary statistics file contains the sample size column "N".
* `n2`: integer; sample size for `ss2`; default is `NULL`. No need to specify if summary statistics file contains the sample size column "N".
* `n3`: integer; sample size for `ss3`; default is `NULL`. No need to specify if summary statistics file contains the sample size column "N".
* `alpha`: numeric or data.frame; default is 0. correlation between spousal genotype (i.e., Corr(Gm, Gp)) at each locus.
* `mode`: integer 1, 2, or 3; default is 2.
* `OutDir`: output directory to write the direct and indirect effect summary statistics (.gz files); default is the current directory. If is `NULL`, analysis will run but won't write the results.

To use the function `donuts()` in the R package, the input GWAS summary statistics should be read into R as data.frame, and each must include the following column names:

* CHR: chromosome
* BP: base-pair coordinate
* SNP: variant ID
* A1: effect allele
* A2: non-effect allele
* BETA: effect size for A1 allele
* SE: standard error
* P: p-value

Sample sizes are also needed for computing effective sample sizes for the direct and indirect effects. They can be either included in the input GWAS summary statistics by column "N", or can be specified in `donuts()` function. Note, if sample sizes are specified in `donuts()`, these values will be used.

Depending on the phenotype and/or the research purpose, 2 or 3 GWAS could be used as input. This could be specified by `mode` argument in the `donuts()` function:
* `mode` = 1: 3 inputs are expected: `ss.own` is GWAS-O, `ss2` is GWAS-M, and `ss3` is GWAS-P. This is the most general scenario. In this case, will obtain direct, indirect, indirect maternal, and indirect paternal effects.
* `mode` = 2: 2 inputs are expected: `ss.own` is GWAS-O and `ss2` is GWAS-MP. This is valid only when indirect maternal effect equals to the indirect paternal effect at each locus or when there are equal number of mothers and fathers in GWAS-MP. In this case, we can obtain direct and indirect genetic effects.
* `mode` = 3: 2 inputs are expected: `ss.own` is GWAS-O, and `ss2` is GWAS-M or GWAS-P. If `ss2` is GWAS-M (GWAS-P), you're assuming indirect paternal (maternal) effect is 0. In this case, we'll have direct, indirect, and indirect maternal (paternal) genetic effects.

Assortative mating: argument `alpha` = Corr(Gm, Gp) is the correlation between spousal genotypes at each locus. Default is 0 (i.e., random mating). You can also specify `alpha` for each SNP using a data.frame with two columns: "SNP" for variant ID and "alpha" for the correlation.

Some other features:
* Will first check whether there are duplicated SNPs (by their IDs) and SNPs with duplicated IDs will be removed.
* Will keep only the SNPs that are present in all the input GWAS (and `alpha` if it's a data.frame).
* Will compare `A1` and `A2` among the input GWAS and multiply `BETA` by -1 when necessary so that it refers to the effect of a same allele.
* Will report a very detailed log.


## Getting started

### Install `DONUTS` R package

Open `R` and run the following codes to install `DONUTS` package:

```
library(devtools)
install_github("qlu-lab/DONUTS")

library(DONUTS)  # after successful installation, load the package
?donuts  # check the manual for the donuts() function
```

Alternatively, you can also directly install the pre-packed package from [Releases](https://github.com/qlu-lab/DONUTS/releases).


### Example codes

Below are several example `R` codes. You may need to modify the codes to meet your needs. The first step is to run [LD score regression](https://github.com/bulik/ldsc) to get the genetic covariance intercept. I skipped this step in the examples. You may need to check their  [website](https://github.com/bulik/ldsc) for LDSC tutorials.

Example 1. Assuming we have GWAS-O and GWAS-M summary statistics, already got their genetic covariance intercept, and we are willing to assume that the indirect paternal effect is 0 (therefore, `mode = 3`), then we can use the following R codes to run the analysis to get the direct, indirect, and indirect maternal effects:

```
library(tidyverse)
library(data.table)
library(DONUTS)
options(stringsAsFactors = F)

df.own <- as.data.frame(fread("example/ss1.sumstats.gz"))  # GWAS-O sumstats
df.mat <- as.data.frame(fread("example/ss2.sumstats.gz"))  # GWAS-M sumstats

str(df.own)
str(df.mat)  # check that the summary statistics have correct column names

l12 <- 0.1379  # genetic covariance intercept reported by LDSC

ss.out <- donuts(ss.own = df.own, ss2 = df.mat, mode = 3, 
                 l12 = l12, OutDir = "results")  

# output to "results" folder; in the "results" folder, you'll find the following files:
# direct_effect.sumstats.gz: direct effect summary statistics
# indirect_effect.sumstat.gz: indirect effect summary statistics
# indirect_ss2_effect.sumstats.gz: since the mode is 3 and ss2 is GWAS-M, we can also have the indirect maternal effect
# all_aligned.sumstats.gz: contains both the input sumstats and output results; SNP information are aligned to ss.own 

```

Example 2. Using user defined `alpha`. The default is random mating for each SNP (i.e., `alpha = 0`). However, user could also prepare a text file containing two columns "SNP" (variant ID) and "alpha" (Corr(Gm, Gp) for each SNP), read it into R as a data.frame and pass it to `donuts()`.

```
library(tidyverse)
library(data.table)
library(DONUTS)
options(stringsAsFactors = F)

df.own <- as.data.frame(fread("example/ss1.sumstats.gz"))  # GWAS-O sumstats
df.mat <- as.data.frame(fread("example/ss2.sumstats.gz"))  # GWAS-M sumstats

str(df.own)
str(df.mat)  # check that the summary statistics have correct column names

l12 <- 0.1379  # genetic covariance intercept reported by LDSC

alpha <- as.data.frame(fread("example/SNP_alpha.txt.gz"))  # assortative mating
str(alpha)  # must include two columns: "SNP" and "alpha"

ss.out2 <- donuts(ss.own = df.own, ss2 = df.mat, mode = 3, 
                  l12 = l12, alpha = alpha, OutDir = NULL)  # run analysis, but won't write the results

```

Note, in the example data, to save space and for the purpose of illustration, only a few SNPs are kept in the GWAS summary statistics and thus LDSC won't run with so few SNPs if you use our example data.

## Citation

## License
All rights reserved for [Lu-Laboratory](https://qlu-lab.org/).

