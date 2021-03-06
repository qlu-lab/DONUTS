% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f.donuts.R
\name{donuts}
\alias{donuts}
\title{DONUTS: Decomposing nature and nurture using GWAS summary statistics}
\usage{
donuts(
  ss.own,
  ss2,
  ss3 = NULL,
  l12 = 0,
  l13 = 0,
  l23 = 0,
  n1 = NULL,
  n2 = NULL,
  n3 = NULL,
  alpha = 0,
  mode = 2,
  OutDir = getwd()
)
}
\arguments{
\item{ss.own}{a data.frame; GWAS-O summary statistics}

\item{ss2}{a data.frame; 2nd input GWAS summary statistics.
Depending on the \code{mode}, this could be GWAS-M, GWAS-P, or GWAS-MP. see Details.}

\item{ss3}{a data.frame; default is NULL. 3rd input GWAS summary statistics; This is GWAS-P when \code{mode} == 1. see Details.}

\item{l12}{numeric; \code{ss.own} and \code{ss2} LDSC genetic covariance intercept; default is 0.}

\item{l13}{numeric; \code{ss.own} and \code{ss3} LDSC genetic covariance intercept; default is 0.}

\item{l23}{numeric; \code{ss2} and \code{ss3} LDSC genetic covariance intercept; default is 0.}

\item{n1}{integer; Sample size for ss.own; default is NULL. Only needs to be specified when sample size is not included in ss.own summary statistics.
If specified, this number will be used in the analysis.}

\item{n2}{integer; Sample size for ss2; default is NULL. Only needs to be specified when sample size is not included in ss2 summary statistics.
If specified, this number will be used in the analysis.}

\item{n3}{integer; Sample size for ss3; default is NULL. Only needs to be specified when sample size is not included in ss2 summary statistics.
If specified, this number will be used in the analysis.}

\item{alpha}{numeric or a data.frame; correlation between spousal genotypes (i.e., Corr(Gm, Gp)) at each locus; default is 0.
This value measures the degree of assortative mating. \code{alpha} can also take a data.frame where user can specify the spousal correlation at each locus
The data.frame must contain two columns "SNP" and "alpha", where "SNP" is the column name for SNP ID (rs#) and
"alpha" is the column name for the SNP-level spousal correlation.
When specified as a data.frame, will keep only the overlapping SNPs (by their rs#) among input GWAS and alpha.}

\item{mode}{integer 1, 2, or 3; default is 2; specify analysis scenario -- see Details.}

\item{OutDir}{Output directory to write the direct and indirect effect summary statistics files.
Default is the current directory. If is NULL, the output files won't be written (but will still return the results as a data.frame).}
}
\value{
Returns a data.frame containing both the input and output summary statistics. The basic information about the SNPs are copied from \code{ss.own}.
The contents of this data.frame will be different depending on the \code{mode}.

\code{CHR}: chromosome

\code{BP}: base-pair coordinate

\code{SNP}: variant IDs

\code{A1}: effect allele

\code{A2}: non-effect allele

\code{alpha}: Corr(Gm, Gp) at each locus

\verb{beta.\{own, ss2, ss3, dir, ind, ind.mat, ind.pat, ind.ss2\}}: effect sizes in the input GWAS summary statistics and for the direct and indirect effects.

\verb{se.\{own, ss2, ss3, dir, ind, ind.mat, ind.pat, ind.ss2\}}: standard errors in the input GWAS summary statistics and for the direct and indirect effects.

\verb{p.\{own, ss2, ss3, dir, ind, ind.mat, ind.pat, ind.ss2\}}: p-values in the input GWAS summary statistics and for the direct and indirect effects.

\verb{n.\{own, ss2, ss3, dir, ind, ind.mat, ind.pat, ind.ss2\}}: sample sizes in the input GWAS summary statistics and the effective sample sizes for the direct and indirect effects.

Tutorials and examples can be found at: https://github.com/qlu-lab/DONUTS
}
\description{
This function takes 2 or 3 GWAS summary statistics as input and output summary statistics for the
direct and indirect genetic effects.
}
\details{
This function will first check whether there are duplicated SNPs using variant IDs. SNPs with duplicated IDs will be removed.
Then, it will take intersection of the SNPs among all the inputs (and also with SNPs in \code{alpha} if it's a data.frame) and only the overlapping SNPs will be kept in the output.
The first input summary statistics' \code{A1} and \code{A2} will be used. That is, the other input's \code{BETA} will multiply by -1 if A1 and A2 are flipped,
or will be re-coded as \code{NA} if the alleles cannot be matched.

GWAS-O: standard GWAS of own phenotype ~ own genotype

GWAS-M: offspring phenotype ~ mother's genotype

GWAS-P: offspring phenotype ~ father's genotype

GWAS-MP: offspring phenotype ~ parental genotype, where we pool together mothers and fathers from different families to run the GWAS

The input GWAS summary statistics must contain the following columns with exactly the following column names (they can contain additional columns, but those will not be used):

CHR: chromosome

BP: base-pair coordinate

SNP: variant IDs

A1: effect allele

A2: non-effect allele

BETA: effect size

SE: standard error

P: p-value

They can also contain "N" column for the sample size at each locus.
If the summary statistics does not contain "N", they can be specified by \code{n1}, \code{n2}, or \code{n3} for the 3 input, respectively.
Note, if the sample size is specified by \code{n1}, \code{n2}, or \code{n3}, these values will be used even if the input summary statistics contains "N" column.

The default value for \code{alpha} is 0. You can also specify the spousal correlation at each SNP using a data.frame for \code{alpha}.
If you want to do so, the data.frame \code{alpha} must contain 2 columns: "SNP" column for the variant ID and "alpha" column for the spousal correlation.
When \code{alpha} is specified as a data.frame, only the overlapping SNPs with those in the input summary statistics will be kept in the output.

When \code{mode} == 1, 3 inputs are expected: \code{ss.own} is GWAS-O, \code{ss2} is GWAS-M, and \code{ss3} is GWAS-P.
The returned data.frame will contain the input summary statistics and the direct, indirect, indirect maternal, and indirect paternal effects.
If \code{OutDir} is not NULL, will write summary statistics for the direct effect (\code{direct_effect.sumstats.gz}),
indirect effect (\code{indirect_effect.sumstats.gz}), indirect maternal effect (\code{indirect_maternal_effect.sumstats.gz}),
indirect paternal effect (\code{direct_effect.sumstats.gz}), and a file containing everythings (\code{all_aligned.sumstats.gz}).

When \code{mode} == 2, 2 inputs are expected: \code{ss.own} is GWAS-O, \code{ss2} is GWAS-MP.
The returned data.frame will contain the input summary statistics and the direct and indirect effects.
If \code{OutDir} is not NULL, will write summary statistics for the direct effect (\code{direct_effect.sumstats.gz}) and
indirect effect (\code{indirect_effect.sumstats.gz}), and a file containing everythings (\code{all_aligned.sumstats.gz}).

When \code{mode} == 3, 2 inputs are expected: \code{ss.own} is GWAS-O, \code{ss2} is GWAS-M or GWAS-P.
If \code{ss2} is GWAS-M, you're assuming the indirect paternal effect is 0. If \code{ss2} is GWAS-P, you're assuming the indirect maternal effect is 0.
If \code{OutDir} is not NULL, will write summary statistics for the direct effect (\code{direct_effect.sumstats.gz}),
indirect effect (\code{indirect_effect.sumstats.gz}), indirect maternal (if \code{ss2} is GWAS-M) or indirect paternal (if \code{ss2} is GWAS-P) effect (\code{indirect_ss2_effect.sumstats.gz}),
and a file containing everythings (\code{all_aligned.sumstats.gz}).

Besides the summary statistics, it is highly recommended to first run LDSC among
any pair of your inputs and use LDSC's genetic covariance intercept to account for possible sample overlap.

Note, since the direct and indirect effects are linear combinations of input GWAS,
it is thus critical that all the input GWAS were done on a same phenotype scale.
}
\author{
Yuchang Wu (ywu423@wisc.edu), University of Wisconsin-Madison
}
