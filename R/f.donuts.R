#' DONUTS: Decomposing nature and nurture using GWAS summary statistics
#'
#' This function takes 2 or 3 GWAS summary statistics as input and output summary statistics for the
#' direct and indirect genetic effects.
#'
#' @param ss.own a data.frame; GWAS-O summary statistics
#' @param ss2 a data.frame; 2nd input GWAS summary statistics.
#' Depending on the `mode`, this could be GWAS-M, GWAS-P, or GWAS-MP. see Details.
#' @param ss3 a data.frame; default is NULL. 3rd input GWAS summary statistics; This is GWAS-P when `mode` == 1. see Details.
#' @param l12 numeric; `ss.own` and `ss2` LDSC genetic covariance intercept; default is 0.
#' @param l13 numeric; `ss.own` and `ss3` LDSC genetic covariance intercept; default is 0.
#' @param l23 numeric; `ss2` and `ss3` LDSC genetic covariance intercept; default is 0.
#' @param n1 integer; Sample size for ss.own; default is NULL. Only needs to be specified when sample size is not included in ss.own summary statistics.
#' If specified, this number will be used in the analysis.
#' @param n2 integer; Sample size for ss2; default is NULL. Only needs to be specified when sample size is not included in ss2 summary statistics.
#' If specified, this number will be used in the analysis.
#' @param n3 integer; Sample size for ss3; default is NULL. Only needs to be specified when sample size is not included in ss2 summary statistics.
#' If specified, this number will be used in the analysis.
#' @param alpha numeric or a data.frame; correlation between spousal genotypes (i.e., Corr(Gm, Gp)) at each locus; default is 0.
#' This value measures the degree of assortative mating. `alpha` can also take a data.frame where user can specify the spousal correlation at each locus
#' The data.frame must contain two columns "SNP" and "alpha", where "SNP" is the column name for SNP ID (rs#) and
#' "alpha" is the column name for the SNP-level spousal correlation.
#' When specified as a data.frame, will keep only the overlapping SNPs (by their rs#) among input GWAS and alpha.
#' @param mode integer 1, 2, or 3; default is 2; specify analysis scenario -- see Details.
#' @param OutDir Output directory to write the direct and indirect effect summary statistics files.
#' Default is the current directory. If is NULL, the output files won't be written (but will still return the results as a data.frame).
#'
#' @return Returns a data.frame containing both the input and output summary statistics. The basic information about the SNPs are copied from `ss.own`.
#' The contents of this data.frame will be different depending on the `mode`.
#'
#' `CHR`: chromosome
#'
#' `BP`: base-pair coordinate
#'
#' `SNP`: variant IDs
#'
#' `A1`: effect allele
#'
#' `A2`: non-effect allele
#'
#' `alpha`: Corr(Gm, Gp) at each locus
#'
#' `beta.{own, ss2, ss3, dir, ind, ind.mat, ind.pat, ind.ss2}`: effect sizes in the input GWAS summary statistics and for the direct and indirect effects.
#'
#' `se.{own, ss2, ss3, dir, ind, ind.mat, ind.pat, ind.ss2}`: standard errors in the input GWAS summary statistics and for the direct and indirect effects.
#'
#' `p.{own, ss2, ss3, dir, ind, ind.mat, ind.pat, ind.ss2}`: p-values in the input GWAS summary statistics and for the direct and indirect effects.
#'
#' `n.{own, ss2, ss3, dir, ind, ind.mat, ind.pat, ind.ss2}`: sample sizes in the input GWAS summary statistics and the effective sample sizes for the direct and indirect effects.
#'
#' Tutorials and examples can be found at: https://github.com/qlu-lab/DONUTS
#'
#' @details This function will first check whether there are duplicated SNPs using variant IDs. SNPs with duplicated IDs will be removed.
#' Then, it will take intersection of the SNPs among all the inputs (and also with SNPs in `alpha` if it's a data.frame) and only the overlapping SNPs will be kept in the output.
#' The first input summary statistics' `A1` and `A2` will be used. That is, the other input's `BETA` will multiply by -1 if A1 and A2 are flipped,
#' or will be re-coded as `NA` if the alleles cannot be matched.
#'
#' GWAS-O: standard GWAS of own phenotype ~ own genotype
#'
#' GWAS-M: offspring phenotype ~ mother's genotype
#'
#' GWAS-P: offspring phenotype ~ father's genotype
#'
#' GWAS-MP: offspring phenotype ~ parental genotype, where we pool together mothers and fathers from different families to run the GWAS
#'
#' The input GWAS summary statistics must contain the following columns with exactly the following column names (they can contain additional columns, but those will not be used):
#'
#' CHR: chromosome
#'
#' BP: base-pair coordinate
#'
#' SNP: variant IDs
#'
#' A1: effect allele
#'
#' A2: non-effect allele
#'
#' BETA: effect size
#'
#' SE: standard error
#'
#' P: p-value
#'
#' They can also contain "N" column for the sample size at each locus.
#' If the summary statistics does not contain "N", they can be specified by `n1`, `n2`, or `n3` for the 3 input, respectively.
#' Note, if the sample size is specified by `n1`, `n2`, or `n3`, these values will be used even if the input summary statistics contains "N" column.
#'
#' The default value for `alpha` is 0. You can also specify the spousal correlation at each SNP using a data.frame for `alpha`.
#' If you want to do so, the data.frame `alpha` must contain 2 columns: "SNP" column for the variant ID and "alpha" column for the spousal correlation.
#' When `alpha` is specified as a data.frame, only the overlapping SNPs with those in the input summary statistics will be kept in the output.
#'
#' When `mode` == 1, 3 inputs are expected: `ss.own` is GWAS-O, `ss2` is GWAS-M, and `ss3` is GWAS-P.
#' The returned data.frame will contain the input summary statistics and the direct, indirect, indirect maternal, and indirect paternal effects.
#' If `OutDir` is not NULL, will write summary statistics for the direct effect (`direct_effect.sumstats.gz`),
#' indirect effect (`indirect_effect.sumstats.gz`), indirect maternal effect (`indirect_maternal_effect.sumstats.gz`),
#' indirect paternal effect (`direct_effect.sumstats.gz`), and a file containing everythings (`all_aligned.sumstats.gz`).
#'
#' When `mode` == 2, 2 inputs are expected: `ss.own` is GWAS-O, `ss2` is GWAS-MP.
#' The returned data.frame will contain the input summary statistics and the direct and indirect effects.
#' If `OutDir` is not NULL, will write summary statistics for the direct effect (`direct_effect.sumstats.gz`) and
#' indirect effect (`indirect_effect.sumstats.gz`), and a file containing everythings (`all_aligned.sumstats.gz`).
#'
#' When `mode` == 3, 2 inputs are expected: `ss.own` is GWAS-O, `ss2` is GWAS-M or GWAS-P.
#' If `ss2` is GWAS-M, you're assuming the indirect paternal effect is 0. If `ss2` is GWAS-P, you're assuming the indirect maternal effect is 0.
#' If `OutDir` is not NULL, will write summary statistics for the direct effect (`direct_effect.sumstats.gz`),
#' indirect effect (`indirect_effect.sumstats.gz`), indirect maternal (if `ss2` is GWAS-M) or indirect paternal (if `ss2` is GWAS-P) effect (`indirect_ss2_effect.sumstats.gz`),
#' and a file containing everythings (`all_aligned.sumstats.gz`).
#'
#' Besides the summary statistics, it is highly recommended to first run LDSC among
#' any pair of your inputs and use LDSC's genetic covariance intercept to account for possible sample overlap.
#'
#' Note, since the direct and indirect effects are linear combinations of input GWAS,
#' it is thus critical that all the input GWAS were done on a same phenotype scale.
#'
#' @author Yuchang Wu (ywu423@wisc.edu), University of Wisconsin-Madison
#'
#' @export
donuts <- function(ss.own, ss2, ss3 = NULL,
                   l12 = 0, l13 = 0, l23 = 0,
                   n1 = NULL, n2 = NULL, n3 = NULL,
                   alpha = 0, mode = 2,
                   OutDir = getwd()){

  # Date:   Oct. 10, 2020
  # Author: Yuchang Wu
  # E-mail: ywu423@wisc.edu
  #
  # this function requires "tidyverse" and "data.table" R packages.
  #
  # this function will take munged sumstats as inputs, find their overlapping SNPs,
  # apply our framework to dissect direct and indirect effects,
  # and ouput new sumstats for the direct and indirect effects
  #
  ### inputs: munged GWAS sumstats, LDSC intercepts, GWAS sample sizes, select mode to meet the study design
  #
  # ss.own: sumstat of GWAS-O (i.e., own phenotype ~ own genotype GWAS)
  # ss2:    the 2nd input GWAS sumstats
  # ss3:    the 3rd input GWAS sumstats
  #
  # l12, l13, l23: intercept from LDSC for genetic correlations
  #           between own vs. ss2, own vs. ss3, and ss2 vs. ss3, respectively
  #           if 3 inputs, ss2 is GWAS-M and ss3 is GWAS-P sumstats
  #
  # alpha: spousal correlation at a locus due to assortative mating
  #        default value set to 0
  #        can take data.frame with two columns: one is SNP ID, and the other is the spousal correlation
  #                                              colunm names must be SNP and alpha, respectively
  #
  # mode: the model for computing the direct and indirect effect;
  #       takes integers 1, 2, or 3
  #       1: 3 inputs (GWAS-O, GWAS-M, GWAS-P)
  #       2: 2 inputs (GWAS-O and GWAS-MP, where either two parents have same indirect effects or same sample size in GWAS-MP)
  #       3: 2 inputs (GWAS-O and GWAS-M or GWAS-P, where only one parent has indirect effect)
  #
  # n1: sample size for the own GWAS; if specified, will overide the sample size column in sumstats
  # n2: sample size for the 2nd input GWAS (ss2); if specified, will overide the sample size column in sumstats
  # n3: sample size for the 3rd input GWAS (ss3); if specified, will overide the sample size column in sumstats
  #
  # plot: logical. If TRUE, will generate the Manhattan and QQ plots to OutDir. Default is TRUE.
  #
  # OutDir: output directory; default is the current directory; if NULL, will not write the results
  #
  ### outputs: a data.frame contaning both input sumstats and the output of direct and indirect effects
  #     each output sumstats contains CHR, BP, SNP, A1, A2, BETA, SE, P, "effective sample size" N,
  #     and sample sizes of each input GWASs



  cat("\n#######################################################################\n")
  cat("# DONUTS: Decomposing nature and nurture using GWAS summary statistics #\n")
  cat("#                     Wu et al.(2021), UW-Madison                      #\n")
  cat("#########################################################################\n")
  cat("\n\tWill take summary statistics from GWASs as inputs. They need to be R's data.frame objects.
        At least 2 sumstats inputs. Can specify LDSC intercept using l12, l13, and l23; defaults are all 0s.
        Can specify sample sizes for each sumstats using n1, n2, and n3, if they are not specified in the input sumstats. Defaults are NULL.
        Can specify alpha if spousal correlation at each loci is known; default is 0 for all SNPs.
        Can specify OutDir for the output directory. Default is current working directory.\n", sep="")


  ### checking: valid inputs? ======================

  # model needs to be an integer:
  if(!mode %in% c(1,2,3)){
    cat("\nError! mode must be 1, 2, or 3!\nmode 1: GWAS-O, GWAS-M, and GWAS-P as inputs.\nmode2: GWAS-O and GWAS-MP as inputs.\nmode3: GWAS-O and GWAS-M or GWAS-P as inputs.\nDefault is mode 2.")
    stopifnot(mode %in% c(1,2,3))
  }

  # LDSC intercept must be numeric:
  if(!is.numeric(l12) & !is.numeric(l13) & !is.numeric(l23)){
    cat("\nError! LDSC intercept must be numeric!\n")
    stopifnot(is.numeric(l12) & is.numeric(l13) & is.numeric(l23))
  }

  # sample sizes must be numeric, if specified:
  if(!is.null(n1)){
    if(!is.numeric(n1)){
      cat("\nError! Sample size for the 1st input sumstats must be numeric!\n")
      stopifnot(is.numeric(n1))
    }
  }

  if(!is.null(n2)){
    if(!is.numeric(n2)){
      cat("\nError! Sample size for the 2nd input sumstats must be numeric!\n")
      stopifnot(is.numeric(n2))
    }
  }

  # alpha must be either a numeric number or a data.frame with at least 2 columns: SNP and alpha
  if(is.numeric(alpha)){
    cat("\nSpousal correaltion at each SNP is set to", alpha, "\n")
  }else if(is.data.frame(alpha)){
    if(all(c("SNP", "alpha") %in% names(alpha))){
      cat("\nSpousal correlation at each SNP will be given by data.frame alpha.\n")
    }else{
      cat("\nError! data frame alpha must contain columns SNP and alpha. Please re-name column names.\n")
      stop()
    }
  }else{
    cat("\nError! alpha must be either a numeric number or a data.frame.\n")
    stop()
  }


  if(mode == 3){
    if(!is.null(n3)){
      if(!is.numeric(n3)){
        cat("\nError! Sample size for the 3rd input sumstats must be numeric!\n")
        stopifnot(is.numeric(n3))
      }
    }
  }

  # if OutDir is not NULL, will write the results to the output directory:
  if(is.null(OutDir)){
    cat("\nOutput directory is set to NULL; will ** NOT ** write the results.\nResults will still be returned as a data.frame.\n")
  }else{
    if(!file.exists(OutDir)){
      cat("\nThe directory \"", OutDir, "\" does NOT exisit! Please check. Quitting now...\n")
      stopifnot(file.exists(OutDir))
    }else{
      cat("\nResults will be written to\n\t", OutDir, "\n", sep = "")
    }
  }


  ### here it begins: ========================

  # report the mode:
  if(mode == 1){
    cat("\nmode == 1; expecting 3 input GWASs:\n")
    if(!is.null(ss.own) & !is.null(ss2) & !is.null(ss3)){
      cat("3 inputs are detected:
  the 1st one will be interpreted as GWAS-O (i.e., own phenotype ~ own genotype GWAS) sumtats,
  the 2nd one will be interpreted as GWAS-M (i.e., offpsring phenotype ~ maternal genotype GWAS) sumstats, and
  the 3rd one will be interpreted as GWAS-P (i.e., offspring phenotype ~ paternal genotype GWAS) sumstats.\n")
    }else{
      cat("\nError! Did not find the 3 input GWASs!\n")
      stop()
    }
  }

  if(mode == 2){
    cat("\nmode == 2 (default); expecting 2 input GWASs:\n")
    if(!is.null(ss.own) & !is.null(ss2)){
      cat("2 inputs are detected:
  the 1st one will be interpreted as GWAS-O (i.e., own phenotype ~ own genotype GWAS) sumtats, and
  the 2nd one will be interpreted as GWAS-MP (i.e., offspring phenotype ~ parental genotype GWAS) sumstats.\n")
      if(!is.null(ss3)){
        cat("\nNote, a 3rd input sumstats is also detected, but mode 2 won't use it.\n")
      }
    }else{
      cat("\nError! Did not find 2 input GWASs!\n")
      stop()
    }
  }

  if(mode == 3){
    cat("\nmode == 3; expecting 2 input GWASs:\n")
    if(!is.null(ss.own) & !is.null(ss2)){
      cat("2 inputs are detected:
  the 1st one will be interpreted as GWAS-O (i.e., own phenotype ~ own genotype GWAS) sumtats, and
  the 2nd one will be interpreted as GWAS-M or GWAS-P (i.e., offspring phenotype ~ maternal/paternal genotype GWAS) sumstats.\n
  Note in this mode, we are assuming only one parent (your 2nd input) has non-zero indirect genetic effect.\n")
      if(!is.null(ss3)){
        cat("\nNote, a 3rd input sumstats is also detected, but mode 3 won't use it.\n")
      }
    }else{
      cat("\nError! Did not find 2 input GWASs!\n")
      stop()
    }
  }

  # report the number of SNPs:
  cat("\n  ", format(nrow(ss.own), big.mark = ","), " SNPs are found in the GWAS-O sumstats. \n  ", format(nrow(ss2), big.mark = ","), " SNPs are found in the 2nd input sumstats.\n", sep = "")

  if(mode == 1){
    cat("  ", format(nrow(ss3), big.mark = ","), " SNPs are found in the 3rd input sumstats.\n", sep = "")
  }

  ## preparations:
  # find intersection SNPs:

  # if there are duplicated SNPs (by SNP ID), remove them
  cat("\nChecking whether there are duplicated SNPs (by their rs# IDs). SNPs having duplicate IDs will be removed.\n")

  if(sum(duplicated(ss.own$SNP)) > 0){
    t = ss.own$SNP[duplicated(ss.own$SNP)]
    ss.own = filter(ss.own, !SNP %in% t)
  }

  if(sum(duplicated(ss2$SNP)) > 0){
    t = ss2$SNP[duplicated(ss2$SNP)]
    ss2 = filter(ss2, !SNP %in% t)
  }

  cond = ss.own$SNP %in% ss2$SNP

  if(mode == 1){  # if the 1st mode, also process the 3rd sumstat
    if(sum(duplicated(ss3$SNP)) > 0){
      t = ss3$SNP[duplicated(ss3$SNP)]
      ss3 = filter(ss3, !SNP %in% t)
    }
    cond2 = ss.own$SNP %in% ss3$SNP
    cond = cond & cond2
  }

  if(is.data.frame(alpha)){
    if(sum(duplicated(alpha$SNP)) > 0){
      t = alpha$SNP[duplicated(alpha$SNP)]
      alpha = filter(alpha, !SNP %in% t)
    }
    cond3 = ss.own$SNP %in% alpha$SNP
    cond = cond & cond3
  }


  snpid = ss.own$SNP[cond]  # snpids of the overlapping SNPs
  snpid = snpid[!is.na(snpid)]  # remove NA snps, if any


  # keep only the overlapping SNPs
  cat("\n", format(length(snpid), big.mark = ","), " overlapping SNPs are found across input sumstats. Will use only the overlapping SNPs.\n", sep = "")

  ss.own = filter(ss.own, SNP %in% snpid)
  ss2 = filter(ss2, SNP %in% snpid)

  if(mode == 1){
    ss3 = filter(ss3, SNP %in% snpid)
  }

  if(is.data.frame(alpha)){
    alpha = filter(alpha, SNP %in% snpid)
  }

  # we'll output new sumstats for direct and indirect effects (input sumstats will also be included):
  # our final data.frame "ss.out" will have the following columns:
  #
  # CHR, BP, SNP, A1, A2 <-- copied from GWAS-O sumstats (the 1st input)
  # alpha <--- spousal correlation
  # beta.{own, ss2, ss3, dir, ind, ind.mat, ind.pat} <-- effect sizes
  # se.{own, ss2, ss3, dir, ind, ind.mat, ind.pat} <-- standard errors
  # p.{own, ss2, ss3, dir, ind, ind.mat, ind.pat} <-- p-values
  # n.{own, ss2, ss3, dir, ind, ind.mat, ind.pat} <-- sample sizes
  #
  ss.out = data.frame(CHR = ss.own$CHR, SNP = ss.own$SNP, BP = ss.own$BP,
                      A1 = ss.own$A1, A2 = ss.own$A2,
                      beta.own = ss.own$BETA, se.own = ss.own$SE, p.own = ss.own$P)

  if(is.null(n1)){ # sample sizes from GWAS-O
    ss.out$n.own = ss.own$N
  }else{
    ss.out$n.own = n1
  }

  # add alpha
  if(is.data.frame(alpha)){
    index = match(ss.out$SNP, alpha$SNP)
    alpha = alpha$alpha[index]
  }
  ss.out$alpha = alpha


  # adding results from other input sumstat:
  # Note, A1 and A2 need to be aligned!
  cat("\nAligning input sumstats according to GWAS-O results so that effect sizes are referring to a same allele among all inputs.\n")

  # add results from the 2nd input gwas:
  index <- match(ss.out$SNP, ss2$SNP)
  ss2 <- ss2[index, ]  # match the 2nd sumstats with the own GWAS by SNP IDs

  cond1 <- ss.out$A1 == ss2$A1 & ss.out$A2 == ss2$A2
  cond2 <- ss.out$A1 == ss2$A2 & ss.out$A2 == ss2$A1

  # make sure GWAS result is for the same allele: (note, could get NA here)
  ss.out$beta.ss2 <- ifelse(cond1, ss2$BETA,
                            ifelse(cond2, -1*ss2$BETA, NA))
  ss.out$se.ss2 <- ss2$SE  # se from the 2nd sumstat
  ss.out$p.ss2 <- ss2$P  # p-value from the 2nd sumstat

  if(is.null(n2)){ # sample sizes from 2nd GWAS
    ss.out$n.ss2 = ss2$N
  }else{
    ss.out$n.ss2 = n2
  }

  if(mode == 1){  # mode 1 which requires 3 inputs, add it also:
    # add results from paternal gwas:
    index = match(ss.out$SNP, ss3$SNP)
    ss3 = ss3[index, ]  # match the rows of 3rd input sumstats with the first 2 inputs

    # check whether A1 and A2 are properly aligned:
    cond1 = ss.out$A1 == ss3$A1 & ss.out$A2 == ss3$A2
    cond2 = ss.out$A1 == ss3$A2 & ss.out$A2 == ss3$A1

    ss.out$beta.ss3 = ifelse(cond1, ss3$BETA,
                             ifelse(cond2, -1*ss3$BETA, NA))
    ss.out$se.ss3 = ss3$SE
    ss.out$p.ss3 = ss3$P

    if(is.null(n3)){ # sample sizes from 3rd GWAS
      ss.out$n.ss3 = ss3$N
    }else{
      ss.out$n.ss3 = n3
    }
  }

  ss.out = na.omit(ss.out)  # A1 A2 may not exactly match, we may have NAs for effect sizes
  ss.out = filter(ss.out, CHR %in% c(1:22))  # use autosome only

  cat("After aligning alleles, ", format(nrow(ss.out), big.mark = ","), " autosome SNPs are left.\n", sep = "")


  # calculate direct and indirect effect sizes using the aligned effect sizes and se:
  # since they are properly aligned, we can stored them in vectors to speed up calculations
  # they're named following:
  # {beta, var, se, p, n}.{own, ss2, ss3, dir, ind}

  beta.own = ss.out$beta.own  # input effect sizes
  beta.ss2 = ss.out$beta.ss2

  se.own = ss.out$se.own  # SE for the input effect size
  se.ss2 = ss.out$se.ss2

  var.own = se.own^2  # variances for the input effect sizes
  var.ss2 = se.ss2^2

  n.own = ss.out$n.own  # input sample sizes
  n.ss2 = ss.out$n.ss2

  alpha = ss.out$alpha  # alpha for the SNP kept in ss.out

  cat("\nApplying the framework at ", format(Sys.time(), "%X, %a. %b. %d, %Y"), sep = "")
  cat("\nCalculating effect sizes, standard errors, p-values, and effective sample sizes for the direct and indirect effects.\n")


  if(mode == 1){  # most general case, 3 inputs: own, maternal, and paternal GWASs

    # the 3rd input sumstats (paternal GWAS):
    beta.ss3 = ss.out$beta.ss3
    se.ss3 = ss.out$se.ss3
    n.ss3 = ss.out$n.ss3

    var.ss3 = se.ss3^2

    # effect sizes for direct and indirect effects:
    beta.dir = (2+alpha)*beta.own - beta.ss2 - beta.ss3
    beta.ind = (2+alpha)*(beta.ss2 + beta.ss3 - (1+alpha)*beta.own)/(2+2*alpha)

    # covariances due to sample overlap and/or assortative mating:
    cov12 = l12*se.own*se.ss2  # covariances among effect sizes due to sample overlap
    cov13 = l13*se.own*se.ss3  # LDSC intercept gives correlation
    cov23 = l23*se.ss2*se.ss3  # correlated due to assortative mating between spouses
                               # note, assortative mating alpha is captured by LDSC intercept l23

    # variances for direct/indirect effects:
    var.dir = var.own*(2+alpha)^2 + var.ss2 + var.ss3 - 2*(2+alpha)*cov12 - 2*(2+alpha)*cov13 + 2*cov23
    var.ind = (var.ss2 + var.ss3 + var.own*(1+alpha)^2 + 2*cov23 - 2*(1+alpha)*cov12 - 2*(1+alpha)*cov13)*((2+alpha)/(2+2*alpha))^2

    # effective sample sizes for direct/indirect effects:
    n.dir = 1/((2+alpha)^2/n.own + 1/n.ss2 + 1/n.ss3 - 2*(2+alpha)*l12*sqrt(1+alpha/2)/(sqrt(n.own)*sqrt(n.ss2)) - 2*(2+alpha)*l13/(sqrt(n.own)*sqrt(n.ss3)) + 2*l23/(sqrt(n.own)*sqrt(n.ss3)))
    n.ind = ((2+2*alpha)/(2+alpha))^2/(1/n.ss2 + 1/n.ss3 + (1+alpha)^2/n.own + 2*l23/(sqrt(n.ss2)*sqrt(n.ss3)) - 2*(1+alpha)*l12/(sqrt(n.own)*sqrt(n.ss2)) - 2*(1+alpha)*l13/(sqrt(n.own)*sqrt(n.ss3)))


    # could also output indirect maternal/paternal/parental effect sumstats
    # effect sizes for indirect maternal and indirect paternal effects:
    beta.ind.mat = (3-alpha^2)*beta.ss2/(2*(1-alpha^2)) + beta.ss3*(1-2*alpha-alpha^2)/(2*(1-alpha^2)) - beta.own*(1+alpha/2)
    beta.ind.pat = (3-alpha^2)*beta.ss3/(2*(1-alpha^2)) + beta.ss2*(1-2*alpha-alpha^2)/(2*(1-alpha^2)) - beta.own*(1+alpha/2)

    ss.out$beta.ind.mat = beta.ind.mat
    ss.out$beta.ind.pat = beta.ind.pat

    # variances for indirect maternal/paternal effects:
    var.ind.mat = var.ss2*((3-alpha^2)/(2-2*alpha^2))^2 + var.ss3*((1-2*alpha-alpha^2)/(2-2*alpha^2))^2 + var.own*(1+alpha/2)^2 + cov23*2*(3-alpha^2)*(1-2*alpha-alpha^2)/(2-2*alpha^2)^2 - cov12*2*(3-alpha^2)*(1+alpha/2)/(2-2*alpha^2) - cov13*2*(1-2*alpha-alpha^2)*(1+alpha/2)/(2-2*alpha^2)
    var.ind.pat = var.ss3*((3-alpha^2)/(2-2*alpha^2))^2 + var.ss2*((1-2*alpha-alpha^2)/(2-2*alpha^2))^2 + var.own*(1+alpha/2)^2 + cov23*2*(3-alpha^2)*(1-2*alpha-alpha^2)/(2-2*alpha^2)^2 - cov13*2*(3-alpha^2)*(1+alpha/2)/(2-2*alpha^2) - cov12*2*(1-2*alpha-alpha^2)*(1+alpha/2)/(2-2*alpha^2)

    # effective sample sizes for indirect maternal and indirect paternal effects:
    n.ind.mat = 1/(((3-alpha^2)/(2-2*alpha^2))^2/n.ss2 + ((1-2*alpha-alpha^2)/(2-2*alpha^2))^2/n.ss3 + (1+alpha/2)/n.own + l23*(3-alpha^2)*(1-2*alpha-alpha^2)/(2*(1-alpha^2)^2*(sqrt(n.ss2)*sqrt(n.ss3))) - l12*(3-alpha^2)*(1+alpha/2)/((1-alpha^2)*sqrt(n.own)*sqrt(n.ss2)) - l13*(1-2*alpha-alpha^2)*(1+alpha/2)/((1-alpha^2)*sqrt(n.own)*sqrt(n.ss3)))
    n.ind.pat = 1/(((3-alpha^2)/(2-2*alpha^2))^2/n.ss3 + ((1-2*alpha-alpha^2)/(2-2*alpha^2))^2/n.ss2 + (1+alpha/2)/n.own + l23*(3-alpha^2)*(1-2*alpha-alpha^2)/(2*(1-alpha^2)^2*(sqrt(n.ss2)*sqrt(n.ss3))) - l13*(3-alpha^2)*(1+alpha/2)/((1-alpha^2)*sqrt(n.own)*sqrt(n.ss3)) - l12*(1-2*alpha-alpha^2)*(1+alpha/2)/((1-alpha^2)*sqrt(n.own)*sqrt(n.ss2)))

    n.ind.mat = ifelse(n.ind.mat >= 0, n.ind.mat, NA)
    n.ind.pat = ifelse(n.ind.pat >= 0, n.ind.pat, NA)

    n.ind.mat = round(n.ind.mat)  # convert to integers
    n.ind.pat = round(n.ind.pat)

    # give warning if varinces are not all positive:
    if(any(var.ind.mat <= 0, na.rm = T)){
      cat("\nWarning:", sum(var.ind.mat <= 0, na.rm = T), "SNPs have non-positive variances for the indirect maternal effect!\n")
    }
    if(any(var.ind.pat <= 0, na.rm = T)){
      cat("\nWarning:", sum(var.ind.pat <= 0, na.rm = T), "SNPs have non-positive variances for the indirect paternal effect!\n")
    }

    var.ind.mat = ifelse(var.ind.mat > 0, var.ind.mat, NA)
    var.ind.pat = ifelse(var.ind.pat > 0, var.ind.pat, NA)

    se.ind.mat = sqrt(var.ind.mat)
    se.ind.pat = sqrt(var.ind.pat)

    ss.out$se.ind.mat = se.ind.mat
    ss.out$se.ind.pat = se.ind.pat

    # p-values for indirect maternal/paternal effects:
    p.ind.mat = pnorm(abs(beta.ind.mat/se.ind.mat), lower.tail = F)*2
    p.ind.pat = pnorm(abs(beta.ind.pat/se.ind.pat), lower.tail = F)*2

    ss.out$p.ind.mat = p.ind.mat
    ss.out$p.ind.pat = p.ind.pat

    # effective sample sizes for indirect maternal/paternal effects:
    ss.out$n.ind.mat = n.ind.mat
    ss.out$n.ind.pat = n.ind.pat

  }else if(mode == 2){  # if 2 inputs: own + parental GWASs

    # effect sizes for direct and indirect effect:
    beta.dir = beta.own*(2+alpha) - 2*beta.ss2
    beta.ind = beta.ss2*(2+alpha)/(1+alpha) - beta.own*(1+alpha/2)

    # variances:
    cov12 = l12*se.own*se.ss2

    var.dir = var.own*(2+alpha)^2 + 4*var.ss2 - 4*(2+alpha)*cov12
    var.ind = var.ss2*((2+alpha)/(1+alpha))^2 + var.own*(1+alpha/2)^2 - cov12*(2+alpha)^2/(1+alpha)

    # effective sample sizes:
    n.dir = 1/((2+alpha)^2/n.own + 4/n.ss2 - 4*(2+alpha)*l12/(sqrt(n.own)*sqrt(n.ss2)))
    n.ind = 1/(((2+alpha)/(1+alpha))^2/n.ss2 + (1+alpha/2)^2/n.own - l12*(2+alpha)^2/((1+alpha)*sqrt(n.own)*sqrt(n.ss2)))

    # note, for this setup, it may not be possible to get indirect parental effect sizes

  }else{  # mode == 3, 2 inputs: own + maternal/paternal GWAS

    # effect sizes for direct and indirect effects:
    beta.dir = (beta.own*(2+alpha) - beta.ss2*(1+alpha))*2/(3-alpha^2)
    beta.ind = (beta.ss2 - beta.own*(1+alpha)/2)*(2+alpha)/(3-alpha^2)

    # variances:
    cov12 = l12*se.own*se.ss2

    var.dir = (var.own*(2+alpha)^2 + var.ss2*(1+alpha)^2 - 2*(1+alpha)*(2+alpha)*cov12)*4/(3-alpha^2)^2
    var.ind = (var.ss2 + var.own*((1+alpha)/2)^2 - (1+alpha)*cov12)*((2+alpha)/(3-alpha^2))^2

    # effective sample sizes:
    n.dir = ((3-alpha^2)/2)^2/((2+alpha)^2/n.own + (1+alpha)^2/n.ss2 - 2*(1+alpha)*(2+alpha)*l12/(sqrt(n.own)*sqrt(n.ss2)))
    n.ind = ((3-alpha^2)/(2+alpha))^2/(1/n.ss2 + ((1+alpha)/2)^2/n.own - (1+alpha)*l12/(sqrt(n.own)*sqrt(n.ss2)))

    # could also output indirect maternal or paternal (depending who's in the 2nd input) effect sumstats:
    beta.ind.ss2 = beta.ind*2

    if(any(var.ind <= 0, na.rm = T)){
      cat("\nWarning:", sum(var.ind <=0, na.rm = T), "SNPs have non-positive variances for the indirect effect!\n")
    }

    var.ind = ifelse(var.ind > 0, var.ind, NA)
    se.ind.ss2 = 2*sqrt(var.ind)
    p.ind.ss2 = pnorm(abs(beta.ind.ss2/se.ind.ss2), lower.tail = F)*2

    n.ind.ss2 = 1/(4*((2+alpha)/(3-alpha^2))^2*(1/n.ss2 + (0.5+alpha/2)^2/(n.own*(1+alpha/2)) - l12*(1+alpha)*sqrt(1+alpha/2)/(sqrt(n.own)*sqrt(n.ss2))))
    n.ind.ss2 = ifelse(n.ind.ss2 >= 0, n.ind.ss2, NA)
    n.ind.ss2 = round(n.ind.ss2)

    ss.out$beta.ind.ss2 = beta.ind.ss2
    ss.out$se.ind.ss2 = se.ind.ss2
    ss.out$p.ind.ss2 = p.ind.ss2
    ss.out$n.ind.ss2 = n.ind.ss2

  }

  # give warning if variances are not all positive:
  if(any(var.dir <= 0, na.rm = T)){
    cat("\nWarning:", sum(var.dir <= 0, na.rm = T), "SNPs have non-positive variances for the direct effect!\n")
  }
  if(any(var.ind <= 0, na.rm = T)){
    cat("\nWarning:", sum(var.ind <=0, na.rm = T), "SNPs have non-positive variances for the indirect effect!\n")
  }
  # give warning if effective sample sizes are not all positive:
  if(any(n.dir <= 0, na.rm = T)){
    cat("\nWarning:", sum(n.dir <= 0, na.rm = T), "SNPs have non-positive effective sample sizes for the direct effect!\n")
  }
  if(any(n.ind <= 0, na.rm = T)){
    cat("\nWarning:", sum(n.ind <= 0, na.rm = T), "SNPs have non-positive effective sample sizes for the indirect effect!\n")
  }


  # standard errors:
  var.dir = ifelse(var.dir > 0, var.dir, NA)
  var.ind = ifelse(var.ind > 0, var.ind, NA)

  se.dir = sqrt(var.dir)
  se.ind = sqrt(var.ind)

  # p-values for direct and indirect effects:
  p.dir = pnorm(abs(beta.dir/se.dir), lower.tail = F)*2
  p.ind = pnorm(abs(beta.ind/se.ind), lower.tail = F)*2

  # effective sample sizes:
  n.dir = ifelse(n.dir >= 0, n.dir, NA)  # in case negative results
  n.dir = round(n.dir) # change to integers

  n.ind = ifelse(n.ind >= 0, n.ind, NA)  # in case negative results
  n.ind = round(n.ind) # change to integers


  # add results to the output data.frame:
  ss.out$beta.dir = beta.dir  # direct effect sumstats:
  ss.out$se.dir = se.dir
  ss.out$p.dir = p.dir
  ss.out$n.dir = n.dir

  ss.out$beta.ind = beta.ind  # indirect effect sumstats:
  ss.out$se.ind = se.ind
  ss.out$p.ind = p.ind
  ss.out$n.ind = n.ind


  ### output the data.frame that contains all the results:

  # direct effect sumstats:
  ss.dir = select(ss.out, CHR, SNP, BP, A1, A2,
                  BETA = beta.dir, SE = se.dir, P = p.dir, N = n.dir)

  # indirect effect sumstats:
  ss.ind = select(ss.out, CHR, SNP, BP, A1, A2,
                  BETA = beta.ind, SE = se.ind, P = p.ind, N = n.ind)

  if(mode == 1){ # if 3 inputs, can also have sumstats for indirect maternal/paternal effects

    # indirect maternal effect sumstats:
    ss.ind.mat = select(ss.out, CHR, SNP, BP, A1, A2,
                        BETA = beta.ind.mat, SE = se.ind.mat, P = p.ind.mat, N = n.ind.mat)

    # indirect paternal effect sumstats:
    ss.ind.pat = select(ss.out, CHR, SNP, BP, A1, A2,
                        BETA = beta.ind.pat, SE = se.ind.pat, P = p.ind.pat, N = n.ind.pat)

  }

  if(mode == 3){  # under mode 3, can also have sumstats for the indirect maternal/paternal effect

    # indirect maternal/paternal effect sumstats
    ss.ind.ss2 = select(ss.out, CHR, SNP, BP, A1, A2,
                        BETA = beta.ind.ss2, SE = se.ind.ss2, P = p.ind.ss2, N = n.ind.ss2)

  }

  cat("Calculations are done at ", format(Sys.time(), "%X, %a. %b. %d, %Y"), sep = "")
  cat("\n\nMost significant SNP for direct effect is ", ss.dir$SNP[which.min(ss.dir$P)], " with P-value = ", format(min(ss.dir$P, na.rm = T), scientific = T, digits = 3), sep = "")
  cat("\nMost significant SNP for indirect effect is ", ss.ind$SNP[which.min(ss.ind$P)], " with P-value = ", format(min(ss.ind$P, na.rm = T), scientific = T, digits = 3), sep = "")


  if(!is.null(OutDir)){
    cat("\n\nWritting the results...")

    # fwrite from data.table package should be much quicker to write:
    # also, let an empty string be written as NA

    fwrite(ss.dir, file = paste0(OutDir, "/direct_effect.sumstats.gz"),
           col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

    fwrite(ss.ind, file = paste0(OutDir, "/indirect_effect.sumstats.gz"),
           col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

    fwrite(ss.out, file = paste0(OutDir, "/all_aligned.sumstats.gz"),
           col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

    if(mode == 1){
      # output the indirect maternal and indirect paternal effects:
      fwrite(ss.ind.mat, file = paste0(OutDir, "/indirect_maternal_effect.sumstats.gz"),
             col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

      fwrite(ss.ind.pat, file = paste0(OutDir, "/indirect_paternal_effect.sumstats.gz"),
             col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

    }

    if(mode == 3){
      # output the sumstats for the indirect maternal or paternal effects (depending on the 2nd input)
      fwrite(ss.ind.ss2, file = paste0(OutDir, "/indirect_ss2_effect.sumstats.gz"),
             col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

    }

    cat("\nDone! The direct, indirect effect sumstats, and also the aligned input sumtats (.gz files) are written to\n\n\t", OutDir, "\n", sep = "")
  }

  cat("\nAlso returns a data.frame that contains both input and output sumstats.\n\n")
  return(ss.out)  # data.frame that contains both input and output sumstats

} # end of the function



