#
# Functions for multivariate QTL scans, including data manipulation.
#
# Author: Kai Fox [kaijordanfox@gmail.com]
# Date: Summer 2019
#



#################################################################
# -- Basic Multivariate Scanning  ----------------------------  #
#################################################################


#' Evaluate fixed-effect genetic model at a single locus.
#'
#' For full explanation of the input data types, see the data
#' structures vignette. \code{trait.table} and \code{geno.table}
#' are as referred to throughout \code{aseQTL}.
#'
#' @param trait Name of the trait in the trait table to be mapped.
#' @param locus Index of the locus in the locus table to be mapped
#' @param trait.table Trait table to be mapped from, containing
#'     response matrices for each sample/subject.
#' @param geno.table Genotype probability table, giving
#'     regression variables for the fixed-effect genetic model.
#' @return Scalar \eqn{p}.
locus.p <- function(trait, locus, trait.table, geno.table) {
  X <- do.call(rbind, lapply(names(geno.table), function(s)
    geno.table[[basename(s)]][locus,]))
  Y <- as.matrix(do.call(rbind, lapply(names(trait.table), function (s)
    trait.table[[s]][trait,])))
  model.geno <- lm(Y ~ A + B + C + D + E + F + G + H, data = data.frame(X))
  model.null <- lm(Y ~ 1, data = data.frame(X))
  return(anova(model.geno, model.null)[-1,"Pr(>F)"])
}


#' Helper interface to \code{locus.p}
#'
#' @param arg A list with two elements, \code{locus} and \code{trait},
#'     that describe respective arguments to \code{locus.p}.
#' @param all.traits,all.genoprobs Corresponding arguments to \code{locus.p}
#' @return The significance measurement from \code{locus.p}.
.scan.helper <- function(arg, all.traits, all.genoprobs, func)
  func(arg$trait, arg$locus, all.traits, all.genoprobs)


#' Perform a multivarate QTL scan
mvscan <- function(all.traits, all.genoprobs, map, cores = 1, loc.cap = NULL, traits = NULL) {
  n.loc <- dim(all.genoprobs[[names(all.genoprobs)[1]]])[1]
  if (!is.null(loc.cap)) n.loc <- min(n.loc, loc.cap)
  if (is.null(traits)) traits <- names(all.traits)

  # Construct iteration argument
  l <- rep(1:n.loc, times = length(traits))
  t <- rep(traits, each = n.loc)
  iter <- lapply(seq_along(l), function(i) list(locus=l[i], trait=t[i]))

  # Run `locus.p` along the iter object via the helper above
  if (cores > 1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("Package 'parallel' required for cores > 1", call = F)
    }
    cl <- parallel::makeCluster(cores)
    all.ps <- parallel::parLapply(
      cl, iter, .scan.helper,
      all.traits = all.traits,
      all.genoprobs = all.genoprobs,
      func = locus.p)
  } else {
    all.ps <- lapply(
      iter, .scan.helper,
      all.traits = all.traits,
      all.genoprobs = all.genoprobs,
      func = locus.p)
  }

  # Reformat the outputs of `locus.p` into a (loci x trait) matrix
  all.ps <- matrix(all.ps, nrow = n.loc, ncol = length(traits))
  rownames(all.ps) <- map[1:n.loc, "marker"]
  colnames(all.ps) <- as.character(traits)
  return(all.ps)
}




#################################################################
# -- Data Manipulation  --------------------------------------  #
#################################################################


#' Load an example dataset of DO/CC allele-specific expression
#'
#' This data is expected to be organized in sample/subject folders
#' with each folder containing a collection of tab-separated-value
#' files which ultimately house the data. This function is not
#' expected to be of use except in testing, until we can integrate
#' with a more standardized data format.
#'
#' @param data.loc Root folder of the data.
#' @param genoprob.fname Within subject folders, name of the '.tsv'
#'     file containing genotype probabilities interpolated out
#'     from markers.
#' @param ase.fname Within subejct folders, name of the file giving
#'     (in tab-separated-value format) the allelic expression of
#'     individuals for a collection of genes. Until more robust
#'     regression methods can be developed this should hold TPM values.
#' @param gene Name of a gene to select from the TPM files.
#' @return List containing values \code{T}, which is the regression
#'     responses and \code{G}, which is another list mapping
#'     sample/subject keys to dataframes of genotype probabilities.
#'     It is unlikely that either \code{Y} or \code{G} will need to be
#'     directly manipulated but rather will be the input to other
#'     \code{aseQTL} functions.
.load.nested.ase <- function(data.loc, genoprob.fname, ase.fname, gene) {
  samps <- list.dirs(data.loc, recursive = F)
  samps <- samps[seq(1, 50)]
  HAPLOS <- c("A","B","C","D","E","F","G","H")
  trait.table <- lapply(samps, function(s) {
    data <- read.csv(file.path(s, ase.fname), sep = '\t')
    rownames(data) <- data[,"locus"]
    return(data[,HAPLOS])
  })
  names(trait.table) <- as.character(lapply(samps, basename))
  geno.table <- lapply(samps, function(s) {
    return(setNames(read.csv(file.path(s, genoprob.fname), sep = '\t'), HAPLOS))
  })
  names(geno.table) <- as.character(lapply(samps, basename))
  return(list(T=trait.table, G=geno.table))
}


