#' segregatr: Segregation Analysis for Identifying Pathogenic Variants
#'
#' An implementation of the full-likelihood Bayes factor (FLB) for evaluating
#' segregation evidence in clinical medical genetics. The method was introduced
#' by Thompson et al. (2003), and further popularised by Bayrak-Toydemir et al.
#' (2008). This implementation allows custom penetrance values and liability
#' classes, and includes specialised pedigree visualisations.
#'
#' @references Thompson D, Easton DF, Goldgar DE. *A full-likelihood method for
#'   the evaluation of causality of sequence variants from family data.* Am J
#'   Hum Genet, 2003. \doi{10.1086/378100}.
#'
#'   Bayrak-Toydemir et al. *Likelihood ratios to assess genetic evidence for
#'   clinical significance of uncertain variants: Hereditary hemorrhagic
#'   telangiectasia as a model.* Exp Mol Pathol, 2008.
#'   \doi{10.1016/j.yexmp.2008.03.006}.
#'
#' @docType package
#' @import pedtools
#' @importFrom pedprobr likelihood
#'
#' @name segregatr
NULL
