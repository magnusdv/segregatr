#' 'Shared' full-likelihood Bayes factor
#'
#' Computes the shared Bayes factor for co-segregation, assuming autosomal
#' dominant inheritance and a single introduction of the variant.
#'
#' @param x A [pedtools::ped()] object.
#' @param carriers A character vector (or coercible to such), containing the ID
#'   labels of pedigree members known to carry one copy of the variant in
#'   question.
#' @param noncarriers A character vector (or coercible to such), containing the
#'   ID labels of pedigree members known *not* to carry the variant in question.
#' @param affected The affected pedigree members.
#' @param unknown Pedigree members with unknown affection status.
#' @param penetrances A numeric vector of length 3 `(f0, f1, f2)`, or a
#'   matrix-like with 3 columns, where row `i` contains the penetrances of
#'   liability class `i`.
#' @param liability A vector of length `pedsize(x)`, containing for each
#'   pedigree member the row number of `penetrances` which should be used for
#'   that individual. If unnamed, it is assumed that the individuals are taken
#'   in order. (If `penetrances` is just a vector, it will be used for all
#'   classes.) If `liability` is NULL (the default), it is set to `1` for all
#'   individuals.
#' @param ... Optional parameters passed on to [segregatr::rareDistr()].
#'
#' @references Ratajska A, Vigeland MD, Wirgenes KV, *et al*. *The use of
#'   segregation analysis in interpretation of sequence variants in SMAD3: A
#'   case report.* Mol Genet Genomic Med, 2023. \doi{10.1002/mgg3.2107}.
#'
#' @return A positive number, the sFLB score.
#'
#' @examples
#'
#' ### Case 1
#'
#' x = halfSibPed(nch1 = 2, type = "maternal")
#'
#' sFLB(x, unknown = 1:3, affected = 4:6, carriers = 4:6,
#'      noncarriers = NULL, penetrances = c(0.1, 0.5, 0.5))
#'
#'
#' ### Ratajska et al. (2023), Family B
#'
#' y = nuclearPed(5, sex = c(2,1,1,1,1)) |>
#'       addDaughter(parents = 3, verbose = FALSE) |>
#'       relabel("asPlot")
#'
#' sFLB(y, unknown = NULL, affected = c(1,4,5,8), carriers = c(1,4,8),
#'      noncarriers = 6:7, penetrances = c(0.01, 0.9, 0.9))
#'
#' @export
sFLB = function(x, carriers = NULL, noncarriers = NULL, affected = NULL,
                unknown = NULL, penetrances = NULL, liability = NULL, ...) {

  # Coerce inputs
  carriers = as.character(carriers)
  noncarriers = as.character(noncarriers)
  affected = as.character(affected)
  unknown = as.character(unknown)

  # Calculate p(g)
  distr = rareDistr(x, typable = x$ID, proband = carriers, ...)

  # Observed genotypes
  observed = rowSums(distr[, carriers, drop = FALSE] == 1) == length(carriers) &
    rowSums(distr[, noncarriers, drop = FALSE] == 0) == length(noncarriers)

  # Assign penetrances
  affvec = setNames(rep(1, pedsize(x)), x$ID)
  affvec[unknown] = NA
  affvec[affected] = 0

  distr_aff = distr[, x$ID, drop = FALSE]

  if (is.vector(penetrances)) {
    for (i in x$ID) {
      distr_aff[[i]] = affvec[i] - penetrances[distr_aff[[i]] + 1]
    }
  } else {
    if (is.null(liability))
      liability = rep(1, pedsize(x))
    names(liability) = x$ID
    for (i in x$ID) {
      distr_aff[[i]] = affvec[i] - penetrances[liability[i], distr_aff[[i]] + 1]
    }
  }

  # Calculate p(g|y)
  lik = apply(distr_aff, 1, function(row) {
    prod(sapply(row, function(i) if (is.na(i)) 1 else i))
  })
  p_H1 = lik * distr$p_H0
  p_H1 = p_H1 / sum(p_H1)

  # BF = p(g|y) / p(g)
  distr_obs = distr[observed, , drop = FALSE]
  p_H1_obs = p_H1[observed]

  bf = abs(sum(p_H1_obs) / sum(distr_obs$p_H0))
  bf
}
