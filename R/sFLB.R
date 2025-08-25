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
#' @import data.table
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
sFLB = function (x, carriers = NULL, noncarriers = NULL, affected = NULL,
                 unknown = NULL, penetrances = NULL, liability = NULL, ...)
{

  # Calculate p(g)
  distr = rareDistr(x, typable = x$ID, proband = carriers, ...)
  observed = apply(distr, 1, function(row) all(row[carriers] == 1) & all(row[noncarriers] == 0))

  # Assign penetrances
  affvec = setNames(rep(1, pedsize(x)), x$ID)
  affvec[unknown] = NA
  affvec[affected] = 0
  if (is.vector(penetrances))
    distr[, (x$ID) := lapply(x$ID, function(i) affvec[i] - penetrances[distr[[i]]+1])]
  else {
    if (is.null(liability))
      liability = rep(1, pedsize(x))
    names(liability) = x$ID
    distr[, (x$ID) := lapply(x$ID, function(i) affvec[i] - penetrances[liability[i], distr[[i]]+1])]
  }

  # Calculate p(g|y)
  distr[, p_H1 := Reduce("*", lapply(.SD, function(i) ifelse(is.na(i), 1, i)))]
  distr[, p_H1 := p_H1 / sum(p_H1)]

  # BF = p(g|y) / p(g)
  distr = distr[observed]
  bf = abs(sum(distr$p_H1) / sum(distr$p_H0))

  bf
}
