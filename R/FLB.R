#' Full-likelihood Bayes factor
#'
#' Computes the Bayes factor for co-segregation, as described by Thompson et al
#' (2003).
#'
#' @param x A `ped` object
#' @param carriers A character vector (or coercible to such), containing the ID
#'   labels of pedigree members known to carry the variant in question.
#' @param noncarriers A character vector (or coercible to such), containing the
#'   ID labels of pedigree members known *not* to carry the variant in question.
#' @param freq A single number strictly between 0 and 1: the population
#'   frequency of the obsverved allele.
#' @param affected The affected pedigree members.
#' @param unknown Pedigree members with unknown affection status.
#' @param proband The ID label of the proband. This person must also be in both
#'   `carriers` and `affected`.
#' @param penetrances Either a numeric vector of length 3 (f0, f1, f2) or a data
#'   frame with 3 or more columns. Each row contains the penetrance values of a
#'   liability class. The first three columns are interpreted as penetrance
#'   values f0, f1, f2 respectively; additional columns are ignored.
#' @param liabilityClasses A vector of length `pedsize(x)`, containing for each
#'   pedigree member the row number of `penetrances` which should be used for
#'   that individual. (If `penetrances` is just a vector, it will be used for
#'   all classes.) If `liabilityClasses` is NULL (the default), it is set to `1`
#'   for all individuals.
#' @param details A logical, indicating if detailed output should be returned.
#' @param plot A logical.
#'
#' @return A positive number. If `details = TRUE`, a list of intermediate
#'   results is returned.
#'
#' @examples
#' FLB(pedtools::nuclearPed(2), carriers = 3:4, aff = 3:4, unknown = 1:2,
#'     freq = 0.0001, penetrances = c(0, 1, 1), proband = 3)
#'
#' @export
FLB = function(x, carriers, noncarriers = NULL, freq,
               affected, unknown = NULL, proband,
               penetrances, liabilityClasses = NULL,
               details = FALSE, plot = FALSE) {

  ### Note to self: Don't mess with the order of input checks.

  if(!is.ped(x))
    stop2("The first argument must be a `ped` object")

  labs = labels(x)

  if(length(proband) == 0 || proband == "")
    stop2("A proband must be specified")
  if(!proband %in% labs)
    stop2("The proband is not recognised as a pedigree member")
  if(!proband %in% affected)
    stop2("The proband must be affected")
  if(!proband %in% carriers)
    stop2("The proband must be a carrier")

  for(ids in list(proband, affected, unknown, carriers, noncarriers)) {
    validID = ids %in% labs
    if(!all(validID))
      stop2("Unknown ID label: ", ids[!validID])
  }

  if(length(err1 <- intersect(affected, unknown)))
    stop2(sprintf("Individual %s cannot be both *affected* and *unknown*", err1[1]))

  if(length(err2 <- intersect(carriers, noncarriers)))
    stop2(sprintf("Individual %s cannot be both a *carrier* and a *non-carrier*", err2[1]))

  if(is.null(freq) || is.na(freq))
    stop2("An allele frequency must be specified")
  if(!is.numeric(freq) || length(freq) != 1 || freq <= 0 || freq >= 1)
    stop2("The allele frequency must be a single number strictly between 0 and 1")


  # Affection status vector, sorted along labs
  aff = logical(pedsize(x))
  aff[internalID(x, affected)] = TRUE
  aff[internalID(x, unknown)] = NA

  # Empty marker (= disease locus)
  dis = m = mProband = marker(x, afreq = c(a = 1 - freq, b = freq))

  # Full marker
  genotype(m, carriers) = c("a", "b")
  genotype(m, noncarriers) = c("a", "a")

  # Proband marker
  genotype(mProband, proband) = c("a", "b")


  if(plot)
    plot(x, m, skip.empty = T, shaded = labs[aff], starred = proband)

  # Utility for setting up likelihood under causative hypothesis
  quickStart = function(locus)
    startdata_causative(x, marker = locus, aff = aff, penetrances = penetrances,
                        liabilityClasses = liabilityClasses)

  # Main Bayes factor
  peelOrder = peelingOrder(x)
  setup1 = list(informativeNucs = peelOrder, startdata = quickStart(m))
  numer1 = likelihood(x, m, setup = setup1)

  setupDis = list(informativeNucs = peelOrder, startdata = quickStart(dis))
  likDis = likelihood(x, dis, setup = setupDis)
  likM = likelihood(x, m)

  denom1 = likDis * likM
  BF1 = numer1/denom1

  # Correction factor
  likMproband = likelihood(x, mProband)
  numer2 = likDis * likMproband

  setup2 = list(informativeNucs = peelOrder, startdata = quickStart(mProband))
  denom2 = likelihood(x, mProband, setup = setup2)

  BF2 = numer2/denom2

  FLB = BF1 * BF2

  if(details)
    return(list(c(FLB = FLB, BF1 = BF1, BF2 = BF2),
                c(numer1 = numer1, denom1 = denom1, numer2 = numer2, denom2 = denom2),
                c(likM = likM, likMproband = likMproband, likDis = likDis)))
  FLB
}
