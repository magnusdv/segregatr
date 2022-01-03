#' Full-likelihood Bayes factor
#'
#' Computes the Bayes factor for co-segregation, as described by Thompson et al.
#' (2003).
#'
#' @param x A [pedtools::ped()] object.
#' @param carriers A character vector (or coercible to such), containing the ID
#'   labels of pedigree members known to carry one copy of the variant in question.
#' @param homozygous A character vector (or coercible to such), containing the ID
#'   labels of pedigree members known to carry two copies of the variant in question.
#' @param noncarriers A character vector (or coercible to such), containing the
#'   ID labels of pedigree members known *not* to carry the variant in question.
#' @param freq A single number strictly between 0 and 1: the population
#'   frequency of the observed allele.
#' @param affected The affected pedigree members.
#' @param unknown Pedigree members with unknown affection status.
#' @param proband The ID label of the proband. This person must also be in both
#'   `carriers` and `affected`.
#' @param penetrances Either a numeric vector of length 3, corresponding to
#'   `(f0, f1, f2)` or a matrix or data frame with 3 columns. Each row contains
#'   the penetrance values of a liability class.
#' @param liability A vector of length `pedsize(x)`, containing for each
#'   pedigree member the row number of `penetrances` which should be used for
#'   that individual. (If `penetrances` is just a vector, it will be used for
#'   all classes.) If `liability` is NULL (the default), it is set to `1` for
#'   all individuals.
#' @param loopBreakers (Relevant only if `x` has loops.) A vector of ID labels
#'   indicating loop breakers. The default value (NULL) initiates automatic loop
#'   breaking, which is recommended in most cases.
#' @param Xchrom A logical, indicating if a model of X-linked inheritance should be applied.
#' @param details A logical, indicating if detailed output should be returned
#'   (for debugging purposes).
#' @param plot A logical.
#' @param ... Optional plot parameters passed on to [pedtools::plot.ped()].
#'
#' @references Thompson D, Easton DF, Goldgar DE. *A full-likelihood method for
#'   the evaluation of causality of sequence variants from family data.* Am J
#'   Hum Genet, 2003. \doi{10.1086/378100}.
#' @return A positive number. If `details = TRUE`, a list of intermediate
#'   results is returned.
#'
#' @examples
#'
#' x = nuclearPed(2)
#'
#' FLB(x, carriers = 3:4, aff = 3:4, unknown = 1:2,
#'     freq = 0.0001, penetrances = c(0, 1, 1), proband = 3)
#'
#' @export
FLB = function(x, carriers = NULL, homozygous = NULL, noncarriers = NULL, freq,
               affected, unknown = NULL, proband,
               penetrances, liability = NULL, loopBreakers = NULL, Xchrom = FALSE,
               details = FALSE, plot = FALSE, ...) {

  ### Note to self: Don't mess with the order of input checks.

  if(!is.ped(x))
    stop2("The first argument must be a `ped` object")

  if(!is.null(x$LOOP_BREAKERS))
    stop2("Pedigrees with pre-broken loops are not allowed as input to `FLB()`")

  labs = labels(x)
  males = which(x$SEX==1)

  if(length(proband) == 0 || proband == "")
    stop2("A proband must be specified")
  if(!proband %in% labs)
    stop2("Unknown proband: ", proband)
  if(!proband %in% affected)
    stop2("The proband must be affected")
  if(!proband %in% c(carriers, homozygous))
    stop2("The proband must be a carrier")

  for(ids in list(proband, affected, unknown, carriers, homozygous, noncarriers)) {
    validID = ids %in% labs
    if(!all(validID))
      stop2("Unknown ID label: ", ids[!validID])
  }

  if(length(err1 <- intersect(affected, unknown)))
    stop2("Individual specified as both affected and unknown: ", err1)

  if(length(err2 <- intersect(c(carriers, homozygous), noncarriers)))
    stop2("Individual specified as both a carrier and a non-carrier: ", err2)

  if(length(err3 <- intersect(carriers, homozygous)))
    stop2("Individual specified as both heterozygous and homozygous carrier: ", err3)

  if(Xchrom & length(err4 <- intersect(males, homozygous)))
    stop2("Male individual specified as a homozygous carrier in an X-linked inheritance model: ", err4)

  if(is.null(freq) || is.na(freq))
    stop2("An allele frequency must be specified")
  if(!is.numeric(freq) || length(freq) != 1 || freq <= 0 || freq >= 1)
    stop2("The allele frequency must be a single number strictly between 0 and 1")

  if(plot) {
    plotSegregation(x, affected, unknown, proband, carriers, homozygous, noncarriers, ...)
  }

  # Empty marker and disease locus)
  dis = m = mProband = marker(x, afreq = c(a = 1 - freq, b = freq), chrom = ifelse(Xchrom, 'X', NA))

  # Full marker
  if (Xchrom) {
    genotype(m, intersect(carriers, males)) = "b"
    genotype(m, setdiff(carriers, males)) = c("a", "b")
    genotype(m, intersect(noncarriers, males)) = "a"
    genotype(m, setdiff(noncarriers, males)) = c("a", "a")
  } else {
    genotype(m, carriers) = c("a", "b")
    genotype(m, noncarriers) = c("a", "a")
  }
  genotype(m, homozygous) = c("b", "b")

  # Proband marker
  genotype(mProband, proband) = genotype(m, proband)

  # Break loops if necessary
  if(hasUnbrokenLoops(x)) {
    x = breakLoops(setMarkers(x, list(dis, m, mProband)), loopBreakers = loopBreakers, verbose = FALSE)

    # Extract updated markers
    dis = x$MARKERS[[1]]
    m = x$MARKERS[[2]]
    mProband = x$MARKERS[[3]]

    # Add duplicates to input vectors, were needed
    lb = x$ID[x$LOOP_BREAKERS[, 'orig']]

    if(length(lbAff <- intersect(affected, lb)))
      affected = unique.default(c(affected, paste0("=", lbAff)))

    if(length(lbUn <- intersect(unknown, lb)))
      unknown = unique.default(c(unknown, paste0("=", lbUn)))

    # TODO: Fix this (easy)
    if(!is.null(liability))
      stop2("Liability classes are not yet implemented when `x` has loops")

    # TODO: Handle proband vs loop breaking better
    if(proband %in% lb)
      stop2("The proband cannot be a loop breaker; please select different loop breaker(s)")
  }

  # Affection status vector, sorted along labs
  aff = logical(pedsize(x))
  aff[internalID(x, affected)] = TRUE
  aff[internalID(x, unknown)] = NA

  # Penetrances and liability classes
  penetMat = fixPenetrances(penetrances)
  liability = liability %||% rep_len(1, pedsize(x))
  if(!all(liability %in% 1:nrow(penetMat)))
    stop2("Illegal liability class: ", setdiff(liability, 1:nrow(penetrances)))


  # Setup for likelihood under causal hypothesis
  peelOrder = peelingOrder(x)
  peelingProcess = pedprobr:::peelingProcess

  if (Xchrom) {
    peeler = function(x, locus) function(dat, sub) pedprobr:::.peel_M_X(dat, sub, SEX = x$SEX)
    startCausal = function(x, locus) startdata_causative_X(x, locus, aff = aff, penetMat = penetMat, liability = liability)
  }
  else {
    peeler = function(x, locus) function(dat, sub) pedprobr:::.peel_M_AUT(dat, sub)
    startCausal = function(x, locus) startdata_causative(x, locus, aff = aff, penetMat = penetMat, liability = liability)
  }

  likelihoodCausal = function(x, locus) {
    peelingProcess(x, locus, startdata = startCausal,
                   peeler = peeler(x, locus), peelOrder = peelOrder)
  }
  # Main Bayes factor
  numer1 = likelihoodCausal(x, m)
  likDis = likelihoodCausal(x, dis)
  likM = likelihood(x, m)

  denom1 = likDis * likM
  BF1 = numer1/denom1

  # Correction factor
  likMproband = likelihood(x, mProband)
  numer2 = likDis * likMproband

  denom2 = likelihoodCausal(x, mProband)

  BF2 = numer2/denom2

  FLB = BF1 * BF2

  if(details)
    return(list(c(FLB = FLB, BF1 = BF1, BF2 = BF2),
                c(numer1 = numer1, denom1 = denom1, numer2 = numer2, denom2 = denom2),
                c(likM = likM, likMproband = likMproband, likDis = likDis)))
  FLB
}
