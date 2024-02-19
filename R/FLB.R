#' Full-likelihood Bayes factor
#'
#' Computes the Bayes factor for co-segregation, as originally described by
#' Thompson et al. (2003).
#'
#' @param x A [pedtools::ped()] object.
#' @param carriers A character vector (or coercible to such), containing the ID
#'   labels of pedigree members known to carry one copy of the variant in
#'   question.
#' @param homozygous A character vector (or coercible to such), containing the
#'   ID labels of pedigree members known to carry two copies of the variant in
#'   question.
#' @param noncarriers A character vector (or coercible to such), containing the
#'   ID labels of pedigree members known *not* to carry the variant in question.
#' @param freq A single number strictly between 0 and 1: the population
#'   frequency of the observed allele.
#' @param affected The affected pedigree members.
#' @param unknown Pedigree members with unknown affection status.
#' @param proband The ID label of the proband. This person must also be in both
#'   `carriers` and `affected`.
#' @param penetrances For autosomal models, a numeric vector of length 3 `(f0,
#'   f1, f2)`, or a matrix-like with 3 columns, where row `i` contains the
#'   penetrances of liability class `i`. For X-linked models, a list of two
#'   vectors named "male" and "female", of lengths 2 `(f0, f1)` and 3 `(f0, f1,
#'   f2)` respectively. Alternatively, each list entry may be matrix-like (with
#'   the same number of columns) where each row represents a liability class.
#' @param liability A vector of length `pedsize(x)`, containing for each
#'   pedigree member the row number of `penetrances` which should be used for
#'   that individual. If unnamed, it is assumed that the individuals are taken
#'   in order. (If `penetrances` is just a vector (or one for each sex in
#'   X-linked models), it will be used for all classes.) If `liability` is NULL
#'   (the default), it is set to `1` for all individuals.
#' @param loopBreakers (Relevant only if `x` has loops.) A vector of ID labels
#'   indicating loop breakers. The default value (NULL) initiates automatic loop
#'   breaking, which is recommended in most cases.
#' @param Xchrom A logical, indicating if a model of X-linked inheritance should
#'   be applied.
#' @param details A logical, indicating if detailed output should be returned
#'   (for debugging purposes).
#' @param plot A logical.
#' @param ... Optional plot parameters passed on to [pedtools::plot.ped()].
#'
#' @references Thompson D, Easton DF, Goldgar DE. *A full-likelihood method for
#'   the evaluation of causality of sequence variants from family data.* Am J
#'   Hum Genet, 2003. \doi{10.1086/378100}.
#'
#' @return A positive number, the FLB score. If `details = TRUE`, a list
#'   including intermediate results.
#'
#' @examples
#'
#' ### Autosomal dominant
#'
#' x = nuclearPed(2)
#'
#' FLB(x, carriers = 3:4, aff = 3:4, unknown = 1:2,
#'     freq = 0.0001, penetrances = c(0, 1, 1), proband = 3)
#'
#'
#' ### Autosomal recessive with phenocopies and reduced penetrance
#'
#' y = nuclearPed(4)
#'
#' FLB(y, carriers = 4:5, homozygous = 3, noncarriers = 6,
#'     aff = 3, unknown = 1:2, freq = 0.0001, proband = 3,
#'     penetrances = c(0.01, 0.01, 0.99), plot = TRUE)
#'
#'
#' ### X-linked recessive
#'
#' z = nuclearPed(3, sex = c(1, 1, 2)) |>
#'   addChildren(mother = 5, nch = 2, sex = 1:2)
#'
#' FLB(z, carriers = c(3, 7), nonc = 4, aff = c(3, 7), unknown = 1:2,
#'     freq = 0.0001, penetrances = list(male = c(0, 1), female = c(0, 0, 1)),
#'     proband = 7, Xchrom = TRUE, plot = TRUE)
#'
#' @export
FLB = function(x, carriers = NULL, homozygous = NULL, noncarriers = NULL, freq = NULL,
               affected = NULL, unknown = NULL, proband = NULL,
               penetrances = NULL, liability = NULL, loopBreakers = NULL, Xchrom = FALSE,
               details = FALSE, plot = FALSE, ...) {

  inputs = checkInput(x, affected = affected, unknown = unknown, proband = proband,
                      carriers = carriers, homozygous = homozygous, noncarriers = noncarriers,
                      freq = freq, Xchrom = Xchrom, requireProband = TRUE, requireFreq = TRUE)
  proband = inputs$proband
  affected = inputs$affected
  unknown = inputs$unknown
  carriers = inputs$carriers
  homozygous = inputs$homozygous
  noncarriers = inputs$noncarriers

  if(!proband %in% affected)
    stop2("The proband must be affected")
  if(!proband %in% c(carriers, homozygous))
    stop2("The proband must be a carrier")
  if(!is.null(x$LOOP_BREAKERS))
    stop2("Pedigrees with pre-broken loops are not allowed")

  if(plot)
    plotSegregation(x, affected, unknown, proband, carriers, homozygous, noncarriers, ...)

  # Ensure liability is a named integer vector of length pedsize
  if(is.null(liability))
    liability = rep_len(1, pedsize(x))
  else if(length(liability) == 1)
    liability = rep_len(liability, pedsize(x))
  else if(length(liability) != pedsize(x))
    stop2(sprintf("Length of `liability` vector must be 1 or %d: ", pedsize(x)), liability)

  if(is.null(names(liability)))
    names(liability) = x$ID
  else if(!setequal(names(liability), x$ID))
    stop2("Unknown ID in names of liability vector: ", setdiff(names(liability), x$ID))

  # Sex of carriers/noncarriers (needed for Xchrom)
  caSex = getSex(x, carriers %||% character(0))    # safeguard against NULL
  ncSex = getSex(x, noncarriers %||% character(0)) # safeguard against NULL

  # Add marker objects: disease locus, main marker and proband genotype
  locAttr = list(afreq = c(a = 1 - freq, b = freq),
                 chrom = if(Xchrom) "X" else NA)
  x = x |>
    addMarker(name = "m", locusAttr = locAttr) |>
    addMarker(name = "mProband", locusAttr = locAttr) |>
    addMarker(name = "dis", locusAttr = locAttr) |>
    setGenotype("m", carriers,    geno = ifelse(Xchrom & caSex == 1, "b", "a/b")) |>
    setGenotype("m", noncarriers, geno = ifelse(Xchrom & ncSex == 1, "a", "a/a")) |>
    setGenotype("m", homozygous, geno = "b/b")

  x = x |>
    setGenotype("mProband", proband, geno = genotype(x, "m", proband))

  # Break loops if necessary
  if(hasUnbrokenLoops(x)) {
    x = breakLoops(x, loopBreakers = loopBreakers, verbose = FALSE)

    # Add duplicates to input vectors, were needed
    lb = x$ID[x$LOOP_BREAKERS[, 'orig']]

    if(length(lbAff <- intersect(affected, lb)))
      affected = unique.default(c(affected, paste0("=", lbAff)))

    if(length(lbUn <- intersect(unknown, lb)))
      unknown = unique.default(c(unknown, paste0("=", lbUn)))

    # Extend liability vector
    # Note: The liab-classes for copies do not actually matter,
    # since their `prob` entry is reset to 1's in peelingProcess()
    liability = c(liability, `names<-`(liability[lb], paste0("=", lb)))

    # TODO: Handle proband vs loop breaking better
    if(proband %in% lb)
      stop2("The proband cannot be a loop breaker; please select different loop breaker(s)")
  }

  # Affection status vector, sorted along labs
  affvec = logical(pedsize(x))
  affvec[internalID(x, affected)] = TRUE
  affvec[internalID(x, unknown)] = NA

  # Penetrances
  penetMat = penet2matrix(penetrances, Xchrom = Xchrom)

  # Check that liability classes are valid row numbers in penetMat
  if(Xchrom) {
    mls = males(x) # internal = FALSE; using names now!
    fem = setdiff(x$ID, mls) # safer than females(x)

    liabM = liability[mls]
    liabF = liability[fem]
    if(!all(liabM %in% 1:nrow(penetMat$male)))
      stop2("Illegal liability class (males): ", setdiff(liabM, 1:nrow(penetMat$male)))
    if(!all(liabF %in% 1:nrow(penetMat$female)))
      stop2("Illegal liability class (females): ", setdiff(liabF, 1:nrow(penetMat$female)))
  }
  else {
    if(!all(liability %in% 1:nrow(penetMat)))
      stop2("Illegal liability class: ", setdiff(liability, 1:nrow(penetMat)))
  }

  # Setup for likelihood under causal hypothesis
  peelOrder = peelingOrder(x)
  peelingProcess = pedprobr:::peelingProcess

  if(Xchrom) {
    peeler = function(x, locus) function(dat, sub) pedprobr:::.peel_M_X(dat, sub, SEX = x$SEX)
    startCausal = function(x, locus) startdata_causative_X(x, locus, aff = affvec, penetMat = penetMat, liability = liability)
  }
  else {
    peeler = function(x, locus) function(dat, sub) pedprobr:::.peel_M_AUT(dat, sub)
    startCausal = function(x, locus) startdata_causative(x, locus, aff = affvec, penetMat = penetMat, liability = liability)
  }

  likelihoodCausal = function(x, locus) {
    peelingProcess(x, locus, startdata = startCausal,
                   peeler = peeler(x, locus), peelOrder = peelOrder)
  }

  dis = getMarkers(x, "dis")[[1]]
  m = getMarkers(x, "m")[[1]]
  mProband = getMarkers(x, "mProband")[[1]]

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
    return(list(c(FLB = FLB, BF1 = BF1, `1/BF2` = BF2),
                c(numer1 = numer1, denom1 = denom1, numer2 = numer2, denom2 = denom2),
                c(likM = likM, likMproband = likMproband, likDis = likDis)))
  FLB
}




checkInput = function(x, proband, affected, unknown, carriers, noncarriers,
                      homozygous, freq = NULL, Xchrom = FALSE,
                      requireProband = FALSE, requireFreq = FALSE) {

  if(!is.ped(x))
    stop2("The first argument must be a connected pedigree")

  # Conversion to character unless NULL
  asChar = function(v) if(!is.null(v)) as.character(v) else NULL

  proband = asChar(proband)
  affected = asChar(affected)
  unknown = asChar(unknown)
  carriers = asChar(carriers)
  noncarriers = asChar(noncarriers)
  homozygous = asChar(homozygous)

  # Proband checks
  if(requireProband && (length(proband) == 0 || proband == ""))
    stop2("A proband must be specified")
  if(length(proband) > 1)
    stop2("At most one proband is permitted")
  if(length(proband) == 1 && !proband %in% x$ID)
    stop2("Unknown proband: ", proband)

  # Other input checks
  allids = c(affected, unknown, carriers, noncarriers, homozygous)
  if(any(!allids %in% x$ID))
    stop2("Unknown ID label: ", setdiff(allids, x$ID))

  if(length(err1 <- intersect(affected, unknown)))
    stop2("Individual specified as both affected and unknown: ", err1)

  if(length(err2 <- intersect(carriers, noncarriers)))
    stop2("Individual specified as both a carrier and a non-carrier: ", err2)

  if(length(err3 <- intersect(carriers, homozygous)))
    stop2("Individual specified as both a (heterozygous) carrier and homozygous: ", err3)

  if(length(err4 <- intersect(noncarriers, homozygous)))
    stop2("Individual specified as both homozygous and a non-carrier: ", err4)

  if(Xchrom && length(err5 <- intersect(males(x), homozygous)))
    stop2("Male individual specified as a homozygous carrier in an X-linked inheritance model: ", err5)

  # Frequency checks
  if(requireFreq && is.null(freq))
      stop2("An allele frequency must be specified")
  if(!is.null(freq) && !isProb(freq, len = 1))
      stop2("The allele frequency must be a single number strictly between 0 and 1")

  list(proband = proband, affected = affected, unknown = unknown, carriers = carriers,
       noncarriers = noncarriers, homozygous = homozygous)
}
