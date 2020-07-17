
startdata_causative = function(x, marker, aff, penetrances, liability = NULL) {

  # Penetrances
  if(is.null(dim(penetrances)))
    penetrances = as.data.frame(as.list(penetrances))

  if(ncol(penetrances) < 3)
    stop2("Wrong format for the penetrance values. Please consult the documentation.")

  penetrances = penetrances[, 1:3, drop = FALSE]

  if(anyNA(penetrances))
    stop2("Some penetrance values are missing")
  if(any(penetrances < 0 | penetrances > 1))
    stop2("All penetrance values must be between 0 and 1, inclusive")

  if(!all(liability %in% 1:nrow(penetrances)))
    stop2("Illegal liability class: ", setdiff(liability, 1:nrow(penetrances)))

  # Build genotype list in internal format
  glist = pedprobr:::.buildGenolist(x, marker, eliminate = 2)

  if (attr(glist, "impossible")) {
    dat = structure(list(), impossible = TRUE)
    return(dat)
  }

  FOU = founders(x, internal = TRUE)

  # Founder inbreeding: A vector of length pedsize(x), with NA's at nonfounders
  # Enables quick look-up e.g. FOU_INB[i].
  FOU_INB = rep(NA_real_, pedsize(x))
  FOU_INB[FOU] = founderInbreeding(x)

  # Allele frequencies (used in HW below)
  afr = afreq(marker)

  # Loop through individuals
  dat = lapply(1:pedsize(x), function(i) {
    h = glist[[i]]

    # Penetrance values
    liab = liability[i]
    penet = as.numeric(penetrances[liab, 1:3])

    # Affection status priors
    nmut = h[1, ] + h[2, ] - 2
    affi = aff[i]
    if(is.na(affi))
      affpriors = rep_len(1, length(nmut))
    else if(affi)
      affpriors = penet[nmut + 1]
    else
      affpriors = 1 - penet[nmut + 1]

    # Total genotype priors
    prob = as.numeric(affpriors)
    if (i %in% FOU)
      prob = prob * pedprobr::HWprob(h[1, ], h[2, ], afr, f = FOU_INB[i])

    # Remove impossible entries
    zer = prob == 0
    if (any(zer)) {
      h = h[, !zer, drop = F]
      prob = prob[!zer]
      if (length(prob) == 0)
        assign("impossible", TRUE, envir = parent.frame())
    }
    list(hap = h, prob = prob)
  })

  # Add impossibility attribute
  impossible = FALSE
  for(dt in dat) if(!length(dt$prob)) {
    impossible = TRUE
    break
  }
  attr(dat, "impossible") = impossible

  dat
}
