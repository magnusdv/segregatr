
startdata_causative_X = function(x, marker, aff, penetMat, liability = NULL) {

  # Founder data
  FOU = founders(x, internal = TRUE)
  FOU_INB = rep(NA_real_, pedsize(x))
  FOU_INB[FOU] = founderInbreeding(x) # enable quick lookup

  # Allele frequencies (used in HW below)
  afr = afreq(marker)

  # Build genotype list in internal format
  glist = pedprobr:::.buildGenolistX(x, marker, eliminate = 2)
  impossible = attr(glist, "impossible")

  # Loop through individuals
  dat = lapply(1:pedsize(x), function(i) {

    # If impossible, finish loop quickly (cannot use break in `apply()`)
    if(impossible)
      return(NULL)

    g = glist[[i]]

    # Liability class
    liab = liability[i]

    if (x$SEX[i] == 1) {
      # Penetrance values (males)
      penet = as.numeric(penetMat$male[liab, ])
      # Affection status priors
      nmut = g$mat - 1
      affi = aff[i]
      if(is.na(affi))
        affpriors = rep_len(1, length(nmut))
      else if(affi)
        affpriors = penet[nmut + 1]
      else
        affpriors = 1 - penet[nmut + 1]
      # Add genotype priors
      prob = as.numeric(affpriors)
      if (i %in% FOU)
        prob = prob * afr[g$mat]
    }
    else {
      # Penetrance values (females)
      penet = as.numeric(penetMat$female[liab, ])
      # Affection status priors
      nmut = g$pat + g$mat - 2
      affi = aff[i]
      if(is.na(affi))
        affpriors = rep_len(1, length(nmut))
      else if(affi)
        affpriors = penet[nmut + 1]
      else
        affpriors = 1 - penet[nmut + 1]
      # Add genotype priors
      prob = as.numeric(affpriors)
      if (i %in% FOU)
        prob = prob * pedprobr::HWprob(g$pat, g$mat, afr, f = FOU_INB[i])
    }

    # Remove impossible entries
    keep = prob > 0
    if(!any(keep)) {
      impossible = TRUE
      return(NULL)
    }

    if(!all(keep)) {
      g$pat = g$pat[keep]
      g$mat = g$mat[keep]
      prob = prob[keep]
    }

    g$prob = prob
    g
  })

  attr(dat, "impossible") = impossible
  dat
}
