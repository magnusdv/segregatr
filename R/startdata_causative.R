
startdata_causative = function(x, marker, aff, penetMat, liability = NULL) {

  nInd = pedsize(x)

  # Founder data
  FOU = founders(x, internal = TRUE)
  FOU_INB = rep(NA_real_, nInd)
  FOU_INB[FOU] = founderInbreeding(x) # enable quick lookup

  # Allele frequencies (used in HW below)
  afr = afreq(marker)

  # Build genotype list in internal format
  glist = pedprobr:::.buildGenolist(x, marker, eliminate = 2)

  if(attr(glist, "impossible"))
    return(glist)

  for(i in 1:nInd) {
    g = glist[[i]]

    # Penetrance values: Liability class = row of penetMat
    lab = x$ID[i]
    liab = liability[[lab]]
    penet = as.numeric(penetMat[liab, ])

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

    # Remove impossible entries
    keep = prob > 0
    if(!all(keep)) {
      g$pat = g$pat[keep]
      g$mat = g$mat[keep]
      prob = prob[keep]
    }
    g$prob = prob

    # Update entry
    glist[[i]] = g

    # If impossible: break
    if(!any(keep)) {
      attr(glist, "impossible") = TRUE
      break
    }
  }

  glist
}
