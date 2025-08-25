#' @importFrom stats aggregate
rareDistr = function(x, typable = NULL, proband = NULL, commonAnc = FALSE) {

  if (is.null(proband)) proband = x$ID
  typable = as.character(typable)
  vFounders = pedtools::founders(x)

  listDistr = lapply(vFounders, function(f) {
    founderUnrel = unrelated(x, f)
    fixed = c(founderUnrel, f)

    distr = data.frame(p_H0 = 1/length(vFounders))
    for(col in fixed) distr[[col]] = 0L
    distr[[f]] = 1L

    if (length(fixed) < pedsize(x)) {
      gens = generations(x, "indiv")[setdiff(x$ID, fixed)]

      for (g in min(gens):max(gens)) {
        currentCols = names(distr)
        cases = names(gens[gens == g])
        cases0 = lapply(cases, function(w) rowSums(distr[, pedtools::parents(x, w), drop = FALSE]) == 0)
        cases01 = rep(list(c(0L,1L)), length(cases))
        names(cases01) = cases

        # Build all new rows as a list
        newDistrList = vector("list", nrow(distr))
        for (i in seq_len(nrow(distr))) {
          baseRow = distr[i, , drop = FALSE]
          casesGeno = cases01
          for (w in seq_along(cases)) {
            if (cases0[[w]][i]) casesGeno[[w]] = 0L
          }
          expanded = expand.grid(casesGeno, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
          expanded = cbind(baseRow[rep(1, nrow(expanded)), currentCols, drop = FALSE], expanded)
          expanded$p_H0 = baseRow$p_H0 / nrow(expanded)
          newDistrList[[i]] = expanded
        }
        distr = do.call(rbind, newDistrList)
      }
    }

    distr = distr[, c(x$ID, "p_H0")]
    distr
  })

  distr = do.call(rbind, listDistr)

  if (commonAnc & length(proband) > 1)
    cols = commonAncestors(x, proband, inclusive = TRUE)
  else
    cols = as.character(proband)

  selection = apply(distr[, x$ID, drop = FALSE], 1, function(row) any(row[cols] == 1))
  distr = distr[selection, ]
  distr$p_H0 = distr$p_H0 / sum(distr$p_H0)

  aggregate(p_H0 ~ ., data = distr[, c(typable, "p_H0")], FUN = sum)
}
