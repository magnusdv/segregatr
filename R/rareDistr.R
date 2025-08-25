#' @import data.table

rareDistr = function (x, typable = NULL, proband = NULL, commonAnc = FALSE) {

  # Full marginal distribution if proband == NULL
  if(is.null(proband))
    proband = x$ID

  typable = as.character(typable)
  vFounders = pedtools::founders(x)

  listDistr = lapply(vFounders, function (f) {
    founderUnrel = unrelated(x, f) # obligate noncarriers
    fixed = c(founderUnrel, f)

    distr = data.table(p_H0 = 1/length(vFounders))
    distr = cbind(distr, t(rep(0L, length(fixed))))
    setnames(distr, c("p_H0", fixed))
    distr[, (f) := 1L]

    if (length(fixed) < pedsize(x)) {
      gens = generations(x, "indiv")[setdiff(x$ID, fixed)]

      for (g in min(gens):max(gens)) {
        currentCols = copy(names(distr))
        cases = names(gens[gens == g])
        cases0 = lapply(cases, function(w) rowSums(distr[, pedtools::parents(x, w), with = FALSE]) == 0)
        cases01 = setNames(rep(list(c(0L, 1L)), length(cases)), cases)
        distr[,
              newCol := list(lapply(
                seq_len(.N), function(i) {
                  casesGeno = cases01
                  for (w in seq_along(cases)) {
                    if (cases0[[w]][i])
                      casesGeno[[w]] = 0L
                  }
                  casesGeno = expand.grid(casesGeno)
                  lapply(casesGeno, unlist)
                }
              ))
        ]
        distr = distr[, p_H0 := p_H0 / sapply(lapply(distr[["newCol"]], "[[", 1), length)]
        distr = distr[, rbindlist(newCol), by = currentCols]
      }
    }
    setcolorder(distr, c(x$ID, "p_H0"))
  })

  distr = rbindlist(listDistr)

  if (commonAnc & length(proband) > 1)
    cols = commonAncestors(x, proband, inclusive = TRUE)
  else
    cols = as.character(proband)
  selection = apply(distr[, -ncol(distr), with = FALSE], 1, function(row) any(row[cols] == 1))
  distr = distr[selection]
  distr[, p_H0 := p_H0 / sum(p_H0)]
  distr[, .(p_H0 = sum(p_H0)), by = typable]

  distr
}
