#' Pedigree plot for segregation analysis
#'
#' Plots a pedigree showing the segregation of a variant.
#'
#' @inheritParams FLB
#' @param cex,margins Arguments passed on to [pedtools::plot.ped()].
#'
#' @examples
#'
#' x = nuclearPed(2)
#' plotSegregation(x, affected = 3:4, unknown = 1:2, proband = 3, carriers = 3:4)
#'
#' y = cousinPed(1, child = TRUE)
#' plotSegregation(y, affected = 9, unknown = 1:2, proband = 9,
#'                 homozygous = 9, deceased = 1:2)
#'
#' @importFrom graphics arrows strheight text
#' @importFrom utils packageVersion
#' @export
plotSegregation = function(x, affected = NULL, unknown = NULL, proband = NULL,
                           carriers = NULL, homozygous = NULL, noncarriers = NULL, cex = 1,
                           margins = 1, pos.geno = "bottom", pos.arrow = "bottomleft", ...) {

  # Input checks
  allids = c(proband, affected, unknown, carriers, noncarriers, homozygous)
  if(any(!allids %in% labels(x)))
    stop2("Unknown ID label: ", setdiff(allids, labels(x)))

  if(length(err1 <- intersect(affected, unknown)))
    stop2("Individual specified as both affected and unknown: ", err1)

  if(length(err2 <- intersect(carriers, noncarriers)))
    stop2("Individual specified as both a carrier and a non-carrier: ", err2)

  if(length(err3 <- intersect(carriers, homozygous)))
    stop2("Individual specified as both a (heterozygous) carrier and homozygous: ", err3)

  if(length(err4 <- intersect(noncarriers, homozygous)))
    stop2("Individual specified as both a non-carrier and homozygous: ", err4)

  hasProband = length(proband) > 0

  if(hasProband && length(proband) > 1)
    stop2("At most one proband is permitted")

  # Invoke automatic margin adjustment unless fully specified
  autoMargins = length(margins) < 4

  # Adjust margins if arrow on left side
  if(hasProband && autoMargins && packageVersion("pedtools") > 2.0) {

    align = .pedAlignment(x, ...)
    idx = match(internalID(x, proband), align$plotord)
    if(align$xall[idx] == 0) {
      margins = rep_len(margins, 4)
      margins[2] = margins[2] + 2.5
    }
  }

  p = plot(x,
           aff = affected,
           textInside = ifelse(1:pedsize(x) %in% unknown, "?", ""),
           cex = cex,
           margins = margins,
           keep.par = TRUE,
           ...)

  # Reorder positions
  p$x = p$x[internalID(x, ids = 1:pedsize(x))]
  p$y = p$y[internalID(x, ids = 1:pedsize(x))]

  # Genotype label position
  lpos = switch(pos.geno,
                "topleft" = c(-3, 2, 0.75),
                "topright" = c(-3, 4, 0.75),
                c(3.25, 3, 0)
  )
  vdist = p$boxh + lpos[1] * abs(strheight("M", cex = cex))  # vertical dist from top of symbol to "+"

  if(!is.null(carriers))
    text(p$x[carriers], p$y[carriers] + vdist, labels = "+", cex = cex*1.5, font = 1, pos = lpos[2], offset = lpos[3])

  if(!is.null(homozygous))
    text(p$x[homozygous], p$y[homozygous] + vdist, labels = "++", cex = cex*1.5, font = 1, pos = lpos[2], offset = lpos[3])

  if(!is.null(noncarriers))
    text(p$x[noncarriers], p$y[noncarriers] + vdist, labels = "-", cex = cex*1.5, font = 1, pos = lpos[2], offset = lpos[3])

  # proband arrow
  if(hasProband) {
    x.mod = ifelse(pos.arrow %in% c("topright", "bottomright"), +1, -1)
    y.mod = ifelse(pos.arrow %in% c("topleft", "topright"), 0, 1)
    corner.x = p$x[proband] + .5*x.mod*p$boxw
    corner.y = p$y[proband] + y.mod*p$boxh
    arrows(corner.x + 1.7*x.mod*p$boxw, corner.y + 0.9*y.mod*p$boxh - 0.9*(1-y.mod)*p$boxh,
           corner.x + .5*x.mod*p$boxw, corner.y,
           lwd = 2, length = .15, xpd = NA)
  }

  # Return plot parameters
  invisible(p)
}
