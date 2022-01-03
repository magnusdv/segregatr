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
#' plotSegregation(x, affected = 3:4, unknown = 1:2, proband = 3,
#'                 carriers = 3:4, margins = c(1,3,1,1))
#'
#' @importFrom graphics arrows strheight text
#' @export
plotSegregation = function(x, affected = NULL, unknown = NULL, proband = NULL,
                           carriers = NULL, homozygous = NULL, noncarriers = NULL, cex = 1,
                           margins = rep(1, 4), ...) {

  # Input checks
  allids = c(proband, affected, unknown, carriers, noncarriers)
  if(any(!allids %in% labels(x)))
    stop2("Unknown ID label: ", setdiff(allids, labels(x)))

  if(length(err1 <- intersect(affected, unknown)))
    stop2("Individual specified as both affected and unknown: ", err1)

  if(length(err2 <- intersect(carriers, noncarriers)))
    stop2("Individual specified as both a carrier and a non-carrier: ", err2)

  p = plot(x,
           aff = affected,
           textInside = ifelse(1:pedsize(x) %in% unknown, "?", ""),
           cex = cex,
           margins = margins,
           keep.par = TRUE,
           ...)

  vdist = p$boxh + 3.25 * abs(strheight("M", cex = cex))  # vertical dist from top of symbol to "+"

  if(!is.null(carriers))
    text(p$x[carriers], p$y[carriers] + vdist, labels = "+", cex = cex*1.5, font = 1, pos = 3, offset = 0)

  if(!is.null(homozygous))
    text(p$x[homozygous], p$y[homozygous] + vdist, labels = "++", cex = cex*1.5, font = 1, pos = 3, offset = 0)

  if(!is.null(noncarriers))
    text(p$x[noncarriers], p$y[noncarriers] + vdist, labels = "-", cex = cex*1.5, font = 1, pos = 3, offset = 0)

  # proband arrow
  if(!is.null(proband)) {
    prob = internalID(x, proband)
    corner.x = p$x[prob] - .5*p$boxw
    corner.y = p$y[prob] + p$boxh
    arrows(corner.x - 1.7*p$boxw, corner.y + 0.9*p$boxh, corner.x - .5*p$boxw, corner.y , lwd = 2, length = .15, xpd = NA)
  }
}
