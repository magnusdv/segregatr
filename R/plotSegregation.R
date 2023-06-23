#' Pedigree plot for segregation analysis
#'
#' Plots a pedigree showing the segregation of a variant.
#'
#' @inheritParams FLB
#' @param cex,margins Arguments passed on to [pedtools::plot.ped()].
#' @param pos.geno Position of genotype labels relative to pedigree symbols;
#'   either "bottom" (default), "topleft" or "topright".
#' @param pos.arrow Position of the proband arrow; either "bottomleft",
#'   "bottomright", "topleft" or "topright".
#'
#' @examples
#'
#' x = nuclearPed(2)
#' plotSegregation(x, proband = 3, carriers = 3:4, noncarriers = 1,
#'                 aff = 3:4, unknown = 1:2)
#'
#' # Same with various options
#' plotSegregation(x, proband = 3, carriers = 3:4, noncarriers = 1,
#'                 aff = 3:4, unknown = 1:2,
#'                 pos.geno = "topright", pos.arrow = "topleft",
#'                 labs = NULL, title = "Family 1", cex.main = 1.5)
#'
#' # Recessive example
#' y = cousinPed(1, child = TRUE)
#' plotSegregation(y, affected = 9, unknown = 1:6, carrier = 7:8,
#'                 homozygous = 9, noncarriers = c(4,6), proband = 9)
#'
#' # Different symbol placements
#' plotSegregation(y, affected = 9, unknown = 1:6, carrier = 7:8,
#'                 homozygous = 9, noncarriers = c(4,6), proband = 9,
#'                 pos.geno = "topleft", pos.arrow = "bottomright")
#'
#' # Incest case
#' y = nuclearPed() |> addChildren(father = 3, mother = 2, nch = 3)
#'
#' plotSegregation(y, proband = 4, aff = 4:6, unknown = 2, carrier = 4:6, deceased = 1,
#'                 pos.geno = "topleft", pos.arrow = "bottomright")
#'
#' @importFrom graphics arrows strheight text
#' @export
plotSegregation = function(x, affected = NULL, unknown = NULL, proband = NULL,
                           carriers = NULL, homozygous = NULL, noncarriers = NULL, cex = 1,
                           margins = 1, pos.geno = "bottom", pos.arrow = "bottomleft", ...) {

  inputs = checkInput(x, affected = affected, unknown = unknown, proband = proband,
                      carriers = carriers, homozygous = homozygous, noncarriers = noncarriers)
  proband = inputs$proband
  affected = inputs$affected
  unknown = inputs$unknown
  carriers = inputs$carriers
  homozygous = inputs$homozygous
  noncarriers = inputs$noncarriers

  hasProband = length(proband) > 0

  pos.geno = match.arg(pos.geno, c("bottom", "topleft", "topright"))
  pos.arrow = match.arg(pos.arrow, c("bottomleft", "bottomright", "topleft", "topright"))

  # Question mark inside those with unknown affection status
  qm = rep_len("?", length(unknown))
  names(qm) = unknown

  # Main alignment
  align = .pedAlignment(x, ...)
  xpos = align$x
  ypos = align$y
  names(xpos) = names(ypos) = x$ID

  # Invoke automatic margin adjustment if not fully specified
  if(autoMargins <- length(margins) < 4) {
    margins = rep_len(margins, 4)
    xr = align$xrange
    yr = align$yrange

    # Extra margins for genotype labels
    typed = c(carriers, noncarriers, homozygous)

    homozL = pos.geno == "topleft" && any(xpos[homozygous] == xr[1])
    carrL  = pos.geno == "topleft" && any(xpos[typed] == xr[1])

    homozR = pos.geno == "topright" && any(xpos[homozygous] == xr[2])
    carrR  = pos.geno == "topright" && any(xpos[typed] == xr[2])

    carrT = pos.geno %in% c("topleft", "topright") && any(ypos[typed] == yr[1])
    carrB = pos.geno == "bottom" && any(ypos[typed] == yr[2])

    # Adjustment for arrow?
    if(hasProband) {
      arrowL = xpos[proband] == xr[1] && pos.arrow %in% c("bottomleft", "topleft")
      arrowR = xpos[proband] == xr[2] && pos.arrow %in% c("bottomright", "topright")
      arrowT = ypos[proband] == yr[1] && pos.arrow %in% c("topleft", "topright")
      arrowB = ypos[proband] == yr[2] && pos.arrow %in% c("bottomleft", "bottomright")
    }
    else
      arrowL = arrowR = arrowT = arrowB = FALSE

    # Combine adjustments (geno + arrow)
    extraMargins = c(if(carrB) 0.75 else if(arrowB) 0.5 else 0,
                     if(arrowL) 2.5  else if(homozL) 1.5 else if(carrL) 0.75 else 0,
                     if(arrowT) 1.25 else if(carrT) 1 else 0,
                     if(arrowR) 2.5  else if(homozR) 1.5 else if(carrR) 0.75 else 0)

    margins = margins + extraMargins
  }

  # Plot pedigree
  p = plot(x, aff = affected, textInside = qm, cex = cex,
           margins = margins, keep.par = TRUE, ...)

  boxw = p$scaling$boxw
  boxh = p$scaling$boxh

  # Genotype label position
  lpos = switch(pos.geno,
                bottom = c(3.25, 3, 0),
                topleft = c(-3, 2, 0.75),
                topright = c(-3, 4, 0.75))

  vdist = boxh + lpos[1] * abs(strheight("M", cex = cex))  # vertical dist from top of symbol to "+"

  if(!is.null(carriers))
    text(xpos[carriers], ypos[carriers] + vdist, labels = "+", cex = cex*1.5, font = 1, pos = lpos[2], offset = lpos[3])

  if(!is.null(homozygous))
    text(xpos[homozygous], ypos[homozygous] + vdist, labels = "++", cex = cex*1.5, font = 1, pos = lpos[2], offset = lpos[3])

  if(!is.null(noncarriers))
    text(xpos[noncarriers], ypos[noncarriers] + vdist, labels = "-", cex = cex*1.5, font = 1, pos = lpos[2], offset = lpos[3])

  # Proband arrow
  if(hasProband) {
    mod = switch(pos.arrow,
           bottomleft = list(x = -1, y = 1),
           bottomright = list(x = 1, y = 1),
           topleft = list(x = -1, y = 0),
           topright = list(x = 1, y = 0))

    corner.x = xpos[proband] + .5*mod$x * boxw
    corner.y = ypos[proband] +    mod$y * boxh
    arrows(x0 = corner.x + 1.7*mod$x * boxw,
           y0 = corner.y + 0.9*mod$y * boxh - 0.9*(1-mod$y) * boxh,
           x1 = corner.x + 0.5*mod$x * boxw,
           y1 = corner.y,
           lwd = 2*cex, length = .15, xpd = NA)
  }

  # Return plot parameters
  invisible(p)
}
