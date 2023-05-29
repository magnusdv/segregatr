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
#' x = cousinPed(1)
#' plotSegregation(x, affected = c(3,6,8), unknown = 1, carrier = 2:5,
#'                 homozygous = 8, noncarriers = 6:7, proband = 8)
#'
#' # Different placement of genotypes and proband arrow
#' plotSegregation(x, affected = c(3,6,8), unknown = 1, carrier = 2:5,
#'                 homozygous = 8, noncarriers = 6:7, proband = 8,
#'                 pos.geno = "topleft", pos.arrow = "bottomright")
#'
#' @importFrom graphics arrows strheight text
#' @importFrom utils packageVersion
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
                bottom = c(3.25, 3, 0),
                topleft = c(-3, 2, 0.75),
                topright = c(-3, 4, 0.75),
                stop2('`pos.geno` must be either "bottom" , "topleft" or "topright": ', pos.geno))

  vdist = p$boxh + lpos[1] * abs(strheight("M", cex = cex))  # vertical dist from top of symbol to "+"

  if(!is.null(carriers))
    text(p$x[carriers], p$y[carriers] + vdist, labels = "+", cex = cex*1.5, font = 1, pos = lpos[2], offset = lpos[3])

  if(!is.null(homozygous))
    text(p$x[homozygous], p$y[homozygous] + vdist, labels = "++", cex = cex*1.5, font = 1, pos = lpos[2], offset = lpos[3])

  if(!is.null(noncarriers))
    text(p$x[noncarriers], p$y[noncarriers] + vdist, labels = "-", cex = cex*1.5, font = 1, pos = lpos[2], offset = lpos[3])

  # proband arrow
  if(hasProband) {
    mod = switch(pos.arrow,
           bottomleft = list(x = -1, y = 1),
           bottomright = list(x = 1, y = 1),
           topleft = list(x = -1, y = 0),
           topright = list(x = 1, y = 0),
           stop2('`pos.arrow` must be either "bottomleft", "bottomright", "topleft" or "topright": ', pos.arrow))

    corner.x = p$x[proband] + .5*mod$x*p$boxw
    corner.y = p$y[proband] + mod$y*p$boxh
    arrows(x0 = corner.x + 1.7*mod$x*p$boxw,
           y0 = corner.y + 0.9*mod$y*p$boxh - 0.9*(1-mod$y)*p$boxh,
           x1 = corner.x + .5*mod$x*p$boxw,
           y1 = corner.y,
           lwd = 2, length = .15, xpd = NA)
  }

  # Return plot parameters
  invisible(p)
}
