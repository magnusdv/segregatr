### Internal functions for preparing and checking penetrance values

# Convert penetrance vectors (or keyword) to liability-style matrices.
penet2matrix = function(penetrances, Xchrom = FALSE) {

  if(is.null(p <- penetrances))
    stop2("No penetrance values given")

  # If keyword, convert to vector form
  if(length(p) == 1) {
    p = switch(p,
               AD = c(0,1,1),
               AR = c(0,0,1),
               XD = list(male = c(0,1), female = c(0,1,1)),
               XR = list(male = c(0,1), female = c(0,0,1)))
  }

  # X-style list?
  Xpen = is.list(p) && setequal(names(p), c("male", "female"))

  # Check for discordance
  if(Xchrom && !Xpen)
    stop2("For X-linked models, `penetrances` must be a list with elements `male` and `female`")
  if(!Xchrom && Xpen)
    stop2("`penetrances` is a list with elements `male` and `female`. Did you forget to specify `Xchrom = TRUE`?")

  # If X model handle males/females separately
  if(Xpen)
    list(male =   .fixPenetrances(p$male,   maleX = TRUE),
         female = .fixPenetrances(p$female, maleX = FALSE))
  else
    .fixPenetrances(p)
}


.fixPenetrances = function(p, maleX = FALSE) {
  orig = p
  if(is.vector(p) && !is.list(p)) {
    if(length(p) != 3 - maleX)
      stop2("`penetrances` vector must have length 3 (or 2 for males on X): ", orig)
    dim(p) = c(1, length(p))
  }
  else if(is.data.frame(p))
    p = as.matrix(p)
  else if(!is.matrix(p))
    stop2("`penetrances` must be either a vector or matrix-like")

  if(ncol(p) != 3 - maleX)
    stop2("Illegal number of columns in `penetrances`: ", ncol(p))

  mode(p) = "numeric"
  bad = is.na(p) | p < 0 | p > 1
  if(any(bad))
    stop2("Illegal penetrance value: ", orig[bad])

  # If colnames, sort. Otherwise set colnames
  nms = colnames(p)
  if(!is.null(nms) && all(nms %in% c("f0", "f1", "f2")))
    p = if(maleX) p[, c("f0", "f1")] else p[, c("f0", "f1", "f2")]
  else
    colnames(p) = if(maleX) c("f0", "f1") else c("f0", "f1", "f2")

  rownames(p) = 1:nrow(p)

  p
}
