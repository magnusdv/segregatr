stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

# Test that input is a positive (or similar) integer.
isCount = function(x, minimum = 1) {
  isTRUE(length(x) == 1 &&
           (is.integer(x) || (is.numeric(x) && x == as.integer(x))) &&
           x >= minimum)
}

.mysetdiff = function(x, y) {
  unique.default(x[match(x, y, 0L) == 0L])
}

# Fast intersection. NB: assumes no duplicates!
.myintersect = function (x, y) {
  y[match(x, y, 0L)]
}


# Same as in paramlink2
fixPenetrances = function(p, maleX = FALSE) {
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
