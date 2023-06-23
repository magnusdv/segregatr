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

isProb = function(x, len = NULL)
  is.numeric(x) && (is.null(len) || length(x) == len) &&
    all(!is.na(x) & x >= 0 & x <= 1)

`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

.mysetdiff = function(x, y) {
  unique.default(x[match(x, y, 0L) == 0L])
}

# Fast intersection. NB: assumes no duplicates!
.myintersect = function (x, y) {
  y[match(x, y, 0L)]
}


