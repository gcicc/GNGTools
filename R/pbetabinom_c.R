#' Replacement for pbetanimom
#'
#' @param q vector of quarantines
#' @param size vector of totals
#' @param m vector of probabilities of success
#' @param s vector of over dispersion parameters
#'
#' @return This is a modification of the pbetanimom function from the rmutil package.
#'  The original function would throw an error for non-nonsensical values rather than returning 0 or 1.
#' @export

pbetabinom_c <- function (q, size, m, s) {
  # if (any(q < 0))
  #     stop("q must contain non-negative values")
  if (any(size < 0))
    stop("size must contain non-negative values")
  if (any(m <= 0) || any(m >= 1))
    stop("m must lie between 0 and 1")
  if (any(s <= 0))
    stop("s must be positive")
  len <- max(length(q), length(m), length(s), length(size))
  if (length(q) != len) {
    if (length(q) == 1)
      q <- rep(q, len)
    else stop("length of q incorrect")
  }
  if (length(size) != len) {
    if (length(size) == 1)
      size <- rep(size, len)
    else stop("size must be the same length as q")
  }
  if (length(m) != len) {
    if (length(m) == 1)
      m <- rep(m, len)
    else stop("m and q must have the same length")
  }
  if (length(s) != len) {
    if (length(s) == 1)
      s <- rep(s, len)
    else stop("s and q must have the same length")
  }
  if (any(q > size)){
    # Updating to correctly deal with this siduation
    # stop("q must be <= size")
    if(all(q>size)) return(rep(1,len))
    Val <- q<= size
    out <- rep(1,len)
    out[Val] <- pbetabinom_c(q[Val],size[Val],m[Val],s[Val])
    return(out)
  }
  if (any(q < 0)){
    # stop("q must contain non-negative values")
    if(all(q < 0)) return(rep(0,len))
    Val <- q>=0
    out <- rep(0,len)
    out[Val] <- pbetabinom_c(q[Val],size[Val],m[Val],s[Val])
    return(out)
  }
  t <- s * m
  u <- s * (1 - m)
  res <- vector("numeric", length(q))
  for (i in 1:length(q)) {
    qq <- 0:q[i]
    res[i] <- sum(exp(lbeta(qq + t[i], size[i] - qq + u[i]) -
                        lbeta(t[i], u[i]) + lchoose(size[i], qq)))
  }
  res
}
