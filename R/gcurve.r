#' @title gcurve
#' @param expr expression to be plotted
#' @param from lower bound of x-values
#' @param to upper bound of x-values
#' @param n number of points to plot
#' @param add logical; if TRUE add to an already existing plot; if NA start a new plot taking the defaults for the limits and log-scaling of the x-axis from the previous plot. Taken as FALSE (with a warning if a different value is supplied) if no graphics device is open.
#' @param type plot type: see plot.default.
#' @param xname character string giving the name to be used for the x axis
#' @param xlab labels and graphical parameters can also be specified as arguments.
#' @param ylab For the "function" method of plot, ... can include any of the other arguments of curve, except expr
#' @param log For the "function" method of plot, ... can include any of the other arguments of curve, except expr
#' @param xlim NULL or a numeric vector of length 2; if non-NULL it provides the defaults for c(from, to) and, unless add = TRUE, selects the x-limits of the plot – see plot.window.
#' @param ... additional items passed to curve
#' @param category optional text column appends to data.frame returned
#' @return A data.frame is returned with x and y-values and an optional column called category
#' @examples
#' my.gcurve <- gcurve(expr = dnorm(x, mean=0, sd=1),from=-4, to = 4, n= 1001,
#' category= "Standard Normal")
#' head(my.gcurve)
#' @description Returns a data.frame associated with a call to base::curve.
gcurve <- function (expr, from = NULL, to = NULL, n = 101, add = FALSE,
                    type = "l", xname = "x", xlab = xname, ylab = NULL,
                    log = NULL, xlim = NULL, category = NULL, ...)
{
  sexpr <- substitute(expr)
  if (is.name(sexpr)) {
    expr <- call(as.character(sexpr), as.name(xname))
  }
  else {
    if (!((is.call(sexpr) || is.expression(sexpr)) && xname %in%
          all.vars(sexpr)))
      stop(gettextf("'expr' must be a function, or a call or an expression containing '%s'",
                    xname), domain = NA)
    expr <- sexpr
  }
  if (dev.cur() == 1L && !identical(add, FALSE)) {
    warning("'add' will be ignored as there is no existing plot")
    add <- FALSE
  }
  addF <- identical(add, FALSE)
  if (is.null(ylab))
    ylab <- deparse(expr)
  if (is.null(from) || is.null(to)) {
    xl <- if (!is.null(xlim))
      xlim
    else if (!addF) {
      pu <- par("usr")[1L:2L]
      if (par("xaxs") == "r")
        pu <- extendrange(pu, f = -1 / 27)
      if (par("xlog"))
        10  ^  pu
      else pu
    }
    else c(0, 1)
    if (is.null(from))
      from <- xl[1L]
    if (is.null(to))
      to <- xl[2L]
  }
  lg <- if (length(log))
    log
  else if (!addF && par("xlog"))
    "x"
  else ""
  if (length(lg) == 0)
    lg <- ""
  if (grepl("x", lg, fixed = TRUE)) {
    if (from <= 0 || to <= 0)
      stop("'from' and 'to' must be > 0 with log=\"x\"")
    x <- exp(seq.int(log(from), log(to), length.out = n))
  }
  else x <- seq.int(from, to, length.out = n)
  ll <- list(x = x)
  names(ll) <- xname
  y <- eval(expr, envir = ll, enclos = parent.frame())
  for.return <- data.frame(x = x, y = y)
  if (is.null(category) == FALSE)
    for.return$category <- factor(category)
  return(for.return)
}


# gcurve(expr = x^2, from=-3, to = 3, category="My parabola")
