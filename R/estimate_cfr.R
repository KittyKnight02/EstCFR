#' Estimate the case fatality ratio (CFR) of an outbreak.
#'
#' Takes a dataset and a method and returns an estimate of the CFR together with
#' a lower and upper confidence limit, optionally recalculated at several time
#' points.
#'
#' @section Data formats:
#' `estimate_cfr()` never modifies your data: if the data are not in the right
#' format it stops and asks you to run `prepare_dataset()`. Two dataset types are
#' supported, selected with `dataset_type`:
#' \describe{
#'   \item{`"count"`}{One row per reporting date, in order of increasing date, with
#'   *cumulative* (non-decreasing) counts. Required columns: `date` (class `Date`),
#'   `cases_total`, `deaths_total`. The `"resolved"` method also needs
#'   `recoveries_total`.}
#'   \item{`"linelist"`}{One row per case. Required columns: `start_date`
#'   (class `Date`; onset or first report), `outcome_date` (class `Date`) and
#'   `event` (character: `"death"`, `"recovery"`, or `NA` for a case with no outcome
#'   yet). `event` and `outcome_date` must be missing together.}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{`"naive"`}{deaths / cases. Clopper-Pearson (default) or Jeffreys CI. Count
#'   or linelist.}
#'   \item{`"resolved"`}{deaths / (deaths + recoveries). Clopper-Pearson (default) or
#'   Jeffreys CI. Count or linelist.}
#'   \item{`"delay_adjusted"`}{Nishiura et al. (2009) delay-adjusted CFR; cases are
#'   weighted by the probability that their outcome would be known by now.
#'   Profile-likelihood CI. Count data need `delay_dist`; for a linelist the
#'   onset-to-death delay is fitted from the data if `delay_dist` is `NULL`.}
#'   \item{`"competing_risks"`}{Aalen-Johansen cumulative incidence of death with
#'   recovery as a competing event. Greenwood-type variance; Wald CI (logit by
#'   default) or bootstrap. Linelist only.}
#'   \item{`"adapted_kaplan_meier"`}{Kaplan-Meier style estimator of Ghani et al.
#'   (2005). Delta-method variance; Wald CI (logit by default). Linelist only.}
#'   \item{`"parametric_mixture"`}{Mixture model of Ghani et al. (2005): death
#'   probability with separate parametric time-to-death and time-to-recovery
#'   distributions, fitted by maximum likelihood. Likelihood-ratio CI. Linelist
#'   only.}
#' }
#'
#' @param data Dataset of the form outputted by the function prepare_dataset; see
#'   the "Data formats" section.
#' @param method One of "naive", "resolved", "delay_adjusted", "competing_risks",
#'   "adapted_kaplan_meier", "parametric_mixture"
#' @param dataset_type One of "count" (cumulative counts of cases and deaths) or "linelist"
#' @param delay_dist Delay distribution used for "delay_adjusted": a function
#'   returning the cumulative distribution function of the onset-to-death delay
#'   at a vector of (whole) days, e.g. `function(x) pgamma(x, shape = 2, scale = 5)`.
#'   If "estimate" and the data are a linelist, a distribution of family
#'   `parametric_delay_type` is fitted to the observed onset-to-death delays.
#' @param parametric_delay_type Family of the parametric delay distribution,
#'   one of "gamma", "weibull" or "lognormal". Used by "parametric_mixture", and
#'   by "delay_adjusted" when the delay is fitted from a linelist.
#' @param CI_method Method used to calculate confidence intervals. `NULL` (the
#'   default) uses the first option listed for the method:
#'   \itemize{
#'     \item `"naive"`, `"resolved"`: `"clopper_pearson"` (exact) or `"jeffreys"`
#'       (Bayesian interval from a Beta(0.5, 0.5) prior; narrower).
#'     \item `"delay_adjusted"`: `"profile_likelihood"`.
#'     \item `"competing_risks"`: `"logit"` (Wald on the logit scale), `"wald"`
#'       (untransformed Wald, clipped to 0-1), `"loglog"` (Wald on the
#'       log(-log(1 - CFR)) scale) or `"bootstrap"` (percentile interval from 1000
#'       resamples of cases; set a seed for reproducibility).
#'     \item `"adapted_kaplan_meier"`, as in the Stata module casefat: `"logit"`,
#'       `"wald"` (module option `untrans`), `"greenwood_logit"` and
#'       `"greenwood_wald"` (module option `greenwood`, which swaps the covariance
#'       of the survival curve for the Greenwood one).
#'     \item `"parametric_mixture"`: `"likelihood_ratio"`.
#'   }
#' @param time_points Time point at which to calculate CFR. Either input as vector of
#'   dates of class Date (see as.Date()) or "final" to mean last date in dataset.
#'   Only data available on or before each date are used.
#' @param alpha Significance level; 0.05 gives 95% confidence intervals.
#' @return A data frame with one row per time point and columns `time_point`,
#'   `method`, `cfr`, `ci_lower` and `ci_upper`. Values are `NA` (with a warning)
#'   when the method cannot be applied to the data, for example when no
#'   recoveries have been observed.
#' @importFrom stats pgamma pweibull plnorm dgamma dweibull dlnorm dbinom dpois
#' @importFrom stats qbeta qnorm qchisq plogis qlogis optim uniroot
#' @export
#' @examples
#' counts <- data.frame(
#'   date = as.Date("2014-06-01") + 0:4,
#'   cases_total = c(10, 25, 40, 60, 80),
#'   deaths_total = c(2, 8, 15, 25, 33),
#'   recoveries_total = c(1, 3, 6, 10, 15)
#' )
#' estimate_cfr(counts, method = "naive")
#' estimate_cfr(counts, method = "resolved")
estimate_cfr <- function(data,
                         method,
                         dataset_type = "count",
                         delay_dist = NULL,
                         parametric_delay_type = "gamma",
                         CI_method = NULL,
                         time_points = "final",
                         alpha = 0.05
                         ){
  methods <- c("naive", "resolved", "delay_adjusted", "competing_risks",
               "adapted_kaplan_meier", "parametric_mixture")
  if (!is.character(method) || length(method) != 1 || !method %in% methods) {
    stop("`method` must be one of: ", paste0('"', methods, '"', collapse = ", "), ".",
         call. = FALSE)
  }
  if (!is.character(dataset_type) || length(dataset_type) != 1 ||
      !dataset_type %in% c("count", "linelist")) {
    stop('`dataset_type` must be either "count" or "linelist".', call. = FALSE)
  }
  if (!is.character(parametric_delay_type) || length(parametric_delay_type) != 1 ||
      !parametric_delay_type %in% c("gamma", "weibull", "lognormal")) {
    stop('`parametric_delay_type` must be one of "gamma", "weibull" or "lognormal".',
         call. = FALSE)
  }
  binomial_cis <- c("clopper_pearson", "jeffreys")
  ci_options <- list(
    naive = binomial_cis,
    resolved = binomial_cis,
    delay_adjusted = "profile_likelihood",
    competing_risks = c("logit", "wald", "loglog", "bootstrap"),
    adapted_kaplan_meier = c("logit", "wald", "greenwood_logit", "greenwood_wald"),
    parametric_mixture = "likelihood_ratio"
  )
  allowed <- ci_options[[method]]
  if (is.null(CI_method)) {
    CI_method <- allowed[1]    # the first option is each method's default
  } else if (!is.character(CI_method) || length(CI_method) != 1 || !CI_method %in% allowed) {
    stop('For method "', method, '", `CI_method` must be NULL or one of: ',
         paste0('"', allowed, '"', collapse = ", "), ".", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  if (!is.null(delay_dist) && !is.function(delay_dist)) {
    stop("`delay_dist` must be NULL or a function returning the delay CDF.",
         call. = FALSE)
  }

  check_dataset(data, dataset_type, method, delay_dist)
  time_points <- resolve_time_points(time_points, data, dataset_type)

  rows <- lapply(time_points, function(tp) {
    d <- data_at_time(data, dataset_type, tp)
    res <- switch(method,
      naive = cfr_naive(d, dataset_type, CI_method, alpha),
      resolved = cfr_resolved(d, dataset_type, CI_method, alpha),
      delay_adjusted = cfr_delay_adjusted(d, dataset_type, CI_method, alpha,
                                          delay_dist, parametric_delay_type),
      competing_risks = cfr_competing_risks(d, CI_method, alpha),
      adapted_kaplan_meier = cfr_adapted_kaplan_meier(d, CI_method, alpha),
      parametric_mixture = cfr_parametric_mixture(d, CI_method, alpha,
                                                  parametric_delay_type)
    )
    data.frame(time_point = tp, method = method, cfr = res$CFR,
               ci_lower = res$ci_lower, ci_upper = res$ci_upper)
  })
  do.call(rbind, rows)
}


# ---------------------------------------------------------------------------
# Input checking
# ---------------------------------------------------------------------------

# Stop with a message that sends the user to prepare_dataset().
bad_format <- function(...) {
  stop(..., "\nUse prepare_dataset() to put your data in the format estimate_cfr() needs.",
       call. = FALSE)
}

# Check that the dataset is in the format the chosen method needs. Changes nothing:
# either returns `data` untouched (invisibly) or stops with an error.
#
# Count data: columns `date` (Date, strictly increasing), cumulative `cases_total`
# and `deaths_total` (plus `recoveries_total` for "resolved"), all non-decreasing.
# Linelist: columns `start_date` and `outcome_date` (Date) and `event`
# ("death", "recovery" or NA for no outcome yet).
check_dataset <- function(data, dataset_type, method, delay_dist = NULL) {
  if (!is.data.frame(data)) {
    bad_format("`data` must be a data frame.")
  }
  if (nrow(data) == 0) {
    bad_format("`data` has no rows.")
  }

  required <- if (dataset_type == "count") {
    c("date", "cases_total", "deaths_total", if (method == "resolved") "recoveries_total")
  } else {
    c("start_date", "outcome_date", "event")
  }
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    bad_format("`data` is missing required column(s) for dataset_type = \"", dataset_type,
               "\": ", paste0("`", missing, "`", collapse = ", "), ".")
  }

  if (dataset_type == "count") {
    # (not format problems: these methods cannot use count data at all)
    if (method %in% c("competing_risks", "adapted_kaplan_meier", "parametric_mixture")) {
      stop('Method "', method, '" needs individual-level data: use a linelist ',
           'with dataset_type = "linelist".', call. = FALSE)
    }
    if (method == "delay_adjusted" && is.null(delay_dist)) {
      stop('Method "delay_adjusted" with count data needs `delay_dist`.', call. = FALSE)
    }

    if (!inherits(data$date, "Date")) bad_format("Column `date` must be of class Date.")
    if (anyNA(data$date)) bad_format("Column `date` contains missing values.")
    if (anyDuplicated(data$date)) bad_format("Column `date` contains duplicated dates.")
    if (is.unsorted(data$date)) bad_format("Rows must be in order of increasing `date`.")

    # cumulative counts: every count column present must be valid and non-decreasing
    count_cols <- intersect(c("cases_total", "deaths_total", "recoveries_total"), names(data))
    for (col in count_cols) {
      x <- data[[col]]
      if (!is.numeric(x)) bad_format("Column `", col, "` must be numeric.")
      if (anyNA(x)) bad_format("Column `", col, "` contains missing values.")
      if (any(x < 0)) bad_format("Column `", col, "` contains negative values.")
      if (any(diff(x) < 0)) {
        bad_format("Column `", col, "` must be cumulative (non-decreasing over time).")
      }
    }
    if (any(data$deaths_total > data$cases_total)) {
      bad_format("`deaths_total` exceeds `cases_total` on at least one date.")
    }
    if (method == "resolved" &&
        any(data$deaths_total + data$recoveries_total > data$cases_total)) {
      bad_format("`deaths_total` + `recoveries_total` exceeds `cases_total` on at least one date.")
    }
    return(invisible(data))
  }

  # linelist
  for (col in c("start_date", "outcome_date")) {
    if (!inherits(data[[col]], "Date")) bad_format("Column `", col, "` must be of class Date.")
  }
  if (anyNA(data$start_date)) bad_format("Column `start_date` contains missing values.")
  if (!is.character(data$event)) {
    bad_format('Column `event` must be a character column with values "death", "recovery" or NA.')
  }
  other <- setdiff(unique(data$event[!is.na(data$event)]), c("death", "recovery"))
  if (length(other) > 0) {
    bad_format('Column `event` must only contain "death", "recovery" or NA; found: ',
               paste0('"', other[seq_len(min(5, length(other)))], '"', collapse = ", "), ".")
  }
  if (!any(data$event %in% c("death", "recovery"))) {
    bad_format('Column `event` contains no "death" or "recovery" values.')
  }
  if (any(is.na(data$event) != is.na(data$outcome_date))) {
    bad_format("`event` and `outcome_date` must both be given or both be missing for every case.")
  }
  if (any(data$outcome_date < data$start_date, na.rm = TRUE)) {
    bad_format("`outcome_date` is before `start_date` for at least one case.")
  }
  invisible(data)
}

# The dates at which to calculate the CFR.
resolve_time_points <- function(time_points, data, dataset_type) {
  if (dataset_type == "count") {
    first <- min(data$date)
    last <- max(data$date)
  } else {
    first <- min(data$start_date)
    last <- if (all(is.na(data$outcome_date))) max(data$start_date) else
      max(data$outcome_date, na.rm = TRUE)
  }
  if (is.character(time_points) && length(time_points) == 1 && time_points == "final") {
    return(last)
  }
  if (!inherits(time_points, "Date") || length(time_points) == 0 || anyNA(time_points)) {
    stop('`time_points` must be "final" or a vector of dates of class Date.', call. = FALSE)
  }
  if (any(time_points < first)) {
    stop("`time_points` contains a date before the start of the data (", first, ").",
         call. = FALSE)
  }
  time_points
}

# Restrict the data to what was known at time point `tp`.
#  count:    rows up to tp.
#  linelist: cases starting by tp; cases with no outcome by tp (event NA, or outcome after
#            tp) get event "censored" in this internal copy;
#            adds `time` = days from start to outcome (or to tp if censored).
data_at_time <- function(data, dataset_type, tp) {
  if (dataset_type == "count") {
    return(data[data$date <= tp, , drop = FALSE])
  }
  d <- data[data$start_date <= tp, , drop = FALSE]
  observed <- !is.na(d$event) & !is.na(d$outcome_date) & d$outcome_date <= tp
  d$event[!observed] <- "censored"
  d$outcome_date[!observed] <- as.Date(NA)
  end <- d$outcome_date
  end[!observed] <- tp
  d$time <- as.numeric(end - d$start_date)
  attr(d, "analysis_date") <- tp
  d
}

# Final deaths / cases / recoveries from either kind of dataset.
tally <- function(d, dataset_type) {
  if (dataset_type == "count") {
    last <- d[which.max(d$date), ]
    list(cases = last$cases_total, deaths = last$deaths_total,
         recoveries = if ("recoveries_total" %in% names(d)) last$recoveries_total else NA_real_)
  } else {
    list(cases = nrow(d), deaths = sum(d$event == "death"),
         recoveries = sum(d$event == "recovery"))
  }
}

na_result <- function(msg = NULL, keep = NA_real_) {
  if (!is.null(msg)) warning(msg, call. = FALSE)
  list(CFR = keep, ci_lower = NA_real_, ci_upper = NA_real_)
}


# ---------------------------------------------------------------------------
# Estimators. Each takes the data as known at one time point and returns
# list(CFR, ci_lower, ci_upper).
# ---------------------------------------------------------------------------

# CI for a binomial proportion x / n.
#  "clopper_pearson": exact interval.
#  "jeffreys": equal-tailed interval from the Beta(x + 0.5, n - x + 0.5) posterior
#              (as Stata's `cii, jeffreys`); limit set to 0 if x = 0 and 1 if x = n.
binomial_ci <- function(x, n, alpha, method = "clopper_pearson") {
  x <- round(x)
  n <- round(n)
  if (method == "jeffreys") {
    lower <- if (x == 0) 0 else qbeta(alpha / 2, x + 0.5, n - x + 0.5)
    upper <- if (x == n) 1 else qbeta(1 - alpha / 2, x + 0.5, n - x + 0.5)
  } else {
    lower <- if (x == 0) 0 else qbeta(alpha / 2, x, n - x + 1)
    upper <- if (x == n) 1 else qbeta(1 - alpha / 2, x + 1, n - x)
  }
  c(lower, upper)
}

cfr_naive <- function(data, dataset_type, CI_method, alpha){
  tl <- tally(data, dataset_type)
  if (tl$cases <= 0) return(na_result("No cases observed; naive CFR is undefined."))
  ci <- binomial_ci(tl$deaths, tl$cases, alpha, CI_method)
  list(CFR = tl$deaths / tl$cases, ci_lower = ci[1], ci_upper = ci[2])
}

cfr_resolved <- function(data, dataset_type, CI_method, alpha){
  tl <- tally(data, dataset_type)
  resolved <- tl$deaths + tl$recoveries
  if (resolved <= 0) return(na_result("No resolved cases (deaths + recoveries); resolved CFR is undefined."))
  ci <- binomial_ci(tl$deaths, resolved, alpha, CI_method)
  list(CFR = tl$deaths / resolved, ci_lower = ci[1], ci_upper = ci[2])
}

# Nishiura et al. delay-adjusted CFR. Daily case counts are weighted by the
# probability F(age) that the outcome of a case of that age would be known; the
# death count is then treated as binomial (or Poisson, for many cases and a
# small CFR) in the expected number of known outcomes.
cfr_delay_adjusted <- function(data, dataset_type, CI_method, alpha,
                               delay_dist, parametric_delay_type){
  if (dataset_type == "count") {
    # fill in unreported days by carrying cumulative counts forward
    dates <- seq(min(data$date), max(data$date), by = "day")
    idx <- findInterval(as.numeric(dates), as.numeric(data$date))
    cases_daily <- diff(c(0, data$cases_total[idx]))
    total_deaths <- data$deaths_total[nrow(data)]
  } else {
    tp <- attr(data, "analysis_date")
    dates <- seq(min(data$start_date), tp, by = "day")
    cases_daily <- as.numeric(table(factor(as.character(data$start_date),
                                           levels = as.character(dates))))
    total_deaths <- sum(data$event == "death")
    if (is.null(delay_dist)) {
      delays <- as.numeric(data$outcome_date - data$start_date)[data$event == "death"]
      delay_dist <- tryCatch(fit_delay_cdf(delays, parametric_delay_type),
                             error = function(e) e)
      if (inherits(delay_dist, "error")) {
        return(na_result(paste0("Could not fit a delay distribution: ",
                                conditionMessage(delay_dist))))
      }
    }
  }

  ages <- rev(seq_along(cases_daily) - 1)   # days since each daily cohort started
  known_outcome_prob <- pmin(pmax(as.numeric(delay_dist(ages)), 0), 1)
  if (anyNA(known_outcome_prob)) {
    stop("`delay_dist` returned missing values.", call. = FALSE)
  }
  total_cases <- sum(cases_daily)
  total_outcomes <- sum(cases_daily * known_outcome_prob)

  p_mid <- if (round(total_outcomes) > 0) total_deaths / round(total_outcomes) else NA_real_
  out <- severity_profile_likelihood(total_cases, total_deaths, total_outcomes, p_mid, alpha)
  list(CFR = out[1], ci_lower = out[2], ci_upper = out[3])
}

# Profile-likelihood estimate and CI of the severity p for `total_deaths`
# deaths among `total_outcomes` (expected) known outcomes. c(estimate, lower, upper).
severity_profile_likelihood <- function(total_cases, total_deaths, total_outcomes,
                                        p_mid, alpha, poisson_threshold = 1000) {
  na3 <- rep(NA_real_, 3)
  if (sum(c(total_cases, total_deaths, total_outcomes) == 0) >= 2) return(na3)
  if (total_outcomes < total_deaths) {
    warning("Total deaths (", total_deaths, ") exceed expected known outcomes (",
            round(total_outcomes), "); returning NA.", call. = FALSE)
    return(na3)
  }
  n <- round(total_outcomes)
  k <- round(total_deaths)
  p_grid <- seq(1e-4, 1, by = 1e-4)
  if (total_cases < poisson_threshold || p_mid >= 0.05) {
    loglik <- dbinom(k, n, p_grid, log = TRUE)
  } else {
    warning("Total cases = ", total_cases, " and p = ", signif(p_mid, 3),
            ": using Poisson approximation to binomial likelihood.", call. = FALSE)
    loglik <- dpois(k, p_grid * n, log = TRUE)
  }
  if (!any(is.finite(loglik))) return(na3)
  keep <- which(loglik >= max(loglik[is.finite(loglik)]) - qchisq(1 - alpha, 1) / 2)
  c(p_grid[which.max(loglik)], p_grid[min(keep)], p_grid[max(keep)])
}

# Maximum-likelihood fit of a delay distribution; returns its CDF as a function.
# Delays are whole days, so a recorded 0 means "under a day" and is set to 0.5.
fit_delay_cdf <- function(delays, family) {
  x <- delays[is.finite(delays) & delays >= 0]
  if (length(x) < 2) stop("fewer than 2 usable delays")
  x <- ifelse(x == 0, 0.5, x)
  if (family == "lognormal") {
    mu <- mean(log(x))
    sdlog <- sqrt(mean((log(x) - mu)^2))
    if (!is.finite(sdlog) || sdlog <= 0) stop("delays have no variability")
    return(function(q) plnorm(q, mu, sdlog))
  }
  if (stats::var(x) == 0) stop("delays have no variability")
  if (family == "gamma") {
    start <- log(c(mean(x)^2 / stats::var(x), stats::var(x) / mean(x)))   # log(shape, scale)
    nll <- function(par) -sum(dgamma(x, shape = exp(par[1]), scale = exp(par[2]), log = TRUE))
    fit <- optim(start, nll, method = "Nelder-Mead", control = list(maxit = 2000))
    shape <- exp(fit$par[1]); scale <- exp(fit$par[2])
    return(function(q) pgamma(q, shape = shape, scale = scale))
  }
  start <- log(c(1, mean(x)))   # log(shape, scale)
  nll <- function(par) -sum(dweibull(x, shape = exp(par[1]), scale = exp(par[2]), log = TRUE))
  fit <- optim(start, nll, method = "Nelder-Mead", control = list(maxit = 2000))
  shape <- exp(fit$par[1]); scale <- exp(fit$par[2])
  function(q) pweibull(q, shape = shape, scale = scale)
}

# Aalen-Johansen cumulative incidence of death with recovery as a competing
# event. Returns the estimate and, if `variance`, its delta-method (Greenwood-type)
# variance, propagated through the state probabilities (survival, death CIF,
# recovery CIF).
aalen_johansen_death <- function(times, events, variance = TRUE) {
  event_times <- sort(unique(times[events %in% c("death", "recovery")]))
  if (length(event_times) == 0) return(c(estimate = NA_real_, variance = NA_real_))

  x <- c(1, 0, 0)               # survival, death CIF, recovery CIF
  Sigma <- matrix(0, 3, 3)
  recovery_seen <- FALSE
  for (t in event_times) {
    at_risk <- sum(times >= t)
    d_death <- sum(times == t & events == "death")
    d_rec <- sum(times == t & events == "recovery")
    if (at_risk <= 0 || d_death + d_rec <= 0) next
    if (d_rec > 0) recovery_seen <- TRUE

    a <- d_death / at_risk
    b <- d_rec / at_risk
    cc <- max(0, 1 - a - b)
    S_prev <- x[1]

    if (variance) {
      A <- matrix(c(cc, 0, 0,
                    a,  1, 0,
                    b,  0, 1), nrow = 3, byrow = TRUE)
      B <- matrix(c(-S_prev, -S_prev,
                     S_prev,  0,
                     0,       S_prev), nrow = 3, byrow = TRUE)
      V <- matrix(c(a * (1 - a), -a * b,
                    -a * b,      b * (1 - b)), nrow = 2, byrow = TRUE) / at_risk
      Sigma <- A %*% Sigma %*% t(A) + B %*% V %*% t(B)
      # everyone left died and no recovery was seen: the state is known exactly
      if (!recovery_seen && a == 1) Sigma <- matrix(0, 3, 3)
    }
    x <- c(S_prev * cc, x[2] + S_prev * a, x[3] + S_prev * b)
  }
  c(estimate = x[2], variance = if (variance) max(Sigma[2, 2], 0) else NA_real_)
}

# CI_method: "logit" (default), "wald", "loglog" or "bootstrap".
cfr_competing_risks <- function(data, CI_method, alpha, n_boot = 1000){
  times <- data$time
  events <- data$event
  fit <- aalen_johansen_death(times, events)
  if (is.na(fit["estimate"])) {
    return(na_result("No deaths or recoveries observed; competing-risks CFR is undefined."))
  }
  estimate <- unname(fit["estimate"])
  se <- sqrt(unname(fit["variance"]))
  z <- qnorm(1 - alpha / 2)
  eps <- 1e-12
  p <- min(max(estimate, eps), 1 - eps)

  if (CI_method == "bootstrap") {
    n <- length(times)
    boot <- vapply(seq_len(n_boot), function(i) {
      idx <- sample.int(n, n, replace = TRUE)
      unname(aalen_johansen_death(times[idx], events[idx], variance = FALSE)["estimate"])
    }, numeric(1))
    ci <- unname(stats::quantile(boot, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE))
  } else if (CI_method == "wald") {
    ci <- pmin(pmax(estimate + c(-1, 1) * z * se, 0), 1)
  } else if (CI_method == "loglog") {
    # g(p) = log(-log(1 - p)) is increasing, so the limits map straight back
    se_g <- se / ((1 - p) * -log(1 - p))
    g <- log(-log(1 - p))
    ci <- 1 - exp(-exp(g + c(-1, 1) * z * se_g))
  } else {
    se_logit <- se / (p * (1 - p))
    ci <- plogis(qlogis(p) + c(-1, 1) * z * se_logit)
  }
  list(CFR = estimate, ci_lower = ci[1], ci_upper = ci[2])
}

# Ghani et al. (2005) Kaplan-Meier style estimator: CFR = theta0 / (theta0 + theta1),
# where theta0 and theta1 are the cumulative incidences of death and recovery from
# the composite (death-or-recovery) survival curve. Follows the Stata module
# casefat (Griffin & Ghani).
# CI_method:
#   "logit" (default): covariance of the survival curve as in casefat's default;
#                      Wald CI on the logit scale.
#   "wald":            same covariance; untransformed Wald CI (casefat `untrans`).
#   "greenwood_logit": Greenwood covariance (casefat `greenwood`); logit Wald CI.
#   "greenwood_wald":  Greenwood covariance; untransformed Wald CI.
cfr_adapted_kaplan_meier <- function(data, CI_method, alpha){
  greenwood <- CI_method %in% c("greenwood_logit", "greenwood_wald")
  untrans <- CI_method %in% c("wald", "greenwood_wald")

  times <- data$time
  events <- data$event
  event_times <- sort(unique(times[events %in% c("death", "recovery")]))
  if (length(event_times) == 0) {
    return(na_result("No deaths or recoveries observed; Kaplan-Meier CFR is undefined."))
  }
  y <- vapply(event_times, function(t) sum(times >= t), numeric(1))
  d_death <- vapply(event_times, function(t) sum(times == t & events == "death"), numeric(1))
  d_rec <- vapply(event_times, function(t) sum(times == t & events == "recovery"), numeric(1))
  d_all <- d_death + d_rec

  if (sum(d_death) == 0 || sum(d_rec) == 0) {
    # with only one outcome observed the variance collapses to 0 and the CI would
    # be falsely precise, so withhold the estimate as well
    return(na_result(paste0("Need both outcomes for the Kaplan-Meier estimator (got ",
                            sum(d_death), " deaths, ", sum(d_rec), " recoveries).")))
  }

  h_dead <- d_death / y
  h_rec <- d_rec / y
  S <- cumprod(1 - d_all / y)             # composite survival, THETA in casefat
  S_prev <- c(1, S[-length(S)])
  theta0 <- sum(S_prev * h_dead)
  theta1 <- sum(S_prev * h_rec)
  cfr <- theta0 / (theta0 + theta1)

  # OM: covariance matrix of the composite survival curve at the event times
  m <- length(event_times)
  later <- pmax(row(diag(m)), col(diag(m)))
  earlier <- pmin(row(diag(m)), col(diag(m)))
  if (greenwood) {
    # Greenwood variance of S (undefined where everyone at risk fails at once);
    # OM[j, k] = var[k] * S[j] / S[k] for k < j, var[j] on the diagonal
    terms <- ifelse(y > d_all, d_all / (y * (y - d_all)), NA_real_)
    varS <- S^2 * cumsum(terms)
    if (anyNA(varS)) {
      return(na_result(paste0(
        "Greenwood variance is undefined: everyone still at risk had an outcome at the ",
        "last event time (no censored case outlasts it). Use CI_method = \"logit\" or \"wald\"."),
        keep = cfr))
    }
    OM <- varS[earlier] * S[later] / S[earlier]
    dim(OM) <- c(m, m)
    diag(OM) <- varS
  } else {
    # OM[j, k] = T[later] * (1 - T[earlier]) / nstar, nstar = (cases + events) / 2
    nstar <- (length(times) + sum(d_all)) / 2
    OM <- S[later] * (1 - S[earlier]) / nstar
    dim(OM) <- c(m, m)
  }

  var_dead <- sum(S^2 * h_dead / y) + drop(t(h_dead) %*% OM %*% h_dead)
  var_rec <- sum(S^2 * h_rec / y) + drop(t(h_rec) %*% OM %*% h_rec)
  cov01 <- drop(t(h_dead) %*% OM %*% h_rec)
  z <- qnorm(1 - alpha / 2)

  if (untrans) {
    var_cfr <- (theta1^2 * var_dead + theta0^2 * var_rec - 2 * theta0 * theta1 * cov01) /
      (theta0 + theta1)^4
    if (!is.finite(var_cfr) || var_cfr < 0) {
      return(list(CFR = cfr, ci_lower = NA_real_, ci_upper = NA_real_))
    }
    se <- sqrt(var_cfr)
    return(list(CFR = cfr, ci_lower = cfr - z * se, ci_upper = cfr + z * se))
  }
  var_logit <- var_dead / theta0^2 + var_rec / theta1^2 - 2 * cov01 / (theta0 * theta1)
  if (!is.finite(var_logit) || var_logit < 0) {
    return(list(CFR = cfr, ci_lower = NA_real_, ci_upper = NA_real_))
  }
  p <- min(max(cfr, 1e-12), 1 - 1e-12)
  list(CFR = cfr,
       ci_lower = plogis(qlogis(p) - z * sqrt(var_logit)),
       ci_upper = plogis(qlogis(p) + z * sqrt(var_logit)))
}

# Ghani et al. (2005) mixture model. A case dies with probability p; given death
# (recovery) the time to outcome follows a parametric death (recovery)
# distribution. Cases with no outcome yet contribute p*S_d(t) + (1-p)*S_r(t).
# Parameters on the optimiser scale: c(logit(p), then four log-scale or location
# parameters: death shape/scale, recovery shape/scale).
cfr_parametric_mixture <- function(data, CI_method, alpha, parametric_delay_type){
  family <- parametric_delay_type
  events <- data$event
  # whole-day delays: an outcome at day 0 means "under a day", so use 0.5
  times <- ifelse(data$time == 0 & events != "censored", 0.5, data$time)
  times <- pmax(times, 1e-12)

  n_dead <- sum(events == "death")
  n_rec <- sum(events == "recovery")
  if (n_dead == 0 || n_rec == 0) {
    # each outcome needs its own time distribution, which is unidentifiable with none
    return(na_result(paste0("Need both outcomes to fit the mixture model (got ",
                            n_dead, " deaths, ", n_rec, " recoveries).")))
  }

  p0 <- min(max(mean(events == "death"), 1e-4), 1 - 1e-4)
  bounds <- mixture_nuisance_bounds(family)

  # The likelihood can have local optima, so fit from two starting points (scale
  # starting at the median of all times, and of resolved times only) and keep the best.
  fits <- lapply(list(times, times[events != "censored"]), function(tt) {
    med <- max(stats::median(tt), 1e-3)
    start <- if (family == "lognormal") {
      c(qlogis(p0), log(0.5), log(med), log(0.5), log(med))
    } else {
      c(qlogis(p0), 0, log(med), 0, log(med))
    }
    optim(start, function(par) mixture_negloglik(par, times, events, family),
          method = "L-BFGS-B",
          lower = c(-10, bounds$lower), upper = c(10, bounds$upper),
          control = list(maxit = 5000))
  })
  fit <- fits[[which.min(vapply(fits, function(f) f$value, numeric(1)))]]
  if (fit$convergence != 0) {
    warning("Mixture model optimisation did not converge: ", fit$message, call. = FALSE)
  }
  ci <- mixture_lr_ci(times, events, family, fit, alpha)
  list(CFR = unname(plogis(fit$par[1])), ci_lower = ci[1], ci_upper = ci[2])
}

mixture_nuisance_bounds <- function(family) {
  if (family == "lognormal") {
    list(lower = c(-8, -10, -8, -10), upper = c(8, 10, 8, 10))
  } else {
    list(lower = rep(-8, 4), upper = rep(8, 4))
  }
}

# Negative log-likelihood of the mixture model.
mixture_negloglik <- function(par, times, events, family) {
  p <- plogis(par[1])
  nu <- par[-1]
  if (family == "gamma") {
    pdf_d <- function(x) dgamma(x, shape = exp(nu[1]), scale = exp(nu[2]))
    sf_d <- function(x) pgamma(x, shape = exp(nu[1]), scale = exp(nu[2]), lower.tail = FALSE)
    pdf_r <- function(x) dgamma(x, shape = exp(nu[3]), scale = exp(nu[4]))
    sf_r <- function(x) pgamma(x, shape = exp(nu[3]), scale = exp(nu[4]), lower.tail = FALSE)
  } else if (family == "weibull") {
    pdf_d <- function(x) dweibull(x, shape = exp(nu[1]), scale = exp(nu[2]))
    sf_d <- function(x) pweibull(x, shape = exp(nu[1]), scale = exp(nu[2]), lower.tail = FALSE)
    pdf_r <- function(x) dweibull(x, shape = exp(nu[3]), scale = exp(nu[4]))
    sf_r <- function(x) pweibull(x, shape = exp(nu[3]), scale = exp(nu[4]), lower.tail = FALSE)
  } else {   # lognormal: nu = (log sdlog_d, meanlog_d, log sdlog_r, meanlog_r)
    pdf_d <- function(x) dlnorm(x, meanlog = nu[2], sdlog = exp(nu[1]))
    sf_d <- function(x) plnorm(x, meanlog = nu[2], sdlog = exp(nu[1]), lower.tail = FALSE)
    pdf_r <- function(x) dlnorm(x, meanlog = nu[4], sdlog = exp(nu[3]))
    sf_r <- function(x) plnorm(x, meanlog = nu[4], sdlog = exp(nu[3]), lower.tail = FALSE)
  }
  # extreme parameter values during optimisation can give NaN/Inf densities;
  # treat them as 0 (a very unlikely observation) rather than breaking the optimiser
  finite0 <- function(v) { v[!is.finite(v)] <- 0; v }
  eps <- 1e-300
  dead <- events == "death"
  rec <- events == "recovery"
  cens <- !(dead | rec)
  ll <- suppressWarnings(
    sum(log(p + eps) + log(finite0(pdf_d(times[dead])) + eps)) +
      sum(log(1 - p + eps) + log(finite0(pdf_r(times[rec])) + eps)) +
      sum(log(p * finite0(sf_d(times[cens])) + (1 - p) * finite0(sf_r(times[cens])) + eps))
  )
  -ll
}

# Profile log-likelihood at a fixed CFR: maximise over the timing parameters only.
mixture_profile_loglik <- function(p_fixed, times, events, family, nuisance_start) {
  p_fixed <- min(max(p_fixed, 1e-12), 1 - 1e-12)
  bounds <- mixture_nuisance_bounds(family)
  fit <- optim(nuisance_start,
               function(nu) mixture_negloglik(c(qlogis(p_fixed), nu), times, events, family),
               method = "L-BFGS-B", lower = bounds$lower, upper = bounds$upper,
               control = list(maxit = 5000))
  -fit$value
}

# Likelihood-ratio CI for the CFR: the values of p whose profile log-likelihood
# is within qchisq(1 - alpha, 1) / 2 of the maximum, found by a grid search for a
# bracket followed by root finding.
mixture_lr_ci <- function(times, events, family, fit, alpha) {
  p_hat <- plogis(fit$par[1])
  target <- -fit$value - qchisq(1 - alpha, 1) / 2
  nuisance_start <- fit$par[-1]
  cache <- new.env()
  g <- function(p) {
    p <- min(max(p, 1e-12), 1 - 1e-12)
    key <- format(round(p, 12), digits = 15)
    if (is.null(cache[[key]])) {
      cache[[key]] <- mixture_profile_loglik(p, times, events, family, nuisance_start)
    }
    cache[[key]] - target
  }
  if (g(p_hat) < 0) return(c(NA_real_, NA_real_))

  eps <- 1e-8
  lower <- NA_real_
  if (p_hat > eps) {
    grid <- seq(eps, p_hat, length.out = 25)
    vals <- vapply(grid, g, numeric(1))
    i <- which(vals[-25] < 0 & vals[-1] >= 0)[1]
    lower <- if (!is.na(i)) {
      uniroot(g, c(grid[i], grid[i + 1]), tol = 1e-10)$root
    } else if (vals[1] >= 0) eps else NA_real_
  }
  upper <- NA_real_
  if (p_hat < 1 - eps) {
    grid <- seq(p_hat, 1 - eps, length.out = 25)
    vals <- vapply(grid, g, numeric(1))
    i <- which(vals[-25] >= 0 & vals[-1] < 0)[1]
    upper <- if (!is.na(i)) {
      uniroot(g, c(grid[i], grid[i + 1]), tol = 1e-10)$root
    } else if (vals[25] >= 0) 1 - eps else NA_real_
  }
  c(lower, upper)
}
