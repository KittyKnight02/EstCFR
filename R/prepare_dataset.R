#' Put a dataset into the format needed by estimate_cfr
#'
#' Takes a dataset with your own column names and date format and returns a data
#' frame in the format `estimate_cfr()` requires (see "Data formats" in
#' [estimate_cfr()]). Your data are not modified; a new data frame is returned.
#'
#' For count data, columns are renamed, dates parsed, rows sorted by date and, if
#' the counts are per-day (`cumulative = FALSE`), converted to cumulative counts.
#'
#' For a linelist, columns are renamed, dates parsed and outcomes recoded to
#' `"death"`, `"recovery"` or `NA` (no outcome yet). Outcome values other than the
#' death and recovery labels you give are treated as "no outcome yet" and reported in
#' a message. Cases with a missing start date or an outcome date before the start
#' date are dropped with a warning. Note that delay distributions used must reflect the
#' start date chosen (e.g. use a first report to death delay for death if start
#' date = first report). Where possible the delay distribution estimate should be from
#' the same outbreak as the data.
#'
#' Problems that cannot be fixed automatically, such as missing, negative or
#' non-numeric counts or duplicated dates, stop with an error that says where the
#' problem is. Cumulative counts that go down also stop with an error, unless
#' `correct_decreases = TRUE`.
#'
#' @param data Your dataset, a data frame.
#' @param dataset_type One of "count" (one row per date) or "linelist" (one row per
#'   case).
#' @param date_format Format of the dates if they are stored as text, as used by
#'   [strptime()], e.g. `"%d/%m/%Y"` for 31/12/2014 or `"%Y-%m-%d"` for 2014-12-31.
#'   Not needed for columns that are already of class Date.
#' @param cumulative Count data only: `TRUE` if `cases_total`, `deaths_total` (and
#'   `recoveries_total`) are already cumulative totals, `FALSE` if they are new cases,
#'   deaths and recoveries per day, which are then converted to cumulative totals.
#' @param correct_decreases Count data with `cumulative = TRUE` only. If `FALSE` (the
#'   default), cumulative counts that go down are an error. If `TRUE`, each downward
#'   revision is corrected by scaling all earlier values in that column by
#'   `current / previous` (the new value over the value just before it) and rounding to
#'   the nearest whole number (halves to even, as numpy's `rint`). Rows are processed
#'   in date order and each column separately, and later revisions act on the already
#'   corrected earlier values. The result is a series that never goes down.
#' @param date,cases_total,deaths_total,recoveries_total Count data only: the names of
#'   the columns in `data` holding the date and the cases, deaths and recoveries.
#'   `recoveries_total` is optional (only the "resolved" method needs it): leave it
#'   out if your data have no recoveries.
#' @param start_date,outcome_date,event Linelist only: the names of the columns in
#'   `data` holding the date the case started (onset or first report), the date of its
#'   outcome, and what the outcome was.
#' @param death_label,recovery_label Linelist only: the values in the `event` column
#'   that mean death and recovery respectively. Can be vectors, e.g. `c("died", "dead")`. Matching
#'   ignores case and surrounding spaces.
#' @return A data frame ready for [estimate_cfr()]. For count data: columns `date`,
#'   `cases_total`, `deaths_total` and, if given, `recoveries_total`. For a linelist:
#'   `start_date`, `outcome_date` and `event`.
#' @export
#' @examples
#' # per-day counts with day/month/year dates and different column names
#' daily <- data.frame(
#'   day = c("01/06/2014", "02/06/2014", "03/06/2014"),
#'   new_cases = c(5, 4, 11),
#'   new_deaths = c(1, 2, 5)
#' )
#' prepare_dataset(daily, dataset_type = "count", date_format = "%d/%m/%Y",
#'                 cumulative = FALSE, date = "day", cases_total = "new_cases",
#'                 deaths_total = "new_deaths", recoveries_total = NULL)
#'
#' # linelist
#' cases <- data.frame(
#'   onset = c("2014-06-01", "2014-06-02", "2014-06-02"),
#'   outcome_on = c("2014-06-08", "2014-06-12", NA),
#'   status = c("Died", "Recovered", "Unknown")
#' )
#' prepare_dataset(cases, dataset_type = "linelist", date_format = "%Y-%m-%d",
#'                 start_date = "onset", outcome_date = "outcome_on", event = "status",
#'                 death_label = "died", recovery_label = "recovered")
prepare_dataset <- function(data,
                            dataset_type,
                            date_format = NULL,
                            cumulative,
                            correct_decreases = FALSE,
                            date = "date",
                            cases_total = "cases_total",
                            deaths_total = "deaths_total",
                            recoveries_total = "recoveries_total",
                            start_date = "start_date",
                            outcome_date = "outcome_date",
                            event = "event",
                            death_label = "death",
                            recovery_label = "recovery"
                            ){
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (nrow(data) == 0) stop("`data` has no rows.", call. = FALSE)
  if (missing(dataset_type) || !is.character(dataset_type) || length(dataset_type) != 1 ||
      !dataset_type %in% c("count", "linelist")) {
    stop('`dataset_type` must be given, as either "count" or "linelist".', call. = FALSE)
  }
  if (!is.null(date_format) && (!is.character(date_format) || length(date_format) != 1)) {
    stop('`date_format` must be a single text string such as "%d/%m/%Y".', call. = FALSE)
  }

  if (dataset_type == "count") {
    if (missing(cumulative) || !is.logical(cumulative) || length(cumulative) != 1 ||
        is.na(cumulative)) {
      stop("For count data, `cumulative` must be given as TRUE (counts are already ",
           "cumulative totals) or FALSE (counts are per day).", call. = FALSE)
    }
    if (!is.logical(correct_decreases) || length(correct_decreases) != 1 ||
        is.na(correct_decreases)) {
      stop("`correct_decreases` must be TRUE or FALSE.", call. = FALSE)
    }
    out <- prepare_counts(data, date_format, cumulative, correct_decreases, date, cases_total, deaths_total,
                          recoveries_total, recoveries_given = !missing(recoveries_total))
  } else {
    out <- prepare_linelist(data, date_format, start_date, outcome_date, event,
                            death_label, recovery_label)
  }

  # safety net: what we return must pass the checks estimate_cfr() will run
  tryCatch(check_dataset(out, dataset_type),
           error = function(e) stop("prepare_dataset() could not produce a valid dataset: ",
                                    conditionMessage(e), call. = FALSE))
  out
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Pull one column out of `data` by name.
get_column <- function(data, col, arg) {
  if (!is.character(col) || length(col) != 1 || is.na(col)) {
    stop("`", arg, "` must be the name of a column in `data`, as a single text string.",
         call. = FALSE)
  }
  if (!col %in% names(data)) {
    stop("Column \"", col, "\" (given as `", arg, "`) is not in `data`. Columns found: ",
         paste0('"', names(data), '"', collapse = ", "), ".", call. = FALSE)
  }
  data[[col]]
}

# Convert a column to class Date. Text dates are read with `date_format`.
parse_dates <- function(x, date_format, arg) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) {
    stop("Column given as `", arg, "` must hold dates, or dates written as text.",
         call. = FALSE)
  }
  if (is.null(date_format)) {
    stop("Column given as `", arg, "` holds dates as text, so `date_format` is needed, ",
         'e.g. "%d/%m/%Y" for 31/12/2014 or "%Y-%m-%d" for 2014-12-31.', call. = FALSE)
  }
  x <- trimws(x)
  x[x == ""] <- NA
  out <- as.Date(x, format = date_format)
  bad <- is.na(out) & !is.na(x)
  if (any(bad)) {
    examples <- unique(x[bad])
    stop(sum(bad), " value(s) in the column given as `", arg, "` do not match ",
         'date_format "', date_format, '", e.g. ',
         paste0('"', examples[seq_len(min(3, length(examples)))], '"', collapse = ", "), ".",
         call. = FALSE)
  }
  out
}

# Convert a count column to numeric, refusing text that is not a number.
parse_counts <- function(x, arg) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    num <- suppressWarnings(as.numeric(trimws(x)))
    bad <- is.na(num) & !is.na(x) & trimws(x) != ""
    if (any(bad)) {
      stop("Column given as `", arg, "` has values that are not numbers, e.g. \"",
           x[bad][1], '".', call. = FALSE)
    }
    return(num)
  }
  if (!is.numeric(x)) stop("Column given as `", arg, "` must be numeric.", call. = FALSE)
  as.numeric(x)
}

# Scale earlier values back in time at each downward revision: when
# x[i] < x[i-1], all of x[1..i-1] become round(x[1..i-1] * x[i] / x[i-1]).
# R's round() and numpy's rint both round halves to even.
correct_downward_revisions <- function(x) {
  for (i in seq_along(x)[-1]) {
    if (x[i] < x[i - 1]) {
      x[seq_len(i - 1)] <- round(x[seq_len(i - 1)] * (x[i] / x[i - 1]))
    }
  }
  x
}

prepare_counts <- function(data, date_format, cumulative, correct_decreases, date, cases_total, deaths_total,
                           recoveries_total, recoveries_given) {
  # recoveries are optional: if the argument was left alone and the default column is
  # not there, there simply are none
  if (is.null(recoveries_total) ||
      (!recoveries_given && !recoveries_total %in% names(data))) {
    recoveries_total <- NULL
  }

  out <- data.frame(date = parse_dates(get_column(data, date, "date"), date_format, "date"),
                    cases_total = parse_counts(get_column(data, cases_total, "cases_total"),
                                               "cases_total"),
                    deaths_total = parse_counts(get_column(data, deaths_total, "deaths_total"),
                                                "deaths_total"))
  if (!is.null(recoveries_total)) {
    out$recoveries_total <- parse_counts(get_column(data, recoveries_total, "recoveries_total"),
                                         "recoveries_total")
  }
  count_cols <- setdiff(names(out), "date")

  if (anyNA(out$date)) {
    stop(sum(is.na(out$date)), " row(s) have a missing date; every row needs a date.",
         call. = FALSE)
  }
  if (anyDuplicated(out$date)) {
    dup <- unique(out$date[duplicated(out$date)])
    stop("Dates must be unique (one row per date). Duplicated: ",
         paste(format(dup[seq_len(min(3, length(dup)))]), collapse = ", "), ".", call. = FALSE)
  }
  for (col in count_cols) {
    if (anyNA(out[[col]])) {
      stop("Column `", col, "` has ", sum(is.na(out[[col]])), " missing value(s). Fill them ",
           "in (e.g. with 0 if nothing was reported) before calling prepare_dataset().",
           call. = FALSE)
    }
    if (any(out[[col]] < 0)) {
      stop("Column `", col, "` has negative values.", call. = FALSE)
    }
  }

  out <- out[order(out$date), , drop = FALSE]
  rownames(out) <- NULL

  if (!cumulative) {
    out[count_cols] <- lapply(out[count_cols], cumsum)
  } else if (correct_decreases) {
    too_many <- function(d) {
      bad <- d$deaths_total > d$cases_total
      if (!is.null(d$recoveries_total)) bad <- bad | d$deaths_total + d$recoveries_total > d$cases_total
      bad
    }
    already_bad <- any(too_many(out))
    out[count_cols] <- lapply(out[count_cols], correct_downward_revisions)
    if (!already_bad && any(too_many(out))) {
      stop("Correcting the downward revisions (`correct_decreases = TRUE`) made deaths",
           if (!is.null(out$recoveries_total)) " (or deaths + recoveries)" else "",
           " exceed cases on ", format(out$date[which(too_many(out))[1]]),
           ", which was not the case in your original data. Each column is scaled and ",
           "rounded separately, so the corrected series can become inconsistent. Check the ",
           "revisions in your data around that date or apply a different correction.", call. = FALSE)
    }
  } else {
    for (col in count_cols) {
      drop_at <- which(diff(out[[col]]) < 0)
      if (length(drop_at) > 0) {
        stop("`", col, "` is not cumulative: it goes down on ", format(out$date[drop_at[1] + 1]),
             ". If these are per-day counts, use cumulative = FALSE; to correct downward ",
           "revisions instead, use correct_decreases = TRUE.", call. = FALSE)
      }
    }
  }

  if (any(out$deaths_total > out$cases_total)) {
    stop("`deaths_total` exceeds `cases_total` on ",
         format(out$date[which(out$deaths_total > out$cases_total)[1]]), ".", call. = FALSE)
  }
  if (!is.null(recoveries_total) &&
      any(out$deaths_total + out$recoveries_total > out$cases_total)) {
    stop("`deaths_total` + `recoveries_total` exceeds `cases_total` on ",
         format(out$date[which(out$deaths_total + out$recoveries_total > out$cases_total)[1]]),
         ".", call. = FALSE)
  }
  out
}

prepare_linelist <- function(data, date_format, start_date, outcome_date, event,
                             death_label, recovery_label) {
  out <- data.frame(
    start_date = parse_dates(get_column(data, start_date, "start_date"), date_format, "start_date"),
    outcome_date = parse_dates(get_column(data, outcome_date, "outcome_date"), date_format,
                               "outcome_date"),
    event = NA_character_,
    stringsAsFactors = FALSE
  )

  # recode outcomes
  if (!is.character(death_label) || !is.character(recovery_label) ||
      length(death_label) == 0 || length(recovery_label) == 0) {
    stop("`death_label` and `recovery_label` must be text values (or vectors of them).",
         call. = FALSE)
  }
  clean <- function(v) tolower(trimws(as.character(v)))
  if (length(intersect(clean(death_label), clean(recovery_label))) > 0) {
    stop("`death_label` and `recovery_label` must not share any values.", call. = FALSE)
  }
  raw <- get_column(data, event, "event")
  lc <- clean(raw)
  out$event <- ifelse(lc %in% clean(death_label), "death",
                      ifelse(lc %in% clean(recovery_label), "recovery", NA_character_))
  unrecognised <- unique(as.character(raw)[is.na(out$event) & !is.na(raw) & trimws(lc) != ""])
  if (length(unrecognised) > 0) {
    message("Treated as 'no outcome yet': ",
            paste0('"', unrecognised[seq_len(min(5, length(unrecognised)))], '"', collapse = ", "),
            if (length(unrecognised) > 5) ", ..." else "",
            ". (Recognised: death = ", paste0('"', death_label, '"', collapse = "/"),
            ", recovery = ", paste0('"', recovery_label, '"', collapse = "/"), ".)")
  }

  # an outcome needs both its type and its date
  no_date <- !is.na(out$event) & is.na(out$outcome_date)
  if (any(no_date)) {
    warning(sum(no_date), " case(s) have an outcome but no outcome date; treated as no ",
            "outcome yet.", call. = FALSE)
    out$event[no_date] <- NA_character_
  }
  no_event <- is.na(out$event) & !is.na(out$outcome_date)
  if (any(no_event)) {
    warning(sum(no_event), " case(s) have an outcome date but no recognised death or ",
            "recovery; the outcome date was ignored.", call. = FALSE)
    out$outcome_date[no_event] <- as.Date(NA)
  }

  # impossible records
  no_start <- is.na(out$start_date)
  if (any(no_start)) {
    warning("Dropped ", sum(no_start), " case(s) with a missing start date.", call. = FALSE)
  }
  before <- !no_start & !is.na(out$outcome_date) & out$outcome_date < out$start_date
  if (any(before)) {
    warning("Dropped ", sum(before), " case(s) whose outcome date is before the start date.",
            call. = FALSE)
  }
  out <- out[!(no_start | before), , drop = FALSE]
  rownames(out) <- NULL

  if (!any(out$event %in% c("death", "recovery"))) {
    stop("No deaths or recoveries found in the `event` column (using death = ",
         paste0('"', death_label, '"', collapse = "/"), ", recovery = ",
         paste0('"', recovery_label, '"', collapse = "/"), "). Check `death_label` and ",
         "`recovery_label`.", call. = FALSE)
  }
  out
}
