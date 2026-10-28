# Period / population restriction of FAERS objects.
#
# `faers_restrict()` narrows a standardized (and optionally de-duplicated)
# FAERS object to a sub-population defined by report-level conditions:
# reporting period, date window, sex, age range, age group and country.
# All conditions are combined by intersection: a report must satisfy every
# enabled condition to survive.  The implementation converges every condition
# onto the demo table's `primaryid` and delegates the actual row filtering to
# `faers_keep()`, so both backends (in-memory data.tables and duckdb) are
# supported for free and the `@standardization` / `@deduplication` flags are
# preserved untouched.

# Human-readable labels for the aggregated drop warning.
RESTRICT_LABELS <- c(
    period = "period",
    date_range = "date range",
    sex = "sex",
    age = "age range",
    age_grp = "age group",
    country = "country"
)

#' Restrict a FAERS object to a sub-population
#'
#' @description `faers_restrict()` narrows a standardized FAERS object to
#' reports matching period, date-window and demographic conditions. All
#' enabled conditions are intersected: a report survives only if it satisfies
#' every one of them. Reports excluded by any enabled condition — including
#' those dropped because their value is missing or unparseable (e.g. partial
#' dates) — are counted per condition and reported in a single aggregated
#' warning, mirroring the package's philosophy of never silently truncating
#' data.
#'
#' @param .object A [FAERSascii] object. Should be standardized (and
#'   de-duplicated for downstream signal scans).
#' @param .period Character vector of period keys, e.g. `c("2019q4",
#'   "2020q1")`. Case-insensitive. `NULL` (default) keeps all periods present
#'   in the object.
#' @param .date_range Length-2 vector `c(start, end)` of dates. Accepts
#'   `Date` objects or strings in `YYYY-MM-DD` or `YYYYMMDD` form; both ends
#'   are inclusive. `NULL` (default) applies no date window.
#' @param .date_field The demo date column to compare against `.date_range`.
#'   `"fda_dt"` (FDA receipt date, essentially complete) is the default
#'   because it is the most complete; `"event_dt"` and `"rept_dt"` are
#'   available. Rows whose date is missing or partial (e.g. year-only `1999`)
#'   cannot be placed in the window and are excluded (counted in the
#'   aggregated warning) unless `.date_partial = "keep"`.
#' @param .sex Character vector of sexes to keep, e.g. `c("Female", "Male")`.
#'   Valid values are `Female` and `Male`; unmatched non-NA values are an
#'   error. Reports with missing sex are excluded whenever `.sex` is enabled.
#' @param .age Length-2 numeric `c(lower, upper)` closed interval on
#'   `age_in_years` (the standardized age in years; all age codes are
#'   converted during standardization). `NULL` applies no age filter. Reports
#'   with missing age are excluded whenever `.age` is enabled.
#' @param .age_grp Character vector of age-group codes to keep. The recoded
#'   `age_grp` column ("A" adult, "E" elderly, "I" infant, "N" neonate) is
#'   missing for older quarters, so this filter is only reliable for recent
#'   periods; missing values are excluded whenever enabled. `NULL` applies no
#'   age-group filter.
#' @param .country Character vector of occurrence countries to keep
#'   (`occr_country`). `NULL` applies no country filter; missing countries are
#'   excluded whenever enabled.
#' @param .invert If `TRUE`, the semantics invert: reports matching **all**
#'   enabled conditions are removed instead of kept.
#' @param .date_partial How to handle reports whose date cannot be parsed
#'   (missing or partial dates like `1999` or `200301`) when a date window is
#'   enabled: `"exclude"` (default) removes them and counts them in the
#'   aggregated warning; `"keep"` retains them unfiltered.
#'
#' @return A [FAERSascii] object of the same backend (in-memory or duckdb)
#'   with the same `@standardization` / `@deduplication` flags as the input,
#'   restricted to reports satisfying all enabled conditions. The input
#'   object is never modified.
#' @seealso [faers_keep()], [faers_filter()]
#' @examples
#' \donttest{
#' # The bundled sample is already standardized.
#' data <- readRDS(system.file("extdata", "standardized_data.rds", package = "faers"))
#'
#' # keep female reports aged 18-85
#' sub <- faers_restrict(data, .sex = "Female", .age = c(18, 85))
#'
#' # period keys
#' sub2 <- faers_restrict(data, .period = c("2004q1", "2017q2"))
#'
#' # date window on the FDA receipt date
#' sub3 <- faers_restrict(data, .date_range = c("2004-01-01", "2004-12-31"))
#' }
#' @export
faers_restrict <- function(
    .object, .period = NULL, .date_range = NULL,
    .date_field = c("fda_dt", "event_dt", "rept_dt"), .sex = NULL,
    .age = NULL, .age_grp = NULL, .country = NULL, .invert = FALSE,
    .date_partial = c("exclude", "keep")) {
    .date_field <- match.arg(.date_field)
    .date_partial <- match.arg(.date_partial)
    .period <- restrict_check_period(.period)
    .date_range <- restrict_check_dates(.date_range)
    .sex <- restrict_check_values(.sex, c("Female", "Male"), ".sex")
    .age <- restrict_check_age(.age)
    .age_grp <- restrict_check_values(
        .age_grp, c("A", "E", "I", "N"), ".age_grp"
    )
    .country <- restrict_check_values_free(.country, ".country")

    demo <- faers_get(.object, "demo")
    enabled <- restrict_enabled_conditions(
        .period, .date_range, .sex, .age, .age_grp, .country
    )
    if (!length(enabled)) {
        cli::cli_warn(
            "No restriction condition is enabled: object returned unchanged"
        )
        return(.object)
    }

    # ---- per-condition keep masks on demo rows (TRUE = keep) ----
    masks <- vector("list", length(enabled))
    names(masks) <- enabled
    n_missing <- integer(length(enabled))
    names(n_missing) <- enabled
    for (cond in enabled) {
        res <- switch(cond,
            period = restrict_mask_period(demo, .period),
            date_range = restrict_mask_date(
                demo, .date_range, .date_field, .date_partial
            ),
            sex = restrict_mask_values(demo$sex, .sex),
            age = restrict_mask_age(demo$age_in_years, .age),
            age_grp = restrict_mask_values(demo$age_grp, .age_grp),
            country = restrict_mask_values(demo$occr_country, .country)
        )
        masks[[cond]] <- res$mask
        n_missing[[cond]] <- res$n_missing
    }

    keep_mask <- Reduce(`&`, masks)
    if (isTRUE(.invert)) {
        keep_mask <- !keep_mask
    }
    ids <- unique(as.character(demo$primaryid[keep_mask]))

    # ---- aggregated drop report (never silently truncate) ----
    n_removed <- sum(!keep_mask)
    if (n_removed > 0L && !isTRUE(.invert)) {
        msgs <- vapply(enabled, function(cond) {
            sprintf(
                "%s: %d removed (%d missing/unparseable)",
                RESTRICT_LABELS[[cond]], sum(!masks[[cond]]), n_missing[[cond]]
            )
        }, character(1L), USE.NAMES = FALSE)
        cli::cli_warn(c(
            "{.fn faers_restrict} removed {n_removed} report{?s} to build the requested sub-population", # nolint
            stats::setNames(msgs, rep("x", length(msgs)))
        ))
    } else if (n_removed > 0L) {
        cli::cli_warn(
            "{.fn faers_restrict} removed {n_removed} report{?s} matching all enabled conditions" # nolint
        )
    }

    # delegate to faers_keep so both backends + flags are handled identically
    faers_keep(.object, ids, invert = FALSE)
}

# ---------- validation helpers ----------

restrict_check_period <- function(period) {
    if (is.null(period)) {
        return(NULL)
    }
    if (!is.character(period) || !length(period)) {
        cli::cli_abort(
            "{.arg .period} must be a character vector like {.val {c(\"2019q4\", \"2020q1\")}}"
        )
    }
    period <- tolower(trimws(period))
    bad <- !grepl("^[0-9]{4}q[1-4]$", period)
    if (any(bad)) {
        cli::cli_abort(c(
            "{.arg .period} entries must look like {.val YYYYqN} with N in 1-4",
            x = "Invalid: {.val {period[bad]}}",
            i = "Examples: {.val {c(\"2020q1\", \"2020q2\")}}"
        ))
    }
    unique(period)
}

restrict_check_dates <- function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    if (!is.atomic(x) || length(x) != 2L || anyNA(x)) {
        cli::cli_abort(
            "{.arg .date_range} must be a length-2 vector c(start, end)"
        )
    }
    d <- vapply(x, restrict_parse_one_date, numeric(1L))
    d <- as.Date(d, origin = "1970-01-01")
    if (anyNA(d)) {
        cli::cli_abort(c(
            "{.arg .date_range} must be two parseable dates",
            i = "Accepted forms: `Date`, {.val YYYY-MM-DD} or {.val YYYYMMDD}"
        ))
    }
    if (d[[1L]] > d[[2L]]) {
        cli::cli_abort("{.arg .date_range} start must not be after the end")
    }
    d
}

restrict_check_values <- function(x, allowed, arg) {
    if (is.null(x)) {
        return(NULL)
    }
    if (!is.character(x) || !length(x)) {
        cli::cli_abort("{.arg {arg}} must be a non-empty character vector")
    }
    x <- trimws(x)
    bad <- unique(x[!x %chin% allowed])
    if (length(bad)) {
        cli::cli_abort(c(
            "{.arg {arg}} only accepts {.val {allowed}}",
            x = "Invalid value{?s}: {.val {bad}}"
        ))
    }
    unique(x)
}

# Free-form values (country): only shape is validated.
restrict_check_values_free <- function(x, arg) {
    if (is.null(x)) {
        return(NULL)
    }
    if (!is.character(x) || !length(x)) {
        cli::cli_abort("{.arg {arg}} must be `NULL` or a non-empty character vector")
    }
    unique(trimws(x))
}

restrict_check_age <- function(x) {
    if (is.null(x)) {
        return(NULL)
    }
    ok <- is.numeric(x) && length(x) == 2L && !anyNA(x) &&
        all(is.finite(x)) && x[[1L]] <= x[[2L]]
    if (!ok) {
        cli::cli_abort(c(
            "{.arg .age} must be `NULL` or a length-2 numeric c(lower, upper) with lower <= upper", # nolint
            i = "Ages are compared in years against the standardized {.field age_in_years} column" # nolint
        ))
    }
    as.numeric(x)
}

# ---------- mask builders ----------
# Each returns list(mask = logical, n_missing = integer).  Reports whose
# condition value is missing/unparseable get FALSE when the condition is
# enabled (dates: kept when .date_partial = "keep").

restrict_mask_values <- function(v, keep_vals) {
    v <- as.character(v)
    v_trim <- trimws(v)
    keep <- !is.na(v) & v_trim %chin% keep_vals
    list(mask = keep, n_missing = sum(is.na(v) | (!is.na(v) & v_trim == "")))
}

restrict_mask_age <- function(v, range) {
    keep <- !is.na(v) & v >= range[[1L]] & v <= range[[2L]]
    list(mask = keep, n_missing = sum(is.na(v)))
}

restrict_mask_period <- function(demo, period) {
    keys <- restrict_period_keys(demo$year, demo$quarter)
    list(mask = keys %chin% period, n_missing = 0L)
}

restrict_period_keys <- function(year, quarter) {
    q <- tolower(trimws(as.character(quarter)))
    q <- sub("^q", "", q) # tolerate both "q1" and "1" forms
    paste(as.character(year), q, sep = "q")
}

restrict_mask_date <- function(demo, date_range, date_field, date_partial) {
    raw <- demo[[date_field]]
    dates <- restrict_parse_dates(raw)
    ok <- !is.na(dates)
    keep <- ok & dates >= date_range[[1L]] & dates <= date_range[[2L]]
    if (identical(date_partial, "keep")) {
        keep[!ok] <- TRUE
    }
    list(mask = keep, n_missing = sum(!ok))
}

# Parse demo date columns (YYYYMMDD numerics, or character forms) into Date;
# partial dates (year-only or year-month) stay NA so they never match a
# window silently.
restrict_parse_dates <- function(x) {
    if (inherits(x, "Date")) {
        return(x)
    }
    if (is.numeric(x)) {
        s <- sprintf("%08.0f", x)
    } else {
        s <- trimws(as.character(x))
    }
    out <- as.Date(rep(NA_character_, length(s)))
    ok <- !is.na(s) & nchar(s) == 8L & !grepl("[^0-9]", s)
    if (any(ok)) {
        out[ok] <- as.Date(s[ok], format = "%Y%m%d")
    }
    # also accept ISO strings if a user hand-edited a table
    rest <- which(!ok & !is.na(s) & nchar(s) == 10L)
    if (length(rest)) {
        iso <- restrict_parse_one_date(s[rest])
        hit <- !is.na(iso)
        out[rest[hit]] <- as.Date(iso[hit], origin = "1970-01-01")
    }
    out
}

# Parse a single date value -> days since epoch (numeric) or NA_real_.
restrict_parse_one_date <- function(x) {
    if (inherits(x, "Date")) {
        return(as.numeric(x))
    }
    if (is.numeric(x)) {
        x <- sprintf("%08.0f", x)
    } else {
        x <- trimws(as.character(x))
    }
    d <- suppressWarnings(as.Date(x, format = "%Y%m%d"))
    if (is.na(d)) {
        d <- tryCatch(
            suppressWarnings(as.Date(x)), # ISO YYYY-MM-DD
            error = function(e) as.Date(NA)
        )
    }
    if (is.na(d)) {
        NA_real_
    } else {
        as.numeric(d)
    }
}

restrict_enabled_conditions <- function(
    period, date_range, sex, age, age_grp, country) {
    args <- list(
        period = period, date_range = date_range, sex = sex,
        age = age, age_grp = age_grp, country = country
    )
    names(args)[!vapply(args, is.null, logical(1L))]
}
