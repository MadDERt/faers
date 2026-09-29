# Stratified whole-database pharmacovigilance signal scanning.
#
# Extends the scan engine of scan.R / scan-db.R: the counting key grows from
# (drug, event) to (stratum, drug, event), where strata are categorical
# columns of the standardized `demo` field.  Every report carries exactly one
# stratum value, so the report atomicity -- and thus the exactness of the
# chunked memory backend -- is preserved.  Each combination gets its own 2x2
# contingency table whose background is the stratum itself (a + b + c + d ==
# n_stratum).
#
# Missing strata values are grouped into an explicit "Missing" stratum by
# default so no report is silently dropped; `.na_stratum = "drop"` excludes
# them instead.  Columns whose missing share exceeds STRATA_MAX_NA_SHARE
# trigger a warning before any counting happens.

STRATA_MISSING_LABEL <- "Missing"
STRATA_MAX_NA_SHARE <- 0.3
# Names that strata columns must not collide with in the output.
STRATA_RESERVED <- c(
    "drug", "event", "a", "b", "c", "d", "n_drug", "n_event",
    "n_stratum", "strat_id", "expected"
)

#' Stratified whole-database pharmacovigilance signal scanning
#'
#' @description `faers_phv_scan_stratified()` runs the whole-database scan of
#' [faers_phv_scan] within strata defined by columns of the standardized
#' `demo` field. Every combination of `(stratum, drug, event)` gets its own
#' 2x2 contingency table, so the disproportionality analysis contrasts drug
#' exposure against the **stratum-internal** background rather than the whole
#' database. This surfaces subgroup-specific signals (e.g. a drug-event pair
#' standing out only among elderly female reports) that a whole-database scan
#' averages away.
#'
#' @details
#' Reports with missing values in a stratification column are grouped into an
#' explicit `"Missing"` stratum by default, so no report is silently dropped
#' and stratum sizes are always visible through the returned `n_stratum`
#' column. The 2x2 table of the `"Missing"` stratum is internally valid (it
#' contrasts exposure against non-exposure within that stratum), but it must
#' not be interpreted as a demographic subgroup. When a stratification column
#' is missing for more than 30% of the reports a warning is raised, because
#' the `"Missing"` stratum then tends to dominate the results; choose another
#' column or set `.na_stratum = "drop"` in that case. Dropping is opt-in
#' because it silently shrinks every stratum's background and missingness
#' itself can be informative (e.g. older quarterly formats lack `age_grp`).
#'
#' The drug and event sides are normalized exactly as in [faers_phv_scan]
#' (`lower(trim(.))`, empty/`unknown`/`?` drug names dropped, empty events
#' dropped), and drug names are aggregated on their exact normalized value,
#' so brand-name variants are not merged.
#'
#' Pairs below `.min_a` are evaluated **per stratum**, which is what keeps
#' stratified output - and any expensive method - tractable: stratifying
#' multiplies the number of candidate pairs while shrinking every stratum's
#' background, so a larger `.min_a` is often appropriate.
#'
#' Both backends produce identical results: the memory backend chunks over
#' `primaryid` blocks (exact, because a report never spans two chunks) and
#' the `database = "duckdb"` backend runs the whole stratified aggregation as
#' a single out-of-core SQL statement.
#'
#' @param .object A [FAERSascii] object, standardized with [faers_standardize]
#' and de-duplicated with [faers_dedup].
#' @param .strata A character vector of column names of the standardized
#' `demo` field defining the strata; multiple columns are crossed (a report
#' belongs to the stratum given by its combination of values). All columns
#' must be categorical (character) columns of the standardized `demo` data.
#' Good candidates: `"sex"` (default), `"age_grp"`, `occp_cod`,
#' `occr_country`, `reporter_country`.
#' @param .events A string, the event column of the standardized `reac`
#' field; any column accepted by [faers_phv_scan] (raw `pt`, `meddra_code`,
#' `meddra_pt`, or a MedDRA hierarchy column). Defaults to `"pt"`.
#' @param .drug_field A string, the column of the standardized `drug` field
#' used to define drugs: `"drugname"` (default) or `"prod_ai"`.
#' @param .drug_pattern An optional regular expression; only drug names
#' matching it (case-insensitively) are kept in the result.
#' @param .min_a A single integer, the minimum number of reports co-listing a
#' drug-event pair **within a stratum** for it to be kept. Defaults to `3L`.
#' @param .na_stratum Either `"keep"` (default), which groups reports with
#' missing strata values into an explicit `"Missing"` stratum, or `"drop"`,
#' which excludes them from the stratified scan.
#' @param .methods An atomic character, the disproportionality methods passed
#' to [phv_signal]. Defaults to `c("ror", "prr")`.
#' @param .chunk_size A single integer, the number of unique `primaryid`s
#' processed per chunk by the memory backend. Only relevant for
#' `database = "memory"`.
#' @param .phv_signal_params Other arguments passed to [phv_signal].
#' @param BPPARAM A [BiocParallel::BiocParallelParam-class] object.
#' @param ... Unused arguments, included for S4 generic/method consistency.
#' @return A [data.table][data.table::data.table] with one row per kept
#' stratum-drug-event combination, sorted by the lower bound of the first
#' method's confidence interval in descending order: the strata columns
#' (named after `.strata`), the drug column (named after `.drug_field`), the
#' event column (named after `.events`), the contingency table columns `a`,
#' `b`, `c`, `d`, the stratum size `n_stratum`, and the columns of
#' [phv_signal]. The four cells always sum to `n_stratum`, the number of
#' distinct reports of the stratum, and `"Missing"` strata mark reports with
#' missing strata values (they are internally valid comparisons but must not
#' be interpreted demographically).
#' @examples
#' # the sample data below is standardized but not de-duplicated; real usage
#' # requires faers_standardize() + faers_dedup() before scanning
#' std_data <- readRDS(system.file("extdata", "standardized_data.rds",
#'     package = "faers"
#' ))
#' std_data@deduplication <- TRUE
#' \dontrun{
#' # scan drug-event pairs separately for each reported sex
#' res <- faers_phv_scan_stratified(std_data, .strata = "sex", .min_a = 1L)
#' head(res)
#'
#' # cross two demographic columns and restrict to a drug class
#' faers_phv_scan_stratified(data, .strata = c("sex", "age_grp"),
#'     .drug_pattern = "insulin|humulin", .min_a = 3L
#' )
#' }
#' @seealso [faers_phv_scan], [phv_signal]
#' @export
#' @aliases faers_phv_scan_stratified
#' @name faers_phv_scan_stratified
methods::setGeneric("faers_phv_scan_stratified", function(.object, ...) {
    standardGeneric("faers_phv_scan_stratified")
})

#' @rdname faers_phv_scan_stratified
#' @export
#' @method faers_phv_scan_stratified FAERSascii
methods::setMethod("faers_phv_scan_stratified", "FAERSascii", function(
    .object, .strata = "sex", .events = "pt", .drug_field = "drugname",
    .drug_pattern = NULL, .min_a = 3L, .methods = c("ror", "prr"),
    .na_stratum = c("keep", "drop"), .chunk_size = 1e6L,
    .phv_signal_params = list(), BPPARAM = BiocParallel::SerialParam()
) {
    .na_stratum <- match.arg(.na_stratum)
    assert_string(.drug_field, allow_empty = FALSE)
    assert_string(.events, allow_empty = FALSE)
    assert_string(.drug_pattern, allow_null = TRUE, allow_empty = FALSE)
    assert_number_whole(.min_a, min = 1)
    assert_number_whole(.chunk_size, min = 1)
    assert_(.phv_signal_params, is.list, "a list")
    if (!is.character(.strata) || !length(.strata) || anyDuplicated(.strata)) {
        cli::cli_abort(c(
            "{.arg .strata} must be a non-empty character vector of distinct column names",
            i = "Good candidates: {.val {c(\"sex\", \"age_grp\", \"occp_cod\", \"occr_country\")}}"
        ))
    }
    strata_bad <- intersect(.strata, STRATA_RESERVED)
    if (length(strata_bad) > 0L) {
        cli::cli_abort(
            "{.arg .strata} must not use the reserved column{?s} {.val {strata_bad}}"
        )
    }
    if (!.object@standardization) {
        cli::cli_abort("{.arg .object} must be standardized using {.fn faers_standardize}")
    }
    if (!.object@deduplication) {
        cli::cli_abort(c(
            "{.arg .object} must be de-duplicated using {.fn faers_dedup}",
            i = "Scanning raw duplicates would heavily bias the drug-event counts"
        ))
    }

    # ---- resolve the strata columns (categorical demo columns) ----
    demo_cols <- if (is.null(.object@db)) {
        names(.object@data$demo)
    } else {
        DBI::dbListFields(.object@db@con, db_field_table("demo"))
    }
    strata_bad <- setdiff(.strata, demo_cols)
    if (length(strata_bad) > 0L) {
        cli::cli_abort(c(
            "{.val {strata_bad}} {?is/are} not column{?s} of the standardized {.field demo} data",
            i = "Available columns: {.val {demo_cols}}"
        ))
    }
    .strata_check_types(.object, .strata)
    strata_na_shares <- if (is.null(.object@db)) {
        .strata_missing_shares(.object@data$demo, .strata)
    } else {
        .strata_missing_shares_db(.object@db@con, .strata)
    }
    for (col in names(strata_na_shares)) {
        if (strata_na_shares[[col]] > STRATA_MAX_NA_SHARE) {
            cli::cli_warn(c(
                "Stratification column {.field {col}} is missing for {round(100 * strata_na_shares[[col]])}% of the reports",
                i = "The {.val {STRATA_MISSING_LABEL}} stratum will dominate the results",
                i = "Consider {.code .na_stratum = \"drop\"} or another stratification column"
            ))
        }
    }

    # ---- resolve the drug column (same rules as faers_phv_scan) ----
    drug_cols <- if (is.null(.object@db)) {
        names(.object@data$drug)
    } else {
        DBI::dbListFields(.object@db@con, db_field_table("drug"))
    }
    if (!.drug_field %chin% drug_cols) {
        cli::cli_abort(c(
            "{.arg .drug_field} must be a column of the {.field drug} data",
            i = "Available columns: {.val {drug_cols}}"
        ))
    }

    # ---- resolve the event column (same whitelist as faers_phv_scan) ----
    event_allowed <- c(
        "pt", "meddra_code", "meddra_pt",
        names(faers_meddra(.object, use = "hierarchy"))
    )
    if (!.events %chin% event_allowed) {
        cli::cli_abort(c(
            "{.arg .events} must be a column of the standardized {.field reac} data",
            i = "Available columns: {.val {event_allowed}}"
        ))
    }

    # ---- count all stratum x drug x event combinations ----
    if (!is.null(.object@db)) {
        counts <- scan_strata_counts_db(
            .object, strata_cols = .strata, drug_col = .drug_field,
            event_col = .events, min_a = .min_a, na_stratum = .na_stratum
        )
    } else {
        counts <- scan_strata_counts_memory(
            .object, strata_cols = .strata, drug_col = .drug_field,
            event_col = .events, min_a = .min_a, chunk_size = .chunk_size,
            na_stratum = .na_stratum
        )
    }

    # ---- build one 2x2 table per (stratum, pair); background = n_stratum ----
    counts[, b := n_drug - a] # nolint
    counts[, c := n_event - a] # nolint
    counts[, d := n_stratum - (n_drug + n_event - a)] # nolint
    counts[, c("n_drug", "n_event") := NULL] # nolint
    data.table::setcolorder(
        counts, c(.strata, "drug", "event", "a", "b", "c", "d", "n_stratum")
    )

    # ---- optional drug whitelist ----
    if (!is.null(.drug_pattern)) {
        keep <- grepl(.drug_pattern, counts$drug, ignore.case = TRUE)
        if (!any(keep)) {
            pattern <- .drug_pattern
            cli::cli_warn("No drug name matches the pattern {.val {pattern}}")
        }
        counts <- counts[keep]
    }

    # ---- warn early about expensive methods on big candidate sets ----
    used_methods <- .methods %||% SCAN_ALL_METHODS
    if (length(intersect(used_methods, SCAN_EXPENSIVE_METHODS)) &&
        nrow(counts) > SCAN_EXPENSIVE_THRESHOLD) {
        cli::cli_warn(c(
            "{nrow(counts)} candidate pairs will be analyzed with the expensive method{?s} {.val {intersect(used_methods, SCAN_EXPENSIVE_METHODS)}}",
            i = "Consider a larger {.arg .min_a} or dropping these methods"
        ))
    }

    # ---- disproportionality analysis on the whole candidate table ----
    signal <- do.call(
        phv_signal,
        c(
            counts[, c("a", "b", "c", "d")],
            list(methods = .methods, BPPARAM = BPPARAM),
            .phv_signal_params
        )
    )
    counts[, names(signal) := signal]

    # ---- finalize: rename and sort by the first method's CI lower bound ----
    data.table::setnames(counts, c("drug", "event"), c(.drug_field, .events))
    ci_col <- grep("_ci_low$", names(counts), value = TRUE)[1L]
    if (!is.na(ci_col)) {
        data.table::setorderv(
            counts, c(ci_col, .strata, .drug_field, .events),
            order = c(-1L, rep(1L, length(.strata) + 2L)), na.last = TRUE
        )
    } else {
        data.table::setorderv(
            counts, c("a", .strata, .drug_field, .events),
            order = c(-1L, rep(1L, length(.strata) + 1L)), na.last = TRUE
        )
    }
    counts[]
})

# Validate that every strata column is categorical (character on the memory
# backend, VARCHAR on the duckdb backend) so both backends group on the same
# values and produce identical output.
.strata_check_types <- function(object, strata_cols) {
    if (is.null(object@db)) {
        demo <- object@data$demo
        missing_cols <- setdiff(strata_cols, names(demo))
        if (length(missing_cols) > 0L) {
            cli::cli_abort(c(
                "{.val {missing_cols}} {?is/are} not column{?s} of the standardized {.field demo} data",
                i = "Available columns: {.val {setdiff(names(demo), c(\"primaryid\", \"caseid\"))}}"
            ))
        }
        for (col in strata_cols) {
            if (!is.character(demo[[col]])) {
                cli::cli_abort(c(
                    "Stratification column {.field {col}} must be a categorical (character) column",
                    i = "Numeric columns such as {.val age} must be grouped into a categorical column first (e.g. {.val age_grp})"
                ))
            }
        }
        NULL
    } else {
        demo_types <- DBI::dbGetQuery(object@db@con, sprintf(
            "DESCRIBE %s", db_field_table("demo")
        ))
        bad <- setdiff(strata_cols, demo_types$column_name)
        if (length(bad) > 0L) {
            cli::cli_abort(c(
                "{.val {bad}} {?is/are} not column{?s} of the standardized {.field demo} data",
                i = "Available columns: {.val {demo_types$column_name}}"
            ))
        }
        non_char <- strata_cols[!grepl("VARCHAR|CHAR|TEXT",
            demo_types$column_type[match(strata_cols, demo_types$column_name)],
            ignore.case = TRUE
        )]
        if (length(non_char) > 0L) {
            cli::cli_abort(c(
                "Stratification column{?s} {.val {non_char}} must be categorical (VARCHAR) columns",
                i = "Numeric columns must be grouped into a categorical column first"
            ))
        }
        NULL
    }
}

# Share of missing (NA or empty after trimming) values per strata column.
.strata_missing_shares <- function(demo, strata_cols) {
    vapply(strata_cols, function(col) {
        v <- demo[[col]]
        mean(is.na(v) | trimws(v) == "")
    }, numeric(1L))
}

# Missing shares via one SQL aggregate (the demo table is never collected).
.strata_missing_shares_db <- function(con, strata_cols) {
    selects <- vapply(seq_along(strata_cols), function(i) {
        q <- db_qident(strata_cols[i])
        sprintf(
            "AVG(CASE WHEN %1$s IS NULL OR trim(%1$s) = '' THEN 1.0 ELSE 0.0 END) AS miss_%2$d",
            q, i
        )
    }, character(1L))
    row <- DBI::dbGetQuery(
        con,
        sprintf(
            "SELECT COUNT(*) AS n, %s FROM %s",
            paste(selects, collapse = ", "), db_field_table("demo")
        )
    )
    stats::setNames(
        as.numeric(unlist(row[paste0("miss_", seq_along(strata_cols))])),
        strata_cols
    )
}

# Memory-backend stratified counting: same normalization, integer coding and
# chunking as `scan_counts_memory()`, with the stratum key carried through
# every step.  Reports are atomic in their stratum (one demo row per
# primaryid), so per-chunk counts sum to the exact database-wide counts.
scan_strata_counts_memory <- function(object, strata_cols, drug_col, event_col,
                                      min_a, chunk_size, na_stratum) {
    demo <- faers_get(object, "demo")
    strata_dt <- unique(
        demo[, c("primaryid", strata_cols), with = FALSE],
        by = "primaryid"
    )
    for (col in strata_cols) {
        data.table::set(
            strata_dt, j = col, value = as.character(strata_dt[[col]])
        )
    }
    if (identical(na_stratum, "keep")) {
        for (col in strata_cols) {
            v <- strata_dt[[col]]
            data.table::set(
                strata_dt, j = col,
                value = data.table::fifelse(
                    is.na(v) | trimws(v) == "", STRATA_MISSING_LABEL, v
                )
            )
        }
    } else {
        keep <- !Reduce(`|`, lapply(
            strata_dt[, strata_cols, with = FALSE],
            function(v) is.na(v) | trimws(v) == ""
        ))
        strata_dt <- strata_dt[keep]
    }
    strata_dt[, strat_id := .GRP, by = strata_cols] # nolint
    strat_map <- unique(strata_dt[, c("strat_id", strata_cols), with = FALSE])
    strat_tot <- strata_dt[, .(n_stratum = .N), by = strat_id] # nolint

    drug <- faers_get(object, "drug")
    reac <- faers_get(object, "reac")

    # ---- normalize + deduplicate both sides (identical to the scan) ----
    drug_u <- drug[, .(primaryid, drug = tolower(trimws(get(drug_col))))]
    drug_u <- drug_u[!is.na(drug) & !drug %chin% SCAN_UNKNOWN_DRUGS]
    drug_u <- unique(drug_u, by = c("primaryid", "drug"), cols = character())
    reac_u <- reac[, .(primaryid, event = get(event_col))]
    reac_u <- reac_u[!is.na(event) & event != ""]
    reac_u <- unique(reac_u, by = c("primaryid", "event"), cols = character())

    if (!nrow(drug_u) || !nrow(reac_u) || !nrow(strata_dt)) {
        return(scan_strata_empty_counts(strata_cols))
    }

    # ---- restrict both sides to reports with a stratum, add the key ----
    drug_u <- drug_u[strata_dt, on = "primaryid", nomatch = NULL]
    reac_u <- reac_u[strata_dt, on = "primaryid", nomatch = NULL]

    drug_u[, drugid := .GRP, by = drug] # nolint
    reac_u[, eventid := .GRP, by = event] # nolint
    drug_ids <- unique(drug_u[, .(drugid, drug)])
    event_ids <- unique(reac_u[, .(eventid, event)])
    d_min <- drug_u[, .(primaryid, drugid, strat_id)]
    r_min <- reac_u[, .(primaryid, eventid, strat_id)]
    data.table::setkey(d_min, primaryid)
    data.table::setkey(r_min, primaryid)

    # per-stratum totals: both sides are unique per report + value
    drug_tot <- drug_u[, .(n_drug = .N), by = .(drugid, strat_id)]
    event_tot <- reac_u[, .(n_event = .N), by = .(eventid, strat_id)]

    # chunk over reports; reports are atomic so chunk counts add up exactly
    ids <- unique(d_min$primaryid) # sorted: the table is keyed
    starts <- seq.int(1L, length(ids), by = chunk_size)
    ends <- pmin(starts + chunk_size - 1L, length(ids))
    n_chunks <- length(starts)
    bar_id <- cli::cli_progress_bar(
        "Scanning stratum x drug x event combinations",
        type = "iterator", total = n_chunks,
        format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}",
        format_done = "Scanned {.val {n_chunks}} chunk{?s} of stratified drug-event pairs in {cli::pb_elapsed}",
        clear = FALSE
    )
    parts <- lapply(seq_len(n_chunks), function(k) {
        block <- ids[starts[k]:ends[k]]
        d_b <- d_min[J(block), on = "primaryid", nomatch = NULL]
        r_b <- r_min[J(block), on = "primaryid", nomatch = NULL]
        tri <- d_b[r_b,
            on = "primaryid", nomatch = NULL, allow.cartesian = TRUE
        ]
        cli::cli_progress_update(id = bar_id)
        tri[, .(a = .N), by = .(strat_id, drugid, eventid)]
    })
    counts <- data.table::rbindlist(parts, use.names = TRUE)
    counts <- counts[, .(a = sum(a)), by = .(strat_id, drugid, eventid)]
    counts <- counts[a >= min_a]

    # map integer ids back to names and attach per-stratum totals
    counts <- strat_map[counts, on = "strat_id"] # strat columns first
    counts[drug_ids, on = "drugid", drug := i.drug] # nolint
    counts[event_ids, on = "eventid", event := i.event] # nolint
    counts[drug_tot, on = c("drugid", "strat_id"), n_drug := i.n_drug] # nolint
    counts[event_tot,
        on = .(eventid, strat_id),
        n_event := i.n_event # nolint
    ]
    counts[strat_tot, on = "strat_id", n_stratum := i.n_stratum] # nolint
    counts[, c("drugid", "eventid", "strat_id") := NULL]
    data.table::setcolorder(
        counts,
        c(strata_cols, "drug", "event", "a", "n_drug", "n_event", "n_stratum")
    )
    counts[]
}

scan_strata_empty_counts <- function(strata_cols) {
    out <- data.table::data.table(
        drug = character(), event = character(), a = integer(),
        n_drug = integer(), n_event = integer(), n_stratum = integer()
    )
    for (col in strata_cols) data.table::set(out, j = col, value = character())
    data.table::setcolorder(out, c(strata_cols, names(out)))
    out
}

utils::globalVariables(c(
    "a", "b", "c", "d", "drug", "event", "drugid", "eventid", "strat_id",
    "n_drug", "n_event", "n_stratum", "n", "primaryid", "J",
    "i.drug", "i.event", "i.n_drug", "i.n_event", "i.n_stratum"
))
