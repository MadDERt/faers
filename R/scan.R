#' Whole-database pharmacovigilance signal scanning
#'
#' @description `faers_phv_scan()` enumerates all `drug x event` pairs present
#' in a standardized, de-duplicated [FAERSascii] object, builds the 2x2
#' contingency table for every pair in one aggregated pass, and runs
#' disproportionality analysis ([phv_signal]) on the whole table. In contrast
#' to [faers_phv_signal], which answers "*which* events are associated with
#' *this* drug?", the scan answers "*which* drug-event pairs stand out across
#' the whole database?", enabling hypothesis-generating pharmacovigilance.
#'
#' @details
#' The background distribution is the supplied object itself: for every pair
#' `(drug, event)`, `a` is the number of distinct patients reported with both,
#' `n1.` the patients exposed to the drug, `n.1` the patients reporting the
#' event, and `n` the total number of distinct patients in `demo`.
#'
#' Before counting, the drug side is normalized (`lower(trim(.))`) and empty,
#' `unknown` or `?` names are dropped; on the event side `NA` and empty values
#' are dropped. Drug names are aggregated on their exact normalized value, so
#' brand-name variants are not merged (a mapping step may be added in the
#' future).
#'
#' Pairs below `.min_a` are never evaluated, which is what keeps the output
#' and the expensive methods tractable. The default `.methods = c("ror",
#' "prr")` are cheap vectorized methods suitable for scanning; `ebgm` and
#' `bcpnn_mcmc` are computationally expensive on large candidate tables and
#' trigger a warning when many pairs are supplied.
#'
#' For whole-database scans (many years of FAERS quarterly data) the
#' `database = "duckdb"` backend of [faers()] is strongly recommended: the
#' drug x event counting runs as a single out-of-core SQL aggregation. The
#' memory backend performs the equivalent computation with integer coding and
#' chunking over `primaryid` blocks (`.chunk_size`), which is exact because a
#' report never spans two chunks, so per-chunk distinct counts sum to the
#' database-wide counts.
#'
#' @param .object A [FAERSascii] object, standardized with [faers_standardize]
#' and de-duplicated with [faers_dedup].
#' @param .events A string, the event column of the standardized `reac` field.
#' Any MedDRA hierarchy column (`pt_name`, `hlt_name`, `hlgt_name`,
#' `soc_name`, ...) added by [faers_standardize] can be used, as well as the
#' raw `pt` column. Defaults to `"pt"`.
#' @param .drug_field A string, the column of the standardized `drug` field
#' used to define drugs: `"drugname"` (default) or `"prod_ai"` (active
#' ingredient).
#' @param .drug_pattern An optional regular expression. If supplied, only drug
#' names matching it (case-insensitively) are kept in the result, enabling
#' scans of drug classes without paying for the whole database.
#' @param .min_a A single integer, the minimum number of reports co-listing a
#' drug-event pair for it to be kept. Defaults to `3L`.
#' @param .methods An atomic character, the disproportionality methods passed
#' to [phv_signal]. Defaults to `c("ror", "prr")`.
#' @param .chunk_size A single integer, the number of unique `primaryid`s
#' processed per chunk by the memory backend. Only relevant for
#' `database = "memory"`; smaller values trade speed for peak memory.
#' @param .phv_signal_params Other arguments passed to [phv_signal].
#' @param BPPARAM A [BiocParallel::BiocParallelParam-class] object.
#' @param ... Unused arguments, included for S4 generic/method consistency.
#' @return A [data.table][data.table::data.table] with one row per kept
#' drug-event pair, sorted by the lower bound of the first method's confidence
#' interval in descending order: the drug column (named after `.drug_field`),
#' the event column (named after `.events`), the contingency table columns
#' `a`, `b`, `c`, `d`, and the columns of [phv_signal].
#' @examples
#' # the sample data below is standardized but not de-duplicated; real usage
#' # requires faers_standardize() + faers_dedup() before scanning
#' std_data <- readRDS(system.file("extdata", "standardized_data.rds",
#'     package = "faers"
#' ))
#' std_data@deduplication <- TRUE
#' \dontrun{
#' # scan the bundled sample at the PT level
#' res <- faers_phv_scan(std_data, .events = "pt", .min_a = 1L)
#' head(res)
#'
#' # for real analyses, run on many quarters with the duckdb backend
#' data <- faers(2015:2024, paste0("q", 1:4),
#'     dir = "faers_data", database = "duckdb"
#' )
#' data <- faers_dedup(faers_standardize(data, meddra_path))
#' res <- faers_phv_scan(data, .events = "pt", .min_a = 3L)
#' }
#' @seealso [phv_signal], [faers_phv_signal]
#' @export
#' @aliases faers_phv_scan
#' @name faers_phv_scan
methods::setGeneric("faers_phv_scan", function(.object, ...) {
    standardGeneric("faers_phv_scan")
})

#' @rdname faers_phv_scan
#' @export
#' @method faers_phv_scan FAERSascii
methods::setMethod("faers_phv_scan", "FAERSascii", function(
    .object, .events = "pt", .drug_field = "drugname", .drug_pattern = NULL,
    .min_a = 3L, .methods = c("ror", "prr"), .chunk_size = 1e6L,
    .phv_signal_params = list(), BPPARAM = BiocParallel::SerialParam()
) {
    assert_string(.drug_field, allow_empty = FALSE)
    assert_string(.events, allow_empty = FALSE)
    assert_string(.drug_pattern, allow_null = TRUE, allow_empty = FALSE)
    assert_number_whole(.min_a, min = 1)
    assert_number_whole(.chunk_size, min = 1)
    assert_(.phv_signal_params, is.list, "a list")
    if (!.object@standardization) {
        cli::cli_abort("{.arg .object} must be standardized using {.fn faers_standardize}")
    }
    if (!.object@deduplication) {
        cli::cli_abort(c(
            "{.arg .object} must be de-duplicated using {.fn faers_dedup}",
            i = "Scanning raw duplicates would heavily bias the drug-event counts"
        ))
    }

    # ---- resolve the drug column ----
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

    # ---- resolve the event column ----
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

    # ---- count all drug x event pairs ----
    if (!is.null(.object@db)) {
        counts <- scan_counts_db(
            .object, drug_col = .drug_field, event_col = .events, min_a = .min_a
        )
    } else {
        counts <- scan_counts_memory(
            .object, drug_col = .drug_field, event_col = .events,
            min_a = .min_a, chunk_size = .chunk_size
        )
    }

    # ---- build one 2x2 table per pair ----
    counts <- scan_build_contingency(counts)

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

    # ---- finalize ----
    data.table::setnames(counts, c("drug", "event"), c(.drug_field, .events))
    ci_col <- grep("_ci_low$", names(counts), value = TRUE)[1L]
    if (!is.na(ci_col)) {
        data.table::setorderv(
            counts, c(ci_col, .drug_field, .events),
            order = c(-1L, 1L, 1L), na.last = TRUE
        )
    } else {
        data.table::setorderv(
            counts, c("a", .drug_field, .events),
            order = c(-1L, 1L, 1L), na.last = TRUE
        )
    }
    counts[]
})

# Methods that are allowed in phv_signal(); duplicated here so the scan can
# warn about expensive ones before running anything.
SCAN_ALL_METHODS <- c(
    "ror", "prr", "chisq", "bcpnn_norm", "bcpnn_mcmc",
    "obsexp_shrink", "fisher", "ebgm"
)
# Methods whose cost scales badly with the number of candidate pairs.
SCAN_EXPENSIVE_METHODS <- c("ebgm", "bcpnn_mcmc")
# Above this many candidate pairs, expensive methods warn.
SCAN_EXPENSIVE_THRESHOLD <- 1e5L
# Drug names dropped before counting (keep in sync with
# DB_SCAN_UNKNOWN_DRUGS in scan-db.R).
SCAN_UNKNOWN_DRUGS <- c("", "unknown", "?")

#' Shared assembly of the per-pair 2x2 contingency table.
#' `counts` must provide `drug`, `event`, `a`, `n_drug`, `n_event`, `n`
#' (total distinct patients).
#' @noRd
scan_build_contingency <- function(counts) {
    counts[, b := n_drug - a] # nolint
    counts[, c := n_event - a] # nolint
    counts[, d := n - (n_drug + n_event - a)] # nolint
    counts[, c("n_drug", "n_event", "n") := NULL]
    data.table::setcolorder(counts, c("drug", "event", "a", "b", "c", "d"))
    counts[]
}

# Memory-backend counting: normalize + deduplicate both sides, integer-code
# them, then count drug x event pairs chunk by chunk.
#
# Chunking is exact because every report's rows share the same `primaryid`, so
# a report never spans two chunks: summing per-chunk COUNT(DISTINCT primaryid)
# yields the database-wide distinct counts.
scan_counts_memory <- function(object, drug_col, event_col, min_a, chunk_size) {
    drug <- faers_get(object, "drug")
    reac <- faers_get(object, "reac")

    # ---- normalize + deduplicate both sides ----
    drug_u <- drug[, .(primaryid, drug = tolower(trimws(get(drug_col))))]
    drug_u <- drug_u[!is.na(drug) & !drug %chin% SCAN_UNKNOWN_DRUGS]
    drug_u <- unique(drug_u, by = c("primaryid", "drug"), cols = character())
    reac_u <- reac[, .(primaryid, event = get(event_col))]
    reac_u <- reac_u[!is.na(event) & event != ""]
    reac_u <- unique(reac_u, by = c("primaryid", "event"), cols = character())

    if (!nrow(drug_u) || !nrow(reac_u)) {
        return(scan_empty_counts())
    }

    n_total <- length(unique(faers_primaryid(object)))

    # ---- integer-code to shrink the join ----
    drug_u[, drugid := .GRP, by = drug] # nolint
    reac_u[, eventid := .GRP, by = event] # nolint
    drug_ids <- unique(drug_u[, .(drugid, drug)])
    event_ids <- unique(reac_u[, .(eventid, event)])
    d_min <- drug_u[, .(primaryid, drugid)]
    r_min <- reac_u[, .(primaryid, eventid)]
    data.table::setkey(d_min, primaryid)
    data.table::setkey(r_min, primaryid)

    # totals per drug / per event (both sides are already unique per pair)
    drug_tot <- drug_u[, .(n_drug = .N), by = .(drugid, drug)]
    event_tot <- reac_u[, .(n_event = .N), by = .(eventid, event)]

    # chunk over reports; reports are atomic so chunk counts add up exactly
    ids <- unique(d_min$primaryid) # sorted: the table is keyed
    starts <- seq.int(1L, length(ids), by = chunk_size)
    ends <- pmin(starts + chunk_size - 1L, length(ids))
    n_chunks <- length(starts)
    bar_id <- cli::cli_progress_bar(
        "Scanning drug x event pairs",
        type = "iterator", total = n_chunks,
        format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}",
        format_done = "Scanned {.val {n_chunks}} chunk{?s} of drug-event pairs in {cli::pb_elapsed}",
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
        tri[, .(a = .N), by = .(drugid, eventid)]
    })
    counts <- data.table::rbindlist(parts, use.names = TRUE)
    counts <- counts[, .(a = sum(a)), by = .(drugid, eventid)]
    counts <- counts[a >= min_a]

    # map integer ids back to names and attach totals
    counts[drug_ids, on = "drugid", drug := i.drug] # nolint
    counts[event_tot,
        on = "eventid",
        `:=`(event = i.event, n_event = i.n_event) # nolint
    ]
    counts[drug_tot, on = "drugid", n_drug := i.n_drug] # nolint
    counts[, c("drugid", "eventid") := NULL]
    counts[, n := n_total]
    data.table::setcolorder(
        counts, c("drug", "event", "a", "n_drug", "n_event", "n")
    )
    counts[]
}

scan_empty_counts <- function() {
    data.table::data.table(
        drug = character(), event = character(), a = integer(),
        n_drug = integer(), n_event = integer(), n = integer()
    )
}

utils::globalVariables(c(
    "a", "b", "c", "d", "drug", "event", "drugid", "eventid",
    "n_drug", "n_event", "n", "primaryid",
    "J", "i.drug", "i.event", "i.n_drug", "i.n_event"
))
