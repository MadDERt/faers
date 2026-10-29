# Tests for faers_phv_scan (whole-database drug x event signal scanning).
#
# The bundled standardized sample (200 reports) is not de-duplicated, so the
# helpers below flip the `deduplication` flag: scan requires both states.  The
# duckdb parity tests reuse the twin-building pattern of test_duckdb.R and are
# skipped when duckdb is not installed.

has_duckdb <- requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("DBI", quietly = TRUE)

scan_sample_object <- function() {
    rds <- readRDS(system.file("extdata", "standardized_data.rds", package = "faers"))
    sl <- list()
    for (sn in methods::slotNames("FAERSascii")) {
        sl[[sn]] <- tryCatch(methods::slot(rds, sn), error = function(e) NULL)
    }
    obj <- do.call(methods::new, c(list(Class = "FAERSascii"), sl))
    obj@deduplication <- TRUE
    obj
}

scan_db_twin <- function(mem_obj) {
    ns <- asNamespace("faers")
    con <- ns$db_connect(":memory:")
    tables <- ns$db_field_tables()
    for (f in names(mem_obj@data)) ns$db_ingest(con, f, mem_obj@data[[f]])
    ns$db_ingest_meddra(con, mem_obj)
    proxies <- stats::setNames(
        lapply(names(mem_obj@data), function(f) {
            methods::new("FAERSdbTbl", con = con, table = tables[[f]])
        }),
        names(mem_obj@data)
    )
    methods::new("FAERSascii",
        data = proxies,
        deletedCases = methods::slot(mem_obj, "deletedCases"),
        year = mem_obj@year, quarter = mem_obj@quarter,
        standardization = mem_obj@standardization,
        deduplication = mem_obj@deduplication,
        format = mem_obj@format,
        meddra = mem_obj@meddra,
        db = methods::new("FAERSdb", con = con, path = ":memory:",
            version = as.character(packageVersion("duckdb")),
            tables = tables)
    )
}

# A tiny hand-crafted standardized object with known drug/event contents:
#   reports: 1, 2, 3, 4
#   drugs:   1 -> ASPIRIN, 1 -> UNKNOWN (dropped), 2 -> "aspirin " (normalized,
#            merges with report 1), 3 -> "?" (dropped), 4 -> IBUPROFEN
#   events:  1 -> HEADACHE, 2 -> HEADACHE (twice, deduped), 3 -> NA (dropped),
#            4 -> NAUSEA
#   pairs:   (aspirin, HEADACHE) a=2;  (ibuprofen, NAUSEA) a=1
scan_synthetic_object <- function() {
    obj <- scan_sample_object()
    obj@data$demo <- data.table::data.table(
        year = 2004L, quarter = "q1", primaryid = as.character(1:4)
    )
    obj@data$drug <- data.table::data.table(
        year = 2004L, quarter = "q1",
        primaryid = c("1", "1", "2", "3", "4"),
        drugname = c("ASPIRIN", "UNKNOWN", "aspirin ", "?", "IBUPROFEN")
    )
    obj@data$reac <- data.table::data.table(
        year = 2004L, quarter = "q1",
        primaryid = c("1", "2", "2", "3", "4"),
        pt = c(
            "HEADACHE", "HEADACHE", "HEADACHE", NA_character_, "NAUSEA"
        ),
        meddra_hierarchy_idx = 1L
    )
    obj
}

testthat::test_that("faers_phv_scan requires standardized and deduplicated data", {
    raw <- faers(c(2004, 2017), c("q1", "q2"), "ascii",
        dir = system.file("extdata", package = "faers"),
        compress_dir = tempdir()
    )
    testthat::expect_error(faers_phv_scan(raw), "standardized")

    std <- readRDS(system.file(
        "extdata", "standardized_data.rds", package = "faers"
    ))
    testthat::expect_error(faers_phv_scan(std), "dedup")
})

testthat::test_that("faers_phv_scan validates its arguments", {
    obj <- scan_sample_object()
    testthat::expect_error(faers_phv_scan(obj, .events = "not_a_column"), "reac")
    testthat::expect_error(
        faers_phv_scan(obj, .drug_field = "not_a_column"), "drug"
    )
    testthat::expect_error(faers_phv_scan(obj, .min_a = 0), "min_a")
    testthat::expect_error(faers_phv_scan(obj, .phv_signal_params = "x"), "list")
})

testthat::test_that("faers_phv_scan counts pairs correctly on hand-crafted data", {
    obj <- scan_synthetic_object()
    out <- faers_phv_scan(obj, .events = "pt", .min_a = 1L)

    # garbage drug names dropped, normalization merged ASPIRIN + "aspirin "
    testthat::expect_setequal(out$drugname, c("aspirin", "ibuprofen"))
    testthat::expect_false(any(out$drugname %chin% c("", "unknown", "?")))

    aspirin <- out[drugname == "aspirin" & pt == "HEADACHE"]
    testthat::expect_equal(aspirin$a, 2L)
    testthat::expect_equal(aspirin$b, 0L) # n_drug = 2
    testthat::expect_equal(aspirin$c, 0L) # n_event(HEADACHE) = 2
    testthat::expect_equal(aspirin$d, 2L) # n = 4

    ibu <- out[drugname == "ibuprofen"]
    testthat::expect_equal(ibu$pt, "NAUSEA")
    testthat::expect_equal(ibu$a, 1L)
    testthat::expect_equal(ibu$b, 0L)
    testthat::expect_equal(ibu$c, 0L) # n_event(NAUSEA) = 1
    testthat::expect_equal(ibu$d, 3L) # 4 - (1 + 1 - 1)
})

testthat::test_that("faers_phv_scan honors .min_a and .drug_pattern", {
    obj <- scan_synthetic_object()

    # below-threshold pairs are dropped entirely
    out <- faers_phv_scan(obj, .events = "pt", .min_a = 2L)
    testthat::expect_equal(nrow(out), 1L)
    testthat::expect_equal(out$drugname, "aspirin")

    # the whitelist keeps only matching drugs
    out_class <- faers_phv_scan(
        obj, .events = "pt", .min_a = 1L, .drug_pattern = "^aspirin$"
    )
    testthat::expect_setequal(out_class$drugname, "aspirin")

    # a pattern matching nothing warns and returns an empty table
    testthat::expect_warning(
        out_empty <- faers_phv_scan(
            obj, .events = "pt", .min_a = 1L, .drug_pattern = "^nomatch$"
        ),
        "No drug name"
    )
    testthat::expect_equal(nrow(out_empty), 0L)
})

testthat::test_that("chunked memory path equals the unchunked one", {
    obj <- scan_sample_object()
    out <- faers_phv_scan(obj, .events = "pt", .min_a = 1L)
    out_chunked <- faers_phv_scan(
        obj, .events = "pt", .min_a = 1L, .chunk_size = 7L
    )
    testthat::expect_identical(out, out_chunked)
})

testthat::test_that("scan pair counts equal directly-computed co-occurrence", {
    obj <- scan_sample_object()
    scan <- faers_phv_scan(
        obj, .events = "pt", .min_a = 1L, .drug_pattern = "^humulin r$"
    )

    drug_tbl <- faers_get(obj, "drug")
    reac_tbl <- faers_get(obj, "reac")
    d_u <- unique(drug_tbl[tolower(trimws(drugname)) == "humulin r",
        .(primaryid, drug = tolower(trimws(drugname)))
    ])
    r_u <- unique(reac_tbl[!is.na(pt) & pt != "", .(primaryid, event = pt)])
    tri <- merge(d_u, r_u,
        by = "primaryid", all = FALSE, allow.cartesian = TRUE
    )
    expected <- tri[, .(a = .N), by = event]

    got <- scan[, .(event = pt, a)]
    testthat::expect_setequal(expected$event, got$event)
    for (i in seq_len(nrow(got))) {
        testthat::expect_equal(
            got$a[i],
            expected[event == got$event[i], a]
        )
    }
})

testthat::test_that("faers_phv_scan output is sorted by signal strength", {
    obj <- scan_sample_object()
    out <- faers_phv_scan(obj, .events = "pt", .min_a = 1L)
    ci <- out$ror_ci_low
    ci <- ci[!is.na(ci)]
    testthat::expect_true(all(diff(ci) <= 0))
})

testthat::test_that("faers_phv_scan matches in db mode (raw + hierarchy events)", {
    testthat::skip_if_not(has_duckdb)
    mem <- scan_sample_object()
    db <- scan_db_twin(mem)
    for (ev in c("pt", "soc_name", "hlgt_name")) {
        testthat::expect_identical(
            suppressWarnings(faers_phv_scan(mem, .events = ev, .min_a = 1L)),
            suppressWarnings(faers_phv_scan(db, .events = ev, .min_a = 1L))
        )
    }
})

testthat::test_that("db_scan_query builds SQL with the right filters", {
    q <- asNamespace("faers")$db_scan_query(
        "drugname", "pt", meddra = FALSE, min_a = 3L
    )
    testthat::expect_match(q, "lower(trim(\"drugname\"))", fixed = TRUE)
    testthat::expect_match(q, "NOT IN ('', 'unknown', '?')", fixed = TRUE)
    testthat::expect_match(q, ">= 3", fixed = TRUE)
    testthat::expect_match(q, "faers_demo", fixed = TRUE)

    q_meddra <- asNamespace("faers")$db_scan_query(
        "drugname", "soc_name", meddra = TRUE, min_a = 3L
    )
    testthat::expect_match(q_meddra, "faers_meddra", fixed = TRUE)
    testthat::expect_match(q_meddra, "meddra_hierarchy_idx", fixed = TRUE)
})

testthat::test_that("faers_phv_scan .role_cod restricts the exposure definition", {
    obj <- scan_sample_object()
    r_all <- suppressMessages(suppressWarnings(
        faers_phv_scan(obj, .min_a = 1L)
    ))
    r_ps <- suppressMessages(suppressWarnings(
        faers_phv_scan(obj, .min_a = 1L, .role_cod = "PS")
    ))
    testthat::expect_gt(nrow(r_all), nrow(r_ps))
    # every PS pair's a is a subset count of the all-roles pair
    m <- merge(
        r_all[, .(drugname, pt, a_all = a)],
        r_ps[, .(drugname, pt, a_ps = a)],
        by = c("drugname", "pt"), all = TRUE
    )
    m[is.na(a_ps), a_ps := 0L]
    testthat::expect_true(all(m$a_ps <= m$a_all))
    # background n is unchanged by the role filter
    # (reconstruct n from one row: n = a + b + c + d)
    if (nrow(r_all) && nrow(r_ps)) {
        n_all <- r_all[1, a + b + c + d]
        n_ps <- r_ps[1, a + b + c + d]
        testthat::expect_equal(n_all, n_ps)
    }
    # multiple roles accepted
    r_pss <- suppressMessages(suppressWarnings(
        faers_phv_scan(obj, .min_a = 1L, .role_cod = c("PS", "SS"))
    ))
    testthat::expect_gte(nrow(r_pss), nrow(r_ps))
    # invalid role code errors with the valid list
    testthat::expect_error(
        suppressMessages(faers_phv_scan(obj, .min_a = 1L, .role_cod = "XX")),
        "PS"
    )
    # NULL keeps everything (default; identical to previous versions)
    r_null <- suppressMessages(suppressWarnings(
        faers_phv_scan(obj, .min_a = 1L, .role_cod = NULL)
    ))
    testthat::expect_identical(r_null, r_all)
})

testthat::test_that("faers_phv_scan_stratified .role_cod filters within strata", {
    obj <- scan_sample_object()
    s_all <- suppressMessages(faers_phv_scan_stratified(
        obj, .strata = "sex", .min_a = 1L
    ))
    s_ps <- suppressMessages(suppressWarnings(
        faers_phv_scan_stratified(
            obj, .strata = "sex", .min_a = 1L, .role_cod = "PS"
        )
    ))
    testthat::expect_gte(nrow(s_all), nrow(s_ps))
    testthat::expect_equal(
        sort(unique(s_all$n_stratum)), sort(unique(s_ps$n_stratum))
    )
    m <- merge(
        s_all[, .(sex, drugname, pt, a_all = a)],
        s_ps[, .(sex, drugname, pt, a_ps = a)],
        by = c("sex", "drugname", "pt"), all = TRUE
    )
    m[is.na(a_ps), a_ps := 0L]
    testthat::expect_true(all(m$a_ps <= m$a_all))
})

testthat::test_that("db_scan_role_filter_sql emits correct SQL", {
    ns <- asNamespace("faers")
    testthat::expect_identical(ns$db_scan_role_filter_sql(NULL), "")
    q <- ns$db_scan_role_filter_sql(c("PS", "SS"))
    testthat::expect_match(q, "role_cod IN ('PS', 'SS')", fixed = TRUE)
    # full query embeds the fragment
    q2 <- asNamespace("faers")$db_scan_query(
        "drugname", "pt", meddra = FALSE, min_a = 3L, role_cod = "PS"
    )
    testthat::expect_match(q2, "role_cod IN ('PS')", fixed = TRUE)
})
