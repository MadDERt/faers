# Tests for faers_phv_scan_stratified (stratified drug x event scanning).
#
# Mirrors test_scan.R: the bundled standardized sample (200 reports) is not
# de-duplicated, so helpers flip the `deduplication` flag; duckdb parity tests
# build an in-memory twin database and are skipped when duckdb is missing.

has_duckdb <- requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("DBI", quietly = TRUE)

strata_sample_object <- function() {
    rds <- readRDS(system.file("extdata", "standardized_data.rds", package = "faers"))
    sl <- list()
    for (sn in methods::slotNames("FAERSascii")) {
        sl[[sn]] <- tryCatch(methods::slot(rds, sn), error = function(e) NULL)
    }
    obj <- do.call(methods::new, c(list(Class = "FAERSascii"), sl))
    obj@deduplication <- TRUE
    obj
}

strata_db_twin <- function(mem_obj) {
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

# A tiny hand-crafted standardized object with known strata contents:
#   reports: 1, 2, 3, 4
#   sex:     1 -> Female, 2 -> Female, 3 -> NA ("Missing"), 4 -> Male
#   drugs:   1 -> ASPIRIN, 1 -> UNKNOWN (dropped), 2 -> "aspirin " (normalized,
#            merges with report 1), 3 -> IBUPROFEN, 4 -> IBUPROFEN
#   events:  1 -> HEADACHE, 2 -> HEADACHE (twice, deduped), 3 -> NAUSEA,
#            4 -> NAUSEA
#   pairs (per stratum):
#     (Female, aspirin, HEADACHE)   a=2, n_stratum=2
#     (Male, ibuprofen, NAUSEA)     a=1, n_stratum=1
#     (Missing, ibuprofen, NAUSEA)  a=1, n_stratum=1
strata_synthetic_object <- function() {
    obj <- strata_sample_object()
    obj@data$demo <- data.table::data.table(
        year = 2004L, quarter = "q1", primaryid = as.character(1:4),
        sex = c("Female", "Female", NA_character_, "Male")
    )
    obj@data$drug <- data.table::data.table(
        year = 2004L, quarter = "q1",
        primaryid = c("1", "1", "2", "3", "4"),
        drugname = c("ASPIRIN", "UNKNOWN", "aspirin ", "IBUPROFEN", "IBUPROFEN")
    )
    obj@data$reac <- data.table::data.table(
        year = 2004L, quarter = "q1",
        primaryid = c("1", "2", "2", "3", "4"),
        pt = c(
            "HEADACHE", "HEADACHE", "HEADACHE", "NAUSEA", "NAUSEA"
        ),
        meddra_hierarchy_idx = 1L
    )
    obj
}

testthat::test_that("faers_phv_scan_stratified requires standardized and deduplicated data", {
    raw <- faers(c(2004, 2017), c("q1", "q2"), "ascii",
        dir = system.file("extdata", package = "faers"),
        compress_dir = tempdir()
    )
    testthat::expect_error(faers_phv_scan_stratified(raw), "standardized")

    std <- readRDS(system.file(
        "extdata", "standardized_data.rds", package = "faers"
    ))
    testthat::expect_error(faers_phv_scan_stratified(std), "dedup")
})

testthat::test_that("faers_phv_scan_stratified validates its arguments", {
    obj <- strata_sample_object()
    testthat::expect_error(
        faers_phv_scan_stratified(obj, .strata = "not_a_column"), "demo"
    )
    testthat::expect_error(
        faers_phv_scan_stratified(obj, .strata = "year"), "character"
    )
    testthat::expect_error(
        faers_phv_scan_stratified(obj, .strata = "a"), "reserved"
    )
    testthat::expect_error(
        faers_phv_scan_stratified(obj, .strata = c("sex", "sex")), "distinct"
    )
    testthat::expect_error(
        faers_phv_scan_stratified(obj, .strata = character()), "character vector"
    )
    testthat::expect_error(
        faers_phv_scan_stratified(obj, .events = "not_a_column"), "reac"
    )
    testthat::expect_error(
        faers_phv_scan_stratified(obj, .drug_field = "not_a_column"), "drug"
    )
    testthat::expect_error(
        faers_phv_scan_stratified(obj, .min_a = 0), "min_a"
    )
})

testthat::test_that("faers_phv_scan_stratified counts per-stratum 2x2 tables on hand-crafted data", {
    obj <- strata_synthetic_object()
    out <- suppressWarnings(faers_phv_scan_stratified(
        obj, .strata = "sex", .events = "pt", .min_a = 1L
    ))

    # one row per (stratum, drug, event); Missing stratum present
    testthat::expect_equal(nrow(out), 3L)
    testthat::expect_setequal(
        out$sex, c("Female", "Male", "Missing")
    )
    testthat::expect_setequal(
        out$drugname, c("aspirin", "ibuprofen")
    )

    aspirin <- out[sex == "Female"]
    testthat::expect_equal(aspirin$drugname, "aspirin")
    testthat::expect_equal(aspirin$pt, "HEADACHE")
    testthat::expect_equal(aspirin$a, 2L)
    testthat::expect_equal(aspirin$b, 0L) # n_drug(F, aspirin) = 2
    testthat::expect_equal(aspirin$c, 0L) # n_event(F, HEADACHE) = 2
    testthat::expect_equal(aspirin$d, 0L)
    testthat::expect_equal(aspirin$n_stratum, 2L)

    ibu_m <- out[sex == "Male"]
    testthat::expect_equal(ibu_m$drugname, "ibuprofen")
    testthat::expect_equal(ibu_m$a, 1L)
    testthat::expect_equal(ibu_m$b, 0L)
    testthat::expect_equal(ibu_m$c, 0L) # n_event(M, NAUSEA) = 1
    testthat::expect_equal(ibu_m$d, 0L)
    testthat::expect_equal(ibu_m$n_stratum, 1L)

    ibu_missing <- out[sex == "Missing"]
    testthat::expect_equal(ibu_missing$a, 1L)
    testthat::expect_equal(ibu_missing$n_stratum, 1L)

    # conservation: a + b + c + d == n_stratum everywhere
    testthat::expect_true(all(out$a + out$b + out$c + out$d == out$n_stratum))
})

testthat::test_that(".min_a applies per stratum, not to the whole database", {
    obj <- strata_synthetic_object()
    # aspirin-HEADACHE co-occurs in reports 1 (Female) and 2 (Male):
    # a = 2 in the whole database but only a = 1 inside each stratum
    obj@data$demo$sex <- c("Female", "Male", NA_character_, "Male")

    out <- suppressWarnings(faers_phv_scan_stratified(
        obj, .strata = "sex", .events = "pt", .min_a = 2L
    ))
    testthat::expect_equal(nrow(out), 0L)

    # the unstratified scan keeps the same pairs (global a = 2)
    out_scan <- suppressWarnings(faers_phv_scan(obj, .events = "pt", .min_a = 2L))
    testthat::expect_equal(nrow(out_scan), 2L)
})

testthat::test_that("missing strata values: keep vs drop", {
    obj <- strata_synthetic_object()

    out_keep <- suppressWarnings(faers_phv_scan_stratified(
        obj, .strata = "sex", .events = "pt", .min_a = 1L
    ))
    testthat::expect_true("Missing" %chin% out_keep$sex)
    testthat::expect_equal(sum(out_keep$n_stratum), 4L) # no report lost

    out_drop <- suppressWarnings(faers_phv_scan_stratified(
        obj, .strata = "sex", .events = "pt", .min_a = 1L,
        .na_stratum = "drop"
    ))
    testthat::expect_false("Missing" %chin% out_drop$sex)
    testthat::expect_equal(sum(out_drop$n_stratum), 3L) # report 3 removed
    # report 3 was the only NAUSEA-in-Missing report; its pair is gone
    testthat::expect_equal(nrow(out_drop), 2L)
})

testthat::test_that("multiple strata columns are crossed", {
    obj <- strata_synthetic_object()
    obj@data$demo$occp_cod <- c("CN", "OT", "CN", NA_character_)
    out <- suppressWarnings(faers_phv_scan_stratified(
        obj, .strata = c("sex", "occp_cod"), .events = "pt", .min_a = 1L
    ))
    testthat::expect_setequal(
        out[, paste(sex, occp_cod, sep = "|")],
        c(
            "Female|CN", # report 1
            "Female|OT", # report 2
            "Missing|CN", # report 3
            "Male|Missing" # report 4
        )
    )
    testthat::expect_true(all(out$n_stratum == 1L))
})

testthat::test_that("a heavily-missing strata column triggers a warning", {
    obj <- strata_synthetic_object()
    obj@data$demo$occr_country <- c("US", NA, NA, "DE")
    testthat::expect_warning(
        faers_phv_scan_stratified(obj, .strata = "occr_country", .min_a = 1L),
        "occr_country"
    )
    # sex is only 9% missing on the sample -> no missing-share warning
    testthat::expect_no_warning(
        faers_phv_scan_stratified(strata_sample_object(), .min_a = 1L),
        message = "missing"
    )
})

testthat::test_that("chunked memory path equals the unchunked one", {
    obj <- strata_sample_object()
    out <- suppressWarnings(faers_phv_scan_stratified(obj, .min_a = 1L))
    out_chunked <- suppressWarnings(faers_phv_scan_stratified(
        obj, .min_a = 1L, .chunk_size = 7L
    ))
    testthat::expect_identical(out, out_chunked)
})

testthat::test_that("strata counts equal directly-computed per-stratum 2x2 tables", {
    obj <- strata_sample_object()
    out <- suppressWarnings(faers_phv_scan_stratified(
        obj, .strata = "sex", .events = "pt", .min_a = 1L,
        .drug_pattern = "^humulin r$"
    ))

    demo <- faers_get(obj, "demo")
    drug_tbl <- faers_get(obj, "drug")
    reac_tbl <- faers_get(obj, "reac")
    d_u <- unique(drug_tbl[tolower(trimws(drugname)) == "humulin r",
        .(primaryid, drug = tolower(trimws(drugname)))
    ])
    r_u <- unique(reac_tbl[!is.na(pt) & pt != "", .(primaryid, event = pt)])
    strata_dt <- unique(demo[, .(primaryid, sex)])

    for (sv in setdiff(unique(strata_dt$sex), NA_character_)) {
        s_ids <- strata_dt[!is.na(sex) & sex == sv, primaryid]
        n_stratum <- length(s_ids)
        d_s <- d_u[primaryid %in% s_ids]
        r_s <- r_u[primaryid %in% s_ids]
        tri <- merge(d_s, r_s,
            by = "primaryid", allow.cartesian = TRUE
        )
        a_by_event <- tri[, .(a = .N), by = event]

        got <- out[sex == sv]
        testthat::expect_setequal(got$pt, a_by_event$event)
        for (i in seq_len(nrow(got))) {
            a <- a_by_event[event == got$pt[i], a]
            n_drug <- nrow(d_s)
            n_event <- sum(r_s$event == got$pt[i])
            testthat::expect_equal(got$a[i], a)
            testthat::expect_equal(got$a[i] + got$b[i], n_drug)
            testthat::expect_equal(got$a[i] + got$c[i], n_event)
            testthat::expect_equal(
                got$d[i], n_stratum - (n_drug + n_event - a)
            )
            testthat::expect_equal(got$n_stratum[i], n_stratum)
        }
    }
})

testthat::test_that("faers_phv_scan_stratified output is sorted by signal strength", {
    obj <- strata_sample_object()
    out <- suppressWarnings(faers_phv_scan_stratified(obj, .min_a = 1L))
    ci <- out$ror_ci_low
    ci <- ci[!is.na(ci)]
    testthat::expect_true(all(diff(ci) <= 0))
})

testthat::test_that("faers_phv_scan_stratified matches in db mode (raw + hierarchy events)", {
    testthat::skip_if_not(has_duckdb)
    mem <- strata_sample_object()
    db <- strata_db_twin(mem)
    for (ev in c("pt", "soc_name", "hlgt_name")) {
        testthat::expect_identical(
            suppressWarnings(faers_phv_scan_stratified(mem, .events = ev, .min_a = 1L)),
            suppressWarnings(faers_phv_scan_stratified(db, .events = ev, .min_a = 1L))
        )
    }
})

testthat::test_that("faers_phv_scan_stratified keep/drop agree in db mode", {
    testthat::skip_if_not(has_duckdb)
    mem <- strata_sample_object()
    db <- strata_db_twin(mem)
    for (na_mode in c("keep", "drop")) {
        testthat::expect_identical(
            suppressWarnings(faers_phv_scan_stratified(
                mem, .strata = "occp_cod", .min_a = 1L, .na_stratum = na_mode
            )),
            suppressWarnings(faers_phv_scan_stratified(
                db, .strata = "occp_cod", .min_a = 1L, .na_stratum = na_mode
            ))
        )
    }
})

testthat::test_that("db_scan_strata_query builds SQL with the right fragments", {
    q <- asNamespace("faers")$db_scan_strata_query(
        "sex", "drugname", "pt", meddra = FALSE, min_a = 3L, na_stratum = "keep"
    )
    testthat::expect_match(q, "faers_demo", fixed = TRUE)
    testthat::expect_match(q, "'Missing'", fixed = TRUE)
    testthat::expect_match(q, "n_stratum", fixed = TRUE)
    testthat::expect_match(q, ">= 3", fixed = TRUE)

    q_drop <- asNamespace("faers")$db_scan_strata_query(
        "sex", "drugname", "pt", meddra = FALSE, min_a = 3L, na_stratum = "drop"
    )
    testthat::expect_match(q_drop, "IS NOT NULL", fixed = TRUE)
})
