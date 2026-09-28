# Tests for faers_phv_composite / faers_phv_signal_composite (signal-eventset.R).
#
# The bundled standardized sample (200 reports) is not de-duplicated; the
# helpers below flip the `deduplication` flag because composite analysis
# requires de-duplicated data.  The duckdb parity tests reuse the twin-building
# pattern of test_duckdb.R / test_scan.R and are skipped when duckdb is not
# installed.

has_duckdb <- requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("DBI", quietly = TRUE)

composite_sample_object <- function() {
    rds <- readRDS(system.file("extdata", "standardized_data.rds", package = "faers"))
    sl <- list()
    for (sn in methods::slotNames("FAERSascii")) {
        sl[[sn]] <- tryCatch(methods::slot(rds, sn), error = function(e) NULL)
    }
    obj <- do.call(methods::new, c(list(Class = "FAERSascii"), sl))
    obj@deduplication <- TRUE
    obj
}

composite_db_twin <- function(mem_obj) {
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

composite_drug_subset <- function(object, pattern) {
    ids <- unique(faers_get(object, "drug")[
        grepl(pattern, drugname, ignore.case = TRUE), primaryid
    ])
    faers_keep(object, primaryid = ids)
}

composite_expected_cells <- function(object, full, event_set, event_type) {
    reac <- faers_get(full, "reac")
    if (is.character(event_set)) {
        event_ids <- unique(reac[get(event_type) %in% event_set, primaryid])
    } else {
        event_ids <- unique(reac[event_set(reac[[event_type]]), primaryid])
    }
    obj_ids <- faers_primaryid(object)
    a <- length(intersect(obj_ids, event_ids))
    data.table::data.table(
        a = a,
        b = length(obj_ids) - a,
        c = length(event_ids) - a,
        d = length(unique(faers_primaryid(full))) -
            (length(obj_ids) + length(event_ids) - a)
    )
}

testthat::test_that("faers_phv_composite builds a consistent 2x2 for pt event sets", {
    full <- composite_sample_object()
    sub <- composite_drug_subset(full, "aspirin")
    out <- faers_phv_composite(sub, c("NAUSEA", "VOMITING"), .event_type = "pt", .full = full)

    expected <- composite_expected_cells(sub, full, c("NAUSEA", "VOMITING"), "pt")
    testthat::expect_identical(out$event, "Composite_Event_Set")
    testthat::expect_identical(out$a, expected$a)
    testthat::expect_identical(out$b, expected$b)
    testthat::expect_identical(out$c, expected$c)
    testthat::expect_identical(out$d, expected$d)
    # conservation: a + b + c + d must equal the full-database report count
    testthat::expect_identical(
        out$a + out$b + out$c + out$d,
        length(unique(faers_primaryid(full)))
    )
})

testthat::test_that("faers_phv_composite supports MedDRA hierarchy event types (regression)", {
    # Regression: .identify_event_set_patients() read the raw `reac` table,
    # which lacks hierarchy columns, so soc_name aborted with
    # "object 'soc_name' not found".
    full <- composite_sample_object()
    sub <- composite_drug_subset(full, "aspirin")

    out_soc <- faers_phv_composite(
        sub, "Psychiatric disorders", .event_type = "soc_name", .full = full
    )
    expected <- composite_expected_cells(sub, full, "Psychiatric disorders", "soc_name")
    testthat::expect_identical(out_soc$a, expected$a)
    testthat::expect_identical(out_soc$b, expected$b)
    testthat::expect_identical(out_soc$c, expected$c)
    testthat::expect_identical(out_soc$d, expected$d)

    reac <- faers_get(full, "reac")
    hlgt_set <- na.omit(unique(reac$hlgt_name))[1]
    out_hlgt <- faers_phv_composite(sub, hlgt_set, .event_type = "hlgt_name", .full = full)
    expected_hlgt <- composite_expected_cells(sub, full, hlgt_set, "hlgt_name")
    testthat::expect_identical(out_hlgt$a, expected_hlgt$a)
    testthat::expect_identical(out_hlgt$b, expected_hlgt$b)
    testthat::expect_identical(out_hlgt$c, expected_hlgt$c)
    testthat::expect_identical(out_hlgt$d, expected_hlgt$d)
})

testthat::test_that("faers_phv_composite accepts function event sets", {
    full <- composite_sample_object()
    sub <- composite_drug_subset(full, "aspirin")
    fn <- function(x) grepl("NAUSEA|VOMITING", x)
    out <- faers_phv_composite(sub, fn, .event_type = "pt", .full = full)
    expected <- composite_expected_cells(sub, full, fn, "pt")
    testthat::expect_identical(out$a, expected$a)
    testthat::expect_identical(out$d, expected$d)
})

testthat::test_that("faers_phv_composite rejects unsupported event types", {
    full <- composite_sample_object()
    testthat::expect_error(
        faers_phv_composite(full, "NAUSEA", .event_type = "not_a_column", .full = full),
        "not_a_column"
    )
})

testthat::test_that("faers_phv_composite requires de-duplicated data", {
    full <- composite_sample_object()
    full@deduplication <- FALSE
    testthat::expect_error(
        faers_phv_composite(full, "NAUSEA", .event_type = "pt", .full = full),
        "de-duplicated"
    )
})

testthat::test_that("faers_phv_signal_composite appends phv_signal columns", {
    full <- composite_sample_object()
    sub <- composite_drug_subset(full, "aspirin")
    out <- faers_phv_signal_composite(
        sub, c("HEADACHE", "DIZZINESS"),
        .event_type = "pt", .methods = c("ror", "prr"), .full = full
    )
    testthat::expect_identical(nrow(out), 1L)
    testthat::expect_true(all(c("ror", "ror_ci_low", "prr") %in% names(out)))
})

testthat::test_that("faers_phv_composite compares two drug subsets (.object2 mode)", {
    full <- composite_sample_object()
    ids_a <- unique(faers_get(full, "drug")[
        grepl("aspirin", drugname, ignore.case = TRUE), primaryid
    ])
    ids_b <- setdiff(unique(faers_get(full, "drug")[
        grepl("metformin", drugname, ignore.case = TRUE), primaryid
    ]), ids_a)
    sub_a <- faers_keep(full, primaryid = ids_a)
    sub_b <- faers_keep(full, primaryid = ids_b)

    out <- suppressWarnings(faers_phv_composite(
        sub_a, "NAUSEA", .event_type = "pt", .object2 = sub_b
    ))
    reac <- faers_get(full, "reac")
    event_ids <- unique(reac[pt == "NAUSEA", primaryid])
    testthat::expect_identical(out$a, length(intersect(ids_a, event_ids)))
    testthat::expect_identical(out$b, length(ids_a) - length(intersect(ids_a, event_ids)))
    testthat::expect_identical(out$c, length(intersect(ids_b, event_ids)))
    testthat::expect_identical(out$d, length(ids_b) - length(intersect(ids_b, event_ids)))
})

testthat::test_that("faers_phv_composite is identical across memory and duckdb backends", {
    testthat::skip_if_not(has_duckdb)
    mem <- composite_sample_object()
    sub_mem <- composite_drug_subset(mem, "aspirin")
    db <- composite_db_twin(mem)
    sub_db <- composite_db_twin(sub_mem)

    for (event_type in c("pt", "soc_name")) {
        out_mem <- faers_phv_composite(
            sub_mem, "Psychiatric disorders", .event_type = event_type, .full = mem
        )
        out_db <- faers_phv_composite(
            sub_db, "Psychiatric disorders", .event_type = event_type, .full = db
        )
        testthat::expect_identical(out_mem, out_db)
    }
})
