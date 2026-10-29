# Tests for faers_restrict (period / population restriction).
#
# Mirrors test_scan-strata.R: the bundled standardized sample is not
# de-duplicated, so helpers flip the `deduplication` flag; duckdb parity tests
# build an in-memory twin database and are skipped when duckdb is missing.

has_duckdb <- requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("DBI", quietly = TRUE)

restrict_sample_object <- function() {
    rds <- readRDS(system.file("extdata", "standardized_data.rds", package = "faers"))
    sl <- list()
    for (sn in methods::slotNames("FAERSascii")) {
        sl[[sn]] <- tryCatch(methods::slot(rds, sn), error = function(e) NULL)
    }
    sl <- Filter(Negate(is.null), sl)
    obj <- do.call(methods::new, c(list(Class = "FAERSascii"), sl))
    obj@deduplication <- TRUE
    obj
}

# in-memory duckdb twin builder (same pattern as test_scan-strata.R)
restrict_db_twin <- function(mem_obj) {
    ns <- asNamespace("faers")
    con <- ns$db_connect(":memory:")
    tables <- ns$db_field_tables()
    for (f in names(mem_obj@data)) ns$db_ingest(con, f, mem_obj@data[[f]])
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
        format = mem_obj@format, meddra = mem_obj@meddra,
        db = methods::new("FAERSdb", con = con, path = ":memory:",
            version = as.character(packageVersion("duckdb")),
            tables = tables
        )
    )
}

testthat::test_that("faers_restrict with no conditions warns and returns input unchanged", {
    obj <- restrict_sample_object()
    testthat::expect_warning(res <- faers_restrict(obj), "No restriction condition")
    testthat::expect_identical(res, obj)
})

testthat::test_that("faers_restrict filters by sex and reports drops in one aggregated warning", {
    obj <- restrict_sample_object()
    demo <- obj@data$demo
    expect_n <- sum(demo$sex == "Female", na.rm = TRUE)
    w <- NULL
    res <- withCallingHandlers(
        faers_restrict(obj, .sex = "Female"),
        warning = function(cnd) {
            w <<- c(w, conditionMessage(cnd))
            invokeRestart("muffleWarning")
        }
    )
    testthat::expect_equal(nrow(res@data$demo), expect_n)
    testthat::expect_length(w, 1L)
    testthat::expect_match(w[[1L]], "removed", fixed = TRUE)
    testthat::expect_match(w[[1L]], "sex", fixed = TRUE)
    testthat::expect_match(w[[1L]], "missing", fixed = TRUE)
    # combined intersection: sex Female AND age 0-120
    res2 <- suppressWarnings(faers_restrict(obj, .sex = "Female", .age = c(0, 120)))
    testthat::expect_lte(nrow(res2@data$demo), expect_n)
})

testthat::test_that("faers_restrict filters by age interval on age_in_years", {
    obj <- restrict_sample_object()
    demo <- obj@data$demo
    expect_n <- sum(
        !is.na(demo$age_in_years) &
            demo$age_in_years >= 18 & demo$age_in_years <= 85
    )
    res <- suppressWarnings(faers_restrict(obj, .age = c(18, 85)))
    testthat::expect_equal(nrow(res@data$demo), expect_n)
    testthat::expect_error(faers_restrict(obj, .age = c(85, 18)), "age")
    testthat::expect_error(faers_restrict(obj, .age = c(18)), "age")
    testthat::expect_error(faers_restrict(obj, .age = c("a", "b")), "age")
})

testthat::test_that("faers_restrict filters by period key", {
    obj <- restrict_sample_object()
    demo <- obj@data$demo
    res <- suppressWarnings(faers_restrict(obj, .period = "2004q1"))
    testthat::expect_equal(nrow(res@data$demo), sum(demo$year == 2004))
    # case-insensitive + whitespace tolerated
    res2 <- suppressWarnings(faers_restrict(obj, .period = " 2004Q1 "))
    testthat::expect_equal(
        sort(faers_primaryid(res)), sort(faers_primaryid(res2))
    )
    testthat::expect_error(faers_restrict(obj, .period = "2020_Q1"), "YYYYqN")
    testthat::expect_error(faers_restrict(obj, .period = "2020q9"), "1-4")
    testthat::expect_error(faers_restrict(obj, .period = 2004), "character")
})

testthat::test_that("faers_restrict applies date windows on the chosen field", {
    obj <- restrict_sample_object()
    demo <- obj@data$demo
    # fda_dt: numeric YYYYMMDD, complete in the sample
    res <- suppressWarnings(
        faers_restrict(obj, .date_range = c("2004-01-01", "2004-12-31"))
    )
    ok <- !is.na(demo$fda_dt)
    d <- as.Date(sprintf("%08.0f", demo$fda_dt[ok]), format = "%Y%m%d")
    testthat::expect_equal(
        nrow(res@data$demo),
        sum(ok & d >= as.Date("2004-01-01") & d <= as.Date("2004-12-31"))
    )
    # numeric YYYYMMDD form is accepted too
    res_num <- suppressWarnings(faers_restrict(
        obj, .date_range = c(20040101, 20040630)
    ))
    expect_num <- sum(
        ok & d >= as.Date("2004-01-01") & d <= as.Date("2004-06-30")
    )
    testthat::expect_equal(nrow(res_num@data$demo), expect_num)
    # event_dt: NA/partial rows (year-only 1999, year-month 201610, etc.)
    # cannot be placed in a window -> excluded by default, kept with
    # .date_partial = "keep". Window is wide; only parseable full dates pass.
    full_ok <- !is.na(demo$event_dt) & nchar(as.character(demo$event_dt)) == 8L
    res_e <- suppressWarnings(faers_restrict(
        obj, .date_range = c("1900-01-01", "2100-01-01"), .date_field = "event_dt"
    ))
    testthat::expect_equal(
        nrow(res_e@data$demo), sum(full_ok)
    )
    res_k <- suppressWarnings(faers_restrict(
        obj, .date_range = c("1900-01-01", "2100-01-01"),
        .date_field = "event_dt", .date_partial = "keep"
    ))
    testthat::expect_equal(nrow(res_k@data$demo), nrow(demo))
    # bad inputs
    testthat::expect_error(
        faers_restrict(obj, .date_range = c("2020-01-01")), "length-2"
    )
    testthat::expect_error(
        faers_restrict(obj, .date_range = c("not-a-date", "2020-01-01")),
        "parseable"
    )
    testthat::expect_error(
        faers_restrict(obj, .date_range = c("2020-12-31", "2020-01-01")),
        "start must not be after"
    )
})

testthat::test_that("faers_restrict combines conditions by intersection", {
    obj <- restrict_sample_object()
    demo <- obj@data$demo
    res <- suppressWarnings(faers_restrict(
        obj, .period = "2004q1", .sex = "Female"
    ))
    expected_ids <- demo[year == 2004 & sex == "Female", unique(primaryid)]
    testthat::expect_setequal(faers_primaryid(res), expected_ids)
    # every downstream field is filtered to the same reports
    testthat::expect_true(all(res@data$drug$primaryid %in% expected_ids))
})

testthat::test_that("faers_restrict keeps missing values of disabled conditions", {
    obj <- restrict_sample_object()
    # age has ~42% NA; a sex-only filter must keep all ages
    res <- suppressWarnings(faers_restrict(obj, .sex = "Female"))
    testthat::expect_true(any(is.na(res@data$demo$age_in_years)))
})

testthat::test_that("faers_restrict .invert removes matching reports", {
    obj <- restrict_sample_object()
    demo <- obj@data$demo
    res <- suppressWarnings(faers_restrict(obj, .sex = "Female", .invert = TRUE))
    testthat::expect_equal(
        nrow(res@data$demo),
        sum(demo$sex != "Female" | is.na(demo$sex))
    )
})

testthat::test_that("faers_restrict preserves flags and never mutates the input", {
    obj <- restrict_sample_object()
    snapshot <- lapply(obj@data, function(x) data.table::copy(x))
    res <- suppressWarnings(faers_restrict(obj, .sex = "Female"))
    testthat::expect_true(res@standardization)
    testthat::expect_true(res@deduplication == obj@deduplication)
    # row contents of the input are untouched (identical() may be FALSE due
    # to data.table's internal selfref pointer, so compare values)
    testthat::expect_true(isTRUE(all.equal(
        obj@data$demo, snapshot$demo, check.attributes = FALSE
    )))
    testthat::expect_true(isTRUE(all.equal(
        obj@data$drug, snapshot$drug, check.attributes = FALSE
    )))
})

testthat::test_that("faers_restrict validates value whitelists", {
    obj <- restrict_sample_object()
    testthat::expect_error(faers_restrict(obj, .sex = "femal"), "Female")
    testthat::expect_error(faers_restrict(obj, .age_grp = "X"), "age group|A|E")
    testthat::expect_error(faers_restrict(obj, .country = 123), "character")
    testthat::expect_error(faers_restrict(obj, .sex = character(0)), "non-empty")
})

testthat::test_that("faers_restrict produces identical id sets on memory and duckdb backends", {
    testthat::skip_if_not(has_duckdb)
    conds <- list(
        list(.sex = "Female"),
        list(.age = c(18, 85)),
        list(.sex = "Female", .age = c(18, 85)),
        list(.period = "2004q1"),
        list(.date_range = c("2004-01-01", "2004-12-31")),
        list(.country = "US")
    )
    for (i in seq_along(conds)) {
        # each condition gets a FRESH pair because the db backend filters by
        # physically deleting rows (faers_keep_db semantics)
        mem2 <- restrict_sample_object()
        db2 <- restrict_db_twin(mem2)
        mem_res <- suppressWarnings(do.call(
            faers_restrict, c(list(.object = mem2), conds[[i]])
        ))
        db_res <- suppressWarnings(do.call(
            faers_restrict, c(list(.object = db2), conds[[i]])
        ))
        testthat::expect_identical(
            sort(faers_primaryid(mem_res)),
            sort(faers_primaryid(db_res)),
            label = sprintf("condition %d", i)
        )
    }
})
