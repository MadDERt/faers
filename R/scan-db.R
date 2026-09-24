# Database-mode faers_phv_scan: push the whole drug x event pair counting
# down to a single DuckDB query.
#
# For every (drug, event) pair the scan needs four aggregates:
#   a    = COUNT(DISTINCT primaryid) among reports containing BOTH
#   n1.  = distinct patients exposed to the drug (any report)
#   n.1  = distinct patients reporting the event
#   n    = distinct patients in demo (the full background)
# Computing them in one SQL statement keeps only the aggregated pair table
# (already filtered by `.min_a`, i.e. small) in memory, so the scan scales to
# the whole FAERS database on the duckdb backend.

# SQL string literal for the drug names that are dropped before counting.
DB_SCAN_UNKNOWN_DRUGS <- c("", "unknown", "?")

# SQL fragment filtering garbage drug names (must stay in sync with
# SCAN_UNKNOWN_DRUGS in scan.R, which the memory path applies in R).
db_scan_drug_filter_sql <- function(drug_col_q) {
    unknowns <- paste(
        vapply(DB_SCAN_UNKNOWN_DRUGS, function(x)
            sprintf("'%s'", gsub("'", "''", x)), character(1L)
        ),
        collapse = ", "
    )
    sprintf(
        "lower(trim(%s)) NOT IN (%s)", drug_col_q, unknowns
    )
}

# Build the full SQL string for the scan aggregation.
db_scan_query <- function(drug_col, event_col, meddra, min_a) {
    drug_tbl <- db_field_table("drug")
    reac_tbl <- db_field_table("reac")
    demo_tbl <- db_field_table("demo")
    drug_col_q <- db_qident(drug_col)

    drug_sel <- sprintf(
        paste(
            "SELECT DISTINCT primaryid, lower(trim(%s)) AS drugname",
            "FROM %s",
            "WHERE %s IS NOT NULL AND %s"
        ),
        drug_col_q, drug_tbl, drug_col_q, db_scan_drug_filter_sql(drug_col_q)
    )
    if (isTRUE(meddra)) {
        reac_sel <- sprintf(
            paste(
                "SELECT DISTINCT f.primaryid AS primaryid, m.%s AS event",
                "FROM %s AS f",
                "LEFT JOIN %s AS m ON f.meddra_hierarchy_idx = m.idx",
                "WHERE m.%s IS NOT NULL AND m.%s <> ''"
            ),
            db_qident(event_col), reac_tbl, DB_MEDDRA_TABLE,
            db_qident(event_col), db_qident(event_col)
        )
    } else {
        reac_sel <- sprintf(
            paste(
                "SELECT DISTINCT primaryid, %s AS event",
                "FROM %s",
                "WHERE %s IS NOT NULL AND %s <> ''"
            ),
            db_qident(event_col), reac_tbl, db_qident(event_col), db_qident(event_col)
        )
    }
    sprintf("
WITH drug_u AS (%s),
     reac_u AS (%s),
     pair AS (
         SELECT g.drugname AS drugname, r.event AS event,
                COUNT(DISTINCT g.primaryid) AS a
         FROM drug_u g
         INNER JOIN reac_u r ON g.primaryid = r.primaryid
         GROUP BY g.drugname, r.event
         HAVING COUNT(DISTINCT g.primaryid) >= %d
     ),
     drug_tot AS (
         SELECT drugname, COUNT(DISTINCT primaryid) AS n_drug
         FROM drug_u GROUP BY drugname
     ),
     event_tot AS (
         SELECT event, COUNT(DISTINCT primaryid) AS n_event
         FROM reac_u GROUP BY event
     )
SELECT p.drugname, p.event, p.a, d.n_drug, e.n_event, x.n
FROM pair p
INNER JOIN drug_tot d ON p.drugname = d.drugname
INNER JOIN event_tot e ON p.event = e.event
CROSS JOIN (SELECT COUNT(DISTINCT primaryid) AS n FROM %s) AS x
",
        drug_sel, reac_sel, as.integer(min_a), demo_tbl
    )
}

# Run the scan aggregation on a db-backed object and return a data.table with
# columns: drug, event, a, n_drug, n_event, n.
scan_counts_db <- function(object, drug_col, event_col, min_a) {
    con <- object@db@con
    meddra <- event_col %chin% names(object@meddra@hierarchy)
    if (meddra && !DB_MEDDRA_TABLE %in% DBI::dbListTables(con)) {
        db_ingest_meddra(con, object)
    }
    out <- DBI::dbGetQuery(
        con, db_scan_query(drug_col, event_col, meddra, min_a)
    )
    out <- data.table::as.data.table(out)
    # COUNT(DISTINCT ...) comes back as BIGINT/double -> coerce to integer,
    # keeping byte-parity with the memory path.
    int_cols <- intersect(c("a", "n_drug", "n_event", "n"), names(out))
    for (col in int_cols) {
        if (is.double(out[[col]])) {
            data.table::set(out, j = col, value = as.integer(out[[col]]))
        }
    }
    data.table::setnames(out, "drugname", "drug")
    data.table::setcolorder(
        out, c("drug", "event", "a", "n_drug", "n_event", "n")
    )
    out
}
