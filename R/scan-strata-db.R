# Database-mode faers_phv_scan_stratified(): push the whole
# stratum x drug x event pair counting down to a single DuckDB query.
#
# The aggregation mirrors scan-db.R, extended with a `strata` CTE built from
# the demo table. For every (stratum, drug, event) combination the query
# returns:
#   a         = COUNT(DISTINCT primaryid) among stratum reports carrying BOTH
#   n_drug    = distinct stratum patients exposed to the drug
#   n_event   = distinct stratum patients reporting the event
#   n_stratum = distinct patients in the stratum (the layer background)
# Only the already-filtered aggregate table enters memory, so stratified
# scans scale to the whole FAERS database on the duckdb backend. Missing
# strata values become 'Missing' (keep mode) or are dropped (drop mode),
# matching the memory backend byte for byte.

# Build the full SQL string for the stratified scan aggregation.
db_scan_strata_query <- function(strata_cols, drug_col, event_col, meddra,
                                 min_a, na_stratum) {
    drug_tbl <- db_field_table("drug")
    reac_tbl <- db_field_table("reac")
    demo_tbl <- db_field_table("demo")
    drug_col_q <- db_qident(drug_col)

    strata_names <- vapply(strata_cols, db_qident, character(1L))
    strata_sel_cols <- paste(strata_names, collapse = ", ")
    strata_case <- vapply(strata_cols, function(col) {
        q <- db_qident(col)
        sprintf(
            "CASE WHEN %1$s IS NULL OR trim(%1$s) = '' THEN '%2$s' ELSE %1$s END AS %1$s",
            q, STRATA_MISSING_LABEL
        )
    }, character(1L))
    strata_sel <- paste(
        sprintf("s.%s", strata_names), collapse = ", "
    )
    strata_grp <- paste(
        sprintf("s.%s", strata_names), collapse = ", "
    )
    strata_keep_where <- if (identical(na_stratum, "drop")) {
        paste(
            sprintf(
                "%1$s IS NOT NULL AND trim(%1$s) <> ''",
                strata_names
            ),
            collapse = " AND "
        )
    } else {
        "TRUE"
    }

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
            db_qident(event_col), reac_tbl,
            db_qident(event_col), db_qident(event_col)
        )
    }
    paste0(
        "WITH drug_u AS (", drug_sel, "),\n",
        "     reac_u AS (", reac_sel, "),\n",
        "     strata AS (\n",
        "         SELECT primaryid, ", paste(strata_case, collapse = ", "), "\n",
        "         FROM ", demo_tbl, "\n",
        "         WHERE ", strata_keep_where, "\n",
        "     ),\n",
        "     pair AS (\n",
        "         SELECT ", strata_sel, ", g.drugname AS drugname,\n",
        "                r.event AS event, COUNT(DISTINCT g.primaryid) AS a\n",
        "         FROM drug_u g\n",
        "         INNER JOIN reac_u r ON g.primaryid = r.primaryid\n",
        "         INNER JOIN strata s ON g.primaryid = s.primaryid\n",
        "         GROUP BY ", strata_grp, ", g.drugname, r.event\n",
        "         HAVING COUNT(DISTINCT g.primaryid) >= ", as.integer(min_a), "\n",
        "     ),\n",
        "     drug_tot AS (\n",
        "         SELECT ", strata_sel, ", g.drugname AS drugname,\n",
        "                COUNT(DISTINCT g.primaryid) AS n_drug\n",
        "         FROM drug_u g\n",
        "         INNER JOIN strata s ON g.primaryid = s.primaryid\n",
        "         GROUP BY ", strata_grp, ", g.drugname\n",
        "     ),\n",
        "     event_tot AS (\n",
        "         SELECT ", strata_sel, ", r.event AS event,\n",
        "                COUNT(DISTINCT r.primaryid) AS n_event\n",
        "         FROM reac_u r\n",
        "         INNER JOIN strata s ON r.primaryid = s.primaryid\n",
        "         GROUP BY ", strata_grp, ", r.event\n",
        "     ),\n",
        "     strat_tot AS (\n",
        "         SELECT ", strata_sel_cols,
        ", COUNT(DISTINCT primaryid) AS n_stratum\n",
        "         FROM strata\n",
        "         GROUP BY ", strata_sel_cols, "\n",
        "     )\n",
        "SELECT ", paste(sprintf("p.%s", strata_names), collapse = ", "),
        ", p.drugname AS drug, p.event, p.a,\n",
        "       d.n_drug, e.n_event, st.n_stratum\n",
        "FROM pair p\n",
        "INNER JOIN drug_tot d ON p.drugname = d.drugname AND ",
        .strata_eq("p", "d", strata_cols), "\n",
        "INNER JOIN event_tot e ON p.event = e.event AND ",
        .strata_eq("p", "e", strata_cols), "\n",
        "INNER JOIN strat_tot st ON ", .strata_eq("p", "st", strata_cols), "\n"
    )
}

# "p.col = d.col AND ..." equality fragment for joining the aggregate tables.
.strata_eq <- function(left, right, strata_cols) {
    paste(
        sprintf(
            "%1$s.%3$s = %2$s.%3$s",
            left, right, vapply(strata_cols, db_qident, character(1L))
        ),
        collapse = " AND "
    )
}

# Run the stratified scan aggregation on a db-backed object and return a
# data.table with columns: <strata cols>, drug, event, a, n_drug, n_event,
# n_stratum.
scan_strata_counts_db <- function(object, strata_cols, drug_col, event_col,
                                  min_a, na_stratum) {
    con <- object@db@con
    meddra <- event_col %chin% names(object@meddra@hierarchy)
    if (meddra && !DB_MEDDRA_TABLE %in% DBI::dbListTables(con)) {
        db_ingest_meddra(con, object)
    }
    out <- DBI::dbGetQuery(
        con,
        db_scan_strata_query(
            strata_cols, drug_col, event_col, meddra, min_a, na_stratum
        )
    )
    out <- data.table::as.data.table(out)
    # COUNT(DISTINCT ...) comes back as BIGINT/double -> coerce to integer,
    # keeping byte-parity with the memory path.
    int_cols <- intersect(
        c("a", "n_drug", "n_event", "n_stratum"), names(out)
    )
    for (col in int_cols) {
        if (is.double(out[[col]])) {
            data.table::set(out, j = col, value = as.integer(out[[col]]))
        }
    }
    data.table::setcolorder(
        out,
        c(strata_cols, "drug", "event", "a", "n_drug", "n_event", "n_stratum")
    )
    out
}
