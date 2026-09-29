# faers 1.5.6

* Added `faers_phv_scan()`, a whole-database signal scanner that enumerates
  all drug x event pairs, builds a contingency table for every pair and runs
  the requested disproportionality methods in one vectorized pass. It ships
  with an optimized backend for `database = "duckdb"` (single out-of-core SQL
  aggregation) and a chunked implementation for the memory backend, and
  supports optional drug whitelisting via `.drug_pattern`.

* Fixed `faers_phv_composite()`: hierarchy event types such as `soc_name`
  no longer abort with "object not found" (event columns are now resolved
  via `faers_get()`, which attaches the MedDRA hierarchy on both the memory
  and the duckdb backends), any standardized `reac` event column is accepted,
  an unused full-database `faers_counts()` call that crashed on duckdb-backed
  objects was removed, de-duplication is now enforced, and the undeclared
  `assertthat` usage was replaced with the internal assertion helpers.

* Added `faers_phv_scan_stratified()`, a stratified variant of the
  whole-database scanner that splits the database into strata defined by
  categorical `demo` columns (`sex`, `age_grp`, `occp_cod`, ...) and builds
  one 2x2 contingency table per `stratum x drug x event` combination against
  the stratum-internal background (`n_stratum`). Reports with missing strata
  values are grouped into an explicit `"Missing"` stratum by default
  (`.na_stratum = "keep"`) or dropped (`"drop"`), with a warning when a
  stratification column is missing for more than 30% of the reports. Both
  backends produce identical results, with the duckdb backend pushing the
  whole stratified aggregation into a single SQL query.

# faers 1.5.5

* Documented the DuckDB backend: `faers()`, `faers_parse()`, `faers_standardize()`
  and the downstream pipeline now accept `database = "duckdb"`, keeping the
  quarterly data out of memory in an on-disk DuckDB database with an identical
  analysis API (`faers_get`, `faers_counts`, `faers_phv_table`, ...).

# faers 1.5.1

* Fixed error in `handle_setopt(h, ...)` caused by unsupported option `multi_timeout`.
  Requires `curl` >= 6.0.0.

# faers 1.1.6

* fda_drugs() now directly use a fixed url to download the data

* faers_meta(internal = TRUE) will always use the cache data in the package. 

* Rename "gndr_cod" into "gender" for periods before 2014q2, and rename "sex" into "gender" for periods after or equal to 2014q2. 

* "sex" was added, which recoded any values other than "F" or "M" as `NA`.

# faers 1.1.4

* fix error when download failed

* meddra augment additional argument `primary_soc` to help save the full meddra data

# faers 0.99.0

* Initial Bioconductor submission.
