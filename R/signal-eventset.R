#' Create contingency table for a specific set of adverse events
#' 
#' @description
#' Build contingency tables for disproportionality analysis focusing on a specific
#' set of adverse events (e.g., irAEs). This function extends `faers_phv_table`
#' to handle composite event definitions.
#' 
#' @param .object A [FAERSascii] object. The unique number of `primaryids` from
#' `.object` will be regarded as `n1.`.
#' @param .object A [FAERSascii] object. It must be standardized and
#' de-duplicated. The unique number of `primaryids` from `.object` will be
#' regarded as `n1.`.
#' @param .event_set A character vector of event terms or a function that
#' defines the set of adverse events of interest.
#' @param .event_type A string specifying the event column of the standardized
#' `reac` data to use, e.g. `"pt"` (Preferred Term), `"meddra_code"`,
#' `"meddra_pt"`, or any MedDRA hierarchy column (`"llt_name"`, `"pt_name"`,
#' `"hlt_name"`, `"hlgt_name"`, `"soc_name"`, ...). Defaults to `"pt"`.
#' @param ... Other arguments passed to specific methods.
#' @return A [data.table][data.table::data.table] object with contingency tables
#' for the specified event set.
#' @export
#' @name faers_phv_composite
methods::setGeneric(
  "faers_phv_composite",
  function(.object, .event_set, ..., .full, .object2) {
    rlang::check_exclusive(.full, .object2)
    standardGeneric("faers_phv_composite")
  }
)

#' @param .full A [FAERSascii] object with data from full database.
#' @inheritParams faers_counts
#' @export
#' @rdname faers_phv_composite
methods::setMethod(
  "faers_phv_composite",
  c(.object = "FAERSascii", .full = "FAERSascii", .object2 = "missing"),
  function(.object, .event_set, .event_type = "pt", ..., .full, .object2) {
    if (!.object@standardization) {
      cli::cli_abort("{.arg .object} must be standardized using {.fn faers_standardize}")
    }
    if (!.object@deduplication) {
      cli::cli_abort("{.arg .object} must be de-duplicated using {.fn faers_dedup}")
    }
    if (!.full@standardization) {
      cli::cli_abort("{.arg .full} must be standardized using {.fn faers_standardize}")
    }
    if (!.full@deduplication) {
      cli::cli_abort("{.arg .full} must be de-duplicated using {.fn faers_dedup}")
    }
    
    full_primaryids <- faers_primaryid(.full)
    interested_primaryids <- faers_primaryid(.object)
    
    if (!all(interested_primaryids %chin% full_primaryids)) {
      cli::cli_abort("Provided {.arg .object} data must be a subset of {.arg .full}")
    }
    
    n <- length(unique(full_primaryids))
    n1. <- length(unique(interested_primaryids))
    
    event_table <- .create_event_set_table(.object, .full, .event_set, .event_type, 
                                           interested_primaryids, full_primaryids, 
                                           n1., n, ...)
    
    return(event_table)
  }
)

#' @param .object2 A [FAERSascii] object with data from another interested drug.
#' @export
#' @rdname faers_phv_composite
methods::setMethod(
  "faers_phv_composite",
  c(.object = "FAERSascii", .full = "missing", .object2 = "FAERSascii"),
  function(.object, .event_set, .event_type = "pt", ..., .full, .object2) {
    # Input validation
    if (!.object@standardization) {
      cli::cli_abort("{.arg .object} must be standardized using {.fn faers_standardize}")
    }
    if (!.object@deduplication) {
      cli::cli_abort("{.arg .object} must be de-duplicated using {.fn faers_dedup}")
    }
    if (!.object2@standardization) {
      cli::cli_abort("{.arg .object2} must be standardized using {.fn faers_standardize}")
    }
    if (!.object2@deduplication) {
      cli::cli_abort("{.arg .object2} must be de-duplicated using {.fn faers_dedup}")
    }
    
    primaryids <- faers_primaryid(.object)
    primaryids2 <- faers_primaryid(.object2)
    
    overlapped_idx <- primaryids %chin% primaryids2
    if (any(overlapped_idx)) {
      cli::cli_warn("{.val {sum(overlapped_idx)}} report{?s} are overlapped between {.arg .object} and {.arg .object2}")
    }
    

    n1. <- length(unique(primaryids))
    n0. <- length(unique(primaryids2))

    event_table <- .create_event_set_table_comparison(.object, .object2, .event_set, 
                                                      .event_type, primaryids, primaryids2,
                                                      n1., n0., ...)
    
    return(event_table)
  }
)

#' Internal function to create contingency table for event set (full database comparison)
#' @keywords internal
.create_event_set_table <- function(object, full, event_set, event_type, 
                                    object_ids, full_ids, n1, n_total, ...) {
  
  event_patients <- .identify_event_set_patients(object, event_set, event_type)
  full_event_patients <- .identify_event_set_patients(full, event_set, event_type)

  composite_event <- data.table::data.table(
    event = "Composite_Event_Set",
    n.1 = length(unique(full_event_patients$primaryid))
  )
  
  # Calculate counts for object
  object_event_count <- length(unique(
    event_patients[primaryid %chin% object_ids, primaryid]
  ))
  
  interested_counts <- data.table::data.table(
    event = "Composite_Event_Set",
    a = object_event_count
  )

  out <- merge(composite_event, interested_counts, by = "event", all = TRUE)
  
  out[, a := data.table::fifelse(is.na(a), 0L, a)]
  out[, b := n1 - a]
  out[, c := n.1 - a]
  out[, d := n_total - (n1 + n.1 - a)]
  out <- out[, !"n.1"]
  
  data.table::setcolorder(out, c("event", "a", "b", "c", "d"))[]
}

#' Internal function to create contingency table for event set (drug comparison)
#' @keywords internal
.create_event_set_table_comparison <- function(object1, object2, event_set, event_type,
                                               object1_ids, object2_ids, n1, n0, ...) {
  
  event_patients1 <- .identify_event_set_patients(object1, event_set, event_type)
  event_patients2 <- .identify_event_set_patients(object2, event_set, event_type)

  a <- length(unique(event_patients1[primaryid %chin% object1_ids, primaryid]))
  c <- length(unique(event_patients2[primaryid %chin% object2_ids, primaryid]))

  out <- data.table::data.table(
    event = "Composite_Event_Set",
    a = a,
    b = n1 - a,
    c = c,
    d = n0 - c
  )
  
  return(out)
}

#' Internal function to identify patients with events in the specified set
#' @keywords internal
.identify_event_set_patients <- function(object, event_set, event_type, ...) {

  # `faers_get()` attaches the MedDRA hierarchy columns on both the memory and
  # the duckdb backends, so any standardized `reac` column (e.g. `pt`,
  # `meddra_code`, `soc_name`) can be used as the event column.
  data_table <- faers_get(object, "reac")
  if (!event_type %chin% names(data_table)) {
    cli::cli_abort(c(
      "{.val {event_type}} is not a column of the standardized {.field reac} data",
      i = "Available columns: {.val {names(data_table)}}"
    ))
  }

  if (is.character(event_set)) {
    event_patients <- data_table[get(event_type) %in% event_set, .(primaryid)]
  } else if (is.logical(event_set) || is.function(event_set)) {
    if (is.function(event_set)) {
      matches <- event_set(data_table[[event_type]])
    } else {
      matches <- event_set
    }
    event_patients <- data_table[matches, .(primaryid)]
  } else {
    cli::cli_abort("Unsupported event_set type. Must be character vector, logical vector, or function.")
  }
  
  return(unique(event_patients))
}

#' Signal detection for specific event sets
#' @export
#' @rdname faers_phv_composite
methods::setGeneric("faers_phv_signal_composite", function(.object, ...) {
  standardGeneric("faers_phv_signal_composite")
})

#' @param .methods Analysis methods to use (passed to [phv_signal]).
#' @param .phv_signal_params Other arguments passed to [phv_signal].
#' @inheritParams phv_signal
#' @seealso [phv_signal]
#' @export
#' @method faers_phv_signal_composite FAERSascii
#' @rdname faers_phv_composite
methods::setMethod(
  "faers_phv_signal_composite", 
  "FAERSascii", 
  function(.object, .event_set, .methods = NULL, ..., 
           .phv_signal_params = list(), BPPARAM = BiocParallel::SerialParam()) {
    
    assert_(.phv_signal_params, is.list, "a list")
    
    out <- faers_phv_composite(.object = .object, .event_set = .event_set, ...)
    
    .__signal__. <- do.call(
      phv_signal, 
      c(
        out[, c("a", "b", "c", "d")],
        list(methods = .methods, BPPARAM = BPPARAM),
        .phv_signal_params
      )
    )
    
    out[, names(.__signal__.) := .__signal__.][]
  }
)

utils::globalVariables(c(".", "a", "b", "c", "d", "n.1", "primaryid", "pt", "soc_name"))
