
<!-- README.md is generated from README.Rmd. Please edit that file -->

# faers

<!-- badges: start -->

[![platform](http://www.bioconductor.org/shields/availability/devel/faers.svg)](https://www.bioconductor.org/packages/devel/bioc/html/faers.html#archives)
[![R-CMD-check](https://github.com/WangLabCSU/faers/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WangLabCSU/faers/actions/workflows/R-CMD-check.yaml)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Ask
DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/WangLabCSU/faers)
<!-- badges: end -->

Here is the polished English version:

------------------------------------------------------------------------

Modern biologics, such as immune checkpoint inhibitors, exhibit complex
toxicity profiles that are often underrepresented in pre-market clinical
trials. While the FAERS database serves as a critical resource for
real-world safety surveillance, its intricate relational structure and
data inconsistencies pose significant barriers to large-scale
epidemiological analyses.

To address these challenges, we developed `faers`, an end-to-end,
reproducible framework for precision pharmacovigilance. The package
streamlines the entire workflow—from raw data acquisition and rigorous
preprocessing to signal detection—empowering researchers to transform
vast spontaneous reporting data into actionable clinical insights.

## Key Features

- 📥 **Data Acquisition**: Automated downloading and parsing of FAERS
  quarterly data (supporting both ASCII and XML formats).

- 🛠️ **Rigorous Preprocessing**: Advanced multi-quarter data merging and
  robust deduplication logic to ensure high data fidelity.

- 🔍 **Terminology Standardization**: Seamless integration with MedDRA,
  RxNorm, and the FDA Drugs API for precise mapping of drugs and adverse
  events.

- 📊 **Advanced Signal Detection**: Comprehensive support for
  disproportionality analysis, including ROR, PRR, BCPNN, and EBGM.

- ⚡ **High-Performance Computing**: Integrated with BiocParallel for
  memory-efficient, parallelized processing of millions of records.

- 🌐 **Knowledge Integration**: Direct support for Athena drug
  vocabularies and Standardised MedDRA Queries (SMQ) for
  mechanism-driven research.

## Installation

To install from Bioconductor, use the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("faers")
```

You can install the development version of `faers` from
[GitHub](https://github.com/WangLabCSU/faers) with:

``` r
if (!requireNamespace("pak")) {
    install.packages("pak",
        repos = sprintf(
            "https://r-lib.github.io/p/pak/devel/%s/%s/%s",
            .Platform$pkgType, R.Version()$os, R.Version()$arch
        )
    )
}
pak::pkg_install("WangLabCSU/faers")
```

## Quick Start

The faers package provides a standardized pipeline that unifies complex
pharmacovigilance workflows. For a comprehensive, step-by-step
demonstration—including data acquisition and a complete Insulin case
study—please refer to our detailed documentation:

👉 **[Full Workflow Tutorial](vignettes/full-workflow.Rmd)**

## sessionInfo

``` r
sessionInfo()
#> R version 4.5.1 (2025-06-13)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Asia/Shanghai
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.5.1    fastmap_1.2.0     cli_3.6.5         tools_4.5.1      
#>  [5] htmltools_0.5.9   rstudioapi_0.17.1 yaml_2.3.12       rmarkdown_2.30   
#>  [9] knitr_1.51        xfun_0.55         digest_0.6.39     rlang_1.1.6      
#> [13] evaluate_1.0.5
```
