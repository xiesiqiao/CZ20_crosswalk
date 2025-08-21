# ==========================================================
# Example: Using CZ20 → {PUMA20, CZ90, CZ2010, CZ2000, PUMA2010} crosswalks
# - Reallocate CZ20 totals to a target geography with population weights
# - Compute population-weighted averages from CZ20 to a target
# ==========================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
})

# --- Paths (adjust to your machine) ---
base_root <- "C:/Users/xiesi/ASU Dropbox/Xie Sijiao/Connor-Kemeny-Storper/Microdata_wealth_estimates/Crosswalk/CZ20_CZ90_xwalk"
cw_dir    <- file.path(base_root, "output", "crosswalks")

# Choose a target geography
cw_file <- file.path(cw_dir, "cz20_to_puma20_popgrid100m_area.csv")   # or: ...cz90..., ...cz2010..., ...cz2000..., ...puma2010...

# Read crosswalk
cw <- readr::read_csv(cw_file, show_col_types = FALSE)

# Inspect columns (all crosswalks have the same schema)
# head(cw)
# names(cw)
# [1] "state" "CZ20" "<TARGET>" "<TARGET>_name" "pop_overlap" "share_of_CZ20" "CZ20_total_pop" "cells" "grid_res_m" "crs_epsg" "method"

# ------------------------------
# EXAMPLE A: Reallocate totals from CZ20 to target (e.g., total jobs)
# ------------------------------
# Pretend we have some CZ20 totals:
cz20_totals <- tibble::tibble(
  CZ20  = c("000001","000002","000003"),
  jobs  = c(1200, 3400, 1800),
  hhlds = c(500,  900,  700)
)

target_col <- setdiff(grep("^PUMA20$|^CZ90$|^CZ2010$|^CZ2000$|^PUMA2010$", names(cw), value = TRUE), "CZ20")
stopifnot(length(target_col) == 1)
tname_col  <- paste0(target_col, "_name")

# Allocate each CZ20 metric by share_of_CZ20 and sum by target
alloc_totals <- cw %>%
  select(CZ20, !!target_col, !!tname_col, share_of_CZ20) %>%
  left_join(cz20_totals, by = "CZ20") %>%
  mutate(across(where(is.numeric) & !c(CZ20 %in% names(.)), ~ .x * share_of_CZ20, .names = "alloc_{.col}")) %>%
  group_by(across(all_of(c(target_col, tname_col)))) %>%
  summarise(across(starts_with("alloc_"), ~ sum(.x, na.rm = TRUE), .names = "{.col}"), .groups = "drop") %>%
  rename_with(~ sub("^alloc_", "", .x), starts_with("alloc_"))

# Result: totals for the target geography
print(alloc_totals)

# ------------------------------
# EXAMPLE B: Population-weighted average from CZ20 to target
# (e.g., average income available at CZ20 with CZ20 population)
# ------------------------------
cz20_avg <- tibble::tibble(
  CZ20            = c("000001","000002","000003"),
  avg_income      = c(52000, 61000, 47000),
  CZ20_population = c(25000, 42000, 18000)  # denominator for weighting
)

weighted_avg <- cw %>%
  select(CZ20, !!target_col, share_of_CZ20, CZ20_total_pop) %>%
  left_join(cz20_avg, by = "CZ20") %>%
  mutate(weight = share_of_CZ20 * CZ20_population,
         numer  = avg_income * weight) %>%
  group_by(across(all_of(target_col))) %>%
  summarise(avg_income_pw = sum(numer, na.rm = TRUE) / sum(weight, na.rm = TRUE),
            pop_used      = sum(weight, na.rm = TRUE),
            .groups = "drop")

print(weighted_avg)

# ------------------------------
# EXAMPLE C: Sanity checks
# ------------------------------
# 1) Shares per CZ20 sum to ~1
share_check <- cw %>%
  group_by(CZ20) %>%
  summarise(sum_share = sum(share_of_CZ20, na.rm = TRUE)) %>%
  filter(abs(sum_share - 1) > 1e-6)
if (nrow(share_check) == 0) message("OK: shares sum to ~1 within each CZ20.")

# 2) Allocation preserves national totals (e.g., jobs)
nat_cz20_jobs <- sum(cz20_totals$jobs)
nat_target_jobs <- sum(alloc_totals$jobs, na.rm = TRUE)
message(sprintf("CZ20 jobs total = %s; Target total = %s (should match)",
                format(nat_cz20_jobs, big.mark=","), format(nat_target_jobs, big.mark=",")))
