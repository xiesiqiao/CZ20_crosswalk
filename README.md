
# 📊 CZ20 Crosswalks & 100m Population Grid (U.S. 2020 census)

This repository provides reproducible code, 100m-resolution population grids (2020 census), and **population-weighted crosswalks** between the new **2020 Commuting Zones (CZ20)** and other standard geographies.

---

## 📂 Contents

* **`/code/`** — R scripts for generating the crosswalks and processing the population grid
* **`/output/population_grid/`** — state-level 100m × 100m population GeoTIFFs (EPSG:5070)
* **`/output/crosswalks/`** — CSV crosswalks from CZ20 to:

  * **CZ90** (Dorn 1990)
  * **CZ2000** (USDA ERS)
  * **CZ2010** (USDA ERS)
  * **PUMA2010** (IPUMS NHGIS)
  * **PUMA2020** (IPUMS NHGIS)

---

## 🧾 Data Sources

* **CZ20 (2020 commuting zones):** Fowler, C.S. (2024). *New Commuting Zone delineation for the U.S. using 2020 data.* [https://doi.org/10.17605/OSF.IO/J256U](https://doi.org/10.17605/OSF.IO/J256U)
* **Population grid (method inspiration):** Zhang, Q. & Zhao, P. (2020). *A fine-scale population distribution dataset for China.* *Big Earth Data, 4(2): 123–139.* [https://doi.org/10.1080/20964471.2020.1776200](https://doi.org/10.1080/20964471.2020.1776200)
* **CZ90:** Dorn, D. (1990). Commuting Zone boundaries.
* **CZ2000, CZ2010:** USDA ERS commuting zones.
* **PUMA2010, PUMA2020:** IPUMS NHGIS shapefiles.

---

## 📐 Crosswalk structure

Each crosswalk is a **CSV file** with the following variables:

| Column           | Description                                                    |
| ---------------- | -------------------------------------------------------------- |
| `state`          | State FIPS code                                                |
| `CZ20`           | CZ20 identifier                                                |
| `CZ20_total_pop` | Total CZ20 population (from grid)                              |
| `share_of_CZ20`  | Share of CZ20 population overlapping with the target geography |
| `<TARGET>`       | Target geography ID (e.g., `PUMA20`, `CZ2010`, etc.)           |
| `<TARGET>_name`  | Target geography name (if available)                           |
| `pop_overlap`    | Population overlapping between CZ20 and the target unit        |
| `cells`          | Number of 100m grid cells overlapping                          |
| `grid_res_m`     | Grid resolution (always 100)                                   |
| `crs_epsg`       | CRS EPSG code (5070)                                           |
| `method`         | `"building_footprint_area"` (weighting method used)            |

---

## 🧮 How to use the crosswalks

* **Totals (e.g., jobs, households):**
  Multiply the CZ20 totals by `share_of_CZ20` and aggregate by the target geography.

* **Rates/means (e.g., average income):**
  Weight by the underlying population. Do **not** simply multiply the mean by the share. Instead, reallocate numerators and denominators separately and recompute.

* **Sanity checks:**
  Within each CZ20, `sum(share_of_CZ20) ≈ 1`.
  National totals are preserved after reallocation.

---

## 🚀 Example Usage (R)

```r
library(readr); library(dplyr)

cw <- read_csv("output/crosswalks/cz20_to_puma20_popgrid100m_area.csv")

# Example: allocate CZ20 job totals to PUMA20
cz20_jobs <- tibble::tibble(CZ20 = c("000001","000002"), jobs = c(1000, 1500))

alloc <- cw %>%
  select(CZ20, PUMA20, PUMA20_name, share_of_CZ20) %>%
  left_join(cz20_jobs, by="CZ20") %>%
  mutate(jobs_alloc = jobs * share_of_CZ20) %>%
  group_by(PUMA20, PUMA20_name) %>%
  summarise(jobs = sum(jobs_alloc, na.rm=TRUE), .groups="drop")

print(alloc)
```

---

## 🔁 Reproducibility

* All processing is scripted in **R** (see `/code/` folder).
* Requires: `sf`, `dplyr`, `data.table`, `terra`, `readr`.
* Duplication instructions are included here:

  ```
  [https://github.com/xiesiqiao/CZ20_crosswalk/tree/main/code]
  ```

---

## Citation

If you use these data, please cite:

* Fowler, C.S. (2024). *New Commuting Zone delineation for the U.S. using 2020 data.* [https://doi.org/10.17605/OSF.IO/J256U](https://doi.org/10.17605/OSF.IO/J256U)
* Zhang, Q. & Zhao, P. (2020). *A fine-scale population distribution dataset for China.* [https://doi.org/10.1080/20964471.2020.1776200](https://doi.org/10.1080/20964471.2020.1776200)

---

## Contact

For questions, please contact:
\[Siqiao Xie] — \[siqiao.xie@colorado.edu]
