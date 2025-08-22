Great — since you’ve now added the new **PUMA2020 → CZ1990** crosswalk, your README should reflect that you have *two types of crosswalk families*:

* **From CZ20 → others** (your original set), and
* **From PUMA20 → CZ90** (your new script/output).

Here’s an updated version of your README with the new info folded in (I kept your structure, just modified the “Contents,” “Crosswalk structure,” and added a short usage note for the new direction).

---

# 📊 CZ20 Crosswalks & 100m Population Grid (U.S. 2020 census)

This repository provides reproducible code, 100m-resolution population grids (2020 census), and **population-weighted crosswalks** linking modern geographies with legacy standards.

Crosswalks are built using **population overlap on a 100m grid**, following Zhang & Zhao (2020).

---

## 📂 Contents

* **`/code/`** — R scripts for generating the crosswalks and processing the population grid
* **`/output/population_grid/`** — state-level 100m × 100m population GeoTIFFs (EPSG:5070)
* **`/output/crosswalks/`** — population-weighted crosswalk CSVs:

  **From CZ20 (2020 commuting zones) to:**

  * **CZ90** (Dorn 1990)
  * **CZ2000** (USDA ERS)
  * **CZ2010** (USDA ERS)
  * **PUMA2010** (IPUMS NHGIS)
  * **PUMA2020** (IPUMS NHGIS)

  **From PUMA2020 to:**

  * **CZ90** (Dorn 1990)

---

## 🧾 Data Sources

* **CZ20 (2020 commuting zones):** Fowler, C.S. (2024). *New Commuting Zone delineation for the U.S. using 2020 data.* [https://doi.org/10.17605/OSF.IO/J256U](https://doi.org/10.17605/OSF.IO/J256U)
* **Population grid (method inspiration):** Zhang, Q. & Zhao, P. (2020). *A fine-scale population distribution dataset for China.* *Big Earth Data, 4(2): 123–139.* [https://doi.org/10.1080/20964471.2020.1776200](https://doi.org/10.1080/20964471.2020.1776200)
* **CZ90:** Dorn, D. (1990). Commuting Zone boundaries.
* **CZ2000, CZ2010:** USDA ERS commuting zones.
* **PUMA2010, PUMA2020:** IPUMS NHGIS shapefiles.

---

## 📐 Crosswalk structure

Each crosswalk is a **CSV file**. The “from” geography is always labeled in the totals/share columns.

**For CZ20 → TARGET:**

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

**For PUMA20 → CZ90:**

| Column             | Description                                         |
| ------------------ | --------------------------------------------------- |
| `state`            | State FIPS code                                     |
| `PUMA20`           | PUMA2020 identifier (GISJOIN)                       |
| `PUMA20_name`      | PUMA2020 name                                       |
| `PUMA20_total_pop` | Total PUMA2020 population (from grid)               |
| `share_of_PUMA20`  | Share of PUMA2020 population overlapping with CZ90  |
| `CZ90`             | CZ1990 identifier                                   |
| `CZ90_name`        | CZ1990 name                                         |
| `pop_overlap`      | Population overlapping between PUMA20 and CZ90      |
| `cells`            | Number of 100m grid cells overlapping               |
| `grid_res_m`       | Grid resolution (always 100)                        |
| `crs_epsg`         | CRS EPSG code (5070)                                |
| `method`           | `"building_footprint_area"` (weighting method used) |

---

## 🧮 How to use the crosswalks

* **Totals (e.g., jobs, households):**
  Multiply the “from” geography totals by `share_of_<FROM>` and aggregate by the target geography.

* **Rates/means (e.g., average income):**
  Weight by the underlying population. Reallocate numerators and denominators separately, then recompute rates.

* **Sanity checks:**
  Within each “from” geography, the shares sum to \~1.
  National totals are preserved after reallocation.

---

## 🚀 Example Usage (R)

### CZ20 → PUMA20 (jobs allocation)

```r
cw <- read_csv("output/crosswalks/cz20_to_puma20_popgrid100m_area.csv")
cz20_jobs <- tibble::tibble(CZ20 = c("000001","000002"), jobs = c(1000, 1500))

alloc <- cw %>%
  select(CZ20, PUMA20, PUMA20_name, share_of_CZ20) %>%
  left_join(cz20_jobs, by="CZ20") %>%
  mutate(jobs_alloc = jobs * share_of_CZ20) %>%
  group_by(PUMA20, PUMA20_name) %>%
  summarise(jobs = sum(jobs_alloc, na.rm=TRUE), .groups="drop")
```

### PUMA20 → CZ90 (population reallocation)

```r
cw <- read_csv("output/crosswalks/puma20_to_cz90_popgrid100m_area.csv")
puma_pop <- tibble::tibble(PUMA20 = c("G0100010","G0100020"), people = c(5000, 8000))

alloc <- cw %>%
  select(PUMA20, CZ90, CZ90_name, share_of_PUMA20) %>%
  left_join(puma_pop, by="PUMA20") %>%
  mutate(people_alloc = people * share_of_PUMA20) %>%
  group_by(CZ90, CZ90_name) %>%
  summarise(people = sum(people_alloc, na.rm=TRUE), .groups="drop")
```

---

## 🔁 Reproducibility

* All processing is scripted in **R** (see `/code/` folder).
* Requires: `sf`, `dplyr`, `data.table`, `terra`, `readr`.
* Duplication instructions are included here:
  [https://github.com/xiesiqiao/CZ20\_crosswalk/tree/main/code](https://github.com/xiesiqiao/CZ20_crosswalk/tree/main/code)

---

## Citation

If you use these data, please cite:



---

## Contact

For questions, please contact:
**Siqiao Xie** — [siqiao.xie@colorado.edu](mailto:siqiao.xie@colorado.edu)

---

Do you want me to also **add a diagram/table** that summarizes both crosswalk “directions” (CZ20→X vs PUMA20→CZ90), so users see at a glance what’s available?
