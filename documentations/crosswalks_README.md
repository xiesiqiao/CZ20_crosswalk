# Population-Weighted Crosswalks: 2020 Commuting Zones, PUMAs, and Historical Geographies (1990–2020)

This repository provides **population-weighted crosswalks** between modern and historical U.S. geographies, including **2020 Commuting Zones (CZ20)**, **2010 & 2020 PUMAs**, and earlier commuting zones (1990, 2000, 2010 ERS).
Crosswalks are constructed using a 100m population raster derived from the 2020 Census.



## Contents

Crosswalk CSV files, one per source–target pair:

**CZ20 → historic commuting zones:**

* `cz20_to_cz90_popgrid100m_area.csv` – CZ20 → 1990 CZ
* `cz20_to_cz2000_popgrid100m_area.csv` – CZ20 → 2000 ERS CZ
* `cz20_to_cz2010_popgrid100m_area.csv` – CZ20 → 2010 ERS CZ
* `cz20_grid_totals.csv` – CZ20 population totals

**PUMA → CZ:**

* `puma20_to_cz20_popgrid100m_area.csv` – PUMA2020 → CZ20
* `puma20_to_cz90_popgrid100m_area.csv` – PUMA2020 → CZ90
* `puma2010_to_cz20_popgrid100m_area.csv` – PUMA2010 → CZ20
* `puma20_grid_totals.csv` – PUMA2020 population totals



## Variables

Each crosswalk contains:

| Variable           | Description                                                           |
| ------------------ | --------------------------------------------------------------------- |
| `<FROM>`           | Source geography ID (e.g., CZ20, PUMA20, PUMA2010)                    |
| `<FROM>_name`      | Source geography name (if available)                                  |
| `<FROM>_total_pop` | Total source geography population (from 100m grid)                    |
| `<TARGET>`         | Target geography ID (CZ90, CZ2000, CZ2010, CZ20, etc.)                |
| `<TARGET>_name`    | Target geography name (if available)                                  |
| `pop_overlap`      | Population overlapping source and target                              |
| `share_of_<FROM>`  | Share of source population in target geography (also named `afactor`) |
| `afactor`          | Duplicate of `share_of_<FROM>` (convenience)                          |
| `cells`            | Number of overlapping 100m grid cells                                 |
| `state`            | State identifier                                                      |
| `grid_res_m`       | Grid resolution (100m)                                                |
| `crs_epsg`         | CRS code (5070)                                                       |
| `method`           | Weighting method used (`building_footprint_area`)                     |



## Methodology

* Source population raster: [Population Grid](../population_grid/) (100m, 2020 Census).
* Overlayed source geography (CZ20 or PUMA) with each target geography.
* Calculated overlap populations and shares.
* Shares are **always relative to the source geography** (CZ20 or PUMA).
* `afactor` is included as a convenience duplicate of the share variable.



## Interpretation

* **share\_of\_CZ20** = fraction of CZ20 population in a target geography.
* **share\_of\_PUMA20 / share\_of\_PUMA2010** = fraction of PUMA population in a target CZ.
* Use shares as weights to reallocate statistics from the source geography to the target.
* For rates/means: reallocate numerators and denominators separately, then recompute.



## Usage Examples

### R (CZ20 → CZ2010)

```r
library(dplyr)
cw <- read.csv("cz20_to_cz2010_popgrid100m_area.csv")
data_cz20 <- data.frame(CZ20=c("12345","67890"), value=c(1000,2000))

weighted <- cw %>% 
  left_join(data_cz20, by="CZ20") %>% 
  mutate(weighted_value=value*afactor) %>% 
  group_by(CZ2010) %>% 
  summarise(value=sum(weighted_value, na.rm=TRUE))
```

### R (PUMA2010 → CZ20)

```r
cw <- read.csv("puma2010_to_cz20_popgrid100m_area.csv")
data_puma10 <- data.frame(PUMA2010=c("G0100010","G0100020"), value=c(4000,7000))

weighted <- cw %>% 
  left_join(data_puma10, by="PUMA2010") %>% 
  mutate(weighted_value=value*afactor) %>% 
  group_by(CZ20) %>% 
  summarise(value=sum(weighted_value, na.rm=TRUE))
```

### Python (PUMA20 → CZ90)

```python
import pandas as pd

cw = pd.read_csv("puma20_to_cz90_popgrid100m_area.csv")
data_puma20 = pd.DataFrame({"PUMA20":["G0100010","G0100020"], "value":[5000,8000]})

weighted = (cw.merge(data_puma20, on="PUMA20")
              .assign(weighted_value=lambda d: d.value*d.afactor)
              .groupby("CZ90", as_index=False)
              .agg(value=("weighted_value","sum")))
```



## Geographic & Temporal Coverage

* **Geographic:** United States
* **Temporal:** 1990, 2000, 2010, 2020



## Data Sources

* **CZ20:** Fowler, C.S. (2024). *New Commuting Zone delineation for the U.S. using 2020 data.* [https://doi.org/10.17605/OSF.IO/J256U](https://doi.org/10.17605/OSF.IO/J256U)
* **CZ90:** Dorn, D. (1990 commuting zones).
* **CZ2000 & CZ2010:** USDA Economic Research Service.
* **PUMAs (2010, 2020):** IPUMS NHGIS.
* **Population raster:** U.S. Census Bureau (2020 population counts).
* **Supplementary:** Zhang, Q. & Zhao, P. (2020). *A fine-scale population distribution dataset for China.* [https://doi.org/10.1080/20964471.2020.1776200](https://doi.org/10.1080/20964471.2020.1776200)



## Citation

Please cite as:

Xie, S. (2025). *Crosswalks between Commuting Zones and Historical Geographies (1990–2020) [Data set]. University of Colorado Boulder. https://doi.org/10.25810/5XP7-YE36*
