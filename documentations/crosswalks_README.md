
# U.S. Commuting Zone (2020) & PUMA Crosswalks to Historical Geographies

This repository provides population-weighted crosswalks between **2020 Commuting Zones (CZ20)**, **2020 PUMAs**, and other geographic units: 1990 commuting zones, 2000 and 2010 ERS commuting zones, and 2010 & 2020 PUMAs. Crosswalks are generated using a 100m population raster derived from the 2020 Census.

## Contents

Crosswalk files (CSV), one per source–target pair:

* `cz20_to_cz90_popgrid100m_area.csv` – CZ20 → 1990 CZ
* `cz20_to_puma20_popgrid100m_area.csv` – CZ20 → 2020 PUMA
* `cz20_to_puma2010_popgrid100m_area.csv` – CZ20 → 2010 PUMA
* `cz20_to_cz2010_popgrid100m_area.csv` – CZ20 → 2010 ERS CZ
* `cz20_to_cz2000_popgrid100m_area.csv` – CZ20 → 2000 ERS CZ
* `puma20_to_cz90_popgrid100m_area.csv` – PUMA20 → 1990 CZ
* `cz20_grid_totals.csv` – CZ20 total population
* `puma20_grid_totals.csv` – PUMA20 total population

## Variables

**For CZ20-based crosswalks (CZ20 → target):**

| Variable         | Description                                  |
| ---------------- | -------------------------------------------- |
| state            | State FIPS or abbreviation                   |
| CZ20             | 2020 Commuting Zone ID                       |
| \[Target]        | Target geography ID (CZ90, PUMA20, etc.)     |
| \[Target]\_name  | Target geography name                        |
| pop\_overlap     | Population overlapping CZ20 and target       |
| share\_of\_CZ20  | Share of CZ20 population in target geography |
| CZ20\_total\_pop | Total population of the CZ20                 |
| cells            | Number of overlapping 100m grid cells        |
| grid\_res\_m     | Grid resolution (100m)                       |
| crs\_epsg        | CRS (5070)                                   |
| method           | "building\_footprint\_area"                  |

**For PUMA20 → CZ90 crosswalk:**

| Variable           | Description                            |
| ------------------ | -------------------------------------- |
| state              | State FIPS or abbreviation             |
| PUMA20             | 2020 PUMA ID (GISJOIN)                 |
| PUMA20\_name       | PUMA name                              |
| CZ90               | 1990 Commuting Zone ID                 |
| CZ90\_name         | 1990 CZ name                           |
| pop\_overlap       | Population overlapping PUMA20 and CZ90 |
| share\_of\_PUMA20  | Share of PUMA20 population in CZ90     |
| PUMA20\_total\_pop | Total population of the PUMA20         |
| cells              | Number of overlapping 100m grid cells  |
| grid\_res\_m       | Grid resolution (100m)                 |
| crs\_epsg          | CRS (5070)                             |
| method             | "building\_footprint\_area"            |

## Methodology

* Source population raster: [Population Grid](../population_grid/).
* Overlayed source geography (CZ20 or PUMA20) with each target geography.
* Calculated overlap population and shares.
* Shares are always relative to the **source geography** (CZ20 or PUMA20).

## Interpretation

* **share\_of\_CZ20** = fraction of CZ20 population assigned to a target geography.
* **share\_of\_PUMA20** = fraction of PUMA20 population assigned to CZ90.
* Use these as weights to reallocate statistics from the source geography to the target geography.
* For rates/means, reallocate numerators and denominators separately, then recompute.

## Usage Examples

### In R (CZ20 → PUMA20)

```r
library(dplyr)
cw <- read.csv("cz20_to_puma20_popgrid100m_area.csv")
data_cz20 <- data.frame(CZ20=c("12345","67890"), value=c(1000,2000))
weighted <- cw %>% 
  left_join(data_cz20, by="CZ20") %>% 
  mutate(weighted_value=value*share_of_CZ20) %>% 
  group_by(PUMA20) %>% 
  summarise(value=sum(weighted_value, na.rm=TRUE))
```

### In R (PUMA20 → CZ90)

```r
cw <- read.csv("puma20_to_cz90_popgrid100m_area.csv")
data_puma20 <- data.frame(PUMA20=c("G0100010","G0100020"), value=c(5000,8000))
weighted <- cw %>% 
  left_join(data_puma20, by="PUMA20") %>% 
  mutate(weighted_value=value*share_of_PUMA20) %>% 
  group_by(CZ90) %>% 
  summarise(value=sum(weighted_value, na.rm=TRUE))
```

### In Python (CZ20 → PUMA20)

```python
import pandas as pd
cw = pd.read_csv("cz20_to_puma20_popgrid100m_area.csv")
data_cz20 = pd.DataFrame({"CZ20":["12345","67890"], "value":[1000,2000]})
weighted = (cw.merge(data_cz20, on="CZ20")
              .assign(weighted_value=lambda d: d.value*d.share_of_CZ20)
              .groupby("PUMA20", as_index=False)
              .agg(value=("weighted_value","sum")))
```

## Data Sources

* **CZ20:** Fowler, C.S. (2024). *New Commuting Zone delineation for the U.S. using 2020 data.* [https://doi.org/10.17605/OSF.IO/J256U](https://doi.org/10.17605/OSF.IO/J256U)
* **CZ90:** Dorn (1990) commuting zones.
* **CZ2000 / CZ2010:** USDA ERS commuting zones.
* **PUMAs:** IPUMS NHGIS (2010, 2020).
* **Population raster:** see [Population Grid](../population_grid/).

## Citations

* Fowler, C.S. (2024). *New Commuting Zone delineation for the U.S. using 2020 data.* [https://doi.org/10.17605/OSF.IO/J256U](https://doi.org/10.17605/OSF.IO/J256U)
* Dorn, D. (2009). "Commuting Zones and Labor Market Areas: A 1990s View."
* USDA Economic Research Service (2000, 2010 commuting zones).
* IPUMS NHGIS, University of Minnesota, [www.nhgis.org](http://www.nhgis.org).
* Zhang, Q. & Zhao, P. (2020). *A fine-scale population distribution dataset for China.* [https://doi.org/10.1080/20964471.2020.1776200](https://doi.org/10.1080/20964471.2020.1776200)