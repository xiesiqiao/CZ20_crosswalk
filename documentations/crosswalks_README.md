# U.S. Commuting Zone (2020) & PUMA Crosswalks to Historical Geographies

This repository provides population-weighted crosswalks between **2020 Commuting Zones (CZ20)**, **2010 & 2020 PUMAs**, and historical commuting zones (1990, 2000, 2010 ERS).  
Crosswalks are generated using a 100m population raster derived from the 2020 Census.

---

## Contents

Crosswalk files (CSV), one per source–target pair:

**CZ20 → historic commuting zones:**
* `cz20_to_cz90_popgrid100m_area.csv` – CZ20 → 1990 CZ  
* `cz20_to_cz2010_popgrid100m_area.csv` – CZ20 → 2010 ERS CZ  
* `cz20_to_cz2000_popgrid100m_area.csv` – CZ20 → 2000 ERS CZ  
* `cz20_grid_totals.csv` – CZ20 total population  

**PUMA → CZ:**
* `puma20_to_cz90_popgrid100m_area.csv` – PUMA2020 → 1990 CZ  
* `puma20_to_cz20_popgrid100m_area.csv` – PUMA2020 → CZ20  
* `puma2010_to_cz20_popgrid100m_area.csv` – PUMA2010 → CZ20  
* `puma20_grid_totals.csv` – PUMA2020 total population  

---

## Variables

````markdown
**For CZ20-based crosswalks (CZ20 → target):**

| Variable         | Description                                  |
| ---------------- | -------------------------------------------- |
| CZ20             | 2020 Commuting Zone ID                       |
| \[Target]        | Target geography ID (CZ90, CZ2000, CZ2010)   |
| \[Target]\_name  | Target geography name                        |
| pop\_overlap     | Population overlapping CZ20 and target       |
| share\_of\_CZ20  | Share of CZ20 population in target geography |
| afactor          | Identical to `share_of_CZ20`                 |
| CZ20\_total\_pop | Total population of the CZ20                 |
| cells            | Number of overlapping 100m grid cells        |
| grid\_res\_m     | Grid resolution (100m)                       |
| crs\_epsg        | CRS (5070)                                   |
| method           | "building\_footprint\_area"                  |

**For PUMA-based crosswalks (PUMA → CZ):**

| Variable           | Description                                        |
| ------------------ | -------------------------------------------------- |
| PUMA20 / PUMA2010  | PUMA ID (GISJOIN)                                  |
| \[PUMA]\_name      | PUMA name                                          |
| \[PUMA]\_total\_pop| Total population of the PUMA                       |
| share\_of\_[PUMA]  | Share of PUMA population in target CZ              |
| afactor            | Identical to `share_of_[PUMA]`                     |
| CZ20 / CZ90        | Target commuting zone ID                           |
| \[Target]\_name    | Target commuting zone name (if available)          |
| pop\_overlap       | Population overlapping between PUMA and CZ target  |
| cells              | Number of overlapping 100m grid cells              |
| grid\_res\_m       | Grid resolution (100m)                             |
| crs\_epsg          | CRS (5070)                                         |
| method             | "building\_footprint\_area"                        |

---

## Methodology

* Source population raster: [Population Grid](../population_grid/).  
* Overlayed source geography (CZ20 or PUMA) with each target geography.  
* Calculated overlap population and shares.  
* Shares are always relative to the **source geography** (CZ20 or PUMA).  
* `afactor` is included as a duplicate of the share variable for convenience.  

---

## Interpretation

* **share\_of\_CZ20** = fraction of CZ20 population assigned to a target geography.  
* **share\_of\_PUMA20 / share\_of\_PUMA2010** = fraction of PUMA population assigned to a CZ.  
* Use these as weights to reallocate statistics from the source geography to the target geography.  
* For rates/means, reallocate numerators and denominators separately, then recompute.  

---

## Usage Examples

### In R (CZ20 → CZ2010)

```r
library(dplyr)
cw <- read.csv("cz20_to_cz2010_popgrid100m_area.csv")
data_cz20 <- data.frame(CZ20=c("12345","67890"), value=c(1000,2000))
weighted <- cw %>% 
  left_join(data_cz20, by="CZ20") %>% 
  mutate(weighted_value=value*afactor) %>% 
  group_by(CZ2010) %>% 
  summarise(value=sum(weighted_value, na.rm=TRUE))
````

### In R (PUMA2010 → CZ20)

```r
cw <- read.csv("puma2010_to_cz20_popgrid100m_area.csv")
data_puma10 <- data.frame(PUMA2010=c("G0100010","G0100020"), value=c(4000,7000))
weighted <- cw %>% 
  left_join(data_puma10, by="PUMA2010") %>% 
  mutate(weighted_value=value*afactor) %>% 
  group_by(CZ20) %>% 
  summarise(value=sum(weighted_value, na.rm=TRUE))
```

### In Python (PUMA20 → CZ90)

```python
import pandas as pd
cw = pd.read_csv("puma20_to_cz90_popgrid100m_area.csv")
data_puma20 = pd.DataFrame({"PUMA20":["G0100010","G0100020"], "value":[5000,8000]})
weighted = (cw.merge(data_puma20, on="PUMA20")
              .assign(weighted_value=lambda d: d.value*d.afactor)
              .groupby("CZ90", as_index=False)
              .agg(value=("weighted_value","sum")))
```

---

## Data Sources

* **CZ20:** Fowler, C.S. (2024). *New Commuting Zone delineation for the U.S. using 2020 data.* [https://doi.org/10.17605/OSF.IO/J256U](https://doi.org/10.17605/OSF.IO/J256U)
* **CZ90:** Dorn (1990) commuting zones.
* **CZ2000 / CZ2010:** USDA ERS commuting zones.
* **PUMAs (2010, 2020):** IPUMS NHGIS.
* **Population raster:** see [Population Grid](../population_grid/).

---

## Citations

* Fowler, C.S. (2024). *New Commuting Zone delineation for the U.S. using 2020 data.*
* Dorn, D. (2009). *Commuting Zones and Labor Market Areas: A 1990s View.*
* USDA Economic Research Service (2000, 2010 commuting zones).
* IPUMS NHGIS, University of Minnesota, [www.nhgis.org](http://www.nhgis.org).
* Zhang, Q. & Zhao, P. (2020). *A fine-scale population distribution dataset for China.*

---

```


