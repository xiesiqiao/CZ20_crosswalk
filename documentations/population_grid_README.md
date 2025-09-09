# U.S. 100m Population Grid (2020)

This repository provides a **high-resolution (100m × 100m) population raster** for the United States, created from 2020 Decennial Census population counts and following the methodological framework of Zhang & Zhao (2020).
The grid serves as the foundation for generating fine-grained spatial crosswalks and supports a wide range of spatial analyses.

---

## Contents

Raster files:

* `*_grid100m_pop_2020_area.tif` – state-level raster files

**Metadata:**

* **CRS:** EPSG:5070 (NAD83 / Conus Albers)
* **Resolution:** 100m × 100m
* **Cell values:** estimated number of residents per cell
* **Coverage:** all 50 states
* **Year:** 2020

---

## Variables

Each raster cell contains:

| Variable   | Description                                 |
| ---------- | ------------------------------------------- |
| Cell value | Population count in that 100m × 100m square |
| CRS        | EPSG:5070 (NAD83 / Conus Albers)            |
| Resolution | 100m × 100m grid cells                      |

---

## Methodology

1. **Source data**: 2020 Decennial Census population counts (U.S. Census Bureau).
2. **Disaggregation**: Counts distributed into 100m grid cells, producing a continuous, uniform surface.
3. **Projection**: Reprojected to EPSG:5070 for compatibility with regional analyses and overlay with census geography.
4. **Framework**: Approach inspired by Zhang & Zhao (2020) on fine-scale population datasets.

---

## Usage

* **Population-weighted crosswalks**: Use as the population surface to construct crosswalks between geographic units.
* **Spatial analysis**: Ideal for uniform population density studies and overlay with geographic boundaries.
* **Aggregation**: Flexible for aggregating to arbitrary units (states, commuting zones, counties, custom boundaries).
* **Clipping**: Can be subset by states, commuting zones, or other regions of interest.

---

## Geographic & Temporal Coverage

* **Geographic:** United States, all 50 states
* **Temporal:** 2020

---

## Data Sources

* **Population counts:** U.S. Census Bureau (2020). *Decennial Census of Population and Housing.*
* **Methodology reference:** Zhang, Q. & Zhao, P. (2020). *A fine-scale population distribution dataset for China.* *Big Earth Data*, 4(2), 123–137. [https://doi.org/10.1080/20964471.2020.1776200](https://doi.org/10.1080/20964471.2020.1776200)
* **Geographic boundaries:** U.S. Census TIGER/Line shapefiles.

---

## Citation

When using this dataset, please cite:

Xie, S. (2025). 100*100m Population Grid (2020 Census) [Data set]. University of Colorado Boulder. https://doi.org/10.25810/PB2D-P578

---

✅ Do you want me to package this into a **ready-to-use `README.md` file** like I did for the crosswalks one?
