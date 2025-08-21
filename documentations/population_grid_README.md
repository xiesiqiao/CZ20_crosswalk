# U.S. 100m Population Grid (2020)

This repository provides a high-resolution (100m × 100m) population raster for the contiguous United States, created using 2020 Census population counts and building on ideas from Zhang & Zhao (2020).  

The grid serves as the foundation for generating fine-grained spatial crosswalks across geographic boundaries.

## Contents
- **Raster files (`*_grid100m_pop_2020_area.tif`)**
  - CRS: EPSG:5070 (NAD83 / Conus Albers)
  - Cell size: 100m × 100m
  - Value: estimated number of residents per cell
  - Coverage: all 50 states

## Methodology
1. **Source population data**: U.S. Census Bureau, 2020 Decennial Census counts.  
2. **Grid construction**: Population counts disaggregated to 100m resolution following the conceptual framework in Zhang & Zhao (2020).  
3. **Projection**: Reprojected to EPSG:5070 for compatibility with regional analyses and overlay with census geography.

## Variables
- Raster cell value: population count in that 100m square.  
- Metadata: CRS (EPSG:5070), grid cell resolution (100m).

## Usage
- Supports population-weighted crosswalk creation across geographies.  
- Suitable for spatial analysis requiring uniform population surface.  
- Can be clipped by states, commuting zones, or custom boundaries.

## Sources & Citations
- U.S. Census Bureau (2020). *Decennial Census of Population and Housing.*  
- Zhang, Q. & Zhao, P. (2020). *A fine-scale population distribution dataset for China.* *Big Earth Data*, 4(2), 123–137. https://doi.org/10.1080/20964471.2020.1776200  
- CRS and geographic boundaries from U.S. Census TIGER/Line shapefiles.

If you use this dataset, please cite:
> Zhang, Q. & Zhao, P. (2020). A fine-scale population distribution dataset for China. *Big Earth Data*, 4(2), 123–137. https://doi.org/10.1080/20964471.2020.1776200
