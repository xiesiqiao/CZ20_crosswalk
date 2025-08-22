# ---- add_puma_ids_by_GISJOIN.R -------------------------------------------

library(sf)
library(dplyr)
library(readr)

# 1) Paths ------------------------------------------------------------------
shp2010 <- "C:/Users/xiesi/ASU Dropbox/Xie Sijiao/Connor-Kemeny-Storper/Microdata_wealth_estimates/Crosswalk/CZ20_CZ90_xwalk/input/processed/boundaries_2010/ipums_puma_2010_tl20.shp"
shp2020 <- "C:/Users/xiesi/ASU Dropbox/Xie Sijiao/Connor-Kemeny-Storper/Microdata_wealth_estimates/Crosswalk/CZ20_CZ90_xwalk/input/processed/boundaries_2020/ipums_puma_2020.shp"

csv_2010 <- "C:/Users/xiesi/ASU Dropbox/Xie Sijiao/Connor-Kemeny-Storper/Microdata_wealth_estimates/Crosswalk/CZ20_CZ90_xwalk/output/crosswalks/cz20_to_puma2010_popgrid100m_area.csv"
csv_2020a <- "C:/Users/xiesi/ASU Dropbox/Xie Sijiao/Connor-Kemeny-Storper/Microdata_wealth_estimates/Crosswalk/CZ20_CZ90_xwalk/output/crosswalks/cz20_to_puma20_popgrid100m_area.csv"
csv_2020b <- "C:/Users/xiesi/ASU Dropbox/Xie Sijiao/Connor-Kemeny-Storper/Microdata_wealth_estimates/Crosswalk/CZ20_CZ90_xwalk/output/crosswalks/puma20_to_cz90_popgrid100m_area.csv"

# 2) Read shapefiles and keep only GISJOIN, PUMA, GEOID ---------------------
p10 <- st_read(shp2010, quiet = TRUE) |>
  st_drop_geometry() |>
  select(GISJOIN, PUMA2010 = PUMA, PUMA10_GEOID = GEOID)

p20 <- st_read(shp2020, quiet = TRUE) |>
  st_drop_geometry() |>
  select(GISJOIN, PUMA2020 = PUMA, PUMA200_GEOID = GEOID)

# 3) Join to CSVs -----------------------------------------------------------
x1 <- read_csv(csv_2010, show_col_types = FALSE) |>
  rename(PUMA2010_GISJOIN=PUMA2010) |>
  left_join(p10, by = c("PUMA2010_GISJOIN"="GISJOIN"))|>
  select(CZ20, PUMA2010,state_fips,afactor,everything())
write_csv(x1, sub("\\.csv$", "", csv_2010))

x2 <- read_csv(csv_2020a, show_col_types = FALSE) |>
  rename(PUMA20_GISJOIN=PUMA20) |>
  left_join(p20, by = c("PUMA20_GISJOIN"="GISJOIN"))|>
  select(CZ20, PUMA2020,state_fips,afactor,everything())
write_csv(x2, sub("\\.csv$", "", csv_2020a))

x3 <- read_csv(csv_2020b, show_col_types = FALSE) |>
  rename(PUMA20_GISJOIN=PUMA20) |>
  left_join(p20, by = c("PUMA20_GISJOIN"="GISJOIN"))|>
  select(PUMA2020,state_fips,CZ90, afactor,everything())
write_csv(x3, sub("\\.csv$", "", csv_2020b))

# --------------------------------------------------------------------------
