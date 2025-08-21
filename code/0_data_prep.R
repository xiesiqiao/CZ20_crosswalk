# Minimal input prep for CZ20_CZ90 crosswalk
# - Keep shapefiles as shapefiles (no reprojection/format change)
# - Keep everything in EPSG:4326 where applicable (buildings already are)
# - Clean NHGIS block CSV down to requested columns

suppressPackageStartupMessages({
  library(fs)
  library(readr)
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ------------------------------------------------------------------
# Root paths (edit if your root changes)
# ------------------------------------------------------------------
root <- "C:/Users/xiesi/ASU Dropbox/Xie Sijiao/Connor-Kemeny-Storper/Microdata_wealth_estimates/Crosswalk/CZ20_CZ90_xwalk"
in_dir <- file.path(root, "input")

in_census <- file.path(in_dir, "2020 census")
in_bldg   <- file.path(in_dir, "building coverage")
in_blocks <- file.path(in_dir, "block boundaries")
in_1990   <- file.path(in_dir, "1990 boundary")
in_2020   <- file.path(in_dir, "2020 boundary")

# outputs (organized but simple)
out_base     <- file.path(in_dir, "processed")
out_bnds1990 <- file.path(out_base, "boundaries_1990")      # copies of your 1990 shapefiles
out_bnds2020 <- file.path(out_base, "boundaries_2020")      # copies of your 2020 shapefiles
out_bldg     <- file.path(out_base, "buildings")            # unzipped per-state GeoJSONs
out_blocks   <- file.path(out_base, "block_boundaries")     # unzipped per-state block shapefiles
out_tables   <- file.path(out_base, "tables")               # cleaned NHGIS CSV

dir_create(c(out_base, out_bnds1990, out_bnds2020, out_bldg, out_blocks, out_tables))

# ------------------------------------------------------------------
# 1) Copy boundary shapefiles AS-IS (no reprojection/format change)
# ------------------------------------------------------------------
copy_shapefile_set <- function(src_dir, pattern_base, dst_dir) {
  # Copies all common shapefile sidecars matching a base name
  parts <- c(".shp", ".shx", ".dbf", ".prj", ".cpg", ".sbn", ".sbx", ".qix")
  for (ext in parts) {
    src <- file.path(src_dir, paste0(pattern_base, ext))
    if (file_exists(src)) file_copy(src, dst_dir, overwrite = TRUE)
  }
}

# 2020: you showed cz20.* and county20.*
copy_shapefile_set(in_2020, "cz20",     out_bnds2020)
copy_shapefile_set(in_2020, "county20", out_bnds2020)

# 1990: you showed cz.* (commuting zones)
copy_shapefile_set(in_1990, "cz", out_bnds1990)

message("✓ Boundaries copied to processed/ (unchanged shapefiles)")

# ------------------------------------------------------------------
# 2) Unzip Microsoft building footprints (keep EPSG:4326 GeoJSON)
#    Input files named like: Alabama.geojson.zip, NewYork.geojson.zip
# ------------------------------------------------------------------
#bldg_zips <- dir_ls(in_bldg, glob = "*.geojson.zip", recurse = FALSE)

#for (z in bldg_zips) {
#  nm <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(z))) # "Alabama"
#  dst_folder <- file.path(out_bldg, nm)
#  dir_create(dst_folder)
#  unzip(z, exdir = dst_folder)
#  message("✓ Buildings: unzipped ", basename(z), " -> ", dst_folder)
#}

# ------------------------------------------------------------------
# 3) Unzip 2020 block boundaries per state (keep shapefiles)
#    Your files look like: nhgis0010_shapefile_tl2020_010_block_2020.zip
#    We'll just expand each into its own subfolder under processed/block_boundaries
# ------------------------------------------------------------------
blk_zips <- dir_ls(in_blocks, glob = "*block_2020.zip", recurse = FALSE)

for (z in blk_zips) {
  # pull 3-digit state code from filename if present
  fips3 <- str_match(basename(z), "tl2020_(\\d{3})_block_2020\\.zip")[,2]
  tag <- ifelse(is.na(fips3), tools::file_path_sans_ext(basename(z)), fips3)
  dst_folder <- file.path(out_blocks, tag)
  dir_create(dst_folder)
  unzip(z, exdir = dst_folder)
  message("✓ Blocks: unzipped ", basename(z), " -> ", dst_folder)
}

# ------------------------------------------------------------------
# 4) Clean NHGIS 2020 block-level population CSV
#    Keep: GISJOIN, STATE, COUNTY, TRACT, U7H001 (rename to pop)
#    If a block column exists (e.g., BLOCKA or BLOCK), we’ll keep it too.
#    Save to processed/tables/block_pop_2020_min.csv
# ------------------------------------------------------------------
# Guess the file based on your screenshot name
cand <- c(
  file.path(in_census, "nhgis0009_ds258_2020_block.csv"),
  file.path(in_census, "nhgis0009_ds258_2020_blck_grp.csv")
)
nhgis_csv <- cand[file_exists(cand)][1]
stopifnot("Can't find NHGIS CSV in '2020 census' folder." = length(nhgis_csv) == 1)

DT <- fread(nhgis_csv, nThread = max(1, parallel::detectCores() - 1))

# normalize column names
cn <- names(DT)

# helper to grab first column matching choices
pick_first <- function(choices, cols) {
  hit <- choices[choices %in% cols]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# most NHGIS block files include: GISJOIN, STATEA, COUNTYA, TRACTA, BLOCKA, U7H001
col_GISJOIN <- pick_first(c("GISJOIN","GISJOIN2"), cn)
col_STATE   <- pick_first(c("STATE","STATEA"), cn)
col_COUNTY  <- pick_first(c("COUNTY","COUNTYA"), cn)
col_TRACT   <- pick_first(c("TRACT","TRACTA"), cn)
col_BLOCK   <- pick_first(c("BLOCK","BLOCKA"), cn)
col_POP     <- pick_first(c("U7H001"), cn)

need <- c(col_GISJOIN, col_STATE, col_COUNTY, col_TRACT, col_POP)
if (any(is.na(need))) {
  stop("Required columns not found. Saw columns: ", paste(cn, collapse = ", "))
}

keep_cols <- c(col_GISJOIN, col_STATE, col_COUNTY, col_TRACT, col_POP)
if (!is.na(col_BLOCK)) keep_cols <- c(keep_cols, col_BLOCK)

DT_min <- DT[, ..keep_cols]

# rename to requested simple names
setnames(DT_min, old = c(col_GISJOIN, col_STATE, col_COUNTY, col_TRACT, col_POP),
         new = c("GISJOIN",   "STATE",   "COUNTY",   "TRACT",   "pop"))
if (!is.na(col_BLOCK)) setnames(DT_min, old = col_BLOCK, new = "BLOCK")

# write cleaned table
out_csv <- file.path(out_tables, "block_pop_2020_min.csv")
fwrite(DT_min, out_csv)
message("✓ Wrote cleaned NHGIS table -> ", out_csv)

cat("\nAll done.\n\nProcessed items:\n",
    "- processed/boundaries_1990/ (copied shapefiles)\n",
    "- processed/boundaries_2020/ (copied shapefiles)\n",
    "- processed/buildings/<State>/*.geojson (EPSG:4326, unzipped)\n",
    "- processed/block_boundaries/<tag>/*.shp (unzipped per state)\n",
    "- processed/tables/block_pop_2020_min.csv (GISJOIN, STATE, COUNTY, TRACT, pop[, BLOCK])\n")
