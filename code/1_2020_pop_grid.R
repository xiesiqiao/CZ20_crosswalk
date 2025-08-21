# ============================================
# Build 100 m Population Grid — All States (2020)
# Output: GeoTIFF only, area-weighted by building footprint
# CRS: EPSG:5070 (USA_Contiguous_Albers)
# Memory-conscious: two-pass over building files + frequent gc()
# ============================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(data.table)
  library(lwgeom)
  library(terra)
})

sf::sf_use_s2(FALSE) # planar ops only (we're in Albers)

# ------------------------------
# Project paths
# ------------------------------
base_input <- "C:/Users/xiesi/ASU Dropbox/Xie Sijiao/Connor-Kemeny-Storper/Microdata_wealth_estimates/Crosswalk/CZ20_CZ90_xwalk/input"
base_root  <- dirname(base_input)
out_dir    <- file.path(base_root, "output", "population_grid")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Inputs
pop_csv      <- file.path(base_input, "processed", "tables", "block_pop_2020_min.csv")
blocks_root  <- file.path(base_input, "processed", "block_boundaries")
bld_root     <- file.path(base_input, "processed", "buildings")

# Settings
crs_target   <- 5070
cell_size_m  <- 100
overwrite_tif <- FALSE       # set TRUE to overwrite existing outputs
keep_tmp     <- FALSE        # set TRUE to keep per-state temp RDS summaries

# ------------------------------
# State lookup (lower-48 + DC)
# ------------------------------
lookup <- tribble(
  ~name,           ~abbr, ~fips,
  "Alabama","AL","01","Arizona","AZ","04","Arkansas","AR","05","California","CA","06",
  "Colorado","CO","08","Connecticut","CT","09","Delaware","DE","10","District of Columbia","DC","11",
  "Florida","FL","12","Georgia","GA","13","Idaho","ID","16","Illinois","IL","17","Indiana","IN","18",
  "Iowa","IA","19","Kansas","KS","20","Kentucky","KY","21","Louisiana","LA","22","Maine","ME","23",
  "Maryland","MD","24","Massachusetts","MA","25","Michigan","MI","26","Minnesota","MN","27","Mississippi","MS","28",
  "Missouri","MO","29","Montana","MT","30","Nebraska","NE","31","Nevada","NV","32","New Hampshire","NH","33",
  "New Jersey","NJ","34","New Mexico","NM","35","New York","NY","36","North Carolina","NC","37","North Dakota","ND","38",
  "Ohio","OH","39","Oklahoma","OK","40","Oregon","OR","41","Pennsylvania","PA","42","Rhode Island","RI","44",
  "South Carolina","SC","45","South Dakota","SD","46","Tennessee","TN","47","Texas","TX","48","Utah","UT","49",
  "Vermont","VT","50","Virginia","VA","51","Washington","WA","53","West Virginia","WV","54","Wisconsin","WI","55",
  "Wyoming","WY","56"
) %>%
  mutate(fips3 = paste0(fips,'0'))

# If you want to restrict which states to run, edit here:
states_to_run <- lookup$name

# Also restrict to states that actually exist on disk (both blocks + buildings)
have_states <- function() {
  bld_dirs   <- list.dirs(bld_root, full.names = FALSE, recursive = FALSE)
  blk_dirs   <- list.dirs(blocks_root, full.names = FALSE, recursive = FALSE)
  # match by name+fips3 presence
  lookup %>%
    filter(name %in% bld_dirs, fips3 %in% blk_dirs) %>%
    pull(name)
}
states_to_run <- intersect(states_to_run, have_states())

# ------------------------------
# Helpers
# ------------------------------
read_blocks_for_state <- function(fips3) {
  shp_dir <- file.path(blocks_root, fips3)
  shps <- list.files(shp_dir, pattern = "\\.shp$", full.names = TRUE)
  if (!length(shps)) stop("No block shapefiles found in: ", shp_dir)
  blocks <- map(shps, \(p) st_read(p, quiet = TRUE)) %>% list_rbind()
  if (!"GISJOIN" %in% names(blocks))
    stop("Block shapefile missing GISJOIN in ", shp_dir)
  blocks %>%
    select(GISJOIN) %>%
    st_make_valid() %>%
    st_transform(crs_target)
}

filter_pop_to_state <- function(pop_all, state_row) {
  # Robust filter:
  # 1) If 'STATE' exists and looks like name/abbr/fips, try those.
  # 2) Otherwise use GISJOIN prefix (G + 0-padded fips).
  fips2 <- state_row$fips
  nm    <- state_row$name
  abbr  <- state_row$abbr
  if ("STATE" %in% names(pop_all)) {
    stcol <- pop_all$STATE
    # normalize character
    st_chr <- if (is.character(stcol)) stcol else as.character(stcol)
    st_chr <- trimws(st_chr)
    keep <- st_chr %in% c(nm, abbr, fips2, as.character(as.integer(fips2)))
    pop_s <- pop_all[keep, c("GISJOIN","pop")]
    if (nrow(pop_s)) return(pop_s)
  }
  # Fallback: GISJOIN prefix
  patt <- paste0("^G0*", as.integer(fips2))
  pop_all %>%
    filter(str_detect(GISJOIN, patt)) %>%
    select(GISJOIN, pop)
}

# Merge-accumulate helper for small keyed tables (sum by key)
accumulate_sum <- function(dt_acc, dt_in, key_col = "GISJOIN", val_col = "area_sum") {
  if (nrow(dt_in) == 0) return(dt_acc)
  if (nrow(dt_acc) == 0) return(copy(dt_in))
  setkeyv(dt_acc, key_col); setkeyv(dt_in, key_col)
  out <- merge(dt_acc, dt_in, by = key_col, all = TRUE, suffixes = c(".x",".y"))
  out[, (val_col) := rowSums(.SD, na.rm = TRUE), .SDcols = paste0(val_col, c(".x",".y"))]
  out[, c(paste0(val_col, c(".x",".y"))) := NULL]
  out[]
}

accumulate_grid <- function(grid_acc, chunk_dt) {
  if (nrow(chunk_dt) == 0) return(grid_acc)
  if (nrow(grid_acc) == 0) return(copy(chunk_dt))
  setkey(grid_acc, ix, iy); setkey(chunk_dt, ix, iy)
  out <- merge(grid_acc, chunk_dt, by = c("ix","iy"), all = TRUE, suffixes = c(".x",".y"))
  out[, pop := rowSums(.SD, na.rm = TRUE), .SDcols = c("pop.x","pop.y")]
  out[, c("pop.x","pop.y") := NULL]
  out[]
}

make_tif_from_grid <- function(grid_sum, x0, y0, cell_size, tif_path, crs_wkt) {
  if (!nrow(grid_sum)) {
    warning("Empty grid; nothing to write at ", tif_path)
    return(invisible(NULL))
  }
  min_ix <- min(grid_sum$ix); max_ix <- max(grid_sum$ix)
  min_iy <- min(grid_sum$iy); max_iy <- max(grid_sum$iy)
  nx <- max_ix - min_ix + 1
  ny <- max_iy - min_iy + 1
  xmin <- x0 + min_ix * cell_size
  xmax <- x0 + (max_ix + 1) * cell_size
  ymin <- y0 + min_iy * cell_size
  ymax <- y0 + (max_iy + 1) * cell_size
  
  r <- terra::rast(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                   ncols = nx, nrows = ny, crs = crs_wkt)
  # compute cell centers and assign values
  cx <- x0 + (grid_sum$ix + 0.5) * cell_size
  cy <- y0 + (grid_sum$iy + 0.5) * cell_size
  cells <- terra::cellFromXY(r, cbind(cx, cy))
  # initialize values lazily, then assign
  r[cells] <- grid_sum$pop
  terra::writeRaster(r, tif_path, overwrite = TRUE,
                     wopt = list(datatype = "FLT4S", gdal = "COMPRESS=LZW"))
  invisible(r)
}

# ------------------------------
# Core per-state builder (two-pass)
# ------------------------------
build_state_tif <- function(state_row, pop_all) {
  state_name  <- state_row$name
  state_abbr  <- state_row$abbr
  fips2       <- state_row$fips
  fips3       <- state_row$fips3
  
  message("\n==============================")
  message("State: ", state_name, " (", state_abbr, ", FIPS ", fips2, ")")
  message("==============================")
  
  tif_out <- file.path(out_dir, paste0(state_abbr, "_grid100m_pop_2020_area.tif"))
  if (!overwrite_tif && file.exists(tif_out)) {
    message("Output exists; skipping: ", tif_out)
    return(invisible(NULL))
  }
  
  # Paths
  buildings_dir <- file.path(bld_root, state_name)
  if (!dir.exists(buildings_dir)) {
    message("No buildings dir: ", buildings_dir, " — skipping.")
    return(invisible(NULL))
  }
  
  # 1) Read population (filtered) + blocks (geometry)
  message("Reading population table…")
  pop_s <- filter_pop_to_state(pop_all, state_row) %>%
    mutate(pop = as.numeric(pop))
  if (!nrow(pop_s)) {
    message("No rows for population in ", state_name, " — skipping.")
    return(invisible(NULL))
  }
  
  message("Reading block boundaries…")
  blocks <- read_blocks_for_state(fips3) # GISJOIN + geometry in 5070
  
  message("Joining blocks + population…")
  blocks_pop <- left_join(blocks, pop_s, by = "GISJOIN")
  rm(blocks); gc()
  
  # Grid origin aligned to 100 m using state blocks bbox
  bb <- st_bbox(blocks_pop)
  x0 <- floor(bb["xmin"] / cell_size_m) * cell_size_m
  y0 <- floor(bb["ymin"] / cell_size_m) * cell_size_m
  
  # List building files (GeoJSON) — process in chunks
  bld_files <- list.files(buildings_dir, pattern = "\\.geojson$", full.names = TRUE)
  if (!length(bld_files)) {
    message("No building geojson files in ", buildings_dir, " — skipping.")
    return(invisible(NULL))
  }
  
  # -------- Pass 1: per-block total footprint area --------
  message("Pass 1/2: summing footprint area per block…")
  dt_area <- data.table(GISJOIN = character(), area_sum = numeric())
  
  for (fp in bld_files) {
    message("  [+] ", basename(fp))
    b <- st_read(fp, quiet = TRUE) %>% st_make_valid() %>% st_transform(crs_target)
    b$footprint_m2 <- as.numeric(st_area(b))
    cent <- st_centroid(b)
    
    # only need GISJOIN in the join target
    j <- st_join(cent, blocks_pop %>% select(GISJOIN), join = st_within, left = FALSE)
    if (!nrow(j)) { rm(b, cent, j); gc(); next }
    
    # align areas from polygons to joined centroids
    idx <- st_nearest_feature(j, b)
    tmp <- data.table(GISJOIN = j$GISJOIN, footprint_m2 = b$footprint_m2[idx])
    tmp <- tmp[!is.na(GISJOIN) & is.finite(footprint_m2) & footprint_m2 > 0]
    tmp_agg <- tmp[, .(area_sum = sum(footprint_m2)), by = GISJOIN]
    
    dt_area <- accumulate_sum(dt_area, tmp_agg, key_col = "GISJOIN", val_col = "area_sum")
    rm(b, cent, j, idx, tmp, tmp_agg); gc()
  }
  
  if (!nrow(dt_area)) {
    message("No buildings landed inside blocks in ", state_name, " — skipping.")
    return(invisible(NULL))
  }
  
  # Attach area_sum to blocks_pop (used in Pass 2)
  blocks_pop2 <- left_join(blocks_pop %>% select(GISJOIN, pop), as.data.frame(dt_area), by = "GISJOIN")
  rm(blocks_pop, dt_area); gc()
  
  # -------- Pass 2: allocate pop to buildings and bin to grid --------
  message("Pass 2/2: allocating and binning to 100 m grid…")
  grid_acc <- data.table(ix = integer(), iy = integer(), pop = numeric())
  
  for (fp in bld_files) {
    message("  [*] ", basename(fp))
    b <- st_read(fp, quiet = TRUE) %>% st_make_valid() %>% st_transform(crs_target)
    b$footprint_m2 <- as.numeric(st_area(b))
    cent <- st_centroid(b)
    
    j <- st_join(cent, blocks_pop2, join = st_within, left = FALSE) %>%
      filter(!is.na(area_sum), area_sum > 0, pop > 0)
    if (!nrow(j)) { rm(b, cent, j); gc(); next }
    
    # align polygon areas to joined centroids
    idx <- st_nearest_feature(j, b)
    areas <- b$footprint_m2[idx]
    
    # area-weighted allocation
    pop_alloc <- (areas / j$area_sum) * j$pop
    keep <- is.finite(pop_alloc) & pop_alloc > 0
    if (!any(keep)) { rm(b, cent, j, idx, areas, pop_alloc); gc(); next }
    
    coords <- st_coordinates(j[keep, ])
    ix <- floor((coords[,1] - x0) / cell_size_m)
    iy <- floor((coords[,2] - y0) / cell_size_m)
    
    chunk <- data.table(ix = ix, iy = iy, pop = pop_alloc[keep])[, .(pop = sum(pop)), by = .(ix, iy)]
    grid_acc <- accumulate_grid(grid_acc, chunk)
    
    rm(b, cent, j, idx, areas, pop_alloc, coords, ix, iy, chunk); gc()
  }
  
  if (!nrow(grid_acc)) {
    message("No allocated population produced in ", state_name, " — skipping.")
    return(invisible(NULL))
  }
  
  # -------- Write GeoTIFF directly from sparse grid --------
  message("Writing GeoTIFF…")
  crs_wkt <- st_crs(crs_target)$wkt
  make_tif_from_grid(grid_acc, x0, y0, cell_size_m, tif_out, crs_wkt)
  
  # -------- Quick closure report --------
  tot_block_pop <- sum(pop_s$pop, na.rm = TRUE)
  tot_alloc_pop <- grid_acc[, sum(pop, na.rm = TRUE)]
  message("Summary for ", state_abbr, ":")
  message("  Blocks total pop:   ", format(round(tot_block_pop), big.mark = ","))
  message("  Raster summed pop:  ", format(round(tot_alloc_pop), big.mark = ","))
  message("  Difference:         ", format(round(tot_block_pop - tot_alloc_pop), big.mark = ","))

}

# ------------------------------
# Load population once; loop states
# ------------------------------
message("Loading population CSV once…")
pop_all <- suppressMessages(readr::read_csv(pop_csv, show_col_types = FALSE))
stopifnot(all(c("GISJOIN","pop") %in% names(pop_all)))

message("States to run (found on disk): ", paste(states_to_run, collapse = ", "))

for (nm in states_to_run) {
  st_row <- lookup %>% filter(name == nm) %>% slice(1)
  try(build_state_tif(st_row, pop_all), silent = FALSE)
  gc()
}

message("\nAll done. TIFs in: ", out_dir)
gc()
