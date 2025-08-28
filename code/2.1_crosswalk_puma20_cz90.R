# ==========================================================
# National Crosswalk: PUMA2020 -> CZ1990  (all states)
# Uses 100 m population GeoTIFFs (EPSG:5070).
# Output: one national CSV; afactor sums ~ 1 per PUMA.
# ==========================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(data.table); library(readr); library(stringr); library(terra); library(tibble)
})
sf::sf_use_s2(FALSE)

# ---- paths ----
project_root <- "C:/Users/xiesi/ASU Dropbox/Xie Sijiao/Connor-Kemeny-Storper/Microdata_wealth_estimates/Crosswalk/CZ20_CZ90_xwalk"
input_dir    <- file.path(project_root, "input")
dir_1990     <- file.path(input_dir, "processed", "boundaries_1990")
dir_2020     <- file.path(input_dir, "processed", "boundaries_2020")

grid_dir <- file.path(project_root, "output", "population_grid")
out_dir  <- file.path(project_root, "output", "crosswalks"); dir.create(out_dir, TRUE, TRUE)

# ---- settings ----
epsg <- 5070
cell_m <- 100
method_tag <- "building_footprint_area"

# ---- helpers ----
read_std <- function(fp, id_col, name_col = NULL) {
  g <- st_read(fp, quiet = TRUE) |> st_make_valid() |> st_transform(epsg)
  g |> mutate(id = as.character(.data[[id_col]]),
              name = if (!is.null(name_col) && name_col %in% names(g)) as.character(.data[[name_col]]) else NA_character_) |>
    select(id, name, geometry)
}
make_code_map <- function(sfobj) data.frame(id = unique(sfobj$id), code = seq_along(unique(sfobj$id)))
bbox_sf <- function(r) sf::st_as_sf(terra::as.polygons(terra::ext(r), crs = terra::crs(r)))
safe_blocks <- function(x){ bs <- try(terra::blocks(x), TRUE); if(inherits(bs,"try-error")||is.null(dim(bs))||nrow(bs)==0L) data.frame(row=1L,nrows=terra::nrow(x)) else bs }
ensure_crs <- function(r, epsg){ if (terra::crs(r, proj=TRUE)=="") terra::crs(r) <- st_crs(epsg)$wkt; r }
safe_read_values <- function(s,row,nrows){ v <- terra::readValues(s,row=row,nrows=nrows); if(is.null(v)) return(NULL); if(is.matrix(v)) return(v); nl <- terra::nlyr(s); if(nl==0L) return(NULL); matrix(v, ncol=nl) }

pairs_for_target <- function(r_pop, from_for_rast, to_for_rast) {
  r_pop <- ensure_crs(r_pop, epsg)
  bb <- bbox_sf(r_pop)
  from_s <- suppressWarnings(st_intersection(from_for_rast, bb)) |> select(code, geometry)
  to_s   <- suppressWarnings(st_intersection(to_for_rast,   bb)) |> select(code, geometry)
  if (nrow(from_s)==0L || nrow(to_s)==0L) return(list(pairs=data.table(), totals=data.table()))
  
  r_from <- terra::rasterize(terra::vect(from_s), r_pop, field="code", background=NA)
  r_to   <- terra::rasterize(terra::vect(to_s),   r_pop, field="code", background=NA)
  
  s123 <- c(r_pop, r_from, r_to); on.exit(try(terra::readStop(s123), silent=TRUE), add=TRUE); terra::readStart(s123)
  bs <- safe_blocks(s123)
  pairs <- data.table(code_from=integer(), code_to=integer(), cells=integer(), pop_overlap=numeric())
  totals<- data.table(code_from=integer(), pop_from=numeric())
  
  for (i in seq_len(nrow(bs))) {
    m <- safe_read_values(s123, bs$row[i], bs$nrows[i]); if (is.null(m)) next
    pop <- m[,1]; cf <- m[,2]; ct <- m[,3]
    
    okf <- !is.na(pop) & !is.na(cf) & pop>0
    if (any(okf)) {
      t1 <- data.table(code_from=cf[okf], pop_from=pop[okf])[, .(pop_from=sum(pop_from)), by=code_from]
      if (nrow(totals)) {
        totals <- merge(totals, t1, by="code_from", all=TRUE, suffixes=c(".x",".y"))
        totals[, pop_from := rowSums(.SD, na.rm=TRUE), .SDcols=c("pop_from.x","pop_from.y")]
        totals[, c("pop_from.x","pop_from.y") := NULL]
      } else totals <- t1
    }
    
    ok <- !is.na(pop) & !is.na(cf) & !is.na(ct) & pop>0
    if (any(ok)) {
      p1 <- data.table(code_from=cf[ok], code_to=ct[ok], cells=1L, pop_overlap=pop[ok])[
        , .(cells=sum(cells), pop_overlap=sum(pop_overlap)), by=.(code_from, code_to)]
      if (nrow(pairs)) {
        pairs <- merge(pairs, p1, by=c("code_from","code_to"), all=TRUE, suffixes=c(".x",".y"))
        pairs[, cells := rowSums(.SD, na.rm=TRUE), .SDcols=c("cells.x","cells.y")]
        pairs[, pop_overlap := rowSums(.SD, na.rm=TRUE), .SDcols=c("pop_overlap.x","pop_overlap.y")]
        pairs[, c("cells.x","cells.y","pop_overlap.x","pop_overlap.y") := NULL]
      } else pairs <- p1
    }
    if (i %% 5 == 0) gc()
  }
  list(pairs=pairs, totals=totals)
}

# Build final table with requested column names for PUMA20 -> CZ90
finalize_p20_to_cz90 <- function(pairs_dt, totals_dt, from_map, to_map, to_names) {
  if (!nrow(pairs_dt)) return(data.frame())
  from_lu <- from_map |> rename(PUMA20 = id, code_from = code)
  to_lu   <- to_map   |> left_join(to_names, by = "id") |>
    rename(CZ90 = id, CZ90_name = name, code_to = code)
  
  pairs_dt |>
    left_join(totals_dt, by="code_from") |>
    left_join(from_lu,   by="code_from") |>
    left_join(to_lu,     by="code_to") |>
    transmute(
      PUMA20,
      PUMA20_name = NA_character_,
      CZ90,
      CZ90_name,
      pop_overlap = pop_overlap,
      share_of_PUMA20 = ifelse(pop_from > 0, pop_overlap / pop_from, 0),
      PUMA20_total_pop = pop_from,
      cells = cells,
      grid_res_m = !!cell_m,
      crs_epsg   = !!epsg,
      method     = !!method_tag
    ) |>
    as.data.frame()
}

# tiny FIPS -> USPS mapper (for convenience in output)
fips_to_usps <- c(
  "01"="AL","02"="AK","04"="AZ","05"="AR","06"="CA","08"="CO","09"="CT","10"="DE","11"="DC",
  "12"="FL","13"="GA","15"="HI","16"="ID","17"="IL","18"="IN","19"="IA","20"="KS","21"="KY",
  "22"="LA","23"="ME","24"="MD","25"="MA","26"="MI","27"="MN","28"="MS","29"="MO","30"="MT",
  "31"="NE","32"="NV","33"="NH","34"="NJ","35"="NM","36"="NY","37"="NC","38"="ND","39"="OH",
  "40"="OK","41"="OR","42"="PA","44"="RI","45"="SC","46"="SD","47"="TN","48"="TX","49"="UT",
  "50"="VT","51"="VA","53"="WA","54"="WV","55"="WI","56"="WY","60"="AS","66"="GU","69"="MP",
  "72"="PR","78"="VI"
)

# ---- read polygons, add integer codes ----
# TO: CZ1990
cz90 <- read_std(file.path(dir_1990, "cz.shp"), "cz", "cz_name")
cz90_map <- make_code_map(cz90); cz90 <- left_join(cz90, cz90_map, by="id")
nm_cz90 <- st_drop_geometry(cz90) |> select(id, name)

# FROM: PUMA2020
puma20 <- read_std(file.path(dir_2020, "ipums_puma_2020.shp"), "GISJOIN", "Name")
puma20_map <- make_code_map(puma20); puma20 <- left_join(puma20, puma20_map, by="id")
nm_p20 <- st_drop_geometry(puma20) |> select(id, name)

# narrow to 'code' for rasterization
cz90_for_rast   <- cz90   |> select(code, geometry)
puma20_for_rast <- puma20 |> select(code, geometry)

# ---- read ALL state population rasters (national run) ----
tifs <- list.files(grid_dir, pattern = "^[A-Z]{2}_grid100m_pop_2020_.*\\.tif$", full.names = TRUE)
stopifnot(length(tifs) > 0)

# national totals for PUMA20 on the grid
acc_totals <- data.table(code_from=integer(), pop_from=numeric())

# accumulator for pairs (no state tagging here)
acc_pairs <- data.table(code_from=integer(), code_to=integer(), cells=integer(), pop_overlap=numeric())

for (ti in tifs) {
  message("\n----- ", basename(ti), " -----")
  r <- terra::rast(ti) |> ensure_crs(epsg)
  
  # Crop FROM to raster bbox; includes any PUMA intersecting this raster
  from_s <- suppressWarnings(st_intersection(puma20_for_rast, bbox_sf(r)))
  if (nrow(from_s) == 0L) { gc(); next }
  
  # totals for this raster (PUMA20)
  t_res <- pairs_for_target(r, from_s, from_s)$totals
  if (nrow(t_res)) {
    if (nrow(acc_totals)) {
      acc_totals <- merge(acc_totals, t_res, by="code_from", all=TRUE, suffixes=c(".x",".y"))
      acc_totals[, pop_from := rowSums(.SD, na.rm=TRUE), .SDcols=c("pop_from.x","pop_from.y")]
      acc_totals[, c("pop_from.x","pop_from.y") := NULL]
    } else acc_totals <- t_res
  }
  
  # pairs to CZ90 (national CZs)
  to_s <- suppressWarnings(st_intersection(cz90_for_rast, bbox_sf(r)))
  if (nrow(to_s) == 0L) { gc(); next }
  res <- pairs_for_target(r, from_s, to_s)$pairs
  if (nrow(res)) {
    acc_pairs <- rbindlist(list(acc_pairs, res), use.names=TRUE, fill=TRUE)
  }
  gc()
}

# collapse duplicates across rasters
if (nrow(acc_pairs)) {
  acc_pairs <- acc_pairs[, .(cells = sum(cells), pop_overlap = sum(pop_overlap)), by = .(code_from, code_to)]
}

# ---- finalize & write (national) ----
if (!nrow(acc_pairs)) {
  message("No rows produced for PUMA20 -> CZ90 (national run)")
} else {
  out <- finalize_p20_to_cz90(
    pairs_dt  = acc_pairs,
    totals_dt = acc_totals,
    from_map  = puma20_map,
    to_map    = cz90_map,
    to_names  = nm_cz90
  )
  
  # drop placeholder name then attach PUMA names (optional)
  out$PUMA20_name <- NULL
  out <- out |>
    left_join(nm_p20 |> rename(PUMA20 = id, PUMA20_name = name), by = "PUMA20")
  
  # --- add state_fips (positions 2–3 of PUMA20) and USPS state ---
  out <- out |>
    mutate(
      state_fips = stringr::str_sub(PUMA20, 2, 3),
      state      = fips_to_usps[ state_fips ] |> as.character()
    )
  
  # convenience: afactor alias
  out$afactor <- out$share_of_PUMA20
  
  # sanity: verify each PUMA sums to ~1 nationally
  sums <- out |>
    group_by(PUMA20) |>
    summarise(afactor_sum = sum(afactor, na.rm=TRUE), .groups = "drop")
  bad <- sums |> filter(abs(afactor_sum - 1) > 1e-6)
  if (nrow(bad)) {
    message("WARNING: ", nrow(bad), " PUMA(s) not summing to 1. If any remain, inspect those PUMA IDs or widen inputs.")
  }
  
  fname <- "puma20_to_cz90_popgrid100m_area.csv"
  readr::write_csv(out, file.path(out_dir, fname))
  message("Wrote: ", fname, " (", nrow(out), " rows)")
}

# Also write national PUMA20 totals on the grid for QA
readr::write_csv(as.data.frame(acc_totals), file.path(out_dir, "puma20_grid_totals_NATIONAL.csv"))

message("\nDone. Crosswalk in: ", out_dir)
