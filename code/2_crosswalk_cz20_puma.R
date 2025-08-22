# ==========================================================
# National crosswalks: CZ20 -> {CZ90, PUMA20, PUMA2010, CZ2010, CZ2000}
# Uses 100 m population GeoTIFFs (EPSG:5070) across ALL states.
# Output: one national CSV per target geography.
# ==========================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(data.table); library(readr); library(stringr); library(terra); library(tibble)
})
sf::sf_use_s2(FALSE)

# ---- paths ----
base_root <- "C:/Users/xiesi/ASU Dropbox/Xie Sijiao/Connor-Kemeny-Storper/Microdata_wealth_estimates/Crosswalk/CZ20_CZ90_xwalk"
input_dir <- file.path(base_root, "input")
dir_1990  <- file.path(input_dir, "processed", "boundaries_1990")
dir_2020  <- file.path(input_dir, "processed", "boundaries_2020")
dir_2010  <- file.path(input_dir, "processed", "boundaries_2010")
dir_cz10  <- file.path(input_dir, "2010 boundaries")
dir_cz00  <- file.path(input_dir, "2000 boundaries")

grid_dir  <- file.path(base_root, "output", "population_grid")
out_dir   <- file.path(base_root, "output", "crosswalks"); dir.create(out_dir, TRUE, TRUE)

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

# ---- final table builder (no 'state' column from raster) ----
finalize_table <- function(pairs_dt, totals_dt, from_map, to_map, to_names, target_label) {
  if (!nrow(pairs_dt)) return(data.frame())
  from_lu <- from_map |> rename(CZ20 = id, code_from = code)
  to_lu   <- to_map   |> left_join(to_names, by = "id") |> 
    rename(!!target_label := id, !!paste0(target_label,"_name") := name, code_to = code)
  
  res <- pairs_dt |>
    left_join(totals_dt, by="code_from") |>
    left_join(from_lu,   by="code_from") |>
    left_join(to_lu,     by="code_to") |>
    transmute(
      CZ20,
      !!target_label := .data[[target_label]],
      !!paste0(target_label,"_name") := .data[[paste0(target_label,"_name")]],
      pop_overlap = pop_overlap,
      share_of_CZ20 = ifelse(pop_from > 0, pop_overlap / pop_from, 0),
      afactor       = share_of_CZ20,
      CZ20_total_pop = pop_from,
      cells = cells,
      grid_res_m = !!cell_m,
      crs_epsg   = !!epsg,
      method     = !!method_tag
    ) |>
    as.data.frame()
  
  # Add PUMA state_fips/state when the target is PUMA*
  if (target_label %in% c("PUMA20","PUMA2010")) {
    fips_to_usps <- c(
      "01"="AL","02"="AK","04"="AZ","05"="AR","06"="CA","08"="CO","09"="CT","10"="DE","11"="DC",
      "12"="FL","13"="GA","15"="HI","16"="ID","17"="IL","18"="IN","19"="IA","20"="KS","21"="KY",
      "22"="LA","23"="ME","24"="MD","25"="MA","26"="MI","27"="MN","28"="MS","29"="MO","30"="MT",
      "31"="NE","32"="NV","33"="NH","34"="NJ","35"="NM","36"="NY","37"="NC","38"="ND","39"="OH",
      "40"="OK","41"="OR","42"="PA","44"="RI","45"="SC","46"="SD","47"="TN","48"="TX","49"="UT",
      "50"="VT","51"="VA","53"="WA","54"="WV","55"="WI","56"="WY","60"="AS","66"="GU","69"="MP",
      "72"="PR","78"="VI"
    )
    tid <- res[[target_label]]
    state_fips <- stringr::str_sub(tid, 2, 3)
    res$state_fips <- state_fips
    res$state      <- unname(fips_to_usps[state_fips])
  }
  res
}

# ---- read polygons, add codes ----
cz20   <- read_std(file.path(dir_2020, "cz20.shp"),                 "CZ20")
cz90   <- read_std(file.path(dir_1990, "cz.shp"),                   "cz",      "cz_name")
puma20 <- read_std(file.path(dir_2020, "ipums_puma_2020.shp"),      "GISJOIN", "Name")
puma10 <- read_std(file.path(dir_2010, "ipums_puma_2010_tl20.shp"), "GISJOIN", "Name")
ers10  <- read_std(file.path(dir_cz10, "ERS10.rep.shp"),            "LM_Code")
ers00  <- read_std(file.path(dir_cz00, "ERS00.shp"),                "LM_Code")

cz20_map   <- make_code_map(cz20);   cz20   <- left_join(cz20,   cz20_map,   by="id")
cz90_map   <- make_code_map(cz90);   cz90   <- left_join(cz90,   cz90_map,   by="id")
puma20_map <- make_code_map(puma20); puma20 <- left_join(puma20, puma20_map, by="id")
puma10_map <- make_code_map(puma10); puma10 <- left_join(puma10, puma10_map, by="id")
ers10_map  <- make_code_map(ers10);  ers10  <- left_join(ers10,  ers10_map,  by="id")
ers00_map  <- make_code_map(ers00);  ers00  <- left_join(ers00,  ers00_map,  by="id")

nm_cz90 <- st_drop_geometry(cz90)   |> select(id, name)
nm_p20  <- st_drop_geometry(puma20) |> select(id, name)
nm_p10  <- st_drop_geometry(puma10) |> select(id, name)
nm_e10  <- st_drop_geometry(ers10)  |> select(id, name)
nm_e00  <- st_drop_geometry(ers00)  |> select(id, name)

cz20_for_rast   <- cz20   |> select(code, geometry)
cz90_for_rast   <- cz90   |> select(code, geometry)
puma20_for_rast <- puma20 |> select(code, geometry)
puma10_for_rast <- puma10 |> select(code, geometry)
ers10_for_rast  <- ers10  |> select(code, geometry)
ers00_for_rast  <- ers00  |> select(code, geometry)

# ---- NATIONAL run over ALL state TIFs ----
tifs <- list.files(grid_dir, pattern = "^[A-Z]{2}_grid100m_pop_2020_.*\\.tif$", full.names = TRUE)
stopifnot(length(tifs) > 0)

acc_totals <- data.table(code_from=integer(), pop_from=numeric())  # national CZ20 totals on grid
acc <- list(
  CZ90     = data.table(),
  PUMA20   = data.table(),
  PUMA2010 = data.table(),
  CZ2010   = data.table(),
  CZ2000   = data.table()
)

for (ti in tifs) {
  message("\n----- ", basename(ti), " -----")
  r <- terra::rast(ti) |> ensure_crs(epsg)
  bb <- bbox_sf(r)
  
  # source = CZ20 inside raster bbox
  from_s <- suppressWarnings(st_intersection(cz20_for_rast, bb))
  if (nrow(from_s) == 0L) { gc(); next }
  
  # totals for this raster
  t_res <- pairs_for_target(r, from_s, from_s)$totals
  if (nrow(t_res)) {
    if (nrow(acc_totals)) {
      acc_totals <- merge(acc_totals, t_res, by="code_from", all=TRUE, suffixes=c(".x",".y"))
      acc_totals[, pop_from := rowSums(.SD, na.rm=TRUE), .SDcols=c("pop_from.x","pop_from.y")]
      acc_totals[, c("pop_from.x","pop_from.y") := NULL]
    } else acc_totals <- t_res
  }
  
  # accumulate pairs for each target
  do_one <- function(target_for_rast, key) {
    t_s <- suppressWarnings(st_intersection(target_for_rast, bb))
    if (nrow(t_s)==0L) return(invisible(NULL))
    res <- pairs_for_target(r, from_s, t_s)$pairs
    if (nrow(res)) acc[[key]] <<- rbindlist(list(acc[[key]], res), use.names=TRUE, fill=TRUE)
    invisible(NULL)
  }
  
  do_one(cz90_for_rast,   "CZ90")
  do_one(puma20_for_rast, "PUMA20")
  do_one(puma10_for_rast, "PUMA2010")
  do_one(ers10_for_rast,  "CZ2010")
  do_one(ers00_for_rast,  "CZ2000")
  
  gc()
}

# collapse duplicates across rasters (national)
collapse_nat <- function(dt) {
  if (!nrow(dt)) return(dt)
  dt[, .(cells = sum(cells), pop_overlap = sum(pop_overlap)), by = .(code_from, code_to)]
}
acc <- lapply(acc, collapse_nat)

# ---- finalize & write (one CSV per geography) ----
final_write <- function(key, target_label, to_map, to_names, filename) {
  if (!nrow(acc[[key]])) { message("No rows for ", key); return(invisible(NULL)) }
  
  out <- finalize_table(
    pairs_dt  = acc[[key]],
    totals_dt = acc_totals,
    from_map  = cz20_map,
    to_map    = to_map,
    to_names  = to_names,
    target_label = target_label
  )
  
  # Optional QA: afactor sums for each CZ20
  qa <- out |>
    group_by(CZ20) |>
    summarise(afactor_sum = sum(afactor, na.rm=TRUE), .groups="drop") |>
    filter(abs(afactor_sum - 1) > 1e-6)
  if (nrow(qa)) message("WARNING (", target_label, "): ", nrow(qa), " CZ20(s) not summing to 1.")
  
  readr::write_csv(out, file.path(out_dir, filename))
  message("Wrote: ", filename, " (", nrow(out), " rows)")
  
}

final_write("CZ90",     "CZ90",     cz90_map,   nm_cz90, "cz20_to_cz90_popgrid100m_area.csv")
final_write("PUMA20",   "PUMA20",   puma20_map, nm_p20,  "cz20_to_puma20_popgrid100m_area.csv")
final_write("PUMA2010", "PUMA2010", puma10_map, nm_p10,  "cz20_to_puma2010_popgrid100m_area.csv")
final_write("CZ2010",   "CZ2010",   ers10_map,  nm_e10,  "cz20_to_cz2010_popgrid100m_area.csv")
final_write("CZ2000",   "CZ2000",   ers00_map,  nm_e00,  "cz20_to_cz2000_popgrid100m_area.csv")

# (Optional) national CZ20 totals on the grid for QA
readr::write_csv(as.data.frame(acc_totals), file.path(out_dir, "cz20_grid_totals_NATIONAL.csv"))

message("\nDone. Crosswalks in: ", out_dir)
