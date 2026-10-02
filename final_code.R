# =====================================================================
# Wolf-moose interaction analysis: iSSF data preparation
# Yukon-Charley Rivers (NPS) moose + wolf GPS collar data
# =====================================================================

# ---------------------------------------------------------------------
# PACKAGES AND PATHS
# ---------------------------------------------------------------------
pkgs = c('tidyverse', 'conflicted', 'geosphere', 'tictoc', 'rlang', 'glmmTMB',
         'survival', 'sf', 'terra', 'AICcmodavg', 'MASS')
# Resolve the dplyr/stats::filter clash up front; assigned into globalenv so it
# outranks anything attached below.
filter <- dplyr::filter
# Loading required packages
sapply(pkgs, require, character = TRUE)

here::here() %>% setwd()

# --- Vegetation class rasters (percent cover from continuous foliar cover layers) ---
salshr = rast("path/salshr.tif")
picmar = rast("path/picmar.tif")
dectre = rast("path/dectre.tif")

# --- Terrain ---
elev = rast("path/Topography/Elevation.tif")
rugg = rast("path/Topography/Ruggedness.tif")
slpe = rast("path/Topography/Slope.tif")
aspt = rast("/path/Topography/Aspect.tif")

# --- Distance-to-nearest-river surface ---
riv_dist=rast("path/NWBR_Rivers_dist.tif") 

# --- ABoVE annual landcover, used for wetland distance ---
ldcv = rast("patha/annual_landcover_2018.tif")

# terra::flip masks other flip() methods
conflicts_prefer(terra::flip)

# format.ssf() and gen.stps() are created in this file
source("SSF_functions.r")

# ---------------------------------------------------------------------
# WOLF PACK COUNTS -> ANNUAL MAX PACK SIZE
# ---------------------------------------------------------------------
# Rows 1-2 are header/spacer junk in the raw sheet.
pck_dat <- read.csv('path/data/wlf_pck_cts_20-25.csv')[-c(1,2),]

# Month columns run May-April to match the biological year.
bio_months <- c("May","Jun","July","Aug","Sept","Oct","Nov","Dec","Jan","Feb","March","April")
month_cols <- lubridate::intersect(bio_months, names(pck_dat))

# Standardize pack names, coerce monthly counts to numeric, and take the annual
# maximum count as the pack-size weight. Bio-year runs 1 May - 30 Apr.
pck_num <- pck_dat %>%
  rename(GroupName = Pack.ID) %>%
  mutate(
    GroupName = case_when(
      GroupName %in% c("Single", "unknown", "unkown") ~ "Lone",
      GroupName %in% c("Lost") ~ "Lost Creek",
      grepl(">", GroupName, fixed = TRUE)             ~ "Lone",  # dispersers / merged records
      TRUE                                            ~ GroupName
    )
  ) %>%
  mutate(across(all_of(month_cols),
                ~ suppressWarnings(as.numeric(na_if(na_if(as.character(.), "<NA>"), ""))))) %>%
  rowwise() %>%
  mutate(
    MaxPack = if (all(is.na(c_across(all_of(month_cols))))) NA_real_
    else max(c_across(all_of(month_cols)), na.rm = TRUE),
    StartDate = make_datetime(year = as.integer(Winter.Year), month = 5, day = 1, tz = "America/Anchorage"),
    EndDate   = {
      y2 <- as.integer(Winter.Year) + 1L
      make_datetime(year = y2, month = 4, day = days_in_month(make_date(y2, 4, 1)), tz = "America/Anchorage")
    }
  ) %>%
  ungroup() %>%
  dplyr::select(GroupName, StartDate, EndDate, MaxPack)


# ---------------------------------------------------------------------
# COLLAR DATA IMPORT AND CLEANING
# ---------------------------------------------------------------------
# Importing wolf data: round to the hour, standardise pack names, then drop lone/unknown
# animals - pack identity is required for the pack-level KDEs and random point selection.
wlf_dat <- read.csv('path/data/yuch_wolf_dat2.csv') %>%
  mutate(Time = mdy_hm(LocalDateTime)) %>%
  dplyr::select(-Year, -OrdinalDate, -UnitCode, -Species, -FixDate, -LocalDateTime) %>%
  rename(Lat = Lat_WGS84, Lon = Lon_WGS84, Sex = Gender) %>%
  mutate(Time = round_date(Time, unit = "hour")) %>%
  mutate(
    GroupName = case_when(
      GroupName %in% c("Single", "unknown", "unkown") ~ "Lone",
      grepl(">", GroupName, fixed = TRUE)             ~ "Lone",
      TRUE                                            ~ GroupName
    )
  ) %>% filter(GroupName != "Great Unknown") %>% filter(GroupName != "Lone")


# Importing moose data: same hourly rounding; animals living outside of wolf collar coverage are excluded
mse_dat = read.csv('path/data/yuch_moose_dat.csv') %>%
  mutate(Time = ymd_hms(LocalDateTime)) %>% dplyr::select(-c(ProjectId, OrdinalDate, UnitCode, Species, Shape, FixDate, GroupName, LocalDateTime)) %>% 
  rename(Lat = Lat_WGS84, Lon = Lon_WGS84) %>% mutate(Time = round_date(Time, unit = "hour")) %>% rename(Sex = Gender) %>% 
  dplyr::filter(!(AnimalId %in% c(2103, 2009, 2019, 2020, 2023, 2021, 2022, 2016, 2108)))

# ---------------------------------------------------------------------
# CALF STATUS INTEGRATION
# ---------------------------------------------------------------------
# Binary CalfStatus indicator associated every moose fix, with any uncertainty
# marked for filtering

# --- Read observation records and flag confirmed calves ---
calf_stat <- read.csv('path/data/calf_status.csv') %>%
  rename(AnimalId = ID) %>%
  mutate(
    Date = mdy(Date),
    BioYear = if_else(month(Date) < 5, year(Date) - 1, year(Date)),
    # Only a positive count counts as a calf; NA and 0 both mean "no calf".
    has_calf = if_else(!is.na(calv_num) & calv_num > 0, 1, 0)
  )

# --- One calf-present interval per moose-year ---
# Interval = 1 May through the last date a calf was confirmed. Moose-years with
# no confirmed sighting are treated as calf-free.
# calf loss is dated to the last confirmed sighting, so the interval is
# right-truncated by survey timing
calf_intervals <- calf_stat %>%
  group_by(AnimalId, BioYear) %>%
  summarize(
    LastConfirmedDate = if(any(has_calf == 1)) max(Date[has_calf == 1]) else as.Date(NA),
    .groups = "drop"
  ) %>%
  filter(!is.na(LastConfirmedDate)) %>%
  mutate(
    CalfStart = make_date(BioYear, 5, 1)
  )

# --- Join onto the GPS fixes ---
mse_dat <- mse_dat %>%
  mutate(
    DateOnly = date(Time),
    BioYear = if_else(month(DateOnly) < 5, year(DateOnly) - 1, year(DateOnly))
  ) %>%
  left_join(calf_intervals, by = c("AnimalId", "BioYear")) %>%
  mutate(
    CalfStatus = case_when(
      # Inside the confirmed interval -> 1; unmatched, pre-May, or post-loss -> 0.
      !is.na(CalfStart) & DateOnly >= CalfStart & DateOnly <= LastConfirmedDate ~ 1,
      TRUE ~ 0
    )
  ) %>%
  dplyr::select(-LastConfirmedDate, -CalfStart, -BioYear, -DateOnly)



# ---------------------------------------------------------------------
# TEMPORAL STRUCTURING AND SEASON SUBSETS
# ---------------------------------------------------------------------
# One wolf location per pack per hour. Prevents packs carrying several collars
# from dominating the KDE.
set.seed(22) # reproducibility
wlf_dat <- wlf_dat %>%
  group_by(Time, GroupName) %>%
  slice_sample(n = 1) %>%   
  ungroup()

# Keep Jan-Feb, May-Sep, Nov-Dec (drops shoulder months Mar, Apr, Oct).
# Bio-year starts 1 May, so Jan-Apr fixes belong to the previous bio-year.
wlf_tme <- wlf_dat %>%
  filter(month(Time) %in% c(1:2, 5:8, 11:12)) %>%
  mutate(
    bio_year = if_else(month(Time) <= 4, year(Time) - 1, year(Time))
  ) %>%
  group_by(bio_year) %>%
  nest() %>% filter(!(bio_year %in% c(2019, 2024, 2025)))   # incomplete collar coverage

mse_tme <- mse_dat %>%
  filter(month(Time) %in% c(1:2, 5:8, 11:12)) %>%
  mutate(
    bio_year = if_else(month(Time) <= 4, year(Time) - 1, year(Time))
  ) %>%
  group_by(bio_year) %>%
  nest() %>% filter(!(bio_year %in% c(2024, 2025)))


# ---------------------------------------------------------------------
# WOLF KDE RISK SURFACES
# ---------------------------------------------------------------------
# Study-area bounding box, WGS84: xmin, xmax, ymin, ymax.
lims <- c(-146.125065, -139.987652, 63.945575, 66.029961)

# Pooled KDE over all collared wolves in a year, normalised to sum to 1.
make_kde_full <- function(df, lims = lims, n = 500) {
  
  kde_out <- kde2d(df$Lon, df$Lat, n = n, lims = lims)
  # kde2d returns z[x, y]; transpose and flip to get raster row order.
  kde_rast <- rast(t(kde_out$z)) %>% 
    flip(direction = "vertical")
  
  ext(kde_rast) <- lims
  crs(kde_rast) <- "EPSG:4326"
  
  kde_rast <- kde_rast / global(kde_rast, "sum", na.rm = TRUE)[[1]]
  
  return(kde_rast)
}

# Annual pooled surface.
wlf_tme <- wlf_tme %>%
  mutate(kde_full = purrr::map(data, make_kde_full, lims = lims, n = 500))

# Season-specific pooled surfaces.
wlf_tme <- wlf_tme %>%
  mutate(
    kde_summer = map(data, ~.x %>% filter(month(Time) %in% c(5:8)) %>% make_kde_full(lims = lims, n = 500)),
    kde_winter = map(data, ~.x %>% filter(month(Time) %in% c(11, 12, 1, 2)) %>% make_kde_full(lims = lims, n = 500))
  )


# removing moose that have moved outside the range of the wolf data for one year
mse_tme$data[[4]] = filter(mse_tme$data[[4]], !(AnimalId %in% c(2304, 2303, 2302)))


# Common 500 x 500 template that all pack KDEs are projected onto.
tot_kde <- rast(
  xmin = lims[1], xmax = lims[2],
  ymin = lims[3], ymax = lims[4],
  nrows = 500, ncols = 500,
  crs = "EPSG:4326"
)

# Split each year's fixes by pack.
wlf_tme <- wlf_tme %>%
  mutate(pack_data = map(data, function(df) {
    df %>%
      group_by(GroupName) %>%
      nest() %>%
      ungroup() %>%
      rename(locations = data)
  }))

# Pack that is outside of the contiguous wolf coverage zone
wlf_tme$pack_data[[4]] = filter(wlf_tme$pack_data[[4]], GroupName != "Tangle Lakes")

# Per-pack KDE, snapped to the shared template and NA-filled with 0 so packs
# can be summed. Normalised to sum to 1 within pack.
make_kde_pack <- function(d) {
  kd <- MASS::kde2d(d$Lon, d$Lat, n = 500, lims = lims)
  
  r <- terra::rast(t(kd$z))
  r <- terra::flip(r, direction = "vertical")
  terra::ext(r) <- terra::ext(lims[1], lims[2], lims[3], lims[4])
  terra::crs(r) <- "EPSG:4326"
  r <- terra::project(r, tot_kde, method = "bilinear")
  r <- terra::ifel(is.na(r), 0, r)
  
  r / terra::global(r, "sum", na.rm = TRUE)[[1]]
}

wlf_tme <- wlf_tme %>%
  mutate(pack_data = map(pack_data, function(pack_df) {
    pack_df %>%
      mutate(kde = map(locations, make_kde_pack))
  }))

# Unweighted pack surface: every pack contributes equally regardless of size.
wlf_tme <- wlf_tme %>%
  mutate(annual_kde = map(pack_data, function(pack_df) {
    Reduce(`+`, pack_df$kde)
  }))

# --- Pack-size weighted surface ---
pck_sizes <- pck_num %>%
  mutate(bio_year = year(StartDate)) %>%
  dplyr::select(GroupName, bio_year, MaxPack) %>%
  filter(!is.na(MaxPack))

# Each pack KDE is scaled by its annual maximum count, summed, and renormalised
# to sum to 1. Packs with no count record default to weight = 1.
wlf_tme <- wlf_tme %>%
  mutate(weighted_annual_kde = map2(pack_data, bio_year, function(pd, yr) {
    
    pd_weighted <- pd %>%
      left_join(
        pck_sizes %>% filter(bio_year == yr) %>% dplyr::select(GroupName, MaxPack),
        by = "GroupName"
      ) %>%
      mutate(MaxPack = if_else(is.na(MaxPack), 1, MaxPack))
    
    weighted_kdes <- map2(pd_weighted$kde, pd_weighted$MaxPack, function(k, w) k * w)
    combined <- Reduce(`+`, weighted_kdes)
    
    combined / terra::global(combined, "sum", na.rm = TRUE)[[1]]
  }))

# ---------------------------------------------------------------------
# USED AND AVAILABLE MOOSE STEPS
# ---------------------------------------------------------------------
# Per animal: derive step lengths and turning angles, then draw 15 available
# steps per used step from the fitted step-length / turn-angle distributions.
mse_tme <- mse_tme %>%
  ungroup() %>%
  mutate(animal_data = map(data, function(df) {
    df %>%
      group_by(AnimalId) %>%
      nest() %>%
      ungroup() %>%
      
      mutate(Length = map_int(data, nrow)) %>%
      
      mutate(stps = map(data, format.ssf)) %>%
      
      mutate(r.stps = map(stps, function(d) {
        gen.stps(d, 15)
      }))
  }))

# Flatten back to one row per (used or available) step, keyed by bio_year.
mse_tme <- mse_tme %>%
  mutate(ssf_data = map(animal_data, function(df) {
    df %>%
      dplyr::select(AnimalId, r.stps) %>%
      na.omit() %>%
      unnest(r.stps)
  })) %>%
  dplyr::select(bio_year, ssf_data)

# ---------------------------------------------------------------------
# EXTRACT WOLF DENSITY AT STEP ENDPOINTS
# ---------------------------------------------------------------------
# Extraction happens while coordinates are still WGS84, matching the KDE CRS.
# Endpoint (Lon_2, Lat_2) is the relevant location for an iSSF.
mse_tme <- mse_tme %>%
  ungroup() %>%
  mutate(ssf_data = map2(ssf_data, bio_year, function(df, yr) {
    wolf_full <- wlf_tme$kde_full[wlf_tme$bio_year == yr][[1]]
    wolf_annual <- wlf_tme$annual_kde[wlf_tme$bio_year == yr][[1]]
    wolf_weighted <- wlf_tme$weighted_annual_kde[wlf_tme$bio_year == yr][[1]]
    wolf_summer <- wlf_tme$kde_summer[wlf_tme$bio_year == yr][[1]]
    wolf_winter <- wlf_tme$kde_winter[wlf_tme$bio_year == yr][[1]]
    
    df <- df %>%
      mutate(
        wolf_kde_full = terra::extract(wolf_full, cbind(Lon_2, Lat_2))[[1]],
        wolf_kde_pack = terra::extract(wolf_annual, cbind(Lon_2, Lat_2))[[1]],
        wolf_kde_weighted = terra::extract(wolf_weighted, cbind(Lon_2, Lat_2))[[1]],
        wolf_kde_summer = terra::extract(wolf_summer, cbind(Lon_2, Lat_2))[[1]],
        wolf_kde_winter = terra::extract(wolf_winter, cbind(Lon_2, Lat_2))[[1]],
        
        month = month(Time_2),
        wolf_kde_seasonal = case_when(
          month %in% c(11, 12, 1, 2) ~ wolf_kde_winter,
          month %in% c(5, 6, 7, 8) ~ wolf_kde_summer,
          TRUE ~ NA_real_
        )
      ) %>%
      dplyr::select(-month, -wolf_kde_summer, -wolf_kde_winter)
    
    df
  }))

# ---------------------------------------------------------------------
# PROJECT TO ALASKA ALBERS AND EXTRACT HABITAT COVARIATES
# ---------------------------------------------------------------------
mse_ssf <- mse_tme %>%
  na.omit() %>%
  unnest(ssf_data) %>%
  as.data.frame()

# Reproject step start (Lon_1, Lat_1) to EPSG:3338 in place.
mse_ssf <- mse_ssf %>%
  st_as_sf(coords = c("Lon_1","Lat_1"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = 3338) %>%
  mutate(
    Lon_1 = st_coordinates(.)[, 1],
    Lat_1  = st_coordinates(.)[, 2]) %>%
  st_drop_geometry()

# Same for step end (Lon_2, Lat_2). After this point Lon/Lat hold metres, not
# degrees - all habitat rasters below must also be EPSG:3338, because
# terra::extract() on a coordinate matrix does no CRS conversion.
mse_ssf <- mse_ssf %>%
  st_as_sf(coords = c("Lon_2","Lat_2"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = 3338) %>%
  mutate(
    Lon_2 = st_coordinates(.)[, 1],
    Lat_2  = st_coordinates(.)[, 2]) %>%
  st_drop_geometry()

# Recode zero step lengths to lowest observed step lengths
mse_ssf[which(mse_ssf$SL == 0),]$SL = 3.605903e-12

# --- Distance to wetland, derived from ABoVE landcover ---
pts <- vect(mse_ssf, geom = c("Lon_1", "Lat_1"), crs = "EPSG:3338")

pts_buffer <- terra::buffer(pts, width = 20000)

# Crop to the buffered extent first; distance() over the full ABoVE tile is slow.
r_subset <- crop(ldcv, pts_buffer)

# Classes 15, 11, 13 = Water, Fen, Shallows (see the recode below).
wetland_ids <- c(15, 11, 13)
r_wetland_binary <- ifel(r_subset %in% wetland_ids, 1, NA)

# calculateiong of distance to nearest wetland
r_dist <- distance(r_wetland_binary, unit="m")


# --- Extract all habitat covariates at the step endpoint ---
coords <- as.matrix(mse_ssf[, c("Lon_2", "Lat_2")])

mse_ssf$rugg <- terra::extract(rugg, coords)[,1]
mse_ssf$riv_dst <- terra::extract(riv_dist, coords)[,1]
mse_ssf$salshr <- terra::extract(salshr, coords)[,1]
mse_ssf$picmar <- terra::extract(picmar, coords)[,1]
mse_ssf$dectre <- terra::extract(dectre, coords)[,1]
mse_ssf$elev <- terra::extract(elev, coords)[,1]
mse_ssf$wtlnd <- terra::extract(r_dist, coords)[,1]
mse_ssf$cats <- terra::extract(r_subset, coords)[,1]

# Categorical habitat label. Unlisted codes become NA.
mse_ssf = mse_ssf %>% mutate(habitat = recode(
  cats, '15' = 'Water', '11' = 'Fen', '1' = 'Conif', '2' = 'Decid', '4' = 'Tall Shrub', '5' = 'Low Shrub', '3' = 'Mixed', '6' = 'Tall Shrub',
  '7' = 'Low Shrub', '8' = 'Bare', '9' = 'Bare', '10' = 'Bare',
  '12' = 'Bog', '13' = 'Shallows', '14' = 'Bare'))

# ---------------------------------------------------------------------
# STANDARDISE COVARIATES AND MOVEMENT TERMS
# ---------------------------------------------------------------------
# All predictors centred and scaled across used + available steps together, so
# coefficients are on a common SD scale.
mse_ssf = mse_ssf %>%
  mutate(.,salshr_sc=salshr %>% scale() %>% .[,1]) %>%
  mutate(.,dectre_sc=dectre %>% scale() %>% .[,1]) %>%
  mutate(.,wtlnd_sc=wtlnd %>% scale() %>% .[,1]) %>%
  mutate(.,picmar_sc=picmar %>% scale() %>% .[,1]) %>%
  mutate(.,rugg_sc=rugg %>% scale() %>% .[,1]) %>%
  mutate(.,riv_dst_sc=riv_dst %>% scale() %>% .[,1]) %>%
  # --- Wolf risk covariates: seasonal, pooled, per-pack, pack-size weighted ---
  mutate(.,wlf_kde_sea_sc=wolf_kde_seasonal %>% scale() %>% .[,1]) %>%
  mutate(.,wlf_kde_full_sc=wolf_kde_full %>% scale() %>% .[,1]) %>%
  mutate(.,wolf_kde_pack_sc=wolf_kde_pack %>% scale() %>% .[,1]) %>%
  mutate(.,wlf_kde_wt_sc = wolf_kde_weighted %>% scale() %>% .[,1]) %>%
  # --- Movement terms; SL is converted km -> m ---
  mutate(.,sl_sc=(SL*1000) %>% scale() %>% .[,1]) %>%
  mutate(.,log_sl=log(SL*1000) %>% scale() %>% .[,1]) %>%
  mutate(.,cos_ta=TA %>% cos())

# ---------------------------------------------------------------------
# LOAD CACHED FRAME AND BUILD ANALYSIS SUBSETS
# ---------------------------------------------------------------------
# saveRDS(mse_ssf, file = "data/mse_ssf_full2.rds")

# --- Seasonal x calf-status subsets ---
# Omitting uncertain calf statuses and separating data by season and calf status
mse_ssf_sum_cf = mse_ssf %>% dplyr::filter(month(Time_2) %in% c(5,6,7,8)) %>% na.omit() %>% filter(CalfStatus == 1)
mse_ssf_sum_nc = mse_ssf %>% dplyr::filter(month(Time_2) %in% c(5,6,7,8)) %>% na.omit() %>% filter(CalfStatus == 0)
mse_ssf_win_cf = mse_ssf %>% dplyr::filter(month(Time_2) %in% c(11,12,1,2)) %>% na.omit() %>% filter(CalfStatus == 1)
mse_ssf_win_nc = mse_ssf %>% dplyr::filter(month(Time_2) %in% c(11,12,1,2)) %>% na.omit() %>% filter(CalfStatus == 0)

# Recombined frames for models that estimate a calf-status interaction rather
# than fitting the two groups separately.
mse_ssf_sum_cs = rbind(mse_ssf_sum_cf, mse_ssf_sum_nc)
mse_ssf_win_cs = rbind(mse_ssf_win_cf, mse_ssf_win_nc)

