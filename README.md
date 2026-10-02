# README for 'A state-based landscape of fear: Habitat function, season, and reproductive status define moose response to wolf predation risk'

## Overview
This repository contains the code and data used to conduct an integrated step-selection analysis (iSSA) of female moose habitat selection in response to wolf predation risk in and around Yukon-Charley Rivers National Preserve, interior Alaska (2020–2023 biological years). This analysis models moose movement through a step selection function including habitat covariates, kernel densities of wolf distribution, and calf status. Model selection uses a tiered AIC framework to identify the most supported model structure separately for summer (May–August) and winter (November–February).

Author: Derek A. Arnold (darnold10@alaska.edu)

Moose and wolf GPS collar data were collected by Matthew D. Cameron, Mathew S. Sorum, and Kyle Joly (National Park Service, Yukon-Charley Rivers National Preserve). Code was written by Derek A. Arnold.

## How to Use
The repository is structured as two scripts and a separately sourced function file. `final_code.R` assembles the analysis (steps, available steps, wolf risk surfaces, habitat covariates) and saves it as `data/mse_ssf_full2.rds`. `tiered_model_selection_final.R` reads that cached object and fits the candidate model sets. The file paths noted as `path` in the scripts should be the only things a user needs to modify.

Raster layer paths in `final_code.R` are placeholders (`"path/..."`) and will need to be set to local versions of the layers listed below. 


## Layout

### **`SSF_functions.r`**
Functions sourced by `final_code.R`. 

 * **`format.ssf`**: converts a data frame of relocations for one individual into sequential steps. Returns a data frame with start and end steps as coordinates and times (`Lat_1`, `Lon_1`, `Lat_2`, `Lon_2`, `Time_1`, `Time_2`), step length (`SL`), turning angle (`TA`), `TimeDiff`. `CalfStatus` is carried through and `Type = TRUE` marks observed steps. 
 
 * **`gen.stps`**: generates available steps from each observed step. Fits a gamma distribution to observed step lengths and a von Mises distribution to observed turning angles. It draws a specified number of random steps per observed step, and returns observed and available steps bound together with `Type` (TRUE/FALSE) and a `Step` stratum identifier. 
 
 * **`anglefun`, `bearing.ta`, `bearing.sl`**: geometric helpers called internally by `format.ssf`.


### **`final_code.R`**
Data formatting and covariate construction.

 1. Reads wolf pack counts and derives the annual maximum observed count per pack over the biological year (1 May – 30 April) for use in the a pack-size weighting.
 2. Reads wolf and moose collar datasets, standardizes pack names, and drops lone wolves and dispersers. Moose living outside the wolf collar coverage area are excluded by `AnimalId`.
 3. Adds calf observations to the moose dataset. A moose-year is assigned as `CalfStatus = 1` from 1 May through the last date a calf was confirmed present, and 0 otherwise. Calf loss is therefore dated to the last confirmed sighting to include only known calf status periods.
 4. Thins to one wolf location per pack for each timestep to prevent pseudo-replication of packs density, and subsets both species to the analysis months (May–August, November–February), shoulder months are dropped.
 5. Calculates four wolf density surfaces for each biological year, each is normalized to one. Density surfaces are a full-year KDE (`wolf_kde_full`), season-specific KDEs (combined into `wolf_kde_seasonal`), an unweighted sum of per-pack KDEs (`wolf_kde_pack`), and a pack-size-weighted sum of per-pack KDEs (`wolf_kde_weighted`).
 6. Generates 15 available steps per observed step for each individual and extracts wolf density values at the step endpoint.
 7. Calculates distance-to-wetland surface from the ABoVE landcover classes Water, Fen, and Shallows to match other continuous covariates. Then extracts terrain and vegetation covariates at the step endpoint.
 8. Centers and scales continuous covariates across used and available steps jointly, so coefficients are on a common scale and are comparatively interpretable. Saves the result as `data/mse_ssf_full2.rds`.
 9. Builds the seasonal and calf-status subsets used for model fitting in following script

### **`tiered_model_selection_final.R`**
Candidate model sets and AIC comparisons. All models are conditional logistic regressions fit through the Poisson reformulation. Summer models are fit to `mse_ssf_sum_cs` and winter models to `mse_ssf_win_cs`.


### **`data`**
The underlying raw GPS relocations for moose and wolves and associated code are not yet provided. Release of precise locations of individually marked animals in a population subject to harvest is restricted under federal and Alaska state law, and the National Park Service withholds such information. Data structures are described below for those used in the analysis.

### **`model_outputs`**
Destination for fitted model objects. Created by the user before running `tiered_model_selection_final.R`; the script writes `t1a_*.rds`, `t1b_*.rds`, `t1c_*.rds`, `t2_*.rds`, and `t3_*.rds` here, each suffixed `_sum` or `_win`.

## Data objects

### **`yuch_wolf_dat2.csv`**
Wolf GPS relocations, 253,263 records from 78 collared wolves, June 2019 – 2025, with 35 distinct group labels (including dispersers and unknown-affiliation records that are filtered out in `final_code.R`).

AnimalId - individual wolf identifier
FixDate - date of the fix in UTC (m/d/Y H:M)
LocalDateTime - fix timestamp in AK local time (m/d/Y H:M); this is the field used in the analysis
Year - calendar year of the fix
OrdinalDate - day of calendar year
UnitCode - NPS unit code 
Species - Wolf
Gender - Female or Male
GroupName - pack affiliation. Values containing ">" denote dispersal from one pack to another; "Single", "unknown", "unkown", and dispersal records are recoded to "Lone" and dropped
Lat_WGS84, Lon_WGS84 - latitude and longitude in decimal degrees (EPSG:4326)

### **`yuch_moose_dat.csv`** 
`final_code.R` drops nine individuals that ranged outside the wolf collar coverage area (2009, 2016, 2019, 2020, 2021, 2022, 2023, 2103, 2108), and drops three more (2302, 2303, 2304) from a single biological year for where they dispersed outside the range. Further subsetting to the analysis seasons and to fixes with known calf status yields the 28 individuals reported in the manuscript.

FixId - unique record identifier
ProjectId - project code (YUCH_Moose)
AnimalId - individual moose identifier, matching `ID` in `calf_status.csv`
FixDate - fix timestamp in UTC (Y-m-d H:M:S)
LocalDateTime - fix timestamp in AK local time (Y-m-d H:M:S)
Year - calendar year of the fix
OrdinalDate - day of calendar year
UnitCode - NPS unit code
Species - Moose
Gender - Female for all records
Lat_WGS84, Lon_WGS84 - latitude and longitude in decimal degrees (EPSG:4326)


### **`calf_status.csv`**
Calf observation records from aerial and ground surveys, 163 observations of 35 collared cow moose, May 2020 – March 2023.

ID - individual moose identifier, matching `AnimalId` in the moose collar data
Date - observation date (m/d/yy)
calv_num - number of calves observed with the cow (0, 1, 2, or NA). Any positive value marks the cow as accompanied; 0 is treated as no calf

### **`wlf_pck_cts_20-25.csv`**
Highest observed count per pack per month, used to derive pack-size weights. 52 pack-year records across winter years 2020–2024. The first two rows below the header are a field-description row and a blank spacer and are dropped on import.

Pack ID - pack name, recoded on import to match `GroupName` in the collar data
Winter.Year - first calendar year of the biological year (a value of 2020 denotes 1 May 2020 – 30 April 2021)
Bio year - biological year 
PrevYr - binary indicator of whether the pack was collared in the previous year
May, Aug, Sept, Oct, Nov, Dec, Jan, Feb, March, April - highest observed count within the given month. The maximum across these columns within a pack-year is used as the pack-size weight, and packs with no count in a year default to a weight of 1


## Spatial covariate layers
The following rasters are read from local paths in `final_code.R` and are not included here.

 * Continuous foliar cover, percent cover: `salshr` (willow), `picgla` (white spruce), `picmar` (black spruce), `dectre` (deciduous trees), `bettre`, `alnus`, `sphagn`, `wetsed`. Only willow, white spruce, and deciduous cover enter the reported models
 * Topography: elevation, ruggedness, slope, aspect. Only ruggedness (with a quadratic term) is retained in the reported models
 * Distance to nearest river: `NWBR_Rivers_dist.tif`.
 * ABoVE annual landcover, 2018: used to derive distance to wetland from the Water, Fen, and Shallows classes. 

## Software
Analyses were run in R with `tidyverse`, `glmmTMB`, `AICcmodavg`, `terra`, `sf`, `MASS`, `CircStats`, `lubridate`


