# get required libraries
library(foreign)
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(exactextractr)

# read functions
source("functions.R")

# get some data

# read in woody vegetation extent data for 2011
woody <- rast(paste(getwd(), "/input/woody_cover/woody_nsw_2011.tif", sep = "")) %>% round()
woody_mask <- woody %>% classify(cbind(c(0, 1), c(NA, 1)))

# read in properties
props <- vect("input/properties/Property_EPSG4283.gdb", layer = "property") %>% project("EPSG:3308") %>% as_sf()

# do some data preprocessing

# read in nvr data, reclassify for treatment (category 1 and category 2 regulated) and control (category 2 - vulnerable and sensitive) then aggregate, reproject, and snap to woody cover layer
nvr_incl <- rast(paste(getwd(), "/input/nvr_mapping/naluma_nsw_2017_abel0_c20221212_u9.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(NA, 1), c(1, NA))) %>% project(woody)
gc()
# get treatment layer (1 = treatment, 0 = control)
nvr_treat <- rast(paste(getwd(), "/input/nvr_mapping/naluma_nsw_2017_abkl0_c20221212_u9.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(NA, 3, 4, 5, 6), c(1, 0, 0, NA, 0))) %>% project(woody)
gc()
# get control layer (1 = control, 0 = treatment)
nvr_contr <- rast(paste(getwd(), "/input/nvr_mapping/naluma_nsw_2017_abkl0_c20221212_u9.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(NA, 3, 4, 5, 6), c(0, 1, 1, NA, 1))) %>% project(woody)
gc()

# generate include layer
# multiply by the woody layer as we only care about raster cells that were woody in 2011
include <- nvr_incl * woody_mask
# generate treatment and control layers
treatment <- include * nvr_treat
treatment <- treatment %>% round()
control <- include * nvr_contr
control <- control %>% round()

# free up memory
rm(nvr_incl, nvr_treat, nvr_contr)
gc()

# calculate whether properties have at least one raster cell that is subject to the LLS Act
props_lls <- include %>% exact_extract(props, "max") %>% as_tibble() %>% rename(include = value)

# bind zonal values to properties
props <- cbind(props, props_lls)

# remove properties that are not inside the LLS region and don't contain woody vegetation
props_incl <- props %>% filter(!is.na(include))

# remove duplicates (i.e., those RIDs with the same propid)
# only keep the RID field
duplicates <- duplicated(props_incl$propid)
props_incl <- props_incl %>% filter(!duplicates) %>% select(RID, propid)

# save properties to shapefile
st_write(props_incl, "input/analysis_data/props_incl.shp", append = FALSE)

# free up memory
rm(props_lls)
rm(props)
gc()

# read in koala habitat data and reclassify so 1 = koala habitat, 0 = non-habitat then aggregate, reproject, and snap to woody cover layer
# then aggregate, reproject, and snap to woody cover layer
khab <- rast(paste(getwd(), "/input/koala_habitat/KoalaHabitatSuitabilityModelClasses.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(1, 2, 3, 4, 5, 6), c(0, 0, 0, 1, 1, 1))) %>% project(woody)

# set up tibble to hold results of zonal analysis
zonal_woody_treat <- as_tibble(props_incl) %>% select(RID)
zonal_woody_contr <- as_tibble(props_incl) %>% select(RID)
zonal_koala_treat <- as_tibble(props_incl) %>% select(RID)
zonal_koala_contr <- as_tibble(props_incl) %>% select(RID)

# get zonal statistics for baseline woody vegetation extent and koala habitat
woody_treat <- woody * treatment
woody_contr <- woody * control
khab_treat <- khab * woody_treat
khab_contr <- khab * woody_contr
zonal_woody_treat <- bind_cols(zonal_woody_treat, woody_treat %>% exact_extract(props_incl, "sum") %>% ceiling() %>% as_tibble())
zonal_woody_contr <- bind_cols(zonal_woody_contr, woody_contr %>% exact_extract(props_incl, "sum") %>% ceiling() %>% as_tibble())
zonal_koala_treat <- bind_cols(zonal_koala_treat, khab_treat %>% exact_extract(props_incl, "sum") %>% ceiling() %>% as_tibble())
zonal_koala_contr <- bind_cols(zonal_koala_contr, khab_contr %>% exact_extract(props_incl, "sum") %>% ceiling() %>% as_tibble()) 

# create clearing file names
clearingfiles <- c("slats_2010_2011.tif", "slats_2011_2012.tif", "slats_2012_2013.tif", "slats_2013_2014.tif", "slats_2014_2015.tif", "slats_2015_2016.tif", "slats_2016_2017.tif", "slats_2017_2018.tif", "slats_2018_2019.tif", "slats_2019_2020.tif", "slats_2020_2021.tif")

for (i in clearingfiles) {
	clear <- rast(paste(getwd(), "/input/slats/", i, sep = "")) %>% ceiling() %>% classify(cbind(c(1, 2, 3, 4), c(0, 1, 0, 0))) %>% crop(woody)
	wtreat <- clear * woody_treat
	wcontr <- clear * woody_contr
	ktreat <- clear * khab_treat
	kcontr <- clear * khab_contr
	zonal_woody_treat <- bind_cols(zonal_woody_treat, wtreat %>% exact_extract(props_incl, "sum") %>% ceiling() %>% as_tibble())
	zonal_woody_contr <- bind_cols(zonal_woody_contr, wcontr %>% exact_extract(props_incl, "sum") %>% ceiling() %>% as_tibble())
	zonal_koala_treat <- bind_cols(zonal_koala_treat, ktreat %>% exact_extract(props_incl, "sum") %>% ceiling() %>% as_tibble())
	zonal_koala_contr <- bind_cols(zonal_koala_contr, kcontr %>% exact_extract(props_incl, "sum") %>% ceiling() %>% as_tibble()) 

	rm(clear, wtreat, wcontr, ktreat, kcontr)
	gc()
}

# set column names
names(zonal_woody_treat) <- c("RID", "baseline", "clear1011", "clear1112", "clear1213", "clear1314", "clear1415", "clear1516", "clear1617", "clear1718", "clear1819", "clear1920", "clear2021")
names(zonal_woody_contr) <- c("RID", "baseline", "clear1011", "clear1112", "clear1213", "clear1314", "clear1415", "clear1516", "clear1617", "clear1718", "clear1819", "clear1920", "clear2021")
names(zonal_koala_treat) <- c("RID", "baseline", "clear1011", "clear1112", "clear1213", "clear1314", "clear1415", "clear1516", "clear1617", "clear1718", "clear1819", "clear1920", "clear2021")
names(zonal_koala_contr) <- c("RID", "baseline", "clear1011", "clear1112", "clear1213", "clear1314", "clear1415", "clear1516", "clear1617", "clear1718", "clear1819", "clear1920", "clear2021")

# save outputs
write_rds(zonal_woody_treat, "input/analysis_data/woody_treat_zonal.rds")
write_rds(zonal_woody_contr, "input/analysis_data/woody_contr_zonal.rds")
write_rds(zonal_koala_treat, "input/analysis_data/koala_treat_zonal.rds")
write_rds(zonal_koala_contr, "input/analysis_data/koala_contr_zonal.rds")

# load processed data back in if needed - comment out if not needed
zonal_woody_treat <- read_rds("input/analysis_data/woody_treat_zonal.rds")
zonal_woody_contr  <- read_rds("input/analysis_data/woody_contr_zonal.rds")
zonal_koala_treat <- read_rds("input/analysis_data/koala_treat_zonal.rds")
zonal_koala_contr <- read_rds("input/analysis_data/koala_contr_zonal.rds")
props_incl <- vect("input/analysis_data/props_incl.shp")
woody <- rast(paste(getwd(), "/input/woody_cover/woody_nsw_2011.tif", sep = "")) %>% round()
khab <- rast(paste(getwd(), "/input/koala_habitat/KoalaHabitatSuitabilityModelClasses.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(1, 2, 3, 4, 5, 6), c(0, 0, 0, 1, 1, 1))) %>% project(woody)

# get only the koala habitat that is woody vegetation
khab_woody <- khab * woody

# calculate area, join to properties spatial layer, and then convert to a tibble
props_incl_spvec <- props_incl %>% vect()
props_incl_tibble <- props_incl %>% as_tibble() %>% mutate(AREA = expanse(props_incl_spvec, unit = "ha", transform = TRUE))
rm(props_incl_spvec)
gc()

# get land value and join to properties spatial layer
lval <- read.dbf("input/land_value/Prop_value_ha_20230101.dbf") %>% as_tibble()
props_incl_tibble <- props_incl_tibble %>% left_join(lval, by = c("propid" = "PROPERTY_I")) %>% rename(AREA = AREA.x, LVALHA = value_ha)
rm(lval)
gc()

# get slope
dem <- rast(paste(getwd(), "/input/dem/lf_dem1sec_noExclusions_noNegValues.tif",sep = ""))
slope <- terrain(dem, "slope") %>% project(woody)
rm(dem)
gc()
slope_zone <- slope %>% exact_extract(props_incl, "mean")
props_incl_tibble <- bind_cols(props_incl_tibble, as_tibble(slope_zone) %>% rename(SLOPE = value))
rm(slope)
gc()

# get soil fertility
soil_fert <- vect("input/soil/Fertility_NSW_v4_5_211020.shp") %>% project("EPSG:3308") %>% rasterize(woody, field = "Fert_code") %>% round()
soil_fert_zone <- soil_fert %>% exact_extract(props_incl, "mode")
props_incl_tibble <- bind_cols(props_incl_tibble, as_tibble(soil_fert_zone) %>% rename(SFERT = value))
rm(soil_fert)
gc()

# get land-use
land_use <- vect("input/landuse/NSWLanduse_2007shp.shp") %>% project("EPSG:3308") %>% rasterize(woody, field = "LU_NSWMajo")
land_use_zone <- land_use %>% exact_extract(props_incl, "mode")
props_incl_tibble <- bind_cols(props_incl_tibble, as_tibble(land_use_zone) %>% rename(LUSE = value))
rm(land_use)
gc()

# get proportion of property that is woody vegetation and koala habitat
count_woody_zone <- woody %>% exact_extract(props_incl, "count")
count_khab_zone <- khab_woody %>% exact_extract(props_incl, "count")
sum_woody_zone <- woody %>% exact_extract(props_incl, "sum")
sum_khab_zone <- khab_woody %>% exact_extract(props_incl, "sum")
prop_woody_zone <- sum_woody_zone / count_woody_zone
prop_khab_zone <- sum_khab_zone / count_khab_zone
props_incl_tibble <- bind_cols(props_incl_tibble, as_tibble(prop_woody_zone) %>% rename(PWOODY = value))
props_incl_tibble <- bind_cols(props_incl_tibble, as_tibble(prop_khab_zone) %>% rename(PKHAB = value))

props_incl_tibble_final <- props_incl_tibble %>% select(RID, propid, AREA, LVALHA, SLOPE, SFERT, LUSE, PWOODY, PKHAB)

# replace soil fertility values of 98 and 99 with NA values
props_incl_tibble_final <- props_incl_tibble_final %>% mutate(SFERT = ifelse(SFERT == 98 | SFERT == 99, NA, SFERT))

# replace koala habitat proportion values of NaN with NA values
props_incl_tibble_final <- props_incl_tibble_final %>% mutate(PKHAB = ifelse(is.na(PKHAB), NA, PKHAB))

# make SFERT and LUSE factors
props_incl_tibble_final <- props_incl_tibble_final %>% mutate(SFERT = as.factor(SFERT), LUSE = as.factor(LUSE))

# save output
write_rds(props_incl_tibble_final, "input/analysis_data/props_incl_table_covariates.rds")

# now get the data to do the matching
# matched based on property area, slope, land lavue, soil fertility, land used and woody vegetation/koala habitat percentage

# find RIDs needed for the matched samples
# woody
matchingRIDs_woody_treat <- zonal_woody_treat %>% select(RID, baseline) %>% rename(baseline_treat = baseline)
matchingRIDs_woody_contr <- zonal_woody_contr %>% select(RID, baseline) %>% rename(baseline_contr = baseline)
matchingRIDs_woody <- matchingRIDs_woody_treat %>% left_join(matchingRIDs_woody_contr, by = ("RID")) %>% filter((baseline_treat > 0 & baseline_contr == 0) | (baseline_treat == 0 & baseline_contr > 0)) %>% mutate(CI = ifelse(baseline_treat > 0 & baseline_contr == 0, 1, 0))
matchingRIDs_woody <- matchingRIDs_woody %>% mutate(baseline = baseline_treat + baseline_contr) %>% select(-baseline_treat, -baseline_contr)
# koala habitat
matchingRIDs_koala_treat <- zonal_koala_treat %>% select(RID, baseline) %>% rename(baseline_treat = baseline)
matchingRIDs_koala_contr <- zonal_koala_contr %>% select(RID, baseline) %>% rename(baseline_contr = baseline)
matchingRIDs_koala <- matchingRIDs_koala_treat %>% left_join(matchingRIDs_koala_contr, by = ("RID")) %>% filter((baseline_treat > 0 & baseline_contr == 0) | (baseline_treat == 0 & baseline_contr > 0)) %>% mutate(CI = ifelse(baseline_treat > 0 & baseline_contr == 0, 1, 0))
matchingRIDs_koala <- matchingRIDs_koala %>% mutate(baseline = baseline_treat + baseline_contr) %>% select(-baseline_treat, -baseline_contr)

# then join to the RIDs for matching
# and remove any rows with NA values
matchingRIDs_woody <- matchingRIDs_woody %>% left_join(props_incl_tibble_final, by = "RID") %>% select(-baseline, -propid, -PKHAB) %>% mutate(CI = as.integer(CI)) %>% drop_na()
matchingRIDs_koala <- matchingRIDs_koala %>% left_join(props_incl_tibble_final, by = "RID") %>% select(-baseline, -propid, PWOODY) %>% mutate(CI = as.integer(CI)) %>% drop_na()

# save matching pools
write_rds(matchingRIDs_woody, "input/analysis_data/matching_pool_woody.rds")
write_rds(matchingRIDs_koala, "input/analysis_data/matching_pool_koala.rds")
write.dbf(as.data.frame(matchingRIDs_woody), "input/analysis_data/matching_pool_woody.dbf")
write.dbf(as.data.frame(matchingRIDs_koala), "input/analysis_data/matching_pool_koala.dbf")

# BROOKE TO WRITE CODE HERE TO GENERATE THE DATA FOR INPUT INTO THE STATISTICAL MODELS

#These are the files:
# input/analysis_data/matched_mixed_RIDs.rds
# input/analysis_data/matched_mixed_RIDs_old.rds
# These are the RIDs for Courtney’s matched and mixed samples – with the CI variable also included in the match data (use this for the CI variable – already switched the 0s and 1s so not need to do that).
#WOODY MATCHED
woody_matched <- data.frame()
RID_repeated <- rep(matched_mixed_RIDs$woody_matched[,1], each = 11)
RID_rep <- unlist(RID_repeated, use.names = FALSE)
CI_repeated <- rep(matched_mixed_RIDs$woody_matched[,2], each = 11)
CI_rep <- unlist(CI_repeated, use.names = FALSE)

# Add columns RID and CI from matched_mixed_RIDs to woody_matched
woody_matched <- data.frame(
  RID = RID_rep,
  CI = CI_rep, 
  time = NA
)

# Order the dataframe by the RID column
woody_matched <- woody_matched[order(woody_matched$RID), ]
time_repeated <- rep(c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6), length(unique(woody_matched$RID)))
woody_matched$time <- time_repeated
woody_matched$baseline <- NA

not_in_zonal <- c()
# Iterate through each row using a for loop
for (i in 1:nrow(woody_matched)) {
  row <- woody_matched[i, ]  # Access current row
  # Check if RID exists in zonal_woody_treat
  if (row$RID %in% zonal_woody_treat$RID) {
    # Perform actions if RID exist
  # Check if CI column value is 1
  if (row$CI == 1) {
    row_treat <- zonal_woody_treat %>%
    filter(RID == row$RID)
    woody_matched[i, 4] <- row_treat$baseline 
    # Add your actions here
  } else {
    row_contr <- zonal_woody_contr %>%
    filter(RID == row$RID)
    woody_matched[i, 4] <- row_contr$baseline
  }
  # Perform operations on the row
  } else {
    #print(row$RID)
    # Perform actions if RID does not exist (optional)
    not_in_zonal <- c(not_in_zonal, row$RID)
  }
}
# Remove duplicates
unique_not_in_zonal <- unique(not_in_zonal)

# unique_not_in_zonal_new
# unique_not_in_zonal_old 
# # Find common values
# common_values <- intersect(unique_not_in_zonal_new, unique_not_in_zonal_old)
# # Print common values
# print(common_values)

# Remove rows where RID matches values in unique_not_in_zonal
woody_matched <- woody_matched %>%
  filter(!RID %in% unique_not_in_zonal)

#Add loss values
woody_matched$loss <- NA

#x <- 340
for (x in unique(woody_matched$RID)) {
  row <- woody_matched %>% filter(RID == x)
    if (row$CI[1] == 1) {
      subset_zonal_woody <- zonal_woody_treat %>% filter(RID == x)
      # Extract numerical values from the last 11 columns into a vector
      last_11_values <- as.vector(unlist(subset_zonal_woody[, -c(1:2)]))
      woody_matched$loss[woody_matched$RID == x] <- last_11_values
    } else {
      subset_zonal_contr <- zonal_woody_contr %>% filter(RID == x)
      # Extract numerical values from the last 11 columns into a vector
      last_11_values <- as.vector(unlist(subset_zonal_contr[, -c(1:2)]))
      woody_matched$loss[woody_matched$RID == x] <- last_11_values
    }
}

write_rds(woody_matched, "input/joined_data_frames/woody_matched_df.rds")

# #WOODY MIXED --> Courtney's RIDs
# woody_mixed <- data.frame()
# RID_repeated <- rep(matched_mixed_RIDs$woody_mixed[,1], each = 22)
# RID_rep <- unlist(RID_repeated, use.names = FALSE)
# # Create a vector with 11 ones followed by 11 zeros
# vector_ones_zeros <- c(rep(1, 11), rep(0, 11))
# CI_repeated <- rep(vector_ones_zeros, each = nrow(matched_mixed_RIDs$woody_mixed[,1]))
# CI_rep <- unlist(CI_repeated, use.names = FALSE)
# # Add columns RID and CI from matched_mixed_RIDs to woody_matched
# woody_mixed <- data.frame(
#   RID = RID_rep,
#   CI = CI_rep, 
#   time = NA
# )
# # Order the dataframe by the RID column
# woody_mixed <- woody_mixed[order(woody_mixed$RID), ]
# time_repeated <- rep(c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6), ((length(unique(woody_mixed$RID)))*2))
# woody_mixed$time <- time_repeated
# woody_mixed$baseline <- NA
# not_in_zonal <- c()
# # Iterate through each row using a for loop
# for (i in 1:nrow(woody_mixed)) {
#   row <- woody_mixed[i, ]  # Access current row
#   # Check if RID exists in zonal_woody_treat
#   if (row$RID %in% zonal_woody_treat$RID) {
#     # Perform actions if RID exist
#     # Check if CI column value is 1
#     if (row$CI == 1) {
#       row_treat <- zonal_woody_treat %>%
#         filter(RID == row$RID)
#       woody_mixed[i, 4] <- row_treat$baseline 
#       # Add your actions here
#     } else {
#       row_contr <- zonal_woody_contr %>%
#         filter(RID == row$RID)
#       woody_mixed[i, 4] <- row_contr$baseline
#     }
#     # Perform operations on the row
#   } else {
#     #print(row$RID)
#     # Perform actions if RID does not exist (optional)
#     not_in_zonal <- c(not_in_zonal, row$RID)
#   }
# }
# write_rds(woody_mixed, "input/joined_data_frames/woody_mixed_df_baseline_only.rds")
# # Remove duplicates
# unique_not_in_zonal <- unique(not_in_zonal)
# write_rds(unique_not_in_zonal, "input/joined_data_frames/unique_not_in_zonal_woody_mixed.rds")
# # Remove rows where RID matches values in unique_not_in_zonal
# woody_mixed <- woody_mixed %>%
#   filter(!RID %in% unique_not_in_zonal)
# #Add loss values
# woody_mixed$loss <- NA
# #x <- 167
# for (x in unique(woody_mixed$RID)) {
#   row <- woody_mixed %>% filter(RID == x)
#   subset_zonal_woody <- zonal_woody_treat %>% filter(RID == x)
#   # Extract numerical values from the last 11 columns into a vector
#   last_11_values_treat <- as.vector(unlist(subset_zonal_woody[, -c(1:2)]))
#   subset_zonal_contr <- zonal_woody_contr %>% filter(RID == x)
#   # Extract numerical values from the last 11 columns into a vector
#   last_11_values_contr <- as.vector(unlist(subset_zonal_contr[, -c(1:2)]))
#   joined <- c(last_11_values_treat, last_11_values_contr)
#   woody_mixed$loss[woody_mixed$RID == x] <- joined
#   }
# write_rds(woody_mixed, "input/joined_data_frames/woody_mixed_df_cm.rds")

#WOODY MIXED --> Jonathan's RIDs
woody_mixed <- data.frame()
RID_repeated <- rep(zonal_woody_treat$RID, each = 22)
RID_rep <- unlist(RID_repeated, use.names = FALSE)
# Create a vector with 11 ones followed by 11 zeros
vector_ones_zeros <- c(rep(1, 11), rep(0, 11))
CI_repeated <- rep(vector_ones_zeros, length(zonal_woody_treat$RID))
CI_rep <- unlist(CI_repeated, use.names = FALSE)
# Add columns RID and CI from matched_mixed_RIDs to woody_matched
woody_mixed <- data.frame(
  RID = RID_rep,
  CI = CI_rep, 
  time = NA
)
# Order the dataframe by the RID column
woody_mixed <- woody_mixed[order(woody_mixed$RID), ]
time_repeated <- rep(c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6), ((length(unique(woody_mixed$RID)))*2))
woody_mixed$time <- time_repeated
woody_mixed$baseline <- NA
not_in_zonal <- c()
# Iterate through each row using a for loop
for (i in 3000633:nrow(woody_mixed)) {
  row <- woody_mixed[i, ]  # Access current row
  # Check if RID exists in zonal_woody_treat
  if (row$RID %in% zonal_woody_treat$RID) {
    # Perform actions if RID exist
    # Check if CI column value is 1
    if (row$CI == 1) {
      row_treat <- zonal_woody_treat %>%
        filter(RID == row$RID)
      woody_mixed[i, 4] <- row_treat$baseline 
      # Add your actions here
    } else {
      row_contr <- zonal_woody_contr %>%
        filter(RID == row$RID)
      woody_mixed[i, 4] <- row_contr$baseline
    }
    # Perform operations on the row
  } else {
    #print(row$RID)
    # Perform actions if RID does not exist (optional)
    not_in_zonal <- c(not_in_zonal, row$RID)
  }
}
write_rds(woody_mixed, "input/joined_data_frames/woody_mixed_df_baseline_only.rds")
# Remove duplicates
unique_not_in_zonal <- unique(not_in_zonal)
write_rds(unique_not_in_zonal, "input/joined_data_frames/unique_not_in_zonal_woody_mixed.rds")
# Remove rows where RID matches values in unique_not_in_zonal
woody_mixed <- woody_mixed %>%
  filter(!RID %in% unique_not_in_zonal)
#Add loss values
woody_mixed$loss <- NA
#x <- 167
for (x in unique(woody_mixed$RID)) {
  row <- woody_mixed %>% filter(RID == x)
  subset_zonal_woody <- zonal_woody_treat %>% filter(RID == x)
  # Extract numerical values from the last 11 columns into a vector
  last_11_values_treat <- as.vector(unlist(subset_zonal_woody[, -c(1:2)]))
  subset_zonal_contr <- zonal_woody_contr %>% filter(RID == x)
  # Extract numerical values from the last 11 columns into a vector
  last_11_values_contr <- as.vector(unlist(subset_zonal_contr[, -c(1:2)]))
  joined <- c(last_11_values_treat, last_11_values_contr)
  woody_mixed$loss[woody_mixed$RID == x] <- joined
}
write_rds(woody_mixed, "input/joined_data_frames/woody_mixed_df.rds")


#KHAB MATCHED
khab_matched <- data.frame(
  RID = matched_mixed_RIDs$khab_matched[,1],
  CI = matched_mixed_RIDs$khab_matched[,2]
)

khab_matched <- data.frame()

RID_repeated <- rep(matched_mixed_RIDs$khab_matched[,1], each = 11)
RID_rep <- unlist(RID_repeated, use.names = FALSE)
CI_repeated <- rep(matched_mixed_RIDs$khab_matched[,2], each = 11)
CI_rep <- unlist(CI_repeated, use.names = FALSE)

# Add columns RID and CI from matched_mixed_RIDs to woody_matched
khab_matched <- data.frame(
  RID = RID_rep,
  CI = CI_rep, 
  time = NA
)

# Order the dataframe by the RID column
khab_matched <- khab_matched[order(khab_matched$RID), ]
time_repeated <- rep(c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6), length(unique(khab_matched$RID)))
khab_matched$time <- time_repeated
khab_matched$baseline <- NA

not_in_zonal <- c()
# Iterate through each row using a for loop
for (i in 1:nrow(khab_matched)) {
  row <- khab_matched[i, ]  # Access current row
  # Check if RID exists in zonal_woody_treat
  if (row$RID %in% zonal_woody_treat$RID) {
    # Perform actions if RID exist
    # Check if CI column value is 1
    if (row$CI == 1) {
      row_treat <- zonal_woody_treat %>%
        filter(RID == row$RID)
      khab_matched[i, 4] <- row_treat$baseline 
      # Add your actions here
    } else {
      row_contr <- zonal_woody_contr %>%
        filter(RID == row$RID)
      khab_matched[i, 4] <- row_contr$baseline
    }
    # Perform operations on the row
  } else {
    #print(row$RID)
    # Perform actions if RID does not exist (optional)
    not_in_zonal <- c(not_in_zonal, row$RID)
  }
}
# Remove duplicates
unique_not_in_zonal <- unique(not_in_zonal)

# unique_not_in_zonal_new
# unique_not_in_zonal_old 
# # Find common values
# common_values <- intersect(unique_not_in_zonal_new, unique_not_in_zonal_old)
# # Print common values
# print(common_values)

# Remove rows where RID matches values in unique_not_in_zonal
khab_matched <- khab_matched %>%
  filter(!RID %in% unique_not_in_zonal)

#Add loss values
khab_matched$loss <- NA

#x <- 340
for (x in unique(khab_matched$RID)) {
  row <- khab_matched %>% filter(RID == x)
  if (row$CI[1] == 1) {
    subset_zonal_woody <- zonal_woody_treat %>% filter(RID == x)
    # Extract numerical values from the last 11 columns into a vector
    last_11_values <- as.vector(unlist(subset_zonal_woody[, -c(1:2)]))
    khab_matched$loss[khab_matched$RID == x] <- last_11_values
  } else {
    subset_zonal_contr <- zonal_woody_contr %>% filter(RID == x)
    # Extract numerical values from the last 11 columns into a vector
    last_11_values <- as.vector(unlist(subset_zonal_contr[, -c(1:2)]))
    khab_matched$loss[khab_matched$RID == x] <- last_11_values
  }
}

write_rds(khab_matched, "input/joined_data_frames/khab_matched_df.rds")



#KHAB MIXED
khab_mixed <- data.frame(
  RID = matched_mixed_RIDs$khab_mixed[,1]
)



khab_mixed <- data.frame()
RID_repeated <- rep(matched_mixed_RIDs$khab_mixed[,1], each = 22)
RID_rep <- unlist(RID_repeated, use.names = FALSE)
# Create a vector with 11 ones followed by 11 zeros
vector_ones_zeros <- c(rep(1, 11), rep(0, 11))
CI_repeated <- rep(vector_ones_zeros, each = nrow(matched_mixed_RIDs$khab_mixed[,1]))
CI_rep <- unlist(CI_repeated, use.names = FALSE)

# Add columns RID and CI from matched_mixed_RIDs to woody_matched
khab_mixed <- data.frame(
  RID = RID_rep,
  CI = CI_rep, 
  time = NA
)

# Order the dataframe by the RID column
khab_mixed <- khab_mixed[order(khab_mixed$RID), ]
time_repeated <- rep(c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6), ((length(unique(khab_mixed$RID)))*2))
khab_mixed$time <- time_repeated
khab_mixed$baseline <- NA

zonal_woody_treat
zonal_woody_contr

not_in_zonal <- c()
# Iterate through each row using a for loop
for (i in 1600679:nrow(khab_mixed)) {
  row <- khab_mixed[i, ]  # Access current row
  # Check if RID exists in zonal_woody_treat
  if (row$RID %in% zonal_woody_treat$RID) {
    # Perform actions if RID exist
    # Check if CI column value is 1
    if (row$CI == 1) {
      row_treat <- zonal_woody_treat %>%
        filter(RID == row$RID)
      khab_mixed[i, 4] <- row_treat$baseline 
      # Add your actions here
    } else {
      row_contr <- zonal_woody_contr %>%
        filter(RID == row$RID)
      khab_mixed[i, 4] <- row_contr$baseline
    }
    # Perform operations on the row
  } else {
    #print(row$RID)
    # Perform actions if RID does not exist (optional)
    not_in_zonal <- c(not_in_zonal, row$RID)
  }
}
write_rds(khab_mixed, "input/joined_data_frames/khab_mixed_df_baseline_only.rds")

# Remove duplicates
unique_not_in_zonal <- unique(not_in_zonal)
write_rds(unique_not_in_zonal, "input/joined_data_frames/unique_not_in_zonal_khab_mixed.rds")

# Remove rows where RID matches values in unique_not_in_zonal
khab_mixed <- khab_mixed %>%
  filter(!RID %in% unique_not_in_zonal)

#Add loss values
khab_mixed$loss <- NA

#x <- 167
for (x in unique(khab_mixed$RID)) {
  row <- khab_mixed %>% filter(RID == x)
  subset_zonal_treat <- zonal_woody_treat %>% filter(RID == x)
  # Extract numerical values from the last 11 columns into a vector
  last_11_values_treat <- as.vector(unlist(subset_zonal_treat[, -c(1:2)]))
  subset_zonal_contr <- zonal_woody_contr %>% filter(RID == x)
  # Extract numerical values from the last 11 columns into a vector
  last_11_values_contr <- as.vector(unlist(subset_zonal_contr[, -c(1:2)]))
  joined <- c(last_11_values_treat, last_11_values_contr)
  khab_mixed$loss[khab_mixed$RID == x] <- joined
}
write_rds(khab_mixed, "input/joined_data_frames/khab_mixed_df.rds")


