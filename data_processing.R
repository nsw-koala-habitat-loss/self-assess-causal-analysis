# get required libraries
library(foreign)
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(exactextractr)

# read functions
source("functions.R")

# do some data preprocessing

#Set working directory
# Set the working directory
setwd("R:/nsw_habitat_loss/for_brooke")

# read in woody vegetation extent data for 2011
woody <- rast(paste(getwd(), "/input/woody_cover/woody_nsw_2011.tif", sep = "")) %>% round()
woody_mask <- woody %>% classify(cbind(c(0, 1), c(NA, 1)))

# read in nvr data, reclassify for treatment (category 1 and category 2 regulated) and control (category 2 - vulnerable and sensitive) then aggregate, reproject, and snap to woody cover layer
nvr_incl <- rast(paste(getwd(), "/input/nvr_mapping/naluma_nsw_2017_abel0_c20221212_u9.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(NA, 1), c(1, NA))) %>% project(woody)
# get treatment layer (1 = treatment, 0 = control)
nvr_treat <- rast(paste(getwd(), "/input/nvr_mapping/naluma_nsw_2017_abkl0_c20221212_u9.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(NA, 3, 4, 5, 6), c(1, 0, 0, NA, 0))) %>% project(woody)
# get control layer (1 = control, 0 = treatment)
nvr_contr <- rast(paste(getwd(), "/input/nvr_mapping/naluma_nsw_2017_abkl0_c20221212_u9.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(NA, 3, 4, 5, 6), c(0, 1, 1, NA, 1))) %>% project(woody)

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

# read in properties
props <- vect("input/properties/Property_EPSG4283.gdb", layer = "property") %>% project("EPSG:3308") %>% as_sf()

# calculate whether cells have at least one raster cell that is subject to the LLS Act
props_zinclude <- include %>% exact_extract(props, "max") %>% as_tibble() %>% rename(include = value)

# bind zonal values to properties
props <- cbind(props, props_zinclude)

# remove properties that are not inside the LLS region
# only keep the RID field
props_incl <- props %>% filter(!is.na(include)) %>% select(RID)

# save properties to shapefile
st_write(props_incl, "input/analysis_data/props_incl.shp", append = FALSE)

# read in koala habitat data and reclassify so 1 = koala habitat, 0 = non-habitat then aggregate, reproject, and snap to woody cover layer
# then aggregate, reproject, and snap to woody cover layer
khab <- rast(paste(getwd(), "/input/koala_habitat/KoalaHabitatSuitabilityModelClasses.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(1, 2, 3, 4, 5, 6), c(0, 0, 0, 1, 1, 1))) %>% project(woody)
gc()

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

# get the matched and mixed data samples

# current samples
# note here we switch the CI coding as it is the wrong way around - NEED TO CHECK
woody_matched <- as_tibble(read.dbf("input/matched_mixed_properties/match_BACI.dbf")) %>% group_by(RID) %>% summarise(CI = mean(CI)) %>% mutate(CI = ifelse(CI == 0, 1, 0), CI = as.integer(CI)) 
# note here we switch the CI coding as it is the wrong way around - NEED TO CHECK
khab_matched <- as_tibble(read.dbf("input/matched_mixed_properties/khab_match_BACI.dbf")) %>% group_by(RID) %>% summarise(CI = mean(CI)) %>% mutate(CI = ifelse(CI == 0, 1, 0), CI = as.integer(CI)) 
woody_mixed <- unique(as_tibble(read.dbf("input/matched_mixed_properties/mixed_BACI.dbf")) %>% select(RID)) %>% as_tibble() %>% mutate(RID = as.integer(RID)) 
khab_mixed <- unique(as_tibble(read.dbf("input/matched_mixed_properties/khab_mix_BACI.dbf")) %>% select(RID)) %>% as_tibble() %>% mutate(RID = as.integer(RID)) 

# old samples
woody_matched_old <- as_tibble(read.dbf("input/matched_mixed_properties/match_BACI_old.dbf")) %>% group_by(RID) %>% summarise(CI = mean(CI)) %>% mutate(CI = as.integer(CI)) 
khab_matched_old <- as_tibble(read.dbf("input/matched_mixed_properties/khab_match_BACI_old.dbf")) %>% group_by(RID) %>% summarise(CI = mean(CI)) %>% mutate(CI = as.integer(CI)) 
woody_mixed_old <- unique(as_tibble(read.dbf("input/matched_mixed_properties/mixed_BACI_old.dbf")) %>% select(RID)) %>% as_tibble() %>% mutate(RID = as.integer(RID)) 
khab_mixed_old <- unique(as_tibble(read.dbf("input/matched_mixed_properties/khab_mix_BACI_old.dbf")) %>% select(RID)) %>% as_tibble() %>% mutate(RID = as.integer(RID)) 

# compile into a list and save
matched_mixed_RIDs <- list(woody_matched = woody_matched, khab_matched = khab_matched, woody_mixed = woody_mixed, khab_mixed = khab_mixed)
matched_mixed_RIDs_old <- list(woody_matched = woody_matched_old, khab_matched = khab_matched_old, woody_mixed = woody_mixed_old, khab_mixed = khab_mixed_old)
write_rds(matched_mixed_RIDs, "input/analysis_data/matched_mixed_RIDs.rds")
write_rds(matched_mixed_RIDs_old, "input/analysis_data/matched_mixed_RIDs_old.rds")

# load processed data back in if needed
zonal_woody_treat <- read_rds("input/analysis_data/woody_treat_zonal.rds")
zonal_woody_contr  <- read_rds("input/analysis_data/woody_contr_zonal.rds")
zonal_koala_treat <- read_rds("input/analysis_data/koala_treat_zonal.rds")
zonal_koala_contr <- read_rds("input/analysis_data/koala_contr_zonal.rds")
matched_mixed_RIDs <- read_rds("input/analysis_data/matched_mixed_RIDs.rds")
matched_mixed_RIDs_old <- read_rds("input/analysis_data/matched_mixed_RIDs_old.rds")

woody_matched_old <- read.dbf("input/matched_mixed_properties/match_BACI_old.dbf")
khab_matched_old <- read.dbf("input/matched_mixed_properties/khab_match_BACI_old.dbf")
woody_mixed_old <- read.dbf("input/matched_mixed_properties/mixed_BACI_old.dbf")
khab_mixed_old <- read.dbf("input/matched_mixed_properties/khab_mix_BACI_old.dbf")


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
# Transforming any value greater than 0 to 1 in the 'baseline' column of zonal_woody_treat
zonal_woody_treat$test <- ifelse(zonal_woody_treat$baseline > 0, 1, zonal_woody_treat$baseline)
zonal_woody_contr$test <- ifelse(zonal_woody_contr$baseline > 0, 1, zonal_woody_contr$baseline)
rows <- zonal_woody_treat$test+zonal_woody_contr$test
zonal_woody_treat_subset <- zonal_woody_treat[which(rows == 2),]
zonal_woody_contr_subset <- zonal_woody_contr[which(rows == 2),]

# Remove the last column
zonal_woody_treat_subset <- zonal_woody_treat_subset[, -ncol(zonal_woody_treat_subset)]
zonal_woody_contr_subset <- zonal_woody_contr_subset[, -ncol(zonal_woody_contr_subset)]

zonal_woody_treat <- zonal_woody_treat_subset
zonal_woody_contr <- zonal_woody_contr_subset

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
for (i in 338124:nrow(woody_mixed)) {
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
# Transforming any value greater than 0 to 1 in the 'baseline' column of zonal_woody_treat
zonal_koala_treat$test <- ifelse(zonal_koala_treat$baseline > 0, 1, zonal_koala_treat$baseline)
zonal_koala_contr$test <- ifelse(zonal_koala_contr$baseline > 0, 1, zonal_koala_contr$baseline)
rows <- zonal_koala_treat$test+zonal_koala_contr$test
zonal_koala_treat_subset <- zonal_koala_treat[which(rows == 2),]
zonal_koala_contr_subset <- zonal_koala_contr[which(rows == 2),]
# Remove the last column
zonal_koala_treat_subset <- zonal_koala_treat_subset[, -ncol(zonal_koala_treat_subset)]
zonal_koala_contr_subset <- zonal_koala_contr_subset[, -ncol(zonal_koala_contr_subset)]

zonal_koala_treat <- zonal_koala_treat_subset
zonal_koala_contr <- zonal_koala_contr_subset

#Build dataframe
khab_mixed <- data.frame()
RID_repeated <- rep(zonal_koala_treat$RID, each = 22)
RID_rep <- unlist(RID_repeated, use.names = FALSE)
# Create a vector with 11 ones followed by 11 zeros
vector_ones_zeros <- c(rep(1, 11), rep(0, 11))
CI_repeated <- rep(vector_ones_zeros, length(zonal_koala_treat$RID))
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

#Adding baseline values
not_in_zonal <- c()
# Iterate through each row using a for loop
for (i in 1:nrow(khab_mixed)) {
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

#Adding loss values
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


