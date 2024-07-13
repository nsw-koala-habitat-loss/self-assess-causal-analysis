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
#zonal_woody_treat <- read_rds("input/analysis_data/woody_treat_zonal.rds")
#zonal_woody_contr  <- read_rds("input/analysis_data/woody_contr_zonal.rds")
#zonal_koala_treat <- read_rds("input/analysis_data/koala_treat_zonal.rds")
#zonal_koala_contr <- read_rds("input/analysis_data/koala_contr_zonal.rds")
#matched_mixed_RIDs <- read_rds("input/analysis_data/matched_mixed_RIDs.rds")
#matched_mixed_RIDs_old <- read_rds("input/analysis_data/matched_mixed_RIDs_old.rds")

# BROOKE TO WRITE CODE HERE TO GENERATE THE DATA FOR INPUT INTO THE STATISTICAL MODELS
