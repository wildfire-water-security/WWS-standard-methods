## code to get the average dNBR or classified burn severity for a burn stream
  #by Katie Wampler

#step 0: load functions, libraries -------
  library(sf)
  library(terra)
  library(dplyr)
  library(fs)

#load util functions
  source("SOPs/Data Analysis/upstream-basins/util-functions.R")
  source("SOPs/Data Analysis/upstream-basins/burn-severity-utils.R")

#step 1: load data -------
  #load fire data
    #recommended to download from MTBS and use `clean_mtbs` to tidy the files
    #clean_mtbs(file.path(path_home(), "Downloads"), "SOPs/Data Analysis/upstream-basins/example-data/")
    dnbr <- rast(file.path("SOPs/Data Analysis/upstream-basins/example-data/holiday_farm_dnbr.tif"))
    bound <- vect(file.path("SOPs/Data Analysis/upstream-basins/example-data/holiday_farm_burn_bndy.shp"))

  #load basins
   #to get basins from sites, see `upstream-basins.Rmd` otherwise load basins
   basin_files <- list.files("SOPs/Data Analysis/upstream-basins/example-data/example-basins", pattern="gpkg", full.names = TRUE)
   basins <- lapply(basin_files, read_sf) %>% bind_rows()

#step 3: prepare fire data -----
  #only consider dNBR in burned area
   dnbr <- crop(dnbr, bound, mask=TRUE)

  #fill unburned areas with a 0
   ext <- st_sf(st_as_sfc(st_bbox(basins)))
   basin_pj <- vect(st_transform(ext, crs=crs(dnbr)))
   dnbr_expand <- terra::extend(dnbr, basin_pj, fill=0)
   dnbr_fill <- app(dnbr_expand, function(x) ifelse(is.na(x), 0, x))

#step 4: average by basin ------
   dnbr_calc <- calc_basin_metrics(dnbr_fill, "dNBR", basins, calc="mean")

#step 5: optional, classify by burn severity ------
  #uses the dNBR classifications used for the dNBR6 map for the fire
   dnbr_calc <- dnbr_calc %>%
     mutate(burn_sev = case_when(
         dNBR <= 40  ~ "Unburned",
         dNBR <= 320 ~ "Low",
         dNBR <= 660 ~ "Moderate",
         TRUE        ~ "High"))
