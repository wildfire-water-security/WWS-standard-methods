## script to help identify study sites for wildfire work
  #first choice: peat and switch creek

  #load libraries
  library(arcgislayers)
  library(sf)
  library(terra)
  library(dplyr)
  library(hydrogeofetch)
  library(ggplot2)
  library(plotly)
  library(pbapply)
  library(StreamCatTools)
  library(fs)
  library(tidyr)
  library(purrr)
  library(tigris)
  sf_use_s2(FALSE)

#update code here
  stream_size <- c(0.5,3) #km2 size range to aim for, will only look at streams of this size
  ref_dist <- 10 #how many km away would you consider a reference?
  fire_names <- c("Austin", "Grasshopper") #the fire(s) you want to target
  prjname <- "example"
  workdir <- file.path(fs::path_home(), "Documents/Projects/WWS-standard-methods/SOPs/Data Analysis/selecting-fire-sites")

# step 1: download data ------
  #load current fires
    fire_url <- "https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/WFIGS_Interagency_Perimeters_Current/FeatureServer/0"
    curr_fires <- arc_read(fire_url, sf_type = "esriGeometryPolygon")

  #filter to fires of interest
    fires <- curr_fires %>% filter(poly_IncidentName %in% fire_names)
    sf_use_s2(TRUE)
    fires <- st_make_valid(fires)
    sf_use_s2(FALSE)

    fire_file <- file.path(workdir, paste0(prjname, "-fires.gpkg"))
    write_sf(fires, fire_file)


  #get hucs for those burns
    basin_file <- file.path(workdir, paste0(prjname, "-basins.gpkg"))
    if(file.exists(basin_file)){
      basins <- read_sf(basin_file)
    }else{
      basins <- pblapply(1:nrow(fires),function(x){
        nhdplusTools::get_huc(fires[x,], type='huc08')
      }) %>% bind_rows()

      #these take a while to load sometimes, so save them
      write_sf(basins, basin_file)
    }

  #get streams (can take a minute)
    stream_file <- file.path(workdir, paste0(prjname, "-streams.gpkg"))
    if(file.exists(stream_file)){
      streams <- read_sf(stream_file)
    }else{
      streams <- pblapply(1:nrow(basins), function(x){
        nhdplusTools::get_nhdplus(AOI = basins[x,])
      }) %>% bind_rows() %>% mutate(gnis_name = ifelse(gnis_name == " ", comid, gnis_name))

      #these take a while to load sometimes, so save them
      write_sf(streams, stream_file)
    }

  # #get fs roads
  #   road_file <- file.path(workdir, paste0(prjname, "-fsroads.gpkg"))
  #   if(file.exists(road_file)){
  #     fsroads <- read_sf(road_file)
  #   }else{
  #     road_layer <- arc_open("https://apps.fs.usda.gov/arcx/rest/services/EDW/EDW_RoadBasic_01/MapServer/0")
  #     #road_layer <- arc_open("https://gis.odot.state.or.us/arcgis1006/rest/services/transgis/catalog/MapServer/164")
  #     my_bbox <- st_bbox(basins)
  #     roads <- arc_select(road_layer, filter_geom = my_bbox)
  #
  #     roads <- sf::st_transform(roads, crs(basins))
  #     roads <- sf::st_intersection(roads, basins)
  #     fsroads <- st_cast(roads_clip, "MULTILINESTRING")
  #
  #     #these take a while to load sometimes, so save them
  #     write_sf(fsroads, road_file)
  #   }

  #get other roads
    road_file <- file.path(workdir, paste0(prjname, "-roads.gpkg"))
    if(file.exists(road_file)){
      road_lines <- read_sf(road_file)
    }else{
      #get counties we want road data for
        count <- counties(state = "OR")
        count <- sf::st_transform(count, crs(basins))
        count <- sf::st_intersection(count, basins)

      #get roads for those counties
      road_lines <- roads("41", count$COUNTYFP)


      # road_layer <- arc_open("https://gis.odot.state.or.us/arcgis1006/rest/services/transgis/catalog/MapServer/348")
      # my_bbox <- st_bbox(basins)
      # roads <- arc_select(road_layer, filter_geom = my_bbox)

      road_lines <- sf::st_transform(road_lines, crs(basins))
      road_lines <- sf::st_intersection(road_lines, basins)
      road_lines <- st_cast(road_lines, "MULTILINESTRING")

      #these take a while to load sometimes, so save them
      write_sf(road_lines, road_file)
    }


# step 2: down-select to potential sites -----
  #filter to streams of a reasonable size
    stream_pot <- streams %>% filter(totdasqkm >= stream_size[1] & totdasqkm <= stream_size[2])

  #classify as burned or not
    stream_pot$burned <- as.logical(rowSums(st_intersects(stream_pot, fires, sparse = FALSE)))

  #only include unburned sites within 10 km2 of burn
    fires_flat <- st_transform(fires, crs="EPSG:32610")
    burn_buff <- st_buffer(fires_flat, dist = ref_dist*1000)
    burn_buff <- st_transform(burn_buff, crs(basins)) %>% select(poly_IncidentName)
    stream_pot <- sf::st_intersection(stream_pot, burn_buff)
    stream_pot <- st_cast(stream_pot, "MULTILINESTRING")


# step 3: get site landscape information -----
  comid <- stream_pot$comid

  #usgs char
    usgs_char_opts <- nhdplusTools::get_characteristics_metadata() %>% filter(grepl("ACC_", ID))

    usgs_chars <- c("ACC_BFI", "ACC_CONTACT", "ACC_EWT", "ACC_TWI", "ACC_WB5100_ANN",
                    "ACC_TAV7100_ANN", "ACC_CNPY11_BUFF100", "ACC_OM", "ACC_CLAYAVE",
                    "ACC_SANDAVE", "ACC_SILTAVE", "ACC_ROCKDEP", "ACC_WTDEP","ACC_BASIN_SLOPE",
                    "ACC_STREAM_SLOPE")

  suppressWarnings(usgs_landscape <- nhdplusTools::get_catchment_characteristics(varname = usgs_chars, ids=comid) %>%
    select(-"percent_nodata") %>% pivot_wider(names_from = "characteristic_id",
                                              values_from = "characteristic_value"))

  #epa char
    epa_char_opt <- sc_get_params(param='variable_info')

    epa_chars <- c("pctbl2019", "pctconif2019", "pctdecid2019","pctmxfst2019", "pctshrb2019", "bankfulldepth", "bankfullwidth",
                   "elev", "hydrlcond", "pcteolcrs", "pcteolfine", "pctextruvol", "pctcolluvsed")

    epa_landscape <- sc_get_data(comid=comid, metric=epa_chars, aoi=c('ws', 'other'))

  #combine and save
    stream_data <- stream_pot %>% left_join(epa_landscape, by="comid") %>% left_join(usgs_landscape, by="comid")

    write_sf(stream_data, file.path(workdir, paste0(prjname, "-site-opts.gpkg")))

  #determine percentage burned
    #get basins for burned to visualize how much of watershed is burned
    pot_basins <- pblapply(comid, function(x){
      nldi_feature <- list(featureSource = "comid", featureID = x)
      basin <- get_nldi_basin(nldi_feature = nldi_feature)
    }) %>% bind_rows()
    pot_basins <- st_cast(pot_basins, "MULTIPOLYGON")
    pot_basins$comid <- comid
    write_sf(pot_basins, file.path(workdir, paste0(prjname, "-site-opts-basins.gpkg")))



#step 4: view plot -----
  #set up tool tips
  suppressWarnings(p <- ggplot() + geom_sf(data = basins, fill=NA, color="black")  +
                                   geom_sf(data=fires, fill="red", color="darkred", alpha=0.5) +
                                   geom_sf(data=streams, color="darkblue") +
                                   geom_sf(data=stream_pot, aes(color = burned, text=gnis_name)) +
                     scale_color_manual(values = c("blue", "green")) +
                                   geom_sf(data= road_lines, color="gray10", linewidth=0.2) +
                     geom_sf(data=pot_basins, fill=NA, color="orange"))

  ggplotly(p, tooltip = "text")

#step 5: manually subset options, and examine characteristics -----
  sel_sites_burn <- c("Switch Creek", "23810176","23810164", "23810162",
                      "Pan Creek", "Wolf Creek")
  sel_sites_unburn <- c("Butte Creek", "Paste Creek","23810246",
                        "Peat Creek", "23810274", "23810598",
                        "23810650", "Sand Creek", "Pink Creek",
                        "Jack Davis Creek", "23810268", "Bonner Creek",
                        "Wall Creek", "Gingham Creek")

  site_char <- stream_data %>% filter(gnis_name %in% c(sel_sites_burn, sel_sites_unburn)) %>%
    select(comid, gnis_name, lengthkm,totdasqkm,slope,maxelevsmo,minelevsmo,burned:ACC_STREAM_SLOPE)

  write_sf(site_char, file.path(workdir, paste0(prjname, "-selsite-opts.gpkg")))
