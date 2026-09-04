## script to help identify study sites for wildfire work
  #future to do: get prelim burn sev from landsat? plot soil and geology maps, clean up file saving it it pulls from baselayers/writes to

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
  library(dataRetrieval)
  library(lubridate)
  library(boxrdrive)
  library(leaflet)
  library(htmltools)
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

    road_lines <- sf::st_transform(road_lines, crs(basins))
      road_lines <- sf::st_intersection(road_lines, basins)
      road_lines <- st_cast(road_lines, "MULTILINESTRING")

      #these take a while to load sometimes, so save them
      write_sf(road_lines, road_file)
    }

  #get existing WQ/flow data
  sonde_file <- file.path(workdir, paste0(prjname, "-sonde-data.gpkg"))
  if(file.exists(sonde_file)){
    sonde_sites <- read_sf(sonde_file)
  }else{
    #get sonde/continuous data
    usgs_monitor <- read_waterdata_ts_meta(hydrologic_unit_code = basins$huc8)
    sonde_sites <- usgs_monitor %>% filter(parameter_name %in% c("Discharge", "Turbidity, FNU", "fDOM, water, in situ")) %>%
      filter(end >= Sys.Date() - years(2)) %>%
      select(monitoring_location_id, parameter_name) %>% group_by(monitoring_location_id) %>%
      distinct() %>%
      summarise(count=n()) %>% mutate(site_type = ifelse(count == 1, "flow only", "full sonde"))

    #get names
    names <- read_waterdata_monitoring_location(monitoring_location_id = sonde_sites$monitoring_location_id) %>%
      select(monitoring_location_id, monitoring_location_name) %>% st_drop_geometry()

    sonde_sites <- sonde_sites %>% left_join(names)

    write_sf(sonde_sites, sonde_file)

    }

  wq_file <- file.path(workdir, paste0(prjname, "-wq-data.gpkg"))
  if(file.exists(wq_file)){
    wq_sites <- read_sf(wq_file)
  }else{
    #get wq samples
    wqbbox <- st_bbox(basins)
    wqsites <- whatWQPdata(bBox=wqbbox) %>% filter(ResolvedMonitoringLocationTypeName == "Stream")

    wqsites <- lapply(c("Inorganics, Minor, Metals", "Nutrient", "Organics, Other"), function(x){
      readWQPsummary(bBox=wqbbox, siteType = "Stream", characteristicType = x) %>%
        select(MonitoringLocationIdentifier, CharacteristicType,ResultCount,MonitoringLocationLatitude,MonitoringLocationLongitude) %>%
        distinct() %>% group_by(across(-ResultCount)) %>% summarise(nobs = sum(ResultCount), .groups = "drop_last")
    }) %>% bind_rows() %>% filter(nobs > 10)

    wq_sites <- st_as_sf(wqsites, coords=c("MonitoringLocationLongitude","MonitoringLocationLatitude"), crs="EPSG:4326")
    wq_sites <-  sf::st_intersection(wq_sites, basins %>% st_geometry())

    write_sf(wq_sites, wq_file)

  }

  #previous wildfires so we can exclude sites burned previously
    prev_fires <- read_sf(file.path(boxrdrive::box_drive(),
                            "Wildfire_Water_Security/02_Nodes/06_Shared_Data_Projects/01_Spatial_Data/mtbs-data/CONUS-boundaries/mtbs_perims_DD.shp"))
    prev_fires <- st_transform(prev_fires, crs(basins))
    prev_fires <-  sf::st_intersection(prev_fires, basins %>% st_geometry())
    prev_fires <- st_cast(prev_fires, "MULTIPOLYGON")


# step 2: down-select to potential sites -----
  #filter to streams of a reasonable size
    stream_pot <- streams %>% filter(totdasqkm >= stream_size[1] & totdasqkm <= stream_size[2])

  #classify as burned or not
    stream_pot$burned <- as.logical(rowSums(st_intersects(stream_pot, fires, sparse = FALSE)))

    #classify as previously burned (remove these sites)
    stream_pot$prev_burned <- as.logical(rowSums(st_intersects(stream_pot, prev_fires, sparse = FALSE)))
    stream_pot <- stream_pot %>% filter(!prev_burned)

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
    basemap <- leaflet() %>% addProviderTiles(providers$Esri.NatGeoWorldMap) %>%
      addPolygons(data=prev_fires, fill="red", color="darkred", group="Previous Fires") %>%
      addPolygons(data=fires, fill="darksalmon", color="red", group="Fires") %>%
      addPolygons(data=pot_basins, fill=NA, color="orange",opacity = 1, weight=1,
                  group = "Watershed Areas") %>%
    addPolylines(data=streams, color="darkblue", opacity = 1, weight=2,
                 label = ~gnis_name,
                 highlightOptions = highlightOptions(
                   weight = 6,  color = "lightblue",  bringToFront = TRUE)) %>%
      addCircleMarkers(data=sonde_sites, group="Sonde Data", radius =1.5, color = "black",
                       label = ~htmlEscape(monitoring_location_name)) %>%
      addCircleMarkers(data=wq_sites, color="#6a3673", radius =1.5,
                       group="Pre-Fire Data", opacity = 0.7) %>%
      addLayersControl(
        overlayGroups = c("Sonde Data","Pre-Fire Data", "Previous Fires",
                          "Watershed Areas"),
        options = layersControlOptions(collapsed = FALSE))

    basemap

#step 5: manually subset options, and examine characteristics -----
  sel_sites <- c("Switch Creek", "23810176","23810164", "23810162",
                      "Pan Creek", "Wolf Creek")


  #calculate euclidean distance from selected burn sites
    site_char <- stream_data %>%
      select(comid, gnis_name, lengthkm,totdasqkm,slope,maxelevsmo,minelevsmo,burned:ACC_STREAM_SLOPE)
    site_char_scaled <- site_char %>% mutate(across(c(totdasqkm, elevws:ACC_STREAM_SLOPE),~ scale(.x)[,1]))

    get_euclidean_dist <- function(site, site_char_scaled){
      burn <- site_char_scaled %>% filter(gnis_name == site)
      dist <- site_char_scaled %>% mutate(burn_dist =
                                            sqrt(rowSums(across(c(totdasqkm, elevws:ACC_STREAM_SLOPE), ~ (.x - burn[[cur_column()]])^2), na.rm = TRUE))) %>%
        pull(burn_dist)

      return(dist)
    }


  #for each burn site get distance and select five "nearest" sites
   get_ref <- function(name, n = 5){
     best_ref <- site_char
     best_ref$euclid <- get_euclidean_dist(name, site_char_scaled)
     best_ref <- best_ref %>% arrange(euclid) %>% filter(!burned & poly_IncidentName == best_ref$poly_IncidentName[best_ref$gnis_name==name])

     pot_ref <- best_ref$gnis_name[2:(1+n)]
     return(pot_ref)}

#step 6: save geopackage for each selected site ------
  #create base layers we want on every map
   base_gpk <- file.path(workdir, "base-map-layers.gpkg")
     st_write(fires, dsn = base_gpk, layer = "current_fires",append=FALSE)
     st_write(prev_fires, dsn = base_gpk, layer = "previous_fires", append=TRUE)
     st_write(streams, dsn = base_gpk, layer = "streams", append=TRUE)
     st_write(road_lines, dsn = base_gpk, layer = "roads", append=TRUE)
     st_write(sonde_sites, dsn = base_gpk, layer = "sondes", append=TRUE)
     st_write(wq_sites, dsn = base_gpk, layer = "wq_data", append=TRUE)


  save_site_info <- lapply(sel_sites, function(x){
    fire_gpk <- file.path(workdir, paste0("burn-site-", x, ".gpkg"))

    #get streams of interest
    stream_get <- c(x, get_ref(x))
    sub_sites <- site_char %>% filter(gnis_name %in% stream_get)
    sub_scaled <- sub_sites %>% mutate(across(c(totdasqkm, elevws:ACC_STREAM_SLOPE),~ scale(.x)[,1]))
    sub_sites$euclid <- get_euclidean_dist(x, sub_scaled)

    #get basins of interest
    sub_basins <- pot_basins %>% filter(comid %in% sub_sites$comid) %>% distinct()

    st_write(sub_sites, dsn = fire_gpk, layer = "substreams",append=FALSE)
    st_write(sub_basins, dsn = fire_gpk, layer = "subbasins", append=TRUE)

  })

#step 7: compile all potential references ------
  #add any sites to exclude
  exclude <- c("Kansas Creek", "Buckeye Creek", "23810300")


  pot_ref <- unique(as.vector(sapply(sel_sites, get_ref)))
  potential_sites <- setdiff(c(sel_sites, pot_ref), exclude)

  #get lat and long for these sites
  flowline <- nhdplusTools::get_nhdplus(comid = streams$comid[streams$gnis_name %in% potential_sites]) %>%
    mutate(gnis_name = ifelse(gnis_name == " ", comid, gnis_name))
  outlet_point <- nhdplusTools::get_node(flowline, position = "end")
  outlet_point$name <- flowline$gnis_name

  #save as shapefile and csv
  write_sf(outlet_point, file.path(workdir, paste0(prjname, "-proposed-sites.gpkg")))

  points_csv <- outlet_point %>%
    mutate(lon = st_coordinates(.)[, "X"],lat = st_coordinates(.)[, "Y"]) %>%
    st_drop_geometry()

  write.csv(points_csv, file.path(workdir, paste0(prjname, "-proposed-sites.csv")), row.names=FALSE)
