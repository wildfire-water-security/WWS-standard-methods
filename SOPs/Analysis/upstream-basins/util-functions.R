#' Prepare data required to submit to StreamStats
#'
#' Downloads the required datasets and snaps provided sites to the
#' StreamStats streamgrid. The provided `.zip` file can be used to
#' generate the required upstream basin areas via StreamStats batch
#' processing.
#'
#' @param site_csv File path to the `.csv` with the site information.
#' @param huc_code HUC number for the basin encompassing all your sites, this is used for grabbing stream data.
#' @param data_wd File path where you want all the spatial data saved (e.g., `data/spatial-data`).
#' @param study_code Study code associated with the sites, this is used for naming files.
#' @param crs Coordinate reference system associated with the site latitude and longitudes.
#' @param snap_dist Maximum snap distance in map units.
#'
#' @returns a `sf` object with the snapped data points. Also saves the snapped points
#' to the `data_wd` with a message providing the save location.
#' @export
#' @md
prep_streamstats <- function(site_csv, data_wd, study_code, huc_code, crs="EPSG:4326",
                             snap_dist = 1000){
  stopifnot(file.exists(site_csv), dir.exists(dirname(data_wd)), is.character(study_code),
            is.numeric(huc_code) | is.character(huc_code))

  #check directories and create ones for writing data to
    data_wd <- normalizePath(data_wd)
    dir.create(file.path(data_wd, paste0(study_code, "-sites")), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(data_wd, paste0(study_code, "-basins/check-plots")), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(data_wd, "streamstats-output"), recursive = TRUE, showWarnings = FALSE)

  #read in .csv
    data <- read.csv(site_csv)

    #check .csv conforms
    stopifnot(c("siteID","locality","latitude","longitude") %in% colnames(data))

  #turn sites into a gpkg and save
    sites <- data %>% sf::st_as_sf(coords=c("longitude","latitude"), crs=crs)
    sf::write_sf(sites, file.path(data_wd, paste0(study_code, "-sites"), paste0(study_code,"-sites.gpkg")))

    #also save one as _adj for editing
    sf::write_sf(sites, file.path(data_wd, paste0(study_code, "-sites"),
                                  paste0(study_code,"-sites_adj.gpkg")))

  #get stream data
    huc_code <- as.character(huc_code)
    type <- stringr::str_pad(nchar(huc_code), 2, side="left", pad="0")
    stopifnot(type %in% c("02", "04", "06", "08", "10", "12"))

    basin <- nhdplusTools::get_huc(id=huc_code, type=paste0("huc",type))
    nhd_streams <- nhdplusTools::get_nhdplus(AOI = basin)
    sf::write_sf(nhd_streams, file.path(data_wd, paste0(study_code, "-nhd-streams.gpkg")))

  #get the streamstats grid of streams, clip to basin before saving
    state <- unique(basin$states)
    if(length(state) > 1){stop("code currently not set up to handle more than one state, contact Katie Wampler to update")}

    #download grid and save to temp
    url <- paste0("https://streamstats.usgs.gov/streamgrids/", state, "/", state, ".zip")
    if(!file.exists(file.path(tempdir(), paste0("streamgrid-", state, "/streamgrid.tif")))){
      download.file(url = url, destfile = file.path(tempdir(), paste0("streamgrid-", state, ".zip")),
                                                    mode = "wb") # 'wb' for binary
      unzip(file.path(tempdir(), paste0("streamgrid-", state, ".zip")),
            exdir = file.path(tempdir(), paste0("streamgrid-", state)))
    }

    stream_grid <- terra::rast(file.path(tempdir(), paste0("streamgrid-", state, "/streamgrid.tif")))
    basin <- terra::vect(basin)
    basin_pj <- terra::project(basin, terra::crs(stream_grid))
    stream_grid <- terra::crop(stream_grid, basin_pj)
    terra::coltab(stream_grid) <- NULL #remove color table so we don't get warnings
    terra::writeRaster(stream_grid, file.path(data_wd, paste0(study_code, "-streamstats-grid.tif")), overwrite=TRUE)

  #snap points to the grid
    if(file.exists(file.path(data_wd, paste0(study_code, "-sites"),
                             paste0(study_code,"-sites_adj.gpkg")))){
      sites <- sf::read_sf(file.path(data_wd, paste0(study_code, "-sites"),
                                     paste0(study_code,"-sites_adj.gpkg")))
    }

    #ensure same projections
    sites_pj <- sf::st_transform(sites, terra::crs(stream_grid)) %>% dplyr::select(siteID)

    #save layers to use with whitebox -> temp directory
    sf::write_sf(sites_pj, file.path(tempdir(), "sites.shp"))
    stream_grid <- terra::rast(file.path(data_wd, paste0(study_code, "-streamstats-grid.tif")))
    terra::writeRaster(stream_grid, file.path(tempdir(), "streams_crop.tif"), overwrite=TRUE)

    # Extract the streams and snap the pour points to the streams
    whitebox::wbt_init() # initiate whitebox tools
    whitebox::wbt_jenson_snap_pour_points(pour_pts = "sites.shp",
                                streams = "streams_crop.tif",
                                output = paste0("snapped-", study_code,"-sites.shp"),
                                snap_dist = snap_dist,
                                wd = tempdir())

    #load points
    snapped_points <- sf::read_sf(file.path(tempdir(), paste0("snapped-", study_code,"-sites.shp")))

    #plot to visualize, might not look like a lot
    ggplot2::ggplot() + ggplot2::geom_sf(data=nhd_streams, color="darkblue") +
      ggplot2::geom_sf(data=sites, color="black", alpha=0.5) +
      ggplot2::geom_sf(data = snapped_points, color="red")

    #save snapped points to main folder to save for streamstats
    sf::write_sf(snapped_points, file.path(data_wd, paste0(study_code, "-sites"),
                                       paste0("snapped-", study_code,"-sites.gpkg")))

    message("site points have been prepared for StreamStats batch processing: \n",
            normalizePath(file.path(data_wd, paste0(study_code, "-sites"))),
            "\n\nPlease manually check *-sites_adj.gpkg and adjust points to ensure sites are on the correct flowline, then rerun.")

    return(snapped_points)

}


#' Query StreamStats API for upstream basins
#'
#' Uses the coordinates of the snapped points to query the USGS StreamStats API to
#' get the upstream basin for each point
#'
#' @param snapped_points An `sf` object with snapped sites.
#' @param region Two digit state or region code associated with the sites.
#'
#' @returns A `sf` object with a `multipolygon` of the upstream area for each row in `snapped_points`.
#' @export
#'
get_streamstats <- function(snapped_points, region){
  stopifnot(inherits(snapped_points, "sf"), nchar(region) == 2)

  #ensure correct crs
  sites <- sf::st_transform(snapped_points, "EPSG:4326")

  #extract lat and longs
  coords <- sf::st_coordinates(sites)
  sites$lat <- coords[,2]
  sites$long <- coords[,1]
  sites$api_basin <- paste0("https://streamstats.usgs.gov/ss-delineate/v1/delineate/sshydro/", region,
                     "?lat=", sites$lat, "&lon=", sites$long)
  sites$api_char <-  paste0("https://streamstats.usgs.gov/ss-hydro/v1/basin-characteristics/calculate-using-ssdelineate/?region=", region,
                            "&lat=", sites$lat, "&lon=", sites$long)

  #get basins
  download_basin <- function(site, region, lat, long){
    base_url <- "https://streamstats.usgs.gov/"

    #get basin
    resp_del <- httr::GET(paste0(base_url, "/ss-delineate/v1/delineate/sshydro/", region,"?lat=", lat, "&lon=", long))
    httr::stop_for_status(resp_del)
    txt_del <- httr::content(resp_del, as = "text", encoding = "UTF-8")
    del <- jsonlite::fromJSON(txt_del, simplifyVector = FALSE)

    fc_list <- del$bcrequest$wsresp$featurecollection[[1]][[2]]
    basin <- sf::st_read(jsonlite::toJSON(fc_list$feature, auto_unbox = TRUE),quiet = TRUE)

  #get characteristics currently can't get API to work
    # resp_bc <- httr::POST(paste0(base_url, "ss-hydro/v1/basin-characteristics/calculate"),
    #                       config=list(httr::add_headers(accept = "application/json",`Content-Type` = "application/json")),
    #   body = del,
    #   encode = "json")
    #
    # httr::stop_for_status(resp_bc)
    #
    # txt_bc <- httr::content(resp_bc, as = "text", encoding = "UTF-8")
    # bc <- fromJSON(txt_bc, simplifyVector = TRUE)


    return(basin)

  }
  message("downloading basins from streamstats...")
  basins <- pbapply::pblapply(1:nrow(sites),
                              function(x){download_basin(sites$siteID[x], region, sites$lat[x], sites$long[x])})

  #save basins
  basins_df <- basins %>% dplyr::bind_rows() %>% dplyr::mutate(siteID = sites$siteID, .before="GlobalWshd")

  #return to user
  return(basins_df)

}

validate_streamstats <- function(basins, data_wd, study_code){
  stopifnot(inherits(basins, "sf"), dir.exists(data_wd))

  #load data we need for plots
    sites <- sf::read_sf(file.path(data_wd, paste0(study_code, "-sites"),
                               paste0(study_code,"-sites_adj.gpkg")))
    sites <- sf::st_transform(sites, terra::crs(basins))
    nhd_streams <- sf::read_sf(file.path(data_wd, paste0(study_code, "-nhd-streams.gpkg")))
    nhd_streams <- sf::st_transform(nhd_streams, terra::crs(basins))

  #create folder for check plots
    dir.create(file.path(data_wd, paste0(study_code, "-basins/check-plots")), recursive = TRUE, showWarnings = FALSE)

  #loop through, saving each shapefile individually and make a plot to check it's location is correct
    message("saving basin files and creating validation plots...")
    plot_basin <- function(x, sites, data_wd, study_code){
      site <- sites[sites$siteID == x$siteID,]

      #get a buffer around site for the zoom
      site_bf <- sf::st_buffer(site, dist=500)

      #get intersecting stream to site
      sf::st_agr(site_bf) = "constant"       #to prevent warnings
      sf::st_agr(nhd_streams) = "constant"
      streams_named <- unique(sf::st_intersection(site_bf, nhd_streams))
      basin_ext <- sf::st_bbox(site_bf)

      #save basin
      sf::write_sf(x, file.path(data_wd, paste0(study_code, "-basins"), paste0(site$siteID, "-basin.gpkg")))

      #plot zoomed in with stream name that should match
      p1 <- ggplot2::ggplot() + ggplot2::geom_sf(data = x, fill="gray") + ggplot2::geom_sf(data = streams_named) +
        ggplot2::geom_sf(data=site) +  ggplot2::coord_sf(xlim = c(basin_ext[1], basin_ext[3]), ylim = c(basin_ext[2], basin_ext[4])) +
        ggrepel::geom_text_repel(data = streams_named, ggplot2::aes(label = gnis_name, geometry = geom),
                        stat = "sf_coordinates", min.segment.length = 0, alpha=0.5,
                        max.overlaps=20) +
        ggrepel::geom_text_repel(data = site, ggplot2::aes(label = locality, geometry = geom),
                        stat = "sf_coordinates", min.segment.length = 0, color="darkred",
                        max.overlaps=20)

        suppressWarnings(ggplot2::ggsave(file.path(data_wd, paste0(study_code, "-basins/check-plots"), paste0(site$siteID, "_zoom.png")),
             plot=p1, height = 6.13, width = 5.57, units = "in"))

      #plot zoomed out to see full extent
        basin_ext <- sf::st_bbox(x)
        p2 <- ggplot2::ggplot() + ggplot2::geom_sf(data = x, fill="gray") + ggplot2::geom_sf(data = nhd_streams) +
              ggplot2::geom_sf(data=site) +  ggplot2::coord_sf(xlim = c(basin_ext[1], basin_ext[3]), ylim = c(basin_ext[2], basin_ext[4]))


      suppressWarnings(ggplot2::ggsave(file.path(data_wd, paste0(study_code, "-basins/check-plots"), paste0(site$siteID, "_full.png")),
                                       plot=p2, height = 6.13, width = 5.57, units = "in"))


    }

    save <- sapply(1:nrow(basins), function(x){plot_basin(basins[x,], sites, data_wd, study_code)})

    message(paste0("Basins saved here: ", normalizePath(file.path(data_wd, paste0(study_code, "-sites"))), "\nCheck 'check-plots' for accuracy."))
}
