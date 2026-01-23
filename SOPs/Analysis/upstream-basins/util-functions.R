#' Prepare data required to submit to StreamStats
#'
#' Downloads the required datasets and snaps provided sites to the
#' StreamStats streamgrid. The provided `.zip` file can be used to
#' generate the required upstream basin areas via StreamStats batch
#' processing.
#'
#' @param site_csv File path to the `.csv` with the site information or a data.frame.
#' @param huc_code HUC number for the basin encompassing all your sites, this is used for grabbing stream data.
#' @param data_wd File path where you want all the spatial data saved (e.g., `data/spatial-data`).
#' @param study_code Study code associated with the sites, this is used for naming files.
#' @param crs Coordinate reference system associated with the site latitude and longitudes.
#' @param snap_dist Maximum snap distance in map units.
#' @param rewrite Logical, if `TRUE` the `study_code-sites_adj.gpkg` layer will be overwritten.
#'
#' @returns a `sf` object with the snapped data points. Also saves the snapped points
#' to the `data_wd` with a message providing the save location.
#' @export
#' @md
prep_streamstats <- function(site_csv, data_wd, study_code, huc_code, crs="EPSG:4326",
                             snap_dist = 1000, rewrite = FALSE){
  stopifnot(is.data.frame(site_csv) || file.exists(site_csv), dir.exists(dirname(data_wd)), is.character(study_code),
            is.numeric(huc_code) | is.character(huc_code))

  #check directories and create ones for writing data to
    data_wd <- normalizePath(data_wd)
    dir.create(file.path(data_wd, paste0(study_code, "-sites")), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(data_wd, paste0(study_code, "-basins/check-plots")), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(data_wd, "streamstats-output"), recursive = TRUE, showWarnings = FALSE)

  #read in .csv
    if(!is.data.frame(site_csv)){
      data <- read.csv(site_csv)
    }

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
    sf::write_sf(basin, file.path(data_wd, paste0(huc_code, "-basin.gpkg")))

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

    #save snapped points to main folder to save for streamstats
    sf::write_sf(snapped_points, file.path(data_wd, paste0(study_code, "-sites"),
                                       paste0("snapped-", study_code,"-sites.gpkg")))

    message("site points have been prepared for StreamStats batch processing: \n",
            normalizePath(file.path(data_wd, paste0(study_code, "-sites"))),
            "\n\nPlease manually check *-sites_adj.gpkg and adjust points to ensure sites are on the correct flowline, then rerun.")

    return(snapped_points)

}

#' Snap manually adjusted points to the streamgrid
#'
#' Once points have been manually adjusted, use this function to resnap the adjusted points to
#' the streamstats grid.
#'
#' @param data_wd File path where you want all the spatial data saved (e.g., `data/spatial-data`).
#' @param study_code Study code associated with the sites, this is used for naming files.
#' @param snap_dist Maximum snap distance in map units.
#'
#' @returns a `sf` object with the snapped data points. Also saves the snapped points
#' to the `data_wd` with a message providing the save location.
#' @export
#' @md
snap_streamstats <- function(data_wd, study_code, snap_dist = 1000){
  stopifnot(dir.exists(dirname(data_wd)), is.character(study_code),
            is.numeric(snap_dist))


  stream_grid <- terra::rast(file.path(data_wd, paste0(study_code, "-streamstats-grid.tif")))

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

  #save snapped points to main folder to save for streamstats
  sf::write_sf(snapped_points, file.path(data_wd, paste0(study_code, "-sites"),
                                         paste0("snapped-", study_code,"-sites.gpkg")))

  message("site points have been prepared for StreamStats batch processing.")

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
#' @md
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

#' Perform validation checks on the StreamStats basin outputs
#'
#' Saves the basin files to the `data_wd` and generates and saves plots for the
#' user to check to ensure the basins were delineated correctly.
#'
#' @details
#' Two plot types will be saved to a `*-basins/check-plots` folder:
#' - **siteID_full.png:** shows the site, the NHD plus streams, and the basin outline. Check this plot to ensure the basin
#' isn't tiny which suggests the correct stream was not grabbed.
#' - **siteID_zoom.png:** shows the site with the expected `locality` name along with any named
#' NHD plus streams within a buffer around the site. Check this plot to ensure the site locality name matches the name of the
#' stream the site is located on.
#'
#'
#' @param basins An object of class `sf` containing the streamstats basins, likely an output from `get_streamstats`.
#' @param data_wd File path where you want all the spatial data saved (e.g., `data/spatial-data`).
#' @param study_code Study code associated with the sites, this is used for naming files.
#'
#' @returns
#' Three files are save to `data_wd/side_code-basins` for each site:
#' 1. **site_code-basins/siteID-basin.gkg:** the basin polygon
#' 2. **siteID_full.png:** A plot showing the entire delineated basin, see details for more info.
#' 3. **siteID_zoom.png:** A plot showing the streams around the site, see datails for more info.
#'
#' @export
#' @md
#'
#' @examples
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
    plot_basin <- function(x, sites, data_wd, study_code, nhd_streams){
      site <- sites[sites$siteID == x$siteID,]

      site_pj <- sf::st_transform(site, "EPSG:3857")

      #get a buffer around site for the zoom
      site_bf <- sf::st_buffer(site_pj, dist=500)
      site_bf <- sf::st_transform(site_bf, terra::crs(site))

      #get intersecting stream to site
      sf::st_agr(site_bf) = "constant"       #to prevent warnings
      sf::st_agr(nhd_streams) = "constant"
      suppressMessages(streams_named <- unique(sf::st_intersection(site_bf, nhd_streams)))
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

    save <- pbapply::pbsapply(1:nrow(basins), function(x){plot_basin(basins[x,], sites, data_wd, study_code, nhd_streams)})

    message(paste0("Basins saved here: ", normalizePath(file.path(data_wd, paste0(study_code, "-sites"))), "\nCheck 'check-plots' for accuracy."))
}

#' Calculate basin summarized values
#'
#' Takes a data layer and using the provided basin boundaries calculates a metric based on the `calc` value for
#' each basin.
#'
#' @param data Either a file path to a raster (`.tif`, `.tiff`, or `.aig`) or a `spatRaster` object.
#' @param data_name A character description of the data layer, used as the column name for the output `data.frame`.
#' @param basins An object of class `sf` containing the streamstats basins, likely an output from `get_streamstats` or the file path to
#' a directory with the basin files saved as `shp` or `gpkg` files.
#' @param calc Function used to summarise the `data` layer, options include: mean, min, max, sum, isNA, notNA, percent.
#'
#' @returns A `data.frame` with rows equal to the number of basins in `basins` and two columns: the siteID and the `data_name` with the
#' calculated value.
#'
#' @export
#' @md
#'
#' @details
#' If `calc` is `percent` the data layer must be a layer comprised of only 0 and 1, the function will sum up the raster and compare to the
#' size of the overall basin layer to calculate the percent. Percent is reported as a decimal.
#'
#'
#' @examples
calc_basin_metrics <- function(data, data_name, basins, calc){
  stopifnot(inherits(data, "SpatRaster") || file.exists(data),
            inherits(basins, "sf") || dir.exists(basins), calc %in% c("mean", "min", "max", "sum", "isNA", "notNA", "percent"),
            is.character(data_name))

  #get data if only path is provided
  if(!inherits(data, "SpatRaster") & !inherits(data, "sf") & tools::file_ext(data) %in% c("tif", "tiff", "AIG")){data <- terra::rast(data)}

  if(!inherits(basins, "sf")){
    files <- list.files(basins, pattern ="[.]shp|[.]gpkg", full.names = TRUE)
    basins <- lapply(files, sf::read_sf) %>% dplyr::bind_rows()}

  #ensure basins are in the same crs as the raster
  basins <- sf::st_transform(basins, terra::crs(data))

  #get stats for each basin
  basin_vals <- pbapply::pbsapply(seq_len(nrow(basins)), function(x) {
    basin <- terra::vect(basins[x, ])
    data_crop <- terra::crop(data, basins[x,])

    if(calc == "percent"){
      val_sum <- as.numeric(unname(terra::zonal(data_crop, basin, fun = "sum")))
      val_bsn <- as.numeric(unname(terra::zonal(data_crop, basin, fun = "notNA")))

      val <- val_sum / val_bsn #calculate the percentage
    }else{
      val <- as.numeric(unname(terra::zonal(data_crop, basin, fun = calc)))
    }


    return(val)
  })

  basin_stats <- data.frame(siteID = basins$siteID)
  basin_stats[[data_name]] <- basin_vals

  return(basin_stats)
}

#' Interpolate missing MTBS dNBR data
#'
#' Sometimes a raster layer will have areas of "Non-processing" due to data gaps in the Landsat data used to determine the dNBR. While we could just remove
#' these layers, sometimes you don't want to lose that data. This procedure aims to fill those gaps using interpolation using MBA::mba.points() which
#' uses a bivariate scatter of points using multilevel B-splines to interpolate.
#'
#' @param raster Object of class `SpatRaster` to interpolate.
#' @param plots Logical, should plots showing the original raster and the interpolated raster be printed?
#'
#' @returns Object of class `SpatRaster` with the data gaps interpolated.
#' @export
#' @md
#'
#' @examples
fill_mtbs <- function(raster, plots = TRUE){
  stopifnot(inherits(raster, "SpatRaster"))
  library(raster)
  library(pracma)
  library(MBA)
  library(dplyr)

  #interpolate NA values
  # Convert raster to a data frame of coordinates and values
  xyz <- as.data.frame(raster, xy = TRUE)

  #remove NA values (MTBS sets to -32768)
  xy.est <- xyz[xyz[,3] == -32768, 1:2]
  xyz <- xyz[xyz[,3] != -32768,]
  colnames(xyz) <- c("x", "y", "z")

  # Perform interpolation using mba.points
  interp_result <- MBA::mba.points(
    xyz, #existing data
    xy.est, #missing data
  )

  # Convert interpolated results back to a raster
  raster_filled <- rbind(xyz, interp_result$xyz.est)
  raster_filled <- terra::rast(raster_filled, crs=crs(raster))


  if(plots){
    #show comparison
    raster_plot <- raster
    NAflag(raster_plot) <- -32768
    p1 <- ggplot2::ggplot() + tidyterra::geom_spatraster(data=raster_plot)
    p2 <- ggplot2::ggplot() + tidyterra::geom_spatraster(data=raster_filled)

    ggpubr::ggarrange(p1, p2, nrow=1, common.legend = TRUE)
  }

  return(raster_filled)
}

#' Download a zip file via URL or path
#'
#' Downloads the file to the temporary directory, and unzips
#' it to the specified location.
#'
#' @param path URL or file path to a .zip file to download.
#' @param exdir Path to the folder to put unzipped files into.
#' @param junkpaths If `TRUE`, use only the basename of the stored filepath when extracting. The equivalent of unzip -j.
#'
#' @returns Files downloaded from URL to the `exdir`
#' @export
#'
#' @examples
download_zip <- function(path, exdir, junkpaths = FALSE){
  if(file.exists(path)){
    tmp_name <- path
  }else{
    tmp_name <- tempfile(fileext=".zip")
    download.file(url = path, destfile = tmp_name, mode = "wb")
  }

  unzip(tmp_name, exdir = exdir, junkpaths=junkpaths)
}
