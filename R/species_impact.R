#' Polygon rotation
#'
#' \code{poly_rotate} is a helper function that rotates a polygon object around
#' its centroid. Based on code from \href{https://geocompr.robinlovelace.net/geometric-operations.html}{Geocomputation with R}.
#'
#' @param a Numeric value indicating the desired amount of rotation, in degrees.
#'
#' @return A \code{data.frame} with the coordinates of the rotated polygon coordinates.
#'
poly_rotate <- function(a){
  r = a * pi / 180 ## Convert degrees to radians
  matrix(c(cos(r), sin(r), -sin(r), cos(r)), nrow = 2, ncol = 2)
}


#' Gravel pad polygon creation
#'
#' \code{pt_to_pad} is a helper function that converts a point location to a polygon
#' by creating a square buffer of a desired area around the point (which is taken to
#' be the centroid). Used to represent gravel pad area of central processing facility
#' (CPF) and satellite production pads. Used by \link{footprint_generation}.
#'
#' @param x \code{sf POINT} object providing pad centroid location(s).
#' @param area Numeric value indicating the desired area of the final polygon object,
#'   in square meters.
#' @param proj.info Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
#'   common to all spatial objects in the analysis.
#'
#' @return An \code{sf POLYGON} object representing the footprint of CPF or satellite
#'   gravel pads.
#'
pt_to_pad <- function(x, area, proj.info){
  ## Buffer to create a square polygon of the desired area
  buf1 <- sf::st_buffer(x, dist=sqrt(2*(sqrt(area)/2)^2), nQuadSegs=1, endCaptStyle="SQUARE")
  ## Rotate the polygon so that square edges align with the orientation of the raster pixels
  buf.sfc <- sf::st_geometry(buf1)
  buf.rotate <- (buf.sfc - sf::st_centroid(buf.sfc)) * poly_rotate(45) + sf::st_centroid(buf.sfc)
  buf2 <- sf::st_set_geometry(buf1, buf.rotate)
  ## Assign the desired projection
  buf2 <- projection_alignment(buf2, proj.info)
  return(buf2)
}


#' Convert point/line infrastructure to spatial footprints
#'
#' \code{footprint_generation} is a helper function that converts point and line
#' infrastructure locations, created by \code{\link{generate_cpf}} and \code{\link{generate_sat_rd}},
#' to spatial footprints of surface disturbance. This is used for each
#' species-specific impact analysis (\code{\link{impact_caribou}}, \code{\link{impact_shorebird}},
#' and \code{\link{impact_brant}}).
#'
#' @param cpf.sf \code{sf POINT} object indicating central processing
#'   facility (CPF) locations. Created by \code{\link{generate_cpf}}.
#' @param sat.sf \code{sf POINT} object indicating satellite production
#'   facility locations. Created by \code{\link{generate_sat_rd}}.
#' @param rd.sf \code{SF LINESTRING} object indicating gravel road centerlines. Created
#'   by \code{\link{generate_sat_rd}}.
#' @param proj.info Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
#'   common to all spatial objects in the analysis.
#' @param area.cpf Numeric value indicating the average size (sq. m) of a CPF
#'   footprint. Default value equates to 100 acres based on BLM (2019b).
#' @param area.sat Numeric value indicating the average size (sq. m) of a
#'   satellite production facility footprint. Default value equates to 15 acres based
#'   on BLM (2019a) p.B-6.
#' @param road.width Numeric value indicating the estimated road ground
#'   footprint width (m). Default value based on BLM (2019a) p.B-6. See Fullman
#'   et al. (in press) for details.
#'
#' @return List containing three \code{sf POLYGON} objects: \code{cpf_footprints},
#'   \code{sat_footprints}, and \code{surf_disturb}, which contain the CPF footprints,
#'   satellite pad footprints, and combined footprint of CPFs, satellite pads, and roads,
#'   respectively.
#'
#' @references BLM \[Bureau of Land Management\] 2019a. National Petroleum Reserve
#'   in Alaska Draft Integrated Activity Plan and Environmental Impact Statement.
#'   Bureau of Land Management, U.S. Department of the Interior, Anchorage, AK,
#'   USA.
#'
#' BLM 2019b. Willow Master Development Plan Draft Environmental Impact
#'   Statement. Bureau of Land Management, U.S. Department of the Interior.
#'   Anchorage, AK, USA.
#'
#' Fullman TJ, Sullender BK, Cameron MD, Joly K. in press. Simulation modeling
#'   accounts for uncertainty while quantifying ecological effects of development
#'   alternatives. Ecosphere.
footprint_generation <- function(cpf.sf, sat.sf, rd.sf, proj.info, area.cpf = 100*4046.86,
                                 area.sat = 15*4046.86, road.width = 18.85953){

  ## Create a square buffer around CPF locations with size based on CPF area
  cpf.footprints <- pt_to_pad(x=cpf.sf, area=area.cpf, proj.info=proj.info)

  ## Do the same for satellite pads
  sat.footprints <- pt_to_pad(x=sat.sf, area=area.sat, proj.info=proj.info)

  ## Create road footprints by buffering roads by 1/2 the road width since the road extends on either side of the
  ## line
  rd.surface <- sf::st_buffer(x=rd.sf, dist=(road.width/2))

  ## Combine these surface area polygons into a single surface disturbance layer
  surf.disturb <- sf::st_union(rbind(cpf.footprints[,ncol(cpf.footprints)], sat.footprints[,ncol(sat.footprints)], rd.surface[,ncol(rd.surface)]))

  return(list("cpf_footprints"=cpf.footprints, "sat_footprints"=sat.footprints, "surf_disturb"=surf.disturb))
}



#' Calculate quantiles and identify high-quality pixels for species raster
#'
#' \code{calc_highquality} is a helper function that prepares input data for the
#'   caribou and shorebird impact analyses (\code{\link{impact_caribou}} and
#'   \code{\link{impact_shorebird}}) by calculating quantiles of an input raster
#'   and determining the number of high-quality pixels (\emph{sensu} Johnson et
#'   al. 2005). See Fullman et al. (in press) for details.
#'
#' @param x A caribou calving resource selection function (RSF) or shorebird
#'   habitat suitability index (HSI) \code{SpatRaster} object.
#' @param y Western Arctic Herd (WAH) caribou habitat overlap weighting
#'   \code{SpatRaster} object (optional).
#' @param z Numeric value indicating the species-specific threshold value
#'   defining shorebird suitable habitat from habitat suitability values (from
#'   Saalfeld et al. 2013). Optional, only used if \code{sb == TRUE}.
#' @param wah Logical indicator of whether the run is for the WAH. If so, adds
#'   the additional weighting analysis. Defaults to \code{FALSE}.
#' @param sb Logical indicator of whether the run is for a shorebird. If so,
#'   adds a calculation of suitable habitat before calculating quantiles.
#'   Defaults to \code{FALSE}.
#' @param hq.quant Numeric value between 0-1 that indicates the quantile value that
#'   reflects high-quality habitat. Defaults to 0.75, following Johnson et al. (2005).
#'
#' @return \code{data.frame} object containing the raster pixel value denoting
#'   high-quality habitat according to the desired quantile value (\code{quant_hq}),
#'   number of high-quality pixels in the input raster (\code{highquality_orig}),
#'   and (if \code{wah = TRUE}) the weighted high-quality value (\code{highquality_weighted}).
#'
#' @references Fullman TJ, Sullender BK, Cameron MD, Joly K. in press. Simulation modeling
#'   accounts for uncertainty while quantifying ecological effects of development
#'   alternatives. Ecosphere.
#'
#' Johnson CJ, Boyce MS, Ray CL, Cluff DH, Gau RJ, Gunn A, Mulders R. 2005.
#'   Cumulative effects of human developments on arctic wildlife. Wildlife
#'   Monographs 160:1-36.
#'
#' Saalfeld ST, Lanctot RB, Brown SC, Saalfeld DT, Johnson JA, Andres BA, Bart JR.
#'   2013. Predicting breeding shorebird distributions on the Arctic Coastal
#'   Plain of Alaska. Ecosphere 4:16.
calc_highquality <- function(x, y = NULL, z = NULL, wah = FALSE, sb = FALSE, hq.quant = 0.75){
  ## If the run is for a shorebird, calculate suitable habitat using the species-specific HSI threshold.
  if(sb) terra::values(x)[terra::values(x) <= z] <- NA

  ## Calculate the high-quality quantile value for the species raster
  quant <- stats::quantile(terra::values(x), probs=hq.quant, na.rm=TRUE)

  ## Calculate the number of pixels of "high quality" habitat (values in the upper quartile, following
  ## Johnson et al. 2005).
  highquality.orig <- length(stats::na.omit(terra::values(x))[stats::na.omit(terra::values(x)) > quant])
  quant.out <- data.frame(quant_hq = quant, highquality_orig = highquality.orig)

  ## If the run is for the WAH, calculate the weighting information too
  if(wah){
    wah.hq <- x
    terra::values(wah.hq)[terra::values(wah.hq) <= quant.out$quant_hq] <- 0
    terra::values(wah.hq)[terra::values(wah.hq) > quant.out$quant_hq] <- 1
    ## Extract overlap values, weighting values, and cell numbers at upper quartile locations
    wt.df <- data.frame("cell"=terra::cells(wah.hq, terra::ext(wah.hq)), "hq"=as.numeric(terra::values(wah.hq)),
                        "weight"=as.numeric(terra::values(y)))
    ## Constrict this to just the high-quality pixels
    wt.nona <- stats::na.omit(wt.df)
    ## Add the result to quant.out
    quant.out$highquality_weighted <- sum(wt.nona$weight)
  }

  return(quant.out)
}



#' Convert raster values to zero where they overlap a SpatialPolygons object
#'
#' \strong{The \code{infrastructure_overlap_to_zero} is deprecated as of DIA version 0.2.0.
#' Please use \code{\link{raster_to_zero}} instead.}   \code{infrastructure_overlap_to_zero}
#'  is a helper function for the caribou impact analysis (\code{\link{impact_caribou}})
#'  that converts habitat values to 0 where they overlap the footprint of development.
#'
#' @param x Caribou calving resource selection function (RSF) RasterLayer object.
#' @param surf.disturb SpatialPolygonsDataFrame object representing the combined
#'   footprint of central processing facilities (CPFs), satellite pads, and
#'   roads. Created by \code{\link{footprint_generation}}.
#'
#' @return RasterLayer object with values of zero everywhere overlapped by
#'   \code{surf.disturb}.
#'
infrastructure_overlap_to_zero	<- function(x, surf.disturb){
  .Deprecated("raster_to_zero")
  ## Identify the cell numbers of the calving rasters that are overlapped by the physical footprint of
  ## infrastructure (i.e., surf.disturb) and convert those raster pixels to have a value of zero (as
  ## long as they are not NA).
  sd.cells <- raster::extract(x, surf.disturb, small=TRUE, cellnumbers=TRUE)
  for(o in 1:length(sd.cells)){
    cells.tmp <- sd.cells[[o]]
    cells.tmp2 <- stats::na.omit(cells.tmp)
    cells.tmp3 <- raster::rowColFromCell(x, cells.tmp2[,1])
    x[cells.tmp3] <- 0
  }
  return(x)
}



#' Discount raster values based on proximity to development and calculate remaining "high quality" habitat
#'
#' \code{infrastructure_proximity_discounting} is a helper function for the
#'   caribou impact analysis (\code{\link{impact_caribou}}) that discounts
#'   resource selection function (RSF) values that do not directly overlap the
#'   surface footprint of infrastructure based on their proximity to development.
#'   Then calculates the amount of remaining "high quality" habitat (\emph{sensu}
#'   Johnson et al. 2005).
#'
#' @param x \code{SpatRaster} object indicating the caribou calving RSF values,
#'   updated to set values overlapping the development footprint to zero. Created
#'   using \code{\link{raster_to_zero}}.
#' @param surf.disturb \code{sf POLYGON} object representing the combined
#'   footprint of central processing facilities (CPFs), satellite pads, and
#'   roads. Created by \code{\link{footprint_generation}}.
#' @param quant Numerical value identifying the threshold for high-quality pixels of the calving
#'   RSF raster. Defaults to the upper quartile of values. Calculated by \code{\link{calc_highquality}}.
#' @param wah.hq.weight Western Arctic Herd (WAH) caribou habitat overlap
#'   \code{SpatRaster} object (optional, unless analyses are being run for the WAH, in
#'   which case it is required).
#' @param zoi Numeric value indicating the zone of influence (m) in which habitat
#'   quality is discounted. Defaults to 4000 m, following Wilson et al. (2013)
#'   and Cameron et al. (2005).
#' @param wah Logical indicator of whether the run is for the WAH. If so, runs
#'   the additional weighting analysis. Defaults to \code{FALSE}.
#' @param caribou.out Logical indicator of whether the resulting discounted
#'   caribou RSF raster(s) should be written out. Defaults to \code{FALSE}.
#' @param proj.info Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
#'   common to all spatial objects in the analysis.
#' @param wd.loc Character string identifying the base directory containing both
#'   input and output data folders.
#' @param path.out Character string identifying the relative path from the base
#'   directory (\code{wd.loc}) to where output data will be saved. Defaults to
#'   an "Output_Data" folder within the base folder.
#' @param scenario Character vector of strings identifying the scenarios to be run.
#' @param n.iter Integer indicating the desired number of iterations to be run
#'   for each scenario. Defaults to 100.
#' @param z Iterator value used by \code{\link{dia}} to specify the development
#'   scenario being analyzed.
#' @param i Iterator value used by \code{\link{dia}} to run development simulation
#'   and impacts analyses in parallel, indicating the specific iteration being run.
#'
#' @return Integer value specifying the number of high-quality habitat pixels
#'   (weighted value if \code{wah = TRUE}) remaining after discounting.
#'
#' @references Cameron RD, Smith WT, White RG, Griffith B. 2005. Central Arctic
#'   caribou and petroleum development: Distributional, nutritional, and
#'   reproductive implications. Arctic 58:1-9.
#'
#' Johnson CJ, Boyce MS, Ray CL, Cluff DH, Gau RJ, Gunn A, Mulders R. 2005.
#'   Cumulative effects of human developments on arctic wildlife. Wildlife
#'   Monographs 160:1-36.
#'
#' Wilson RR, Liebezeit JR, Loya WM. 2013. Accounting for uncertainty in oil and
#'   gas development impacts to wildlife in Alaska. Conservation Letters
#'   6:350-358.
infrastructure_proximity_discounting <- function(x, surf.disturb, quant, wah.hq.weight = NULL, zoi = 4000,
                                                 wah = FALSE, caribou.out = FALSE, proj.info, wd.loc, path.out=path.out, scenario, n.iter = 100, z = NULL, i = NULL){
  ## Discounting is only done within 4 km of development so calculate distance to development within this
  ## radius
  sd.buf <- sf::st_buffer(x=surf.disturb, dist=zoi)
  sd.mask <- terra::mask(x, mask=terra::vect(sd.buf))
  sd.mask <- projection_alignment(sd.mask, proj.info)
  dist.tmp <- rgeos::gDistance(spgeom1=sf::as_Spatial(surf.disturb), spgeom2=methods::as(raster::raster(sd.mask), "SpatialPoints"), byid=TRUE)
  dist.ras <- sd.mask
  terra::values(dist.ras)[!is.na(terra::values(dist.ras))] <- dist.tmp[,1]

  ## Calculate the habitat discounting adapting the logistic function from Wilson et al. 2013, based on
  ## Cameron et al. 2005 Fig. 3.
  discount.ras <- (exp(-2.9741+.001079*dist.ras)/(1+exp(-2.9741+.001079*dist.ras))) * sd.mask
  calving.discount <- raster::cover(discount.ras, x)

  ## Output the discount raster, if desired
  if(caribou.out) raster::writeRaster(calving.discount, file=paste(wd.loc, "/", path.out, "/", ifelse(wah, "WAH", "TCH"), "_discounted_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), ".tif", sep=""))

  ## Calculate the remaining amount of "high quality" habitat
  highquality.discount <- length(stats::na.omit(terra::values(calving.discount))[stats::na.omit(terra::values(calving.discount)) > quant])

  ## If the run is for the WAH, run the weighting analysis
  if(wah){

    ## Create a raster that gives all remaining WAH high-quality pixels a value of 1 and makes all other
    ## pixels NA
    wah.hq.remain <- calving.discount
    terra::values(wah.hq.remain)[terra::values(wah.hq.remain) <= quant] <- NA
    terra::values(wah.hq.remain)[terra::values(wah.hq.remain) > quant] <- 1

    ## Identify the cell numbers of the remaining high-quality WAH pixels
    wah.remain.all <- data.frame("cell"=terra::cells(wah.hq.remain, terra::ext(wah.hq.remain)),
                                 "hq"=as.numeric(terra::values(wah.hq.remain)), "weight"=as.numeric(terra::values(wah.hq.weight)))

    ## Constrict this down to just the high-quality pixels
    wah.remain.hq <- stats::na.omit(wah.remain.all)

    ## Calculate the weighted sum
    weighted.sum <- sum(wah.remain.hq$weight)

    ## Give it the same name for consistent output
    highquality.discount <- weighted.sum
  }

  return(highquality.discount)
}



#' Caribou impact analysis
#'
#' \code{impact_caribou} runs the caribou impact analysis, discounting caribou
#'   calving habitat values in proximity to infrastructure and calculating the
#'   amount of "high quality" habitat remaining (\emph{sensu} Johnson et al. 2005).
#'
#' @param tch.raster Character string giving the file name of the Teshekpuk
#'   Caribou Herd (TCH) calving raster. Optional, allowing the analysis to
#'   be run for only one herd.
#' @param wah.raster Vector of two character strings, specifying the file names
#'   of the Western Arctic Herd (WAH) calving raster and weighted overlap
#'   raster, respectively. Optional, allowing the analysis to only be run
#'   for one herd.
#' @param tch.calving,wah.calving \code{SpatRaster} object(s) depicting the calving resource
#'   selection function (RSF) for the TCH and WAH, respectively. Optional,
#'   allowing the analysis to be run for only one herd.
#' @param out.df \code{data.frame} object containing the scenario information and
#'   infrastructure summary. Created during the impact analysis by \code{\link{dia}}.
#' @param tch.hq,wah.hq \code{data.frame} object(s) created by \code{\link{calc_highquality}}.
#' @inheritParams infrastructure_proximity_discounting
#'
#' @return \code{out.df data.frame}, updated to contain the caribou impact results.
#'
#' @references Cameron RD, Smith WT, White RG, Griffith B. 2005. Central Arctic
#'   caribou and petroleum development: Distributional, nutritional, and
#'   reproductive implications. Arctic 58:1-9.
#'
#' Johnson CJ, Boyce MS, Ray CL, Cluff DH, Gau RJ, Gunn A, Mulders R. 2005.
#'   Cumulative effects of human developments on arctic wildlife. Wildlife
#'   Monographs 160:1-36.
#'
#' Wilson RR, Liebezeit JR, Loya WM. 2013. Accounting for uncertainty in oil and
#'   gas development impacts to wildlife in Alaska. Conservation Letters
#'   6:350-358.
impact_caribou <- function(surf.disturb, tch.raster = NULL, wah.raster = NULL, tch.calving = NULL,
                           wah.calving = NULL, wah.hq.weight = NULL, out.df, wd.loc, path.out, proj.info, tch.hq = NULL, wah.hq = NULL,
                           zoi = 4000, caribou.out = FALSE, scenario, n.iter = 100, z = NULL, i = NULL){

  ## Run TCH analysis, if desired
  if(!is.null(tch.raster)){

    ## Convert RSF scores to 0 where they overlap footprints of CPFs, satellite pads, or roads
    tch.calving <- raster_to_zero(tch.calving, terra::vect(surf.disturb))

    ## Discount non-overlapping RSF scores based on their proximity to development
    tch.highquality.discount <- infrastructure_proximity_discounting(x=tch.calving, surf.disturb=surf.disturb,
                                                                     quant=tch.hq$quant_hq, wah=FALSE, proj.info=proj.info, zoi=zoi, caribou.out=caribou.out,
                                                                     wd.loc=wd.loc, path.out=path.out, z=z, i=i)

    ## Update the output data.frame
    out.df$tch_discount <- tch.highquality.discount
    out.df$tch_remaining <- tch.highquality.discount / tch.hq$highquality_orig
  }

  ## Run WAH analysis, if desired
  if(!is.null(wah.raster)){

    ## Convert RSF scores to 0 where they overlap footprints of CPFs, satellite pads, or roads
    wah.calving <- raster_to_zero(wah.calving, terra::vect(surf.disturb))

    ## Discount non-overlapping RSF scores based on their proximity to development
    wah.highquality.discount <- infrastructure_proximity_discounting(x=wah.calving, surf.disturb=surf.disturb,
                                                                     quant=wah.hq$quant_hq, wah.hq.weight=wah.hq.weight, wah=TRUE, zoi=zoi, proj.info=proj.info,
                                                                     caribou.out=caribou.out, wd.loc=wd.loc, path.out=path.out, z=z, i=i)

    ## Update the output data.frame
    out.df$wah_discount <- wah.highquality.discount
    out.df$wah_remaining <- wah.highquality.discount / wah.hq$highquality_weighted
  }

  return(out.df)
}


#' Shorebird impact analysis
#'
#' \code{impact_shorebird} runs the shorebird impact analysis, discounting
#'   shorebird suitable habitat values in proximity to infrastructure and
#'   calculating the amount of "high quality" suitable habitat remaining
#'   (\emph{sensu} Johnson et al. 2005).
#'
#' @param sb.stack \code{SpatRaster} object with one layer for each shorebird habitat
#'   suitability index (HSI) raster, based on data from Saalfeld et al. (2013) that
#'   are loaded and combined in the data preparation code of \code{\link{dia}}.
#' @param rd.sf \code{sf LINESTRING} object indicating gravel road centerlines. Created
#'   by \code{\link{generate_sat_rd}}.
#' @param sb.df \code{data.frame} containing the species-specific raster pixel threhsold
#'   values indicating high-quality suitable habitat (defaults to upper quantile of values)
#'   and the number of high-quality suitable pixels in the input data. Created in the
#'   data preparation code of \code{\link{dia}}.
#'
#' @return Temporary \code{data.frame} containing the shorebird impact results.
#'
#' @references Johnson CJ, Boyce MS, Ray CL, Cluff DH, Gau RJ, Gunn A, Mulders R. 2005.
#'   Cumulative effects of human developments on arctic wildlife. Wildlife
#'   Monographs 160:1-36.
#'
#' Saalfeld ST, Lanctot RB, Brown SC, Saalfeld DT, Johnson JA, Andres BA, Bart JR.
#'   2013. Predicting breeding shorebird distributions on the Arctic Coastal
#'   Plain of Alaska. Ecosphere 4:16.
impact_shorebird <- function(sb.stack, rd.sf, sb.df){
  ## Identify impacted pixels using the road centerlines sf object created by the infrastructure
  ## simulation. This is used rather than surf.disturb due to the relatively coarse resolution of the
  ## shorebird data (1 km resolution) compared to the finer scale at which infrastructure data were generated
  ## (120 m resolution).
  sb.disturb.hq <- terra::extract(x=sb.stack, y=terra::vect(sf::as_Spatial(rd.sf)))
  ## Remove any pixels with NAs, as these will have NAs for all species.
  sb.disturb.hq <- stats::na.omit(sb.disturb.hq)
  ## Update sb.df with the number of remaining high-quality suitable pixels by determining the number of
  ## affected high-quality pixels and subtracting this from the original number of high-quality pixels.
  for(ii in 1:nrow(sb.df)){
    sb.df$hq_discount[ii] <- sb.df$hq_orig[ii] - length(sb.disturb.hq[,ii+1][sb.disturb.hq[,ii+1] > sb.df$quant_hq[ii]])
  }
  ## Update sb.df with the proportion of remaining high-quality suitable pixels
  sb.df$hq_remaining <- sb.df$hq_discount / sb.df$hq_orig

  return(sb.df)
}



#' Brant impact analysis
#'
#' \code{impact_brant} runs the brant impact analysis, calculating the area of
#'   expected brant disturbance on molting lakes, based on physical infrastructure
#'   and helicopter overflights, then estimating the habitat area and number of brant
#'   affected.
#'
#' @param brant.lakes \code{sf POLYGON} object containing lake boundaries
#'   with molting brant density data from Schults and Zeller (2019). Created in
#'   the data preparation code of \code{\link{dia}}.
#' @param footprints List of three \code{sf POLYGON} objects: \code{cpf_footprints},
#'   \code{sat_footprints}, and \code{surf_disturb}. Created by \code{\link{footprint_generation}}.
#' @param land.dist Numeric value indicating the distance (m) from central processing
#'   facilities (CPFs) or satellite pads at which helicopters are expected to be under
#'   500 m altitude, based on a descent angle of 10 degrees following FAA guidelines. Used
#'   to determine helicopter impact area around CPFs and satellite pads. Defaults
#'   to 2835.64 m based on the assumption that a helicopter approaching a landing
#'   zone at a 10 degree descent angle reaches a 500 m altitude at \eqn{500m / tan(10 deg) = 2835.64 m}
#'   out from the landing zone.
#' @param heli.disturb Numeric value indicating the distance (m) from helicopters
#'   (under 500 m altitude) at which molting brant are expected to be disturbed.
#'   Defaults to 3570 m based on Jensen (1990) and Miller et al. (1994).
#' @param proximity.effect Numeric value indicating the distance (m) at which
#'   birds are disturbed by infrastructure proximity effects (e.g., dust
#'   deposition). Defaults to 100 m as a conservative estimate.
#' @param proj.info Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
#'   common to all spatial objects in the analysis.
#' @param brant.out Logical indicator of whether a shapefile of the affected lakes should
#'   be written to \code{path.out}. Only written if any lakes are affected. Defaults to
#'   \code{FALSE}.
#' @inheritParams impact_caribou
#'
#' @return Temporary \code{data.frame} containing the brant impact results. If
#'   \code{brant.out =TRUE}, then a shapefile of the affected lakes will also be written
#'   out, if any exist.
#'
#' @references Jensen KC. 1990. Responses of molting Pacific black brant to
#'   experimental aircraft disturbance in the Teshekpuk Lake Special Area, Alaska.
#'   Ph.D. Thesis, Texas A&M University. College Station, TX, USA.
#'
#'   Miller MW, Jensen KC, Grant WE, Weller MW. 1994. A simulation model of
#'   helicopter disturbance of molting Pacific black brant. Ecological Modelling
#'   73:293-309.
#'
#'   Shults BS, Zeller TK. 2019. Abundance and Distribution of Molting
#'   Geese in the Vicinity of Teshekpuk Lake, Alaska, July 2018. Migratory Bird
#'   Management, U.S. Fish and Wildlife Service. Anchorage, AK, USA, 18 pp.
impact_brant <- function(brant.lakes, footprints, land.dist = 2835.64, heli.disturb = 3570, proximity.effect = 100,
                         proj.info, brant.out = FALSE, wd.loc, path.out, scenario,
                         z, i, n.iter){
  ## Calculate the disturbance buffer around CPFs and satellites based on the assumption that
  ## helicopters will stay above 500 m, except when taking off or landing, but may approach CPFs and
  ## satellites from any direction
  cpfBuff <- sf::st_buffer(x=footprints$cpf_footprints, dist=land.dist+heli.disturb)
  satBuff <- sf::st_buffer(x=footprints$sat_footprints, dist=land.dist+heli.disturb)

  ## Account for proximity effects around the physical footprint of infrastructure, adding this to the
  ## helicopter disturbance polygons.
  sd.disturb.buf <- sf::st_buffer(x=footprints$surf_disturb, dist=proximity.effect)
  all.disturb.buf <- sf::st_union(rbind(cpfBuff[,ncol(cpfBuff)], satBuff[,ncol(satBuff)], sd.disturb.buf))

  ## Identify which portions of lakes intersect the disturbance buffer.
  disturbed.lakes.all <- rgeos::gIntersection(sf::as_Spatial(brant.lakes), sf::as_Spatial(all.disturb.buf), byid=TRUE)

  ## If the infrastructure scenario does not overlap the brant lakes, enter zeros for impact,
  ## otherwise calculate the lake area and estimated number of brant disturbed.
  if(class(disturbed.lakes.all) == "NULL"){
    tmp.df <- data.frame("lakes_disturbed_all_ha" = 0, "lakes_disturbed_molt_ha" = 0, "disturbed_geese" = 0)
  } else{
    ## The names of disturbed.lakes.all are the name of the LAKE column in the brant.lakes shapefile and "ID1" for
    ## all.disturb.buf. Use this to relate the brant.lakes attributes to the intersected lakes to allow
    ## calculation of impact.
    tmp.df.disturb <- data.frame(LAKE=as.numeric(substr(names(disturbed.lakes.all), start=1,
                                                        stop=nchar(names(disturbed.lakes.all))-4)),
                                 affected_area_sq_m=sf::st_area(sf::st_as_sf(disturbed.lakes.all)))
    tmp.df2.disturb <- merge(tmp.df.disturb, brant.lakes)

    ## Add the area for all disturbed lakes and those in which brant have been observed molting over
    ## the last 5 years to a tmp.df object. Multiply area by bird density to estimate the number of
    ## disturbed birds. Convert areas to units of hectares, for easier use.
    tmp.df <- data.frame("lakes_disturbed_all_ha"=units::set_units(sum(tmp.df2.disturb$affected_area_sq_m), "ha"),
                         "lakes_disturbed_molt_ha"=units::set_units(sum(tmp.df2.disturb$affected_area_sq_m[tmp.df2.disturb$X5Yr_Av > 0]), "ha"),
                         "disturbed_geese"=units::drop_units(round(sum(tmp.df2.disturb$affected_area_sq_m * tmp.df2.disturb$BLBRpM2),0)))

    ## Save the affected lakes as a shapefile, if desired
    if(brant.out){
      disturbed.shp <- sf::st_as_sf(disturbed.lakes.all)
      disturbed.shp2 <- sf::st_as_sf(cbind(tmp.df.disturb, disturbed.shp))
      sf::st_write(disturbed.shp2, dsn=paste(wd.loc, path.out, sep="/"), layer=paste("Brant_impacted_molting_lakes_shapefile_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), sep=""), driver="ESRI Shapefile")
    }
  }

  return(tmp.df)
}
