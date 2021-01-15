#' Convert point/line infrastructure to spatial footprints
#'
#' \code{footprint_generation} is a helper function that converts point and line
#' infrastructure locations, created by \code{\link{generate_cpf}} and \code{\link{generate_sat_rd}},
#' to spatial footprints of surface disturbance. This is used for each
#' species-specific impact analysis (\code{\link{impact_caribou}}, \code{\link{impact_shorebird}},
#' and \code{\link{impact_brant}}).
#'
#' @param cpf.spdf SpatialPointsDataFrame object indicating central processing
#'   facility (CPF) locations. Created by \code{\link{generate_cpf}}.
#' @param sat.spdf SpatialPointsDataFrame object indicating satellite production
#'   facility locations. Created by \code{\link{generate_sat_rd}}.
#' @param rd.sl SpatialLines object indicating gravel road centerlines. Created
#'   by \code{\link{generate_sat_rd}}.
#' @param proj.info Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
#'   common to all spatial objects in the analysis.
#' @param area.cpf Numeric value indicating the average size (sq. m) of a CPF
#'   footprint. Default value estimated based on BLM (2019b).
#' @param area.sat Numeric value indicating the average size (sq. m) of a
#'   satellite production facility footprint. Default value from BLM (2019a)
#'   p.B-6.
#' @param road.width Numeric value indicating the estimated road ground
#'   footprint width (m). Default value based on BLM (2019a) p.B-6. See Fullman
#'   et al. (in review) for details.
#'
#' @return List containing three SpatialPolygonsDataFrame objects:
#'   \code{cpf_footprints}, \code{sat_footprints}, and \code{surf_disturb}, which
#'   contains the combined footprint of CPFs, satellite pads, and roads.
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
#' Fullman TJ, Sullender BK, Cameron MD, Joly K. in review. Simulation modeling
#'   accounts for uncertainty while quantifying ecological effects of development
#'   alternatives. Ecosphere.
footprint_generation <- function(cpf.spdf, sat.spdf, rd.sl, proj.info, area.cpf = 100*4046.86,
                                 area.sat = 15*4046.86, road.width = 18.85953){

  ## Convert CPF locations to a SpatialPoints object for easier indexing during buffering
  tmp.cpf.sp <- sp::SpatialPoints(cpf.spdf)

  ## Create a square buffer around CPF locations with size based on CPF area
  cpf.bufs <- list()
  for(m in 1:length(tmp.cpf.sp)){
    cpf.bufs[[m]] <- rgeos::gBuffer(tmp.cpf.sp[m], width=sqrt(area.cpf)/2, quadsegs=1, capStyle="SQUARE")
  }
  cpf.footprints <- do.call(raster::bind, cpf.bufs)

  ## Repeat the above process for satellite pads
  tmp.sat.sp <- sp::SpatialPoints(sat.spdf)
  sat.bufs <- list()
  for(n in 1:length(tmp.sat.sp)){
    sat.bufs[[n]] <- rgeos::gBuffer(tmp.sat.sp[n], width=sqrt(area.sat)/2, quadsegs=1, capStyle="SQUARE")
  }
  sat.footprints <- do.call(raster::bind, sat.bufs)

  ## Repeat for roads (Using 1/2 the road width since the road extends on either side of the line)
  rd.surface <- raster::buffer(x=rd.sl, width=(road.width/2))

  ## Combine these surface area polygons into a single surface disturbance layer. To do this I need to
  ## make sure all have identical projection information and alter the row.names of the satellite
  ## records to avoid a conflict.
  cpf.footprints <- projection_alignment(cpf.footprints, proj.info)
  sat.footprints <- projection_alignment(sat.footprints, proj.info)
  row.names(sat.footprints) <- paste(as.character(length(cpf.footprints)+1):(length(cpf.footprints)+length(sat.footprints)))
  surf.disturb <- raster::aggregate(rbind(cpf.footprints, sat.footprints, rd.surface))
  return(list("cpf_footprints"=cpf.footprints, "sat_footprints"=sat.footprints, "surf_disturb"=surf.disturb))
}



#' Calculate quantiles and identify high-quality pixels for species raster
#'
#' \code{calc_highquality} is a helper function that prepares input data for the
#'   caribou and shorebird impact analyses (\code{\link{impact_caribou}} and
#'   \code{\link{impact_shorebird}}) by calculating quantiles of an input raster
#'   and determining the number of high-quality pixels (\emph{sensu} Johnson et
#'   al. 2005). See Fullman et al. (in review) for details.
#'
#' @param x A caribou calving resource selection function (RSF) or shorebird
#'   habitat suitability index (HSI) RasterLayer object.
#' @param y Western Arctic Herd (WAH) caribou habitat overlap weighting
#'   RasterLayer object (optional).
#' @param z Numeric value indicating the species-specific threshold value
#'   defining shorebird suitable habitat from habitat suitability values (from
#'   Saalfeld et al. 2013). Optional, only used if \code{sb == TRUE}.
#' @param wah Logical indicator of whether the run is for the WAH. If so, adds
#'   the additional weighting analysis. Defaults to \code{FALSE}.
#' @param sb Logical indicator of whether the run is for a shorebird. If so,
#'   adds a calculation of suitable habitat before calculating quantiles.
#'   Defaults to \code{FALSE}.
#'
#' @return data.frame object containing the 75% quantile value (\code{quant75}),
#'   number of high-quality pixels in the input raster (\code{highquality.orig}),
#'   and (if \code{wah = TRUE}) the weighted high-quality value (\code{highquality.weighted}).
#'
#' @references Fullman TJ, Sullender BK, Cameron MD, Joly K. in review. Simulation modeling
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
calc_highquality <- function(x, y = NULL, z = NULL, wah = FALSE, sb = FALSE){
  ## If the run is for a shorebird, calculate suitable habitat using the species-specific HSI threshold.
  if(sb) x[raster::values(x) <= z] <- NA

  ## Calculate the quantiles for the species raster
  quant <- stats::quantile(x)
  ## Calculate the number of pixels of "high quality" habitat (values in the upper quartile, following
  ## Johnson et al. 2005).
  highquality.orig <- length(x[x > quant[4]])
  quant.out <- data.frame(quant75 = quant[4], highquality.orig = highquality.orig)

  ## If the run is for the WAH, calculate the weighting information too
  if(wah){
    wah.hq <- x
    wah.hq[wah.hq <= quant.out$quant75] <- 0
    wah.hq[wah.hq > quant.out$quant75] <- 1
    ## Extract overlap values, weighting values, and cell numbers at upper quartile locations
    wt.df <- data.frame(cell=raster::cellsFromExtent(object=wah.hq, extent=wah.hq),
                        hq=raster::values(wah.hq), weight=raster::values(y))
    ## Constrict this to just the high-quality pixels
    wt.nona <- stats::na.omit(wt.df)
    ## Add the result to quant.out
    quant.out$highquality.weighted <- sum(wt.nona$weight)
  }

  return(quant.out)
}



#' Convert raster values to zero where they overlap a SpatialPolygons object
#'
#'\code{infrastructure_overlap_to_zero} is a helper function for the caribou
#'  impact analysis (\code{\link{impact_caribou}}) that converts habitat values
#'  to 0 where they overlap the footprint of development.
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
#' @param x RasterLayer indicating the caribou calving RSF values, updated to
#'   set values overlapping the development footprint to zero. Created by
#'   \code{\link{infrastructure_overlap_to_zero}}.
#' @param surf.disturb SpatialPolygonsDataFrame object representing the combined
#'   footprint of central processing facilities (CPFs), satellite pads, and
#'   roads. Created by \code{\link{footprint_generation}}.
#' @param quant Numerical value identifying the upper quartile of the calving
#'   RSF raster. Calculated by \code{\link{calc_highquality}}.
#' @param wah.hq.weight Western Arctic Herd (WAH) caribou habitat overlap
#'   RasterLayer object (optional, unless analyses are being run for the WAH, in
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
#'   and impacts analyses in parallel.
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
                                                 wah = FALSE, caribou.out = FALSE, proj.info, wd.loc, path.out, scenario, n.iter = 100, z = NULL, i = NULL){
  ## Discounting is only done within 4 km of development so calculate distance to development within this
  ## radius
  sd.buf <- raster::buffer(x=surf.disturb, width=zoi)
  sd.mask <- raster::mask(x, mask=sd.buf)
  sd.mask <- projection_alignment(sd.mask, proj.info)
  dist.tmp <- rgeos::gDistance(spgeom1=surf.disturb, spgeom2=methods::as(sd.mask, "SpatialPoints"), byid=TRUE)
  dist.ras <- sd.mask
  raster::values(dist.ras)[!is.na(raster::values(dist.ras))] <- apply(dist.tmp, 1, min)

  ## Calculate the habitat discounting adapting the logistic function from Wilson et al. (2013), based on
  ## Cameron et al. (2005) Fig. 3.
  discount.ras <- (exp(-2.9741+.001079*dist.ras)/(1+exp(-2.9741+.001079*dist.ras))) * sd.mask
  calving.discount <- raster::cover(discount.ras, x)

  ## Output the discount raster, if desired
  if(caribou.out) raster::writeRaster(calving.discount, file=paste(wd.loc, "/", path.out, "/", ifelse(wah, "WAH", "TCH"), "_discounted_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), ".tif", sep=""))

  ## Calculate the remaining amount of "high quality" habitat
  highquality.discount <- length(calving.discount[calving.discount > quant])

  ## If the run is for the WAH, run the weighting analysis
  if(wah){

    ## Create a raster that gives all remaining WAH high-quality pixels a value of 1 and makes all other
    ## pixels NA
    wah.hq.remain <- calving.discount
    wah.hq.remain[wah.hq.remain <= quant] <- NA
    wah.hq.remain[wah.hq.remain > quant] <- 1

    ## Identify the cell numbers of the remaining high-quality WAH pixels
    wah.remain.all <- data.frame(cell=raster::cellsFromExtent(object=wah.hq.remain, extent=wah.hq.remain),
                                 hq=raster::values(wah.hq.remain), weight=raster::values(wah.hq.weight))

    ## Constrict this down to just the high-quality pixels
    wah.remain.hq <- wah.remain.all[!is.na(wah.remain.all$hq),]

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
#'   Caribou Herd (TCH) calving RasterLayer. Optional, allowing the analysis to
#'   be run for only one herd.
#' @param wah.raster Vector of two character strings, specifying the file names
#'   of the Western Arctic Herd (WAH) calving RasterLayer and weighted overlap
#'   RasterLayer, respectively. Optional, allowing the analysis to only be run
#'   for one herd.
#' @param tch.calving,wah.calving RasterLayer depicting the calving resource
#'   selection function (RSF) for the TCH and WAH, respectively. Optional,
#'   allowing the analysis to be run for only one herd.
#' @param out.df data.frame object containing the scenario information and
#'   infrastructure summary. Created during the impact analysis by \code{\link{dia}}.
#' @param tch.hq,wah.hq data.frame object(s) created by \code{\link{calc_highquality}}.
#' @inheritParams infrastructure_proximity_discounting
#'
#' @return \code{out.df} data.frame, updated to contain the caribou impact results.
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
    tch.calving <- infrastructure_overlap_to_zero(tch.calving, surf.disturb)

    ## Discount non-overlapping RSF scores based on their proximity to development
    tch.highquality.discount <- infrastructure_proximity_discounting(x=tch.calving, surf.disturb,
                                                                     quant=tch.hq$quant75, wah=FALSE, proj.info=proj.info, zoi=zoi, caribou.out=caribou.out,
                                                                     wd.loc=wd.loc, path.out=path.out, z=z, i=i)

    ## Update the output data.frame
    out.df$tch_discount <- tch.highquality.discount
    out.df$tch_remaining <- tch.highquality.discount / tch.hq$highquality.orig
  }

  ## Run WAH analysis, if desired
  if(!is.null(wah.raster)){

    ## Convert RSF scores to 0 where they overlap footprints of CPFs, satellite pads, or roads
    wah.calving <- infrastructure_overlap_to_zero(wah.calving, surf.disturb)

    ## Discount non-overlapping RSF scores based on their proximity to development
    wah.highquality.discount <- infrastructure_proximity_discounting(x=wah.calving, surf.disturb,
                                                                     quant=wah.hq$quant75, wah.hq.weight=wah.hq.weight, wah=TRUE, zoi=zoi, proj.info=proj.info,
                                                                     caribou.out=caribou.out, wd.loc=wd.loc, path.out=path.out, z=z, i=i)

    ## Update the output data.frame
    out.df$wah_discount <- wah.highquality.discount
    out.df$wah_remaining <- wah.highquality.discount / wah.hq$highquality.weighted
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
#' @param sb.stack RasterStack of shorebird habitat suitability index (HSI) rasters,
#'   based on data from Saalfeld et al. (2013) and updated in the data preparation
#'   code of \code{\link{dia}}.
#' @param rd.sl SpatialLines object indicating gravel road centerlines. Created
#'   by \code{\link{generate_sat_rd}}.
#' @param sb.df data.frame containing the upper quantile of suitable high-quality
#'   habitat values and the number of high-quality suitable pixels in the input
#'   data. Created in the data preparation code of \code{\link{dia}}.
#'
#' @return \code{out.df} data.frame, updated to contain the shorebird impact results.
#'
#' @references Johnson CJ, Boyce MS, Ray CL, Cluff DH, Gau RJ, Gunn A, Mulders R. 2005.
#'   Cumulative effects of human developments on arctic wildlife. Wildlife
#'   Monographs 160:1-36.
#'
#' Saalfeld ST, Lanctot RB, Brown SC, Saalfeld DT, Johnson JA, Andres BA, Bart JR.
#'   2013. Predicting breeding shorebird distributions on the Arctic Coastal
#'   Plain of Alaska. Ecosphere 4:16.
impact_shorebird <- function(sb.stack, rd.sl, sb.df){
  ## Identify impacted pixels using the road centerlines SpatialLines object created by the infrastructure
  ## simulation. This is used rather than surf.disturb due to the relatively coarse resolution of the
  ## shorebird data (1 km resolution) compared to the finer scale at which infrastructure data were generated
  ## (120 m resolution).
  sb.disturb.hq <- raster::extract(x=sb.stack, y=rd.sl, df=TRUE)
  ## Remove any pixels with NAs, as these will have NAs for all species.
  sb.disturb.hq <- stats::na.omit(sb.disturb.hq)
  ## Update sb.df with the number of remaining high-quality suitable pixels by determining the number of
  ## affected high-quality pixels and subtracting this from the original number of high-quality pixels.
  for(ii in 1:nrow(sb.df)){
    sb.df$hq_discount[ii] <- sb.df$hq_orig[ii] - length(sb.disturb.hq[,ii+1][sb.disturb.hq[,ii+1] > sb.df$quant75[ii]])
  }
  ## Update sb.df with the proportion of remaining high-quality suitable pixels
  sb.df$hq_remaining <- sb.df$hq_discount / sb.df$hq_orig

  return(sb.df)
}



#' Brant impact analysis
#'
#' \code{impact_brant} runs the brant impact analysis, calculating the area of
#'   expected brant disturbance based on physical infrastructure and helicopter
#'   overflights, then estimating the habitat area and number of brant affected.
#'
#' @param brant.lakes SpatialPolygonsDataFrame object containing lake boundaries
#'   with molting brant density data from Schults and Zeller (2019). Created in
#'   the data preparation code of \code{\link{dia}}.
#' @param footprints List of SpatialPolygonsDataFrame objects, \code{cpf_footprints},
#'   \code{sat_footprints}, and \code{surf_disturb}, created by \code{\link{footprint_generation}}.
#' @param land.dist Numeric value indicating the distance (m) from CPF or
#'   satellite pads at which helicopters are expected to be under 500 m altitude,
#'   based on a descent angle of 10 degrees following FAA guidelines. Used to
#'   determine helicopter impact area around CPFs and satellite pads. Defaults
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
#'
#' @return \code{out.df} data.frame, updated to contain the brant impact results.
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
impact_brant <- function(brant.lakes, footprints, land.dist = 2835.64, heli.disturb = 3570, proximity.effect = 100, proj.info){
  ## Calculate the disturbance buffer around CPFs and satellites based on the assumption that
  ## helicopters will stay above 500 m, except when taking off or landing, but may approach CPFs and
  ## satellites from any direction
  cpfBuff <- raster::buffer(x=footprints$cpf_footprints, width=land.dist+heli.disturb)
  satBuff <- raster::buffer(x=footprints$sat_footprints, width=land.dist+heli.disturb)
  row.names(cpfBuff) <- paste(1:length(cpfBuff))
  row.names(satBuff) <- paste((length(cpfBuff)+1):(length(cpfBuff)+length(satBuff)))

  ## Account for proximity effects around the physical footprint of infrastructure, adding this to the
  ## helicopter disturbance polygons.
  sd.disturb.buf <- raster::buffer(x=footprints$surf_disturb, width=proximity.effect)
  sd.disturb.buf <- projection_alignment(sd.disturb.buf, proj.info)
  all.disturb.buf <- raster::aggregate(rbind(cpfBuff, satBuff, sd.disturb.buf))

  ## Identify which portions of lakes intersect the disturbance buffer.
  disturbed.lakes.all <- rgeos::gIntersection(brant.lakes, all.disturb.buf, byid=TRUE)

  ## If the infrastructure scenario does not overlap the brant lakes, enter zeros for impact,
  ## otherwise calculate the lake area and estimated number of brant disturbed.
  if(class(disturbed.lakes.all) == "NULL"){
    tmp.df <- data.frame("lakes_disturbed_all_ha" = 0, "lakes_disturbed_molt_ha" = 0, "disturbed_geese" = 0)
  } else{
    ## The LAKE column in the brant.lakes shapefile is FID + 1. Use this to relate the brant.lakes
    ## attributes to the intersected lakes to allow calculation of impact.
    tmp.df.disturb <- data.frame(LAKE=as.numeric(substr(names(disturbed.lakes.all), start=1,
                                                        stop=nchar(names(disturbed.lakes.all))-2))+1,
                                 affected_area_sq_m=raster::area(disturbed.lakes.all))
    tmp.df2.disturb <- merge(tmp.df.disturb, as.data.frame(brant.lakes))

    ## Add the area for all disturbed lakes and those in which brant have been observed molting over
    ## the last 5 years to a tmp.df object. Multiply area by bird density to estimate the number of
    ## disturbed birds. Convert areas to units of hectares, for easier use.
    tmp.df <- data.frame("lakes_disturbed_all_ha"=sum(tmp.df2.disturb$affected_area_sq_m) / 10000,
                         "lakes_disturbed_molt_ha"=sum(tmp.df2.disturb$affected_area_sq_m[tmp.df2.disturb$X5Yr_Av > 0]) / 10000,
                         "disturbed_geese"=round(sum(tmp.df2.disturb$affected_area_sq_m * tmp.df2.disturb$BLBRpM2),0))
  }

  return(tmp.df)
}
