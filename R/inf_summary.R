#' Summarize simulated infrastructure
#'
#' Infrastructure summary analysis. This calculates the number of facilities/length
#'  of roads, area of direct surface disturbance, and area within 4 km of
#'  infrastructure for comarison with statistics provided by the Bureau of Land
#'  Management (BLM) in the Integrated Activity Plan Draft Environmental Impact
#'  Statement (IAP DEIS) for the National Petroleum Reserve - Alaska (NPR-A)
#'  (BLM 2019a). These metrics are added to the \code{out.df} data.frame.
#'
#' @param out.df data.frame object containing the current scenario and iteration
#'   information. Created during the impact analysis by \code{\link{dia}}.
#' @param cpf.spdf  SpatialPointsDataFrame object indicating central processing
#'   facility (CPF) locations. Created by \code{\link{generate_cpf}}.
#' @param sat.spdf SpatialPointsDataFrame object indicating satellite production
#'   facility locations. Created by \code{\link{generate_sat_rd}}.
#' @param rd.sl SpatialLines object indicating gravel road centerlines. Created
#'   by \code{\link{generate_sat_rd}}.
#' @param surf.disturb SpatialPolygonsDataFrame object representing the combined
#'   footprint of central processing facilities (CPFs), satellite pads, and
#'   roads. Created by \code{\link{footprint_generation}}.
#' @param npra.file Character string giving the file name of the NPR-A boundary
#'   shapefile, without extension.
#' @param zoi Numeric value indicating the zone of influence (m) in which habitat
#'   quality is discounted. Defaults to 4000 m, as is used in BLM (2019a).
#' @param wd.loc Character string identifying the base directory containing both
#'   input and output data folders.
#' @param path.in Character string identifying the relative path from the base
#'   directory (\code{wd.loc}) to where input data are stored.
#'
#' @return \code{out.df} data.frame, updated to contain the infrastructure
#'   summary results.
#'
#' @references BLM 2019a. National Petroleum Reserve in Alaska Draft Integrated
#'   Activity Plan and Environmental Impact Statement. Bureau of Land Management,
#'   U.S. Department of the Interior, Anchorage, AK, USA.
inf_summary <- function(out.df, cpf.spdf, sat.spdf, rd.sl, surf.disturb, npra.file, zoi = 4000, wd.loc, path.in){
  ## Identify the total number of CPFs (including Willow), satellites, and km of roads. Note that because of
  ## the co-location of the Willow CPF and BT3, Willow is counted both as a satellite and CPF. This code
  ## corrects for that by subtracting 1 from the number of satellite pads.
  out.df$n_cpf <- length(cpf.spdf)
  out.df$n_sat <- length(sat.spdf)-1
  out.df$rd_km <- sum(sp::SpatialLinesLengths(SL=rd.sl, longlat=FALSE)/1000)

  ## Calculate area of surface disturbance based on the surface disturbance footprint created by footprint_generation().
  out.df$surface_disturb_ha <- raster::area(surf.disturb)/10000

  ## Calculate the area within the zoi (default 4 km) of development, first cropping with the NPR-A
  ## boundary polygon to ensure disturbance is only counted within the NPR-A.
  sd.buf <- raster::buffer(x=surf.disturb, width=zoi)
  npra <- rgdal::readOGR(dsn=paste(wd.loc, path.in, sep="/"), layer=npra.file)
  sd.buf.npra <- raster::intersect(x=sd.buf, y=npra)
  out.df$area_within_4km_of_dev_ha <- raster::area(sd.buf.npra)/10000

  ## The DEIS uses units of miles and acres, rather than km and hectares. For comparison, add these too.
  out.df$rd_mi <- out.df$rd_km / 1.609
  out.df$surface_disturb_acres <- out.df$surface_disturb_ha * 2.47105
  out.df$area_within_4km_of_dev_acres <- out.df$area_within_4km_of_dev_ha * 2.47105

  return(out.df)
}
