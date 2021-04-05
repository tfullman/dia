#' Summarize simulated infrastructure
#'
#'  \code{inf_summary} calculates the number of facilities/length
#'  of roads, area of direct surface disturbance, and area within 4 km of
#'  infrastructure for comparison with statistics provided by the Bureau of Land
#'  Management (BLM) in the Integrated Activity Plan Draft Environmental Impact
#'  Statement (IAP DEIS) for the National Petroleum Reserve - Alaska (NPR-A)
#'  (BLM 2019a). These metrics are added to the \code{out.df data.frame}.
#'
#' @param out.df \code{data.frame} object containing the current scenario and iteration
#'   information. Created during the impact analysis by \code{\link{dia}}.
#' @param cpf.sf  \code{sf POINT} object indicating central processing
#'   facility (CPF) locations. Created by \code{\link{generate_cpf}}.
#' @param sat.sf \code{sf POINT} object indicating satellite production
#'   facility locations. Created by \code{\link{generate_sat_rd}}.
#' @param rd.sf \code{sf LINESTRING} object indicating gravel road centerlines. Created
#'   by \code{\link{generate_sat_rd}}.
#' @param surf.disturb \code{sf POLYGON} object representing the combined
#'   footprint of CPFs, satellite pads, and roads. Created by \code{\link{footprint_generation}}.
#' @param npra.file Character string giving the file name of the NPR-A boundary
#'   shapefile, without extension.
#' @param zoi Numeric value indicating the zone of influence (m) in which habitat
#'   quality is discounted. Defaults to 4000 m, as is used in BLM (2019a).
#' @param wd.loc Character string identifying the base directory containing both
#'   input and output data folders.
#' @param path.in Character string identifying the relative path from the base
#'   directory (\code{wd.loc}) to where input data are stored.
#' @param proj.info Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
#'   common to all spatial objects in the analysis.
#'
#' @return \code{out.df data.frame}, updated to contain the infrastructure
#'   summary results.
#'
#' @references BLM 2019a. National Petroleum Reserve in Alaska Draft Integrated
#'   Activity Plan and Environmental Impact Statement. Bureau of Land Management,
#'   U.S. Department of the Interior, Anchorage, AK, USA.
inf_summary <- function(out.df, cpf.sf, sat.sf, rd.sf, surf.disturb, npra.file, zoi = 4000, wd.loc, path.in, proj.info){
  ## Identify the total number of CPFs (including Willow), satellites, and km of roads. Note that because of
  ## the co-location of the Willow CPF and BT3, Willow is counted both as a satellite and CPF. This code
  ## corrects for that by subtracting 1 from the number of satellite pads.
  out.df$n_cpf <- nrow(cpf.sf)
  out.df$n_sat <- nrow(sat.sf)-1
  out.df$rd_km <- units::set_units(sum(sf::st_length(rd.sf)), "km")

  ## Calculate area of surface disturbance based on the surface disturbance footprint created by
  ## footprint_generation(), returning the result in units of hectares.
  out.df$surface_disturb_ha <- units::set_units(sf::st_area(surf.disturb), "ha")

  ## Calculate the area within the zoi (default 4 km) of development, first cropping with the NPR-A
  ## boundary polygon to ensure disturbance is only counted within the NPR-A.
  sd.buf <- sf::st_buffer(x=surf.disturb, dist=zoi)
  npra <- load_spatial(x=npra.file, vec=TRUE, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)
  sd.buf.npra <- sf::st_intersection(x=sd.buf, y=npra)
  out.df$area_within_4km_of_dev_ha <- units::set_units(sf::st_area(sd.buf.npra), "ha")

  ## The DEIS uses units of miles and acres, rather than km and hectares. For comparison, add these too.
  out.df$rd_mi <- udunits2::ud.convert(out.df$rd_km, "km", "miles")
  out.df$surface_disturb_acres <- udunits2::ud.convert(out.df$surface_disturb_ha, "ha", "acres")
  out.df$area_within_4km_of_dev_acres <- udunits2::ud.convert(out.df$area_within_4km_of_dev_ha, "ha", "acres")

  return(out.df)
}
