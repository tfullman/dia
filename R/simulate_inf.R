#' Create infrastructure exclusion zones
#'
#' \code{infrastructure_exclusion_buffer} is a helper function for infrastructure
#'   generation that sets raster values to zero within a certain threshold of a
#'   given feature (like a central processing facility) to ensure desired spacing
#'   of infrastructure is maintained. Used by \code{\link{prep_general_inputs}}
#'   and \code{\link{generate_sat_rd}}.
#'
#' @param x Spatial object identifying the location(s) around which infrastructure
#'   is to be excluded (i.e., around which raster values will be set to zero).
#' @param y \code{SpatRaster} object containing values that will be set to zero.
#' @param threshold Numeric value providing the distance (m) within which raster
#'   values will be set to zero.
#'
#' @return \code{SpatRaster}, \code{y}, updated to have values within \code{threshold} of
#'   infrastructure locations, \code{x}, set to zero.
#'
infrastructure_exclusion_buffer <- function(x, y, threshold){
  ## Convert the infrastructure location(s) to polygons by buffering by the desired threshold.
  x.buf <- sf::st_buffer(x, threshold)
  ## Set raster values to zero within threshold of infrastructure
  y2 <- raster_to_zero(y, terra::vect(x.buf))
  return(y2)
}


#' Simulate locations proportional to suitability and outside unavailable areas
#'
#' \code{infrastructure_spacer} is a helper function for infrastructure simulation
#'   that generates infrastructure locations proportional to a suitability
#'   raster (e.g., relative estimated undiscovered oil) and constrained to not
#'   occur in unavailable areas or within a minimum distance threshold from other
#'   infrastructure. Used by \code{\link{generate_cpf}} and \code{\link{generate_sat_rd}}.
#'
#' @param n.dev Integer value indicating the desired number of simulated
#'   developments (e.g., central processing facilities or satellite pads).
#' @param im Object of class \code{im} that is used by spatstat to randomly
#'   generate locations proportional to suitability and constrained to avoid
#'   restrictions like water or areas closed to development.
#' @param threshold Numeric value indicating the minimum spacing (m) between
#'   simulated infrastructure locations.
#'
#' @return \code{data.frame} containing the xy coordinates of the simulated development.
#'
infrastructure_spacer <- function(n.dev, im, threshold){
  ## Generate infrastructure locations proportional to availability and avoiding constraints
  dev.pts <- spatstat::rpoint(n=n.dev, f=im)
  dev.coords <- as.data.frame(dev.pts)

  ## Ensure a minimum distance between development locations
  if(n.dev > 1){
    dev.dists <- stats::dist(cbind(dev.pts$x, dev.pts$y))
    flag <- ifelse(min(dev.dists) < threshold, 1, 0)
    while(flag == 1){
      dist.mat <- as.matrix(dev.dists)
      dist.mat[upper.tri(dist.mat, diag=TRUE)] <- NA
      dist.df <- data.frame(p1=as.numeric(rownames(dist.mat)), p2=rep(1:ncol(dist.mat), each=ncol(dist.mat)), dist=c(dist.mat))
      dist.df <- stats::na.omit(dist.df)
      cand <- dist.df[dist.df$dist < threshold,]
      to.remove <- unique(cand$p1)
      dev.coords <- dev.coords[-to.remove,]
      dev.pts2 <- spatstat::rpoint(n=n.dev-nrow(dev.coords), f=im)
      dev.coords <- rbind(dev.coords, data.frame(x=dev.pts2$x, y=dev.pts2$y))
      dev.dists <- stats::dist(cbind(dev.coords$x, dev.coords$y))
      flag <- ifelse(min(dev.dists) < threshold, 1, 0)
    }
  }
  return(dev.coords)
}



#' Prepare input data that will be used in all scenarios
#'
#' \code{prep_general_inputs} loads and prepares input data that is common to all
#'   scenarios to be evaluated, like existing infrastructure layers and relative
#'   estimated undiscovered oil.
#'
#' @param rd.exist.file Character string specifying the file name for the
#'   shapefile depicting existing roads as lines.
#' @param pad.exist.file Character string specifying the file name for the
#'   shapefile depicting existing central processing facility (CPF) and satellite
#'   pads as points.
#' @param oil.av.file Character string specifying the file name for the relative
#'   estimated undiscovered oil raster.
#' @param d2cpf Numeric value indicating the minimum allowed distance (m) between
#'   two central processing facilities (CPFs). Defaults to 35,000 m.
#' @param wd.loc Character string identifying the base directory containing both
#'   input and output data folders.
#' @param path.in Character string identifying the relative path from the base
#'   directory (\code{wd.loc}) to where input data are stored. Defaults to an
#'   Input_Data folder within the base folder.
#' @param proj.info Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
#'   common to all spatial objects in the analysis.
#'
#' @return List containing input data used for infrastructure simulation across
#'   all scenarios. This includes the \code{SpatRaster} indicating relative estimated
#'   undiscovered oil (\code{oil_av}), \code{SpatRaster} depicting the infrastructure
#'   exclusion zone around the Willow CPF (\code{willow_cpf_ras}), \code{sf} \code{POINT}
#'   object depicting existing CPF and satellite gravel pads (\code{pad_exist}),
#'   \code{sf} \code{LINESTRING} object depicting existing roads (\code{rd_exist}),
#'   \code{data.frame} with xy coordinates of existing roads (\code{exist_coords}), and
#'   \code{sf} \code{POINT} object depicting the Willow CPF location (\code{willow_spdf}).
#'   Many of these feed into \code{\link{prep_scenario_inputs}}.
#'
prep_general_inputs <- function(rd.exist.file, pad.exist.file, oil.av.file, d2cpf = 35000,
                                wd.loc, path.in, proj.info){

  ## Check that needed files are present
  if(is.null(rd.exist.file) || is.null(pad.exist.file) || is.null(oil.av.file)) stop("Missing data: One or more of rd.exist.file, pad.exist.file, or oil.av.file is missing. These each cannot be NULL if simulate.inf == TRUE.")


  #-------------------------
  ## Existing infrastructure
  #-------------------------

  ## Load the shapefiles representing existing infrastructure. This will include both a point file of CPF
  ## and satellite pad locations, as well as a line file of road locations.
  rd.exist <- load_spatial(x=rd.exist.file, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in, vec=TRUE)
  pad.exist <- load_spatial(x=pad.exist.file, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in, vec=TRUE)

  ## Create an sf object just for Willow's CPF, matching the attribute formatting used below
  willow.cpf <- pad.exist[4,]
  willow.cpf$cpf <- 0
  willow.cpf$x <- sf::st_coordinates(willow.cpf)[,1]
  willow.cpf$y <- sf::st_coordinates(willow.cpf)[,2]

  ## Put the existing infrastructure data attribute table into the format used below
  pad.exist$cpf <- 0
  pad.exist$sat <- 1:nrow(pad.exist)
  pad.exist$x <- sf::st_coordinates(pad.exist)[,1]
  pad.exist$y <- sf::st_coordinates(pad.exist)[,2]

  ## Isolate the coordinates of the existing roads and rearrange to match the format used below
  exist.coords.df <- as.data.frame(sf::st_coordinates(rd.exist))
  exist.coords.df <- exist.coords.df[,c(3:4,1:2)]
  names(exist.coords.df) <- c("cpf", "sat", "x", "y")


  #----------------------------
  ## Estimated undiscovered oil
  #----------------------------

  ## Load the raster indicating the relative probability of of estimated undiscovered oil and make sure it aligns with the desired projection
  oil.av <- load_spatial(x=oil.av.file, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)

  ## Update the estimated undiscovered oil raster so that values around existing infrastructure are zero, preventing
  ## new development from being built too close to existing infrastructure (i.e., the Willow CPF). Create a
  ## raster that identifies the pixels around the existing CPF and gives them a value of 0 and gives all the
  ## other pixels a value of 1.
  willow.cpf.ras <- oil.av
  terra::values(willow.cpf.ras)[!is.na(terra::values(willow.cpf.ras))] <- 1
  willow.cpf.ras <- infrastructure_exclusion_buffer(x=willow.cpf, y=willow.cpf.ras, threshold=d2cpf)


  #----------------------------------------------
  ## Return the outputs needed in other functions
  #----------------------------------------------

  return(list("oil_av"=oil.av, "willow_cpf_ras"=willow.cpf.ras, "pad_exist"=pad.exist, "rd_exist"=rd.exist,
              "exist_coords"=exist.coords.df, "willow_cpf"=willow.cpf))
}



#' Prepare input data needed for an individual scenario
#'
#' Some input data needed for infrastructure simulation is common to all scenarios
#' while other input data is scenario-specific. The former is prepared by
#' \code{\link{prep_general_inputs}}, while \code{prep_scenario_inputs} prepares
#' the scenario-specific data. It thus must be run separately for each scenario.
#'
#' @param scenario Character vector of strings identifying the scenarios to be run.
#' @param oil.av \code{SpatRaster} indicating the relative estimated undiscovered oil,
#'   created by \code{\link{prep_general_inputs}}.
#' @param cost.map.file Character string specifying the file name for the cost map
#'   raster, which has all water areas coded as 0 and all land as 1.
#' @param willow.cpf.ras \code{SpatRaster} that restricts central processing facility
#'   (CPF) generation within \code{d2cpf} of the Willow CPF. Created by
#'   \code{\link{prep_general_inputs}}.
#' @param pad.res.file Character string specifying the file name for the raster
#'   containing scenario-specific restrictions for gravel pad (i.e., CPF and
#'   satellite pad) placement. Restriction values may consist of 0 (no development),
#'   0.1 (reduced likelihood of development), or 1 (unrestricted development).
#' @param road.res.file Character string indicating the file name for the raster
#'   containing scenario-specific restrictions for road placement. Restriction values
#'   may consist of 0 (no development), 0.1 (reduced likelihood of development), or
#'   1 (unrestricted development).
#' @param alt.b.rd.stranded.res.file Character string specifying the file name for
#'   the road development restriction raster for stranded leases under Alternative
#'   B. Optional, allowing the code to be run for scenarios other than Alternative B.
#' @param alt.b.stranded.lease.file Character string specifying the file name for the
#'   raster identifying stranded lease locations under Alternative B. Optional,
#'   allowing the code to be run for scenarios other than Alternative B.
#' @param alt.c.row.file Character string specifying the file name for the raster
#'   identifying right-of-way (ROW) areas under Alternative C where infield roads
#'   are prohibited but connecting infrastructure is allowed. Optional, allowing
#'   the code to be run for scenarios other than Alternative C.
#' @param alt.d.north.file Character string specifying the file name for the
#'   raster identifying discrete areas north of Teshekpuk Lake where roadless
#'   development is possible under Alternative D. Optional, allowing the code to
#'   be run for scenarios other than Alternative D.
#' @inheritParams prep_general_inputs
#' @param z Iterator value used by \code{\link{dia}} to specify the development
#'   scenario being analyzed.
#'
#' @return List containing scenario-specific data for: \code{SpatRaster} indicating
#'   the relative likelihood of pad development as a function of relative
#'   estimated undiscovered oil, development restrictions, and the presence of
#'   waterbodies (\code{pad_prob}), \code{im} version of \code{pad_prob} (\code{pad_im}),
#'   \code{SpatRaster} indicating resistance for least cost path road generation
#'   (\code{road_res}), \code{TransitionLayer} indicating road development cost
#'   (\code{tr_cost_alt}), \code{TransitionLayer} specific to the stranded leases in
#'   Alternative B (\code{tr_cost_B_stranded}), \code{SpatRaster} indicating location
#'   of stranded leases in Alternative B (\code{alt_B_stranded_leases}),
#'   \code{TransitionLayer} specific to the stranded leases in the ROW area south of
#'   Teshekpuk Lake for Alternative C (\code{tr_cost_C_row}), \code{SpatRaster}
#'   identifying discrete areas north of Teshekpuk Lake where roadless
#'   development is possible under Alternative D (\code{alt_D_north}).
#'
#' @references Wilson RR, Liebezeit JR, Loya WM. 2013. Accounting for uncertainty in oil and
#'   gas development impacts to wildlife in Alaska. Conservation Letters
#'   6:350-358.
prep_scenario_inputs <- function(scenario, oil.av, cost.map.file, willow.cpf.ras, pad.res.file,
                                 road.res.file, alt.b.rd.stranded.res.file = NULL, alt.b.stranded.lease.file = NULL,
                                 alt.c.row.file = NULL, alt.d.north.file = NULL, wd.loc, path.in, proj.info, z){

  #-------------------
  ## Cost map creation
  #-------------------

  ## Load the cost.map raster, indicating water and land areas and make sure it aligns with the desired projection
  cost.map <- load_spatial(x=cost.map.file, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)

  ## Create the road-based cost map, which changes water cost to 0.05 per Wilson et al. 2013 due to possibility of building bridges over rivers and waterbodies.
  cost.water <- cost.map
  cost.water[cost.water == 0] <- 0.05


  #-----------------------------------------------------------------------
  ## Create scenario-specific relative likelihood of pad generation raster
  #-----------------------------------------------------------------------

  ## The relative likelihood of CPF and satellite pad generation is based on a combination of estimated
  ## undiscovered oil, scenario-specific development restrictions, and presence of waterbodies. Create this
  ## layer here in the various forms needed for pad simulation.

  ## Load the raster defining areas of restricted development for pads under the current scenario and ensure
  ## it matches the desired projection. Areas where infrastructure is banned have a value of 0, while areas
  ## where these protections generally exist but could be relaxed have a value of 0.1, following Wilson et
  ## al. 2013. Areas without restrictions on infrastructure have a value of 1.
  pad.res <- load_spatial(x=pad.res.file, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)

  ## The estimated undiscovered oil raster is used to simulate the locations of CPFs and satellite pads.
  ## Update the base raster so that water (i.e., anything with a value of zero in the base cost map) gets
  ## a value of zero in the likelihood raster, ensuring that CPFs and satellite pads are not built on
  ## water. Also incorporate the scenario-specific restrictions on development so that areas closed to
  ## development or with a lower probability of development can be incorporated. Multiplying restrictions
  ## and water by estimated undiscovered oil will leave oil values unchanged in unrestricted land areas
  ## (x*1*1), set to zero in water or prohibited areas (x*0*1 or x*1*0), and lowered in areas where
  ## infrastructure is less likely but still possible (x*1*0.1).
  pad.prob <- oil.av * cost.map * pad.res

  ## Ensure that CPFs are not built within d2cpf of the existing Willow CPF by multiplying the
  ## scenario-specific estimated undiscovered oil and restriction raster by the Willow CPF mask.
  pad.prob.inf <- pad.prob * willow.cpf.ras

  ## Convert this to an im.object, as required by spatstat
  pad.im <- maptools::as.im.RasterLayer(pad.prob.inf)


  #--------------------------------------------------------------------------------------------
  ## Create road restriction transition layers to enable road generation using least cost paths
  #--------------------------------------------------------------------------------------------

  ## Load the raster defining areas of restricted development for roads under the current scenario and
  ## ensure it matches the desired projection. Areas where infrastructure is banned have a value of 0, while
  ## areas where these protections generally exist but could be relaxed (e.g., along most river buffers)
  ## have a value of 0.1, following Wilson et al. 2013. Areas without restrictions on infrastructure have a
  ## value of 1.
  road.res <- load_spatial(x=road.res.file, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)

  ## Update this to reduce the probability of development in areas overlapping water.
  road.cost <- road.res * cost.water

  ## Convert the road cost raster to a transition layer that can be used for creating roads using least
  ## cost paths. This has to be in RasterLayer, not SpatRaster, format for the transition() function to work.
  road.cost.ras <- raster::raster(road.cost)
  tr.cost.alt <- gdistance::transition(road.cost.ras, transitionFunction=mean, directions=8)


  #-----------------------------------------------------
  ## Prepare additional data needed in certain scenarios
  #-----------------------------------------------------

  ## Alternatives B-D have extra inputs needed to deal with scenario-specific constraints on development.
  ## Prepare them here.

  ## Alternative B
  if(substr(scenario[z], start=1, stop=5) == "Alt_B"){
    ## Development in the stranded leases west of Teshekpuk Lake under Alternative B will require a
    ## different restriction raster as well as identification of stranded leases (i.e., those completely
    ## surrounded by areas unavailable for leasing and development. Load these and ensure they match the
    ## desired projection.
    alt.b.road.stranded.res <- load_spatial(x=alt.b.rd.stranded.res.file, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)
    alt.b.stranded.leases <- load_spatial(x=alt.b.stranded.lease.file, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)

    ## Calculate road cost for stranded leases, accounting for water
    alt.b.stranded.road.cost <- alt.b.road.stranded.res * cost.water

    ## Create a transition layer for stranded leases
    alt.b.stranded.road.cost.ras <- raster::raster(alt.b.stranded.road.cost)
    tr.cost.alt.b.stranded <- gdistance::transition(alt.b.stranded.road.cost.ras, transitionFunction=mean, directions=8)
  }

  ## Alternative C
  if(substr(scenario[z], start=1, stop=5) == "Alt_C"){
    ## Under Alternative C, infield roads are prohibited but connecting infrastructure in a right of way
    ## (ROW) would be allowed within some of the 50% TCH calving kernel. This means areas that get a value
    ## of 0 for pads and infield roads would have a value of 0.1 for connecting roads. Load the raster
    ## that accounts for this and for existing leases and ensure it has the desired projection.
    alt.c.row.res <- load_spatial(x=alt.c.row.file, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)

    ## Calculate road cost for the ROW areas, accounting for water
    alt.c.row.cost <- alt.c.row.res * cost.water

    ## Create a transition layer for ROW areas
    alt.c.row.cost.ras <- raster::raster(alt.c.row.cost)
    tr.cost.alt.c.row <- gdistance::transition(alt.c.row.cost.ras, transitionFunction=mean, directions=8)
  }

  ## Alternative D
  if(substr(scenario[z], start=1, stop=5) == "Alt_D"){
    ## Alternative D is the only alternative to open areas north of Teshekpuk Lake for leasing and
    ## development. Restrictions to the north and on each side of the lake necessitate that this would be
    ## roadless development, where infield roads could occur linking satellites to a CPF, but not a
    ## connecting road to link the field to other existing infrastructure. There are also some places
    ## where satellites would also have to be roadless. Load the raster that indicates unique values for
    ## each available area to allow for roadless development north of Teshekpuk Lake under this
    ## alternative and ensure it matches the desired projection.
    alt.d.TLnorth <- load_spatial(x=alt.d.north.file, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)

    ## Change values outside of the region north of Teshekpuk Lake set from NA to zero.
    terra::values(alt.d.TLnorth)[is.na(terra::values(alt.d.TLnorth))] <- 0
    alt.d.TLnorth <- terra::mask(alt.d.TLnorth, mask=road.res)
  }


  #------------------------------
  ## Return the needed input data
  #------------------------------

  list.tmp <- list("pad_prob"=pad.prob, "pad_im"=pad.im, "road_res"=road.res, "tr_cost_alt"=tr.cost.alt,
                   "tr_cost_B_stranded"=NULL, "alt_B_stranded_leases"=NULL, "tr_cost_C_row"=NULL, "alt_D_north"=NULL)
  if(substr(scenario[z], start=1, stop=5) == "Alt_B"){
    list.tmp$tr_cost_B_stranded <- tr.cost.alt.b.stranded
    list.tmp$alt_B_stranded_leases <- alt.b.stranded.leases
  }
  if(substr(scenario[z], start=1, stop=5) == "Alt_C") list.tmp$tr_cost_C_row <- tr.cost.alt.c.row
  if(substr(scenario[z], start=1, stop=5) == "Alt_D") list.tmp$alt_D_north <- alt.d.TLnorth

  return(list.tmp)
}



#' Simulate CPF locations
#'
#' \code{generate_cpf} simulates the location of central processing facilities
#' (CPFs) as a function of relative estimated undiscovered oil, development
#' restrictions, and the presence of waterbodies. CPFs are "the operational
#' center for long-term production" of an oilfield (BLM 2019a, p. B-6) so
#' infrastructure simulation starts with CPFs. CPF locations approximate the
#' locations of newly discovered oil accumulations suitable for development
#' (Wilson et al. 2013).
#'
#' @param pad.im Object of type "\code{im}" representing the relative likelihood
#'   of pad development as a function of relative undiscovered oil, development
#'   restrictions, and the presence of waterbodies. Created by \code{\link{prep_scenario_inputs}}.
#' @param n.cpf Indicator of the desired number of CPFs to be generated. This
#'   can either be a single integer value, indicating a fixed number of CPFs, or
#'   a range \code{c(min,max)}, from which a random value will be drawn using a
#'   uniform distribution. Defaults to \code{c(3,7)} to reflect a reasonable
#'   range of potential variability, depending on the number and distribution of
#'   undiscovered oil deposits.
#' @param debug.out Logical indicator of whether intermediate infrastructure .csv
#'   files should be written out for the purpose of debugging code. Defaults to
#'   \code{FALSE}.
#' @param path.out Character string identifying the relative path from the base
#'   directory (\code{wd.loc}) to where output data will be saved. Defaults to
#'   an "Output_Data" folder within the base folder.
#' @param n.iter Integer indicating the desired number of iterations to be run
#'   for each scenario. Defaults to 100.
#' @param scenario Vector of character strings identifying the scenarios being run.
#' @param z Iterator value used by \code{\link{dia}} to specify the development
#'   scenario being analyzed.
#' @param i Iterator value used by \code{\link{dia}} to run development simulation
#'   and impacts analyses in parallel, indicating the specific iteration being run.
#' @inheritParams prep_general_inputs
#'
#'
#' @return List containing a \code{data.frame} of simulated CPF locations (\code{cpf_df}),
#'   \code{sf} \code{POINT} object of CPF locations (\code{cpf_sf}), and integer indicating
#'   the number of simulated CPFs in the current iteration (\code{tmp_n_cpf}).
#'
#' @references BLM \[Bureau of Land Management\] 2019a. National Petroleum Reserve
#'   in Alaska Draft Integrated Activity Plan and Environmental Impact Statement.
#'   Bureau of Land Management, U.S. Department of the Interior, Anchorage, AK,
#'   USA.
#'
#'  Wilson RR, Liebezeit JR, Loya WM. 2013. Accounting for uncertainty in oil and
#'   gas development impacts to wildlife in Alaska. Conservation Letters
#'   6:350-358.
generate_cpf <- function(pad.im, n.cpf = c(3, 7), d2cpf = 35000, debug.out = FALSE, wd.loc, path.out,
                         n.iter=100, proj.info = NULL, scenario = NULL, z = NULL, i = NULL){
  ## Check whether n.cpf is random or fixed. If random, draw a uniform random value.
  tmp.n.cpf <- ifelse(length(n.cpf)==1, n.cpf, sample(n.cpf[1]:n.cpf[2], size=1))

  ## Generate CPF locations proportional to estimated undiscovered oil and constrained to not occur in water
  ## or restricted areas, ensuring a minimum distance (d2cpf) between CPFs.
  cpf.coords <- infrastructure_spacer(n.dev = tmp.n.cpf, im = pad.im, threshold = d2cpf)

  ## Prepare the output data
  cpf.df <- data.frame(cpf=1:tmp.n.cpf, x=cpf.coords$x, y=cpf.coords$y)
  cpf.sf <- sf::st_as_sf(x=cpf.df, coords=c("x", "y"), crs=proj.info)

  ## For debugging purposes, this option will write out intermediate data
  if(debug.out) utils::write.csv(cpf.df, file=paste(wd.loc, "/", path.out, "/debug_CPF_data_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), ".csv", sep=""), row.names=FALSE)

  ## Return the desired outputs
  return(list("cpf_df"=cpf.df, "cpf_sf"=cpf.sf, "tmp_n_cpf"=tmp.n.cpf))
}


#' Road generation using least cost paths
#'
#' \code{gen_lcp_rd} is a helper function that generates a road using least cost
#' paths and outputs the results as both a \code{data.frame} and \code{sf LINESTRING}
#' object. Used by \code{\link{lcp_rds_infield}}, \code{\link{lcp_rds_outfield}} via
#' \code{\link{gen_linkage_rd}}, and \code{\link{generate_sat_rd}}.
#'
#' @param tr Scenario-specific \code{TransitionLayer} object, used by the \code{gdistance}
#'   package for creating least cost paths.
#' @param loc1 Numeric vector of length 2 indicating the x and y coordinates of the first location to be joined with a road.
#' @param loc2 Numeric vector of length 2 indicating the x and y coordinates of the second location to be joined with a road.
#' @param cur.cpf Integer indicating the identity of the central processing facility (CPF) affiliated with the road to be generated.
#' @param cur.sat Integer indicating the identity of the satellite production pad affiliated with the road to be generated.
#' @param proj.info Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
#'   common to all spatial objects in the analysis.
#'
#' @return List of length 2 containing a \code{data.frame} with newly generated road
#'   coordinates (\code{df_tmp}) and \code{sf LINESTRING} object depicting road locations
#'   (\code{sf_tmp}).
#'
gen_lcp_rd <- function(tr, loc1, loc2, cur.cpf, cur.sat, proj.info){
  if(length(loc1) != 2 || class(loc1) != "numeric") stop("loc1 must be a numeric vector of length 2")
  if(length(loc2) != 2 || class(loc2) != "numeric") stop("loc2 must be a numeric vector of length 2")
  ## Calculate least cost path
  rd.tmp <- gdistance::shortestPath(x=tr, origin=loc1, goal=loc2, output="SpatialLines")
  ## Format the results
  rd.tmp.sf <- projection_alignment(sf::st_as_sf(rd.tmp), proj.info)
  rd.tmp.df <- as.data.frame(cbind(cpf=cur.cpf, sat=cur.sat,
                                   x=sf::st_coordinates(rd.tmp.sf)[,1], y=sf::st_coordinates(rd.tmp.sf)[,2]))
  ## Return the desired objects
  return(list("df_tmp"=rd.tmp.df, "sf_tmp"=rd.tmp.sf))
}


#' Connect satellite pads to their CPF via infield roads
#'
#' \code{lcp_rds_infield} is a helper function that uses least cost paths to
#'   connect satellite pads to their parent central processing facility (CPF) via
#'   infield roads. Used by \code{\link{generate_sat_rd}}.
#'
#' @param x \code{data.frame} of satellite pad locations, created within \code{\link{generate_sat_rd}}.
#'   This may contain all satellite pad locations, or a subset, depending on the
#'   scenario being analyzed.
#' @param tr Scenario-specific \code{TransitionLayer} object, used by the \code{gdistance}
#'   package for creating least cost paths.
#' @param road.res \code{SpatRaster} indicating resistance for least cost path road
#'   generation. Created by \code{\link{prep_scenario_inputs}}.
#' @param cpf.df \code{data.frame} of CPF point locations. Created by \code{\link{generate_cpf}}.
#' @param min.dist.exist Proximity order-sorted \code{data.frame} of distance from CPFs
#'   to existing infrastructure. Created within \code{\link{generate_sat_rd}}.
#' @param j Iterator used within \code{\link{generate_sat_rd}} to indicate the
#'   current CPF being analyzed.
#' @inheritParams generate_cpf
#'
#' @return List containing two objects, a temporary \code{data.frame} of generated road
#'   coordinates (\code{tmp_rd_df}) and a temporary \code{sf LINESTRING} object
#'   depicting generated roads (\code{tmp_rd_sf}). These are fed into \code{\link{lcp_rds_outfield}}.
#'
lcp_rds_infield <- function(x, tr, road.res, cpf.df, min.dist.exist, wd.loc, path.out, scenario, n.iter, j,
                            z = NULL, i = NULL, proj.info, debug.out = FALSE){
  ## The closest satellite pad will be connected directly to the CPF using a least cost path, then the rest
  ## will be connected to the closest infrastructure in order of proximity. Determine proximity order.
  ## Iterating through the rest in order of proximity, determine if each is closer to the CPF or an
  ## already connected satellite pad or road and connect to that via least cost paths.
  x$d2cpf <- sqrt((x$x - cpf.df$x[min.dist.exist$cpf[j]])^2 + (x$y - cpf.df$y[min.dist.exist$cpf[j]])^2)
  x <- x[order(x$d2cpf),]

  ## Generate the roads for each satellite pad
  for(m in 1:nrow(x)){
    if(m == 1){
      ## Connect the closest satellite pad to the CPF using least cost paths
      tmp.rd <- gen_lcp_rd(tr=tr, loc1=c(x$x[m], x$y[m]),
                           loc2=c(cpf.df$x[min.dist.exist$cpf[j]], cpf.df$y[min.dist.exist$cpf[j]]),
                           cur.cpf=min.dist.exist$cpf[j], cur.sat=x$sat[m], proj.info=proj.info)
      tmp.rd.sf <- tmp.rd$sf_tmp
      tmp.rd.df <- tmp.rd$df_tmp
    } else{
      ## Continue in proximity order for each remaining satellite pad, determining whether it is closest
      ## to the CPF, a previously connected pad, or a connecting road.
      cand.coords <- rbind(data.frame(dev="cpf", num=min.dist.exist$cpf[j],
                                      x=cpf.df$x[min.dist.exist$cpf[j]], y=cpf.df$y[min.dist.exist$cpf[j]]),
                           data.frame(dev="sat", num=x$sat[1:m-1], x=x$x[1:m-1], y=x$y[1:m-1]),
                           data.frame(dev="rd", num=1:nrow(tmp.rd.df), x=tmp.rd.df$x, y=tmp.rd.df$y))
      cand.coords$dist_sat <- t(terra::distance(x=cbind(x$x[m], x$y[m]), y=cbind(cand.coords$x, cand.coords$y), lonlat=FALSE))
      min.cand <- cand.coords[which.min(cand.coords$dist_sat),]

      ## Connect to the closest feature using least cost paths, as long as the satellite is at least
      ## one pixel away from the feature to which it is being connected. Otherwise, the pad is
      ## already connected and no road is needed.
      if(min.cand$dist_sat > min(terra::res(road.res))){
        tmp.rd2 <- gen_lcp_rd(tr=tr, loc1=c(x$x[m], x$y[m]), loc2=c(min.cand$x, min.cand$y),
                              cur.cpf=min.dist.exist$cpf[j], cur.sat=x$sat[m], proj.info=proj.info)
        tmp.rd.sf <- rbind(tmp.rd.sf, tmp.rd2$sf_tmp)
        tmp.rd.df <- rbind(tmp.rd.df, tmp.rd2$df_tmp)
      }
    }
  }

  ## For debugging purposes, this option will write out intermediate data
  if(debug.out) utils::write.csv(tmp.rd.df, file=paste(wd.loc, "/", path.out, "/debug_road_infield_data_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_cpf_", j, "_", Sys.Date(), ".csv", sep=""), row.names=FALSE)

  ## Return the desired products
  return(list("tmp_rd_df"=tmp.rd.df, "tmp_rd_sf"=tmp.rd.sf))
}


#' Generate linkage road between anchor fields
#'
#' \code{gen_linkage_rd} is a helper function used by \code{\link{lcp_rds_outfield}}
#' to identify the closest point between a new anchor field (a central processing
#' facility and associated satellite production pads and connecting roads) and
#' previously generated/existing infrastructure and then link the two with a road
#' using least cost paths.
#'
#' @param x \code{data.frame} of locations to be connected to existing infrastructure.
#'   For most scenarios this is the full anchor field dataset, \code{tmp.rd.df},
#'   created by \code{\link{lcp_rds_infield}}, but \code{if(roadless.check == 2)}
#'   it might be a subset of this \code{data.frame}, designated as \code{rd.0s}. This
#'   is the same as the object \code{x} input in \code{\link{lcp_rds_outfield}}.
#' @param y \code{data.frame} of all anchor field road locations. This defaults to be
#'   the same as \code{x}, which is valid for most scenarios. However,
#'   \code{if(roadless.check == 2)} it can be specified as the full anchor field
#'   dataset, \code{tmp.rd.df}, created by \code{\link{lcp_rds_infield}}, allowing
#'   this to be passed along to the final dataset, even for roadless portions north
#'   of Teshekpuk Lake that are not connected to existing infrastructure. This
#'   is the same as the object \code{y} input in \code{\link{lcp_rds_outfield}}.
#' @param w \code{data.frame} of coordinates of existing/previously generated
#'   infrastructure to which \code{x} will be connected.
#' @param v Road \code{sf LINESTRING} object to which the newly generated road will
#'   be added. Set to \code{rd.exist} if \code{j == 1} and \code{rd.sf} otherwise.
#' @param road.res \code{SpatRaster} indicating resistance for least cost path road
#'   generation. Created by \code{\link{prep_scenario_inputs}}.
#' @param tmp.rd.sf \code{sf LINESTRING} object depicting roads generated by
#'   \code{\link{lcp_rds_infield}}.
#' @param tr Scenario-specific \code{TransitionLayer} object, used by the \code{gdistance}
#'   package for creating least cost paths.
#' @param min.dist.exist Proximity order-sorted \code{data.frame} of distance from CPFs
#'   to existing infrastructure. Created within \code{\link{generate_sat_rd}}.
#' @inheritParams lcp_rds_outfield
#'
#' @return List containing two objects, a \code{data.frame} of newly generated linkage
#'   road coordinates added to previously generated road coordinates (\code{rd_df})
#'   and a \code{sf LINESTRING} object depicting all these roads (\code{rd_sf}).
#'   These are fed into \code{\link{lcp_rds_outfield}}.
#'
gen_linkage_rd <- function(x, y, w, v, road.res, tmp.rd.sf, tr, min.dist.exist, roadless.check,
                           alt.b.stranded.leases = NULL, alt.d.TLnorth = NULL, scenario, proj.info, z, j){
  ## Determine which part of the newly generated infield roads is closest to some part of the existing
  ## roads, first filtering out any parts of existing roads that overlap areas unavailable for
  ## development (this only applies in Alt B, where some of the Willow road overlaps the area
  ## unavailable for new infrastructure).
  w.tmp <- w
  w.tmp$unav_check <- terra::extract(road.res, cbind(w.tmp$x, w.tmp$y))
  if(substr(scenario[z], start=1, stop=5) == "Alt_B"){
    w.tmp$ntl_check <- terra::extract(alt.b.stranded.leases, cbind(w.tmp$x, w.tmp$y))
  } else if(substr(scenario[z], start=1, stop=5) == "Alt_D"){
    w.tmp$ntl_check <- terra::extract(alt.d.TLnorth, cbind(w.tmp$x, w.tmp$y))
  } else{
    w.tmp$ntl_check <- 0
  }
  if(roadless.check == 1) w.tmp$ntl_check <- 0
  w.tmp <- w.tmp[w.tmp$unav_check != 0 & w.tmp$ntl_check == 0,]
  dist.check <- raster::pointDistance(p1=cbind(x$x, x$y), p2=cbind(w.tmp$x, w.tmp$y), lonlat=FALSE)
  min.elem <- arrayInd(which.min(dist.check), dim(dist.check))

  ## Unless the new roads already overlap existing roads, connect the two using least cost paths
  if(dist.check[min.elem] <= min(terra::res(road.res))){
    ## Add the resulting files to the larger outputs
    rd.df <- rbind(w, y)
    rd.sf <- rbind(v, tmp.rd.sf)
  } else{
    ## Calculate least cost path
    tmp.rd <- gen_lcp_rd(tr=tr, loc1=c(x$x[min.elem[1]], x$y[min.elem[1]]),
                         loc2=c(w.tmp$x[min.elem[2]], w.tmp$y[min.elem[2]]), cur.cpf=min.dist.exist$cpf[j],
                         cur.sat=0, proj.info=proj.info)
    ## Add the resulting files to the larger outputs
    rd.df <- rbind(w, y, tmp.rd$df_tmp)
    rd.sf <- rbind(v, tmp.rd.sf, tmp.rd$sf_tmp)
  }

  ## Return the desired objects
  return(list("rd_df"=rd.df, "rd_sf"=rd.sf))
}



#' Connect anchor field to existing infrastructure
#'
#' \code{lcp_rds_outfield} is a helper function that operates within \code{\link{generate_sat_rd}}
#' to connect an anchor field (a central processing facility, or CPF, and its
#' associated satellites, connected by infield roads created by \code{\link{lcp_rds_infield}})
#' to existing or previously simulated infrastructure.
#'
#' @param x \code{data.frame} of locations to be connected to existing infrastructure.
#'   For most scenarios this is the full anchor field dataset, \code{tmp.rd.df},
#'   created by \code{\link{lcp_rds_infield}}, but \code{if(roadless.check == 2)}
#'   it might be a subset of this \code{data.frame}, designated as \code{rd.0s}.
#' @param y \code{data.frame} of all anchor field road locations. This defaults to be
#'   the same as \code{x}, which is valid for most scenarios. However,
#'   \code{if(roadless.check == 2)} it can be specified as the full anchor field
#'   dataset, \code{tmp.rd.df}, created by \code{\link{lcp_rds_infield}}, allowing
#'   this to be passed along to the final dataset, even for roadless portions north
#'   of Teshekpuk Lake that are not connected to existing infrastructure.
#' @param tmp.rd.sf \code{sf LINESTRING} object depicting anchor field roads. Created by
#'   \code{\link{lcp_rds_infield}}.
#' @param rd.df \code{data.frame} of road coordinates. If this is the first iteration
#'   through \code{\link{generate_sat_rd}} (i.e., \code{j == 1}) then this will
#'   not exist and will be created by the code. Otherwise, the file previously
#'   created/updated by \code{lcp_rds_outfield} will be used.
#' @param rd.sf \code{sf LINESTRING} object depicting road locations. If this is the first
#'   iteration through \code{\link{generate_sat_rd}} (i.e., \code{j == 1}) then
#'   this will not exist and will be created by the code. Otherwise, the file
#'   previously created/updated by \code{lcp_rds_outfield} will be used.
#' @param exist.coords \code{data.frame} containing coordinates of existing infrastructure.
#'   Created by \code{\link{prep_general_inputs}}.
#' @param road.res Scenario-specific road restriction \code{SpatRaster}. Created by
#'   \code{\link{prep_scenario_inputs}}.
#' @param rd.exist \code{sf LINESTRING} object representing existing roads.
#'   Created by \code{\link{prep_general_inputs}}.
#' @param alt.b.stranded.leases \code{SpatRaster} object identifying stranded lease locations
#'   under Alternative B. Created by \code{\link{prep_scenario_inputs}}. Optional,
#'   as this is only needed when running Alternative B.
#' @param tr.cost.alt.c.row \code{TransitionLayer} for right-of-way (ROW) areas under
#'   Alternative C, accounting for connecting roads that are allowed within some
#'   of the 50% Teshekpuk Caribou Herd (TCH) calving kernel under this alternative.
#'   Created by \code{\link{prep_scenario_inputs}}. Optional, as this is only
#'   needed when running Alternative C.
#' @param alt.d.TLnorth \code{SpatRaster} object identifying available areas for development
#'   north of Teshekpuk Lake under Alternative D. Created by \code{\link{prep_scenario_inputs}}.
#'   Optional, as this is only needed when running Alternative D.
#' @param roadless.check Integer value providing a flag that indicates how
#'   roadless development should be dealt with. Created by \code{\link{generate_sat_rd}}.
#' @inheritParams lcp_rds_infield
#'
#' @return List containing two objects, a \code{data.frame} of generated road coordinates
#'   (\code{rd_df}) and an \code{sf LINESTRING} object depicting generated roads (\code{rd_sf}).
#'
lcp_rds_outfield <- function(x, y = x, tmp.rd.sf, rd.df = NULL, rd.sf = NULL, tr, min.dist.exist,
                             exist.coords, road.res, rd.exist, alt.b.stranded.leases = NULL, tr.cost.alt.c.row = NULL,
                             alt.d.TLnorth = NULL, wd.loc, path.out, scenario, n.iter, roadless.check, proj.info, j, z = NULL,
                             i = NULL, debug.out = FALSE){
  ## If this is the first anchor field to be processed, connect it to the existing infrastructure, otherwise
  ## connect it to the closest infrastructure.

  ## If in Alternative C, use the right-of-way (ROW) transition object to allow connecting roads through part
  ## of the TCH 50% calving kernel.
  if(substr(scenario[z], start=1, stop=5) == "Alt_C") tr <- tr.cost.alt.c.row

  ## Connect the anchor field to existing infrastructure
  if(j == 1){
    outfield.rds <- gen_linkage_rd(x=x, y=y, w=exist.coords, v=sf::st_zm(rd.exist["geometry"]),
                                   road.res=road.res, tmp.rd.sf=tmp.rd.sf, tr=tr, min.dist.exist=min.dist.exist,
                                   roadless.check=roadless.check, alt.b.stranded.leases = alt.b.stranded.leases,
                                   alt.d.TLnorth = alt.d.TLnorth, scenario=scenario, proj.info=proj.info, z=z, j=j)
  } else{
    outfield.rds <- gen_linkage_rd(x=x, y=y, w=rd.df, v=rd.sf, road.res=road.res,
                                   tmp.rd.sf=tmp.rd.sf, tr=tr, min.dist.exist=min.dist.exist, roadless.check=roadless.check,
                                   alt.b.stranded.leases = alt.b.stranded.leases, alt.d.TLnorth = alt.d.TLnorth,
                                   scenario=scenario, proj.info=proj.info, z=z, j=j)
  }

  ## For debugging purposes, this option will write out intermediate data
  if(debug.out) utils::write.csv(outfield.rds$rd_df, file=paste(wd.loc, "/", path.out, "/debug_road_outfield_data_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_cpf_", j, "_", Sys.Date(), ".csv", sep=""), row.names=FALSE)

  ## Return the desired objects
  return(outfield.rds)
}



#' Generate satellite pads and roads for each CPF
#'
#' \code{generate_sat_rd} simulates the location of satellite production pads
#' as a function of relative estimated undiscovered oil, development restrictions,
#' and the presence of waterbodies. These are then connected to their parent
#' central processing facility (CPF) with roads using least cost paths, and the
#' entire "anchor facility" is connected via roads to existing or previously
#' simulated infrastructure.
#'
#' @param cpf.sf \code{sf POINT} object depicting CPF point locations. Created by \code{\link{generate_cpf}}.
#' @param pad.exist \code{sf POINT} object of existing CPF and satellite pad locations.
#'   Created by \code{\link{prep_general_inputs}}.
#' @param pad.prob \code{SpatRaster} layer depicting the scenario-specific relative
#'   likelihood of CPF or satellite pad generation, relative to estimated undiscovered
#'   oil, scenario-specific development restrictions, and waterbodies. Created by
#'   \code{\link{prep_scenario_inputs}}.
#' @param maxd2sat Numeric value indicating the maximum allowed distance (m) from
#'   a CPF to an associated satellite production pad. Defaults to 56,327 m based
#'   on BLM (2019c) p.B-9.
#' @param mind2sat Numeric value indicating the minimum allowed distance (m)
#'   between two satellite production pads, or a CPF and associated satellite
#'   production pad. Defaults to 6437.4 m, similar to what is seen in the Willow
#'   Master Development Plan (BLM 2019b).
#' @param n.sat Desired number of satellite pads to be simulated per CPF. This
#'   can either be a fixed integer value or a range, \code{c(min,max)}, from
#'   which a random value will be drawn each iteration using a uniform
#'   distribution. Defaults to a range of \code{c(4, 8)} satellites per CPF.
#' @param tr.cost.alt Scenario-specific \code{TransitionLayer} indicating road
#'   development cost. Created by \code{\link{prep_scenario_inputs}}.
#' @param tr.cost.alt.b.stranded \code{TransitionLayer} specific to the stranded leases
#'   in Alternative B. Created by \code{\link{prep_scenario_inputs}}. Optional, allowing
#'   the code to be run for scenarios other than Alternative B.
#' @inheritParams lcp_rds_infield
#' @inheritParams lcp_rds_outfield
#' @inheritParams prep_general_inputs
#'
#' @return List containing four objects: a \code{data.frame} of satellite pad locations
#'   (\code{sat_df}), \code{sf POINT} object depicting satellite pad locations (\code{sat_sf}),
#'   \code{data.frame} of road coordinates (\code{rd_df}), and \code{sf LINESTRING}
#'   representation of road locations (\code{rd_sf}).
#'
#' @references BLM \[Bureau of Land Management\] 2019b. Willow Master Development
#'   Plan Draft Environmental Impact Statement. Bureau of Land Management, U.S.
#'   Department of the Interior. Anchorage, AK, USA.
#'
#'  BLM 2019c. Arctic National Wildlife Refuge Coastal Plain Oil and Gas Leasing
#'   Program Final Environmental Impact Statement. Bureau of Land Management,
#'   U.S. Department of the Interior, Anchorage, AK, USA.
generate_sat_rd <- function(cpf.df, cpf.sf, pad.exist, exist.coords, pad.prob, maxd2sat = 56327,
                            mind2sat = 6437.4, n.sat = c(4,8), n.iter = 100, alt.b.stranded.leases = NULL, tr.cost.alt.c.row = NULL,
                            alt.d.TLnorth = NULL, tr.cost.alt, road.res, rd.exist, tr.cost.alt.b.stranded = NULL,
                            wd.loc, path.out, proj.info, scenario, z, i, debug.out = FALSE){
  ## Iterate through the CPFs in order of proximity to existing infrastructure so that they can all be
  ## connected to existing infrastructure via roads
  dist.exist <- sf::st_distance(x=cpf.sf, y=rd.exist)
  min.dist.exist <- data.frame(cpf=1:nrow(cpf.df), min_dist_exist=apply(dist.exist, 1, min))
  min.dist.exist <- min.dist.exist[order(min.dist.exist$min_dist_exist),]

  for(j in 1:nrow(min.dist.exist)){

    #----------------------------------------
    ## Establish satellite wells for each CPF
    #----------------------------------------

    ## Refine the pixel image so that points are only drawn from a circle with radius equal to
    ## maxd2sat, but not within mind2sat of the CPF or other previously-generated satellites.
    tmp.disc <- sf::st_buffer(cpf.sf[min.dist.exist$cpf[j],], maxd2sat)
    tmp.win <- terra::mask(pad.prob, terra::vect(tmp.disc))
    ## Current CPF
    tmp.win <- infrastructure_exclusion_buffer(x=cpf.sf[min.dist.exist$cpf[j],], y=tmp.win, threshold=mind2sat)
    ## Existing/previously-generated satellites
    if(j == 1) sat.sf <- pad.exist[,c("cpf", "sat", "geometry")]
    tmp.win <- infrastructure_exclusion_buffer(x=sat.sf, y=tmp.win, threshold=mind2sat)
    tmp.win.im <- maptools::as.im.RasterLayer(tmp.win)

    ## Determine how many satellite wells to generate
    if(length(n.sat)==1) tmp.n.sat <- n.sat else tmp.n.sat <- sample(n.sat[1]:n.sat[2], size=1)

    ## Generate satellite pad locations within the temporary window, ensuring the minimum distance
    ## (mind2sat) between satellite pads
    sat.coords <- infrastructure_spacer(n.dev=tmp.n.sat, im=tmp.win.im, threshold=mind2sat)
    tmp.sat.df <- data.frame(cpf=min.dist.exist$cpf[j], sat=1:tmp.n.sat, x=sat.coords$x, y=sat.coords$y)
    if(j == 1) sat.df <- tmp.sat.df else sat.df <- rbind(sat.df, tmp.sat.df)

    ## For debugging purposes, this option will write out intermediate data
    if(debug.out) utils::write.csv(sat.df, file=paste(wd.loc, "/", path.out, "/debug_satellite_data_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), ".csv", sep=""), row.names=FALSE)

    ## Add these to the previous satellites sp object
    tmp.sat.sf <- sf::st_as_sf(x=tmp.sat.df, coords=c("x", "y"), crs=proj.info)
    sat.sf <- rbind(sat.sf, tmp.sat.sf)


    #------------------------------------------------------------------------------------------------------------
    ## Connect satellites to their CPF via infield roads and connect each anchor field to existing infrastructure
    #------------------------------------------------------------------------------------------------------------

    ## This needs to be handled differently in Alt B for any pads in the stranded leases, which will
    ## require a different road cost raster, and in Alt D for any pads north of Teshekpuk Lake as
    ## these satellites may feature roadless development and thus not be connected to their
    ## corresponding CPF. Check for these situations and handle them accordingly.

    ## If in Alt B, check whether the current CPF and/or its associated satellites are in the
    ## stranded leases. If so assign a roadless.check value of 1. If in Alt D, check whether
    ## the current CPF and/or its associated satellites are within the roadless area N of T Lake. If
    ## so assign a value of 2. If neither, assign roadless.check a value of 0.
    if(substr(scenario[z], start=1, stop=5) == "Alt_B"){
      cpf.val <- terra::extract(alt.b.stranded.leases, terra::vect(cpf.sf[min.dist.exist$cpf[j],]))
      sat.val <- terra::extract(alt.b.stranded.leases, terra::vect(tmp.sat.sf))
      ## Indicate a value of 1 if at least one pad is in the stranded leases, 0 otherwise
      roadless.check <- ifelse(max(rbind(cpf.val, sat.val)[,2])==0, 0, 1)
    } else if(substr(scenario[z], start=1, stop=5) == "Alt_D"){
      cpf.val <- terra::extract(alt.d.TLnorth, terra::vect(cpf.sf[min.dist.exist$cpf[j],]))
      sat.val <- terra::extract(alt.d.TLnorth, terra::vect(tmp.sat.sf))
      ## Indicate a value of 2 if at least one pad is in the roadless area N of T Lake, 0 otherwise
      roadless.check <- ifelse(max(rbind(cpf.val, sat.val)[,2])==0, 0, 2)
    } else{
      roadless.check <- 0
    }

    ## Run things for non-Alt B/D scenarios and Alt B/D iterations in which the CPF and associated
    ## satellites are not in the areas requiring special action.
    if(roadless.check == 0){
      ## Connect satellites to their CPF with infield roads
      rd.data1 <- lcp_rds_infield(x=tmp.sat.df, tr=tr.cost.alt, road.res=road.res, cpf.df=cpf.df,
                                  min.dist.exist=min.dist.exist, debug.out=debug.out, wd.loc=wd.loc, path.out=path.out,
                                  proj.info=proj.info, scenario=scenario, n.iter=n.iter, j=j, z=z, i=i)
      ## Connect anchor field to existing infrastructure
      if(j == 1){
        rd.data2 <- lcp_rds_outfield(x=rd.data1$tmp_rd_df, tmp.rd.sf=rd.data1$tmp_rd_sf,
                                     tr=tr.cost.alt, min.dist.exist=min.dist.exist, exist.coords=exist.coords,
                                     road.res=road.res, rd.exist=rd.exist, alt.b.stranded.leases=alt.b.stranded.leases,
                                     tr.cost.alt.c.row=tr.cost.alt.c.row, alt.d.TLnorth=alt.d.TLnorth, debug.out=debug.out,
                                     wd.loc=wd.loc, path.out=path.out, scenario=scenario, n.iter=n.iter, proj.info=proj.info,
                                     roadless.check=roadless.check, j=j, z=z, i=i)
      } else{
        rd.data2 <- lcp_rds_outfield(x=rd.data1$tmp_rd_df, tmp.rd.sf=rd.data1$tmp_rd_sf,
                                     rd.df=rd.data2$rd_df, rd.sf=rd.data2$rd_sf, tr=tr.cost.alt,
                                     min.dist.exist=min.dist.exist, exist.coords=exist.coords, road.res=road.res,
                                     rd.exist=rd.exist, alt.b.stranded.leases=alt.b.stranded.leases,
                                     tr.cost.alt.c.row=tr.cost.alt.c.row, alt.d.TLnorth=alt.d.TLnorth, debug.out=debug.out,
                                     wd.loc=wd.loc, path.out=path.out, scenario=scenario, n.iter=n.iter, proj.info=proj.info,
                                     roadless.check=roadless.check, j=j, z=z, i=i)
      }
    }

    ## Run things for situations where one or more pads is within the stranded leases under Alt B
    if(roadless.check == 1){
      ## Connect satellites to their CPF with infield roads
      rd.data1 <- lcp_rds_infield(x=tmp.sat.df, tr=tr.cost.alt.b.stranded, road.res=road.res,
                                  cpf.df=cpf.df, min.dist.exist=min.dist.exist, debug.out=debug.out, wd.loc=wd.loc,
                                  path.out=path.out, scenario=scenario, n.iter=n.iter, proj.info=proj.info, j=j, z=z, i=i)
      ## Connect anchor field to existing infrastructure
      if(j == 1){
        rd.data2 <- lcp_rds_outfield(x=rd.data1$tmp_rd_df, tmp.rd.sf=rd.data1$tmp_rd_sf,
                                     tr=tr.cost.alt.b.stranded, min.dist.exist=min.dist.exist, proj.info=proj.info,
                                     exist.coords=exist.coords, road.res=road.res, rd.exist=rd.exist,
                                     alt.b.stranded.leases=alt.b.stranded.leases, tr.cost.alt.c.row=tr.cost.alt.c.row,
                                     alt.d.TLnorth=alt.d.TLnorth, debug.out=debug.out, wd.loc=wd.loc, path.out=path.out,
                                     scenario=scenario, n.iter=n.iter, roadless.check=roadless.check, j=j, z=z, i=i)
      } else{
        rd.data2 <- lcp_rds_outfield(x=rd.data1$tmp_rd_df, tmp.rd.sf=rd.data1$tmp_rd_sf,
                                     rd.df=rd.data2$rd_df, rd.sf=rd.data2$rd_sf, tr=tr.cost.alt.b.stranded,
                                     min.dist.exist=min.dist.exist, exist.coords=exist.coords, road.res=road.res,
                                     rd.exist=rd.exist, alt.b.stranded.leases=alt.b.stranded.leases,
                                     tr.cost.alt.c.row=tr.cost.alt.c.row, alt.d.TLnorth=alt.d.TLnorth, debug.out=debug.out,
                                     wd.loc=wd.loc, path.out=path.out, scenario=scenario, n.iter=n.iter, proj.info=proj.info,
                                     roadless.check=roadless.check, j=j, z=z, i=i)
      }
    }

    ## Run things for situations in which one or more pad is within the roadless area under Alt D
    if(roadless.check == 2){
      ## Identify any satellites within the same region as the CPF. These can be connected with
      ## roads normally.
      sat.same <- sat.val$ID[sat.val[,2] == cpf.val[,2]]
      sat.same.df <- tmp.sat.df[tmp.sat.df$sat %in% sat.same,]

      ## As long as there is at least one record, connect it with roads normally
      if(nrow(sat.same.df) >= 1){
        ## Connect satellites to their CPF with infield roads
        rd.data1 <- lcp_rds_infield(x=sat.same.df, tr=tr.cost.alt, road.res=road.res, cpf.df=cpf.df,
                                    min.dist.exist=min.dist.exist, debug.out=debug.out, wd.loc=wd.loc, path.out=path.out,
                                    proj.info=proj.info, scenario=scenario, n.iter=n.iter, j=j, z=z, i=i)
      }

      ## Connect any satellites that are in the same region but not connected to the CPF
      sat.dif <- sat.val[sat.val[,2] != cpf.val[,2],]
      sat.dif.check <- table(sat.dif[,2])
      sat.dif.check2 <- names(sat.dif.check)[sat.dif.check > 1]
      if(length(sat.dif.check2) > 0){
        tmp.sat.df2 <- tmp.sat.df
        tmp.sat.df2$val <- sat.val[,2]
        sat.dif.df <- tmp.sat.df2[tmp.sat.df2$val %in% as.numeric(sat.dif.check2),]
        ## This needs to be done separately for each region
        for(aa in 1:length(unique(sat.dif.df$val))){
          sat.dif.df2 <- sat.dif.df[sat.dif.df$val == unique(sat.dif.df$val)[aa],]
          ## If there are only two satellites, simply connect them
          if(nrow(sat.dif.df2) == 2){
            tmp.rd <- gen_lcp_rd(tr=tr.cost.alt, loc1=c(sat.dif.df2$x[1], sat.dif.df2$y[1]),
                                 loc2=c(sat.dif.df2$x[2], sat.dif.df2$y[2]), cur.cpf=min.dist.exist$cpf[j],
                                 cur.sat=sat.dif.df2$sat[1], proj.info=proj.info)
            if(aa==1 & nrow(sat.same.df) == 0){
              rd.data1 <- list("tmp_rd_df"=tmp.rd$df_tmp, "tmp_rd_sf"=tmp.rd$sf_tmp)
            } else{
              rd.data1$tmp_rd_df <- rbind(rd.data1$tmp_rd_df, tmp.rd$df_tmp)
              rd.data1$tmp_rd_sf <- rbind(rd.data1$tmp_rd_sf, tmp.rd$sf_tmp)
            }
          } else{
            ## Otherwise, connect them in order of proximity to the satellite furthest west
            ## (note that furthest west is arbitrary, but makes them proceed in a reasonable
            ## order to minimize criscrossing of roads).

            ## First identify the satellite that is farthest west. This will be the one with
            ## the lowest value of x. Isolate this from the other satellites.
            west.sat <- sat.dif.df2[which.min(sat.dif.df2$x),]
            sat.dif.df3 <- sat.dif.df2[sat.dif.df2$sat != west.sat$sat,]
            ## Identify the next closest satellite
            sat.dif.df3$dist_west <- raster::pointDistance(p1=cbind(west.sat$x, west.sat$y),
                                                           p2=cbind(sat.dif.df3$x, sat.dif.df3$y), lonlat=FALSE)
            sat.dif.df3 <- sat.dif.df3[order(sat.dif.df3$dist_west),]
            ## Connect this
            tmp.rd <- gen_lcp_rd(tr=tr.cost.alt, loc1=c(west.sat$x, west.sat$y),
                                 loc2=c(sat.dif.df3$x[1], sat.dif.df3$y[1]), cur.cpf=min.dist.exist$cpf[j],
                                 cur.sat=sat.dif.df3$sat[1], proj.info=proj.info)
            tmp.rd.df.dif <- tmp.rd$df_tmp
            if(aa==1 & nrow(sat.same.df) == 0){
              rd.data1 <- list("tmp_rd_df"=NULL, "tmp_rd_sf"=tmp.rd$sf_tmp)
            } else{
              rd.data1$tmp_rd_sf <- rbind(rd.data1$tmp_rd_sf, tmp.rd$sf_tmp)
            }
            ## Now iterate through the remaining satellites and connect them to the nearest
            ## satellite or road
            for(bb in 2:nrow(sat.dif.df3)){
              ## Identify the candidate coordinates among existing roads and satellites for
              ## closest point
              cand.coords <- rbind(data.frame(dev="sat", num=c(west.sat$sat,
                                                               sat.dif.df3$sat[1:(bb-1)]), x=c(west.sat$x, sat.dif.df3$x[1:(bb-1)]),
                                              y=c(west.sat$y, sat.dif.df3$y[1:(bb-1)])),
                                   data.frame(dev="rd", num=1:nrow(tmp.rd.df.dif), x=tmp.rd.df.dif$x,
                                              y=tmp.rd.df.dif$y))
              cand.coords$dist_sat <- raster::pointDistance(p1=cbind(sat.dif.df3$x[bb],
                                                                     sat.dif.df3$y[bb]),	p2=cbind(cand.coords$x, cand.coords$y), lonlat=FALSE)
              min.cand <- cand.coords[which.min(cand.coords$dist_sat),]

              ## Connect to that feature using least cost paths, as long as the satellite is
              ## at least one pixel away from the feature to which it is being connected.
              ## Otherwise, the pad is already connected and no road is needed.
              if(min.cand$dist_sat > min(raster::res(road.res))){
                tmp.rd <- gen_lcp_rd(tr=tr.cost.alt, loc1=c(sat.dif.df3$x[bb], sat.dif.df3$y[bb]),
                                     loc2=c(min.cand$x, min.cand$y), cur.cpf=min.dist.exist$cpf[j],
                                     cur.sat=sat.dif.df3$sat[bb], proj.info=proj.info)
                tmp.rd.df.dif <- rbind(tmp.rd.df.dif, tmp.rd$df_tmp)
                rd.data1$tmp_rd_sf <- rbind(rd.data1$tmp_rd_sf, tmp.rd$sf_tmp)
              }
            }
            if(aa==1 & nrow(sat.same.df) == 0){
              rd.data1$tmp_rd_df <- tmp.rd.df.dif
            } else{
              rd.data1$tmp_rd_df <- rbind(rd.data1$tmp_rd_df, tmp.rd.df.dif)
            }
          }
        }
      }

      ## For debugging purposes, this option will write out intermediate data
      if(debug.out) utils::write.csv(rd.data1$tmp_rd_df, file=paste(wd.loc, "/", path.out, "/debug_road_infield_data_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_cpf_", j, "_", Sys.Date(), ".csv", sep=""), row.names=FALSE)


      #-------------------------------------------------
      ## Connect anchor field to existing infrastructure
      #-------------------------------------------------

      ## Under roadless development, this will only be done for pads with a value of zero. These
      ## will be connected to the closest infrastructure.

      ## Identify the satellite pads able to be connected to other infrastructure
      sat.0s <- sat.val$ID[sat.val[,2] == 0]

      ## If these exist, pull the roads associated with those pads
      if(length(sat.0s) > 0){
        rd.0s <- rd.data1$tmp_rd_df[rd.data1$tmp_rd_df$sat %in% sat.0s,]
      } else{
        ## Otherwise, check if the CPF has a value of zero (this would be very rare, but could
        ## happen if the CPF is below the roadless area and all satellites are above it. In such
        ## a case the CPF should be connected by road to existing development, even though the
        ## satellites are not. For consistency of code this is called rd.0s too.
        if(cpf.val[,2] == 0){
          rd.0s <- cpf.df[min.dist.exist$cpf[j],]
        } else{
          rd.0s <- data.frame()
        }
      }

      ## If there are roads (or a CPF) to connect, then connect the anchor field to existing
      ## infrastructure. As noted above, this will be done differently depending on whether this is
      ## the first CPF to be processed or not.
      if(nrow(rd.0s) > 0){
        if(j == 1){
          rd.data2 <- lcp_rds_outfield(x=rd.0s, y=rd.data1$tmp_rd_df, tmp.rd.sf=rd.data1$tmp_rd_sf,
                                       tr=tr.cost.alt, min.dist.exist=min.dist.exist, proj.info=proj.info,
                                       exist.coords=exist.coords, road.res=road.res, rd.exist=rd.exist,
                                       alt.d.TLnorth=alt.d.TLnorth, debug.out=debug.out, wd.loc=wd.loc, path.out=path.out,
                                       scenario=scenario, n.iter=n.iter, roadless.check=roadless.check, j=j, z=z, i=i)
        } else{
          rd.data2 <- lcp_rds_outfield(x=rd.0s, y=rd.data1$tmp_rd_df, tmp.rd.sf=rd.data1$tmp_rd_sf,
                                       rd.df=rd.data2$rd_df, rd.sf=rd.data2$rd_sf, tr=tr.cost.alt, proj.info=proj.info,
                                       min.dist.exist=min.dist.exist, exist.coords=exist.coords, road.res=road.res,
                                       rd.exist=rd.exist, alt.d.TLnorth=alt.d.TLnorth, debug.out=debug.out,
                                       wd.loc=wd.loc, path.out=path.out, scenario=scenario, n.iter=n.iter,
                                       roadless.check=roadless.check, j=j, z=z, i=i)
        }
      } else{
        ## This situation should be rare, only being implemented if the CPF and all associated
        ## satellites are north of T Lake and so no part will be connected to existing roads
        ## (i.e., nrow(rd.0s) == 0). This allows an rd.data2 object to be created without a call
        ## to lcp_rds_outfield() so that it can be returned in the output of generate_sat_rd().
        if(j == 1){
          rd.df <- rbind(exist.coords, rd.data1$tmp_rd_df)
          rd.sf <- rbind(rd.exist, rd.data1$tmp_rd_sf)
        } else{
          rd.df <- rbind(rd.data2$rd_df, rd.data1$tmp_rd_df)
          rd.sf <- rbind(rd.data2$rd_sf, rd.data1$tmp_rd_sf)
        }
        rd.data2 <- list("rd_df"=rd.df, "rd_sf"=rd.sf)
      }
    }
  }

  ## Return the final data
  return(list("sat_df"=sat.df, "sat_sf"=sat.sf, "rd_df"=rd.data2$rd_df, "rd_sf"=rd.data2$rd_sf))
}
