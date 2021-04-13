#' Wrapper function to perform the DIA analysis
#'
#' The \code{dia} function is a wrapper that implements the overall Development
#' Impacts Analysis (DIA) model (see Fullman et al. in press for details). The
#' code below is designed to analyze the proposed development alternatives in
#' BLM (2019a), as well as one community-generated proposal (see Fullman et al.
#' in press). It can, however, be adapted to other development scenarios.
#' Infrastructure simulation and species impact analyses can be run individually
#' or jointly.
#'
#' @param wd.loc Character string identifying the base directory containing both
#'   input and output data folders. Defaults to the current working directory,
#'   identified by \code{base::getwd}.
#' @param scenario vector of character strings identifying the names of the
#'   scenarios to be run.
#' @param simulate.inf Logical indicator of whether infrastructure simulation
#'   should be conducted. Defaults to \code{TRUE}.
#' @param n.iter Integer indicating the desired number of iterations to be run
#'   for each scenario. Defaults to 100.
#' @param n.cpf Indicator of the desired number of central processing facilities
#'   (CPFs) to be generated. This can either be a single integer value, indicating
#'   a fixed number of CPFs to be used in all iterations, or a range, \code{c(min,max)},
#'   from which a random value will be drawn using a uniform distribution.
#'   Defaults to \code{c(3,7)} to reflect a reasonable range of potential
#'   variability, depending on the number and distribution of undiscovered oil
#'   deposits.
#' @param d2cpf Numeric value indicating the minimum allowed distance (m) between
#'   two CPFs. Defaults to 35,000 m.
#' @param maxd2sat Numeric value indicating the maximum allowed distance (m) from
#'   a CPF to an associated satellite production pad. Defaults to 56,327 m based
#'   on BLM (2019c) p.B-9.
#' @param mind2sat Numeric value indicating the minimum allowed distance (m) between
#'   two satellite production pads, or a CPF and associated satellite production
#'   pad. Defaults to 6437.4 m, similar to what is seen in BLM (2019b).
#' @param n.sat Desired number of satellite pads per CPF. This can either be fixed
#'   or a range, \code{c(min,max)}, from which a random value will be drawn for
#'   each iteration using a uniform distribution. Defaults to \code{c(4,8)}
#'   satellites per CPF.
#' @param rd.exist.file Character string indicating the file name for the shapefile
#'   depicting existing roads as lines. Optional, allowing the code to be run
#'   without infrastructure generation. However, if infrastructure generation is
#'   desired (i.e., \code{simulate.inf = TRUE}), then this is required.
#' @param pad.exist.file Character string indicating the file name for the
#'   shapefile depicting existing CPF and satellite pads as points. Optional,
#'   allowing the code to be run without infrastructure generation. However, if
#'   \code{simulate.inf = TRUE}, then this is required.
#' @param oil.av.file Character string indicating the file name for the relative
#'   expected undiscovered oil raster. Optional, allowing the code to be run
#'   without infrastructure generation. However, if \code{simulate.inf = TRUE},
#'   then this is required.
#' @param cost.map.file Character string indicating the file name for the cost
#'   map raster, which depicts water areas coded as \code{0} and land as \code{1}.
#'   Optional, allowing the code to be run without infrastructure generation.
#'   However, if \code{simulate.inf = TRUE}, then this is required.
#' @param pad.res.files Character string indicating the file names for the rasters
#'   indicating scenario-specific restrictions for pad (i.e., CPF and satellite
#'   pad) placement. Values may consist of 0, 0.1, or 1. \strong{The length and order of
#'   this object must match that of \code{scenario}}. Optional, allowing the code
#'   to be run without infrastructure generation. However, if \code{simulate.inf = TRUE},
#'   then this is required.
#' @param road.res.files Character string indicating the file names for the rasters
#'   indicating scenario-specific restrictions for road placement. Values may
#'   consist of 0, 0.1, or 1. \strong{The length and order of this object must
#'   match that of \code{scenario}}. Optional, allowing the code to be run without
#'   infrastructure generation. However, if \code{simulate.inf = TRUE}, then this
#'   is required.
#' @param alt.b.rd.stranded.res.file Character string indicating the file name
#'   for the road development restriction raster for stranded leases (i.e., those
#'   completely surrounded by areas closed to leasing) under Alternative B.
#'   Optional, allowing the code to be run without infrastructure generation or
#'   for scenarios other than Alternative B. However, if Alternative B is to be
#'   run this must be specified.
#' @param alt.b.stranded.lease.file Character string indicating the file name for
#'   the raster identifying stranded leases under Alternative B. Optional,
#'   allowing the code to be run without infrastructure generation or for
#'   scenarios other than Alternative B. However, if Alternative B is to be
#'   run this must be specified.
#' @param alt.c.row.file Character string indicating the file name for the raster
#'   identifying right-of-way (ROW) areas under Alternative C where infield roads
#'   are prohibited but connecting infrastructure is allowed (see BLM 2019a for
#'   details). Optional, allowing the code to be run without infrastructure
#'   generation or for scenarios other than Alternative C. However, if Alternative
#'   C is to be run this must be specified.
#' @param alt.d.north.file Character string indicating the file name for the
#'   raster identifying discrete areas north of Teshekpuk Lake where roadless
#'   development is possible under Alternative D (see BLM 2019a for details).
#'   Optional, allowing the code to be run without infrastructure generation or
#'   for scenarios other than Alternative D. However, if Alternative D is to be
#'   run this must be specified.
#' @param npra.file Character string indicating the file name of the NPR-A
#'   boundary shapefile. Used for calculating the infrastructure summary.
#'   Optional, allowing the code to be run without the impact analysis. However,
#'   if the impact analysis is desired this must be specified.
#' @param tch.raster Character string indicating the Teshekpuk Caribou Herd (TCH)
#'   calving raster file name. Optional, allowing the analysis to be run
#'   without the TCH or without any impact analyses. Failure to specify this or
#'   \code{wah.raster} means the caribou impact analysis will not be run.
#' @param wah.raster Vector of two character strings indicating Western Arctic
#'   Herd (WAH) raster file names. The first should specify the calving resource
#'   selection function (RSF) raster, the second the weighting raster. See
#'   Fullman et al. (in press) for details. Optional, allowing the analysis to
#'   be run without the WAH or without any impact analyses. Failure to specify
#'   this or \code{tch.raster} means the caribou impact analysis will not be run.
#' @param shorebird.raster Vector of character strings indicating shorebird
#'   habitat suitability index (HSI) raster file names. Each should be for a different
#'   species. Optional, allowing the analysis to be run without shorebirds or
#'   without any impact analyses. Failure to specify this means the shorebird
#'   impact analysis will not be run.
#' @param shorebird.threshold Vector of numeric values indicating species-specific
#'   thresholds defining what shorebird HSI values equate to suitable versus
#'   unsuitable habitat (Saalfeld et al 2013). Values should range between 0-1.
#'   Optional, only used if a vector is given for \code{shorebird.raster}. Order
#'   of threshold values must be the same as the order of species in \code{shorebird.raster}.
#' @param brant.shp Character string indicating the name of the brant lakes
#'   shapefile, without file extension. Optional, allowing the analysis to be run
#'   without brant analyses or without any impact analyses. Failure to specify
#'   this means the brant impact analysis will not be run.
#' @param hq.quant Numeric value between 0-1 that indicates the pixel value that
#'   reflects high-quality habitat. Defaults to 0.75, following Johnson et al. (2005)'s
#'   use of the upper quartile to represent high-quality habitat.
#' @param zoi Numeric value indicating the zone of influence (m) for calculation
#'   of the surface disturbance buffer and caribou displacement. Defaults to
#'   4000 m, following Cameron et al. (2005) and Wilson et al. (2013).
#' @param area.cpf Numeric value indicating the average size (sq. m) of a CPF
#'   footprint. Default value equates to 100 acres based on BLM (2019b).
#' @param area.sat Numeric value indicating the average size (sq. m) of a
#'   satellite production facility footprint. Default value equates to 15 acres based
#'   on BLM (2019a) p.B-6.
#' @param road.width Numeric value indicating the estimated road ground
#'   footprint width (m). Default value based on BLM (2019a) p.B-6. See Fullman
#'   et al. (in press) for details.
#' @param land.dist Numeric value indicating the distance (m) from CPF or
#'   satellite pads at which helicopters are expected to be under 500 m altitude,
#'   based on a descent angle of 10 degrees following FAA guidelines. Used to
#'   determine helicopter impact areas around CPFs and satellite pads in the brant
#'   impact analysis. Defaults to 2835.64 m based on the assumption that a helicopter
#'   approaching a landing zone at a 10 degree descent angle reaches a 500 m
#'   altitude at \eqn{500m / tan(10 deg) = 2835.64 m} out from the landing zone.
#' @param heli.disturb Numeric value indicating the distance (m) from helicopters
#'   (under 500 m altitude) at which molting brant are expected to be disturbed.
#'   Defaults to 3570 m based on Jensen (1990) and Miller et al. (1994).
#' @param proximity.effect Numeric value indicating the distance (m) at which
#'   birds are disturbed by infrastructure proximity effects (e.g., dust
#'   deposition). Defaults to 100 m as a conservative estimate.
#' @param proj.info Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
#'   common to all spatial objects in the analysis. Defaults to the NAD83 2011 Alaska
#'   Albers projection (\code{EPSG:6393}).
#' @param path.in Character string identifying the relative path from the base
#'   directory (\code{wd.loc}) to where input data are stored. Defaults to an
#'   "Input_Data" folder within the base folder.
#' @param path.out Character string identifying the relative path from the base
#'   directory (\code{wd.loc}) to where output data will be saved. Defaults to
#'   an "Output_Data" folder within the base folder.
#' @param shp.out Logical indicator of whether infrastructure created if
#'   \code{simulate.in = TRUE} should be written out as shapefiles to \code{path.out}.
#'   If so, three shapefiles will be written, one each for CPFs, satellite pads, and
#'   roads. Defaults to \code{TRUE}.
#' @param caribou.out Logical indicator of whether the resulting discounted
#'   caribou RSF raster(s) should be written out. Defaults to \code{FALSE}.
#' @param brant.out Logical indicator of whether a shapefile of the affected lakes should
#'   be written to \code{path.out}. Only written if any lakes are affected. Defaults to
#'   \code{FALSE}.
#' @param suppress.rgdal.proj4.warnings Logical indicator of whether \href{https://cran.r-project.org/package=rgdal}{rgdal}
#'   warnings about discarded datums in coordinate reference system definitions
#'   should be suppressed. These are simply to warn developers about the switch
#'   from \code{proj4string} to \code{WKT CRS} specification, but generate many
#'   unnecessary warnings in the \code{dia} code. Defaults to \code{TRUE},
#'   indicating that warnings will be suppressed.
#' @param run.in.parallel Logical indicator of whether the impact analyses will be
#'   conducted in parallel. Defaults to \code{TRUE}. Used to determine whether
#'   \code{SpatRaster} outputs are packed for distribution among parallel workers.
#'   In future versions of this code, this will also make parallel running optional,
#'   but this is still under development.
#' @param n.cores Integer value, indicating the number of cores to be used for
#'   parallel processing. Defaults to the total number of available cores, as
#'   determined by \code{parallel::detectCores()}. In situations where RAM is limiting,
#'   however, it may be helpful to reduce the number of cores used to run
#'   parallel iterations.
#' @param debug.out Logical indicator of whether intermediate infrastructure .csv
#'   files should be written out for the purpose of debugging code. Defaults to
#'   \code{FALSE}.
#' @param completion.sound Indication of whether the computer should make a sound
#'   when the model run is complete and, if so, which sound to play. Defaults to
#'   no sound (\code{NULL}). Valid inputs are integers corresponding to the sound
#'   argument options for \code{\link[beepr]{beep}}. If \code{!is.null(completion.sound)},
#'   then the \href{https://cran.r-project.org/package=beepr}{beepr} package must be
#'   installed.
#'
#' @details The \code{dia} code assumes there is a \code{path.out/RData}
#'   folder to contain output results. If this does not exist, it will be created by
#'   the \code{dia} function.
#'
#'   \code{dia} is set up to address discrepancies in spatial object
#'   projection, but is not designed to deal with misaligned extent or cell size of
#'   \code{SpatRaster} input data. Prior to using \code{dia}, the user should perform
#'   initial preparation to ensure any raster files have the same extent and
#'   resolution. For tools to address such data preparation, see the \href{https://cran.r-project.org/package=terra}{terra}
#'   package.
#'
#'   To speed up analyses, code is run in parallel. However, infrastructure
#'   simulation can be memory intensive, causing code to fail if all available
#'   cores are used in parallel processing and if \code{simulate.inf = TRUE}.
#'   If this occurs, it is advisable to reduce \code{n.cores}. In tests on our high
#'   performance modeling computer, we found modifying \code{n.cores} so that there
#'   were about 6.5 GB available RAM per core yielded a good compromise, though
#'   further testing is needed and this may vary widely across systems and configurations.
#'
#' @return Nothing is directly returned by \code{dia}. However, if \code{simulate.inf
#'   = TRUE} then .csv files, shapefiles (if \code{shp.out = TRUE}), and an .RData file
#'   containing the results of the infrastructure simulation will be written to
#'   \code{path.out}. If any of the species impact analyses are run, a .csv file
#'   containing the infrastructure summary and species-specific impact results will be
#'   written to \code{path.out}. Regardless, the run times for each iteration are saved
#'   in an .RData file and the overall model run time is saved as a .csv file.
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar%
#' @export
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
#' BLM 2019c. Arctic National Wildlife Refuge Coastal Plain Oil and Gas Leasing
#'   Program Final Environmental Impact Statement. Bureau of Land Management,
#'   U.S. Department of the Interior, Anchorage, AK, USA.
#'
#' Cameron RD, Smith WT, White RG, Griffith B. 2005. Central Arctic
#'   caribou and petroleum development: Distributional, nutritional, and
#'   reproductive implications. Arctic 58:1-9.
#'
#' Fullman TJ, Sullender BK, Cameron MD, Joly K. in press. Simulation modeling
#'   accounts for uncertainty while quantifying ecological effects of development
#'   alternatives. Ecosphere.
#'
#' Jensen KC. 1990. Responses of molting Pacific black brant to
#'   experimental aircraft disturbance in the Teshekpuk Lake Special Area, Alaska.
#'   Ph.D. Thesis, Texas A&M University. College Station, TX, USA.
#'
#' Johnson CJ, Boyce MS, Ray CL, Cluff DH, Gau RJ, Gunn A, Mulders R. 2005.
#'   Cumulative effects of human developments on arctic wildlife. Wildlife
#'   Monographs 160:1-36.
#'
#' Miller MW, Jensen KC, Grant WE, Weller MW. 1994. A simulation model of
#'   helicopter disturbance of molting Pacific black brant. Ecological Modelling
#'   73:293-309.
#'
#' Saalfeld ST, Lanctot RB, Brown SC, Saalfeld DT, Johnson JA, Andres BA, Bart JR.
#'   2013. Predicting breeding shorebird distributions on the Arctic Coastal
#'   Plain of Alaska. Ecosphere 4:16.
#'
#' Wilson RR, Liebezeit JR, Loya WM. 2013. Accounting for uncertainty in oil and
#'   gas development impacts to wildlife in Alaska. Conservation Letters
#'   6:350-358.
#'
#' @examples
#' \dontrun{
#'   dia(wd.loc="C:/DIA", scenario=c("Alt_Aplus", "Alt_A", "Alt_B", "Alt_C", "Alt_D"),
#'     rd.exist.file="NPRA_existing_roads", pad.exist.file="NPRA_existing_pads",
#'     oil.av.file="NPRA_oil_MMBO_Mean_per_pix.tif", cost.map.file="costmap_base.tif",
#'     pad.res.files=c("AltAplus_restrictions_pad.tif", "AltA_restrictions_pad.tif",
#'       "AltB_restrictions_pad.tif", "AltC_restrictions_pad.tif",
#'       "AltD_restrictions_pad.tif"),
#'     road.res.files=c("AltAplus_restrictions_road.tif", "AltA_restrictions_road.tif",
#'       "AltB_restrictions_road.tif", "AltC_restrictions_road.tif",
#'       "AltD_restrictions_road.tif"),
#'     npra.file="NPRA_Boundary", tch.raster="tch_calving.tif",
#'     wah.raster=c("WAH_calving.tif", "WAH_calving_overlap_weight.tif"),
#'     shorebird.raster=c("SESA_HSI_NPRA.tif", "AMGP_HSI_NPRA.tif", "BBPL_HSI_NPRA.tif",
#'       "DUNL_HSI_NPRA.tif", "LBDO_HSI_NPRA.tif", "PESA_HSI_NPRA.tif",
#'       "REPH_HSI_NPRA.tif", "RNPH_HSI_NPRA.tif"),
#'     shorebird.threshold=c(0.06, 0.56, 0.31, 0.05, 0.05, 0.13, 0.05, 0.62),
#'     brant.shp="Lakes_W_Data", n.cores = 30, completion.sound=4)
#' }
dia <- function(wd.loc=getwd(), scenario, simulate.inf = TRUE, n.iter = 100, n.cpf = c(3,7), d2cpf = 35000,
                maxd2sat = 56327, mind2sat = 6437.4, n.sat = c(4,8), rd.exist.file = NULL, pad.exist.file = NULL,
                oil.av.file = NULL, cost.map.file = NULL, pad.res.files = NULL, road.res.files = NULL,
                alt.b.rd.stranded.res.file = NULL, alt.b.stranded.lease.file = NULL, alt.c.row.file = NULL,
                alt.d.north.file = NULL, npra.file = NULL, tch.raster = NULL, wah.raster = NULL, shorebird.raster = NULL,
                shorebird.threshold = NULL, brant.shp = NULL, hq.quant = 0.75, zoi = 4000, area.cpf = 100 * 4046.86,
                area.sat = 15 * 4046.86, road.width = 18.85953, land.dist = 2835.64, heli.disturb = 3570,
                proximity.effect = 100, proj.info = "EPSG:6393", path.in = "Input_Data", path.out = "Output_Data", shp.out = TRUE,
                caribou.out = FALSE, brant.out = FALSE, suppress.rgdal.proj4.warnings = TRUE, run.in.parallel = TRUE,
                n.cores = parallel::detectCores(), debug.out = FALSE, completion.sound = NULL){


  ## Suppress rgdal warnings about dropping datums, if desired
  if(suppress.rgdal.proj4.warnings) withr::local_options("rgdal_show_exportToProj4_warnings"="none")


  ############################################
  ## Prepare input data needed for all runs ##
  ############################################

  #-------------------------------
  ## For infrastructure simulation
  #-------------------------------

  ## Prepare inputs needed across all scenarios for infrastructure simulation, if needed
  if(simulate.inf){
    ## Prepare general inputs
    gen.inputs <- prep_general_inputs(rd.exist.file = rd.exist.file, pad.exist.file = pad.exist.file,
                                      oil.av.file = oil.av.file, d2cpf = d2cpf, wd.loc = wd.loc, path.in = path.in, proj.info = proj.info)
  }


  #-----------------------------
  ## For caribou impact analysis
  #-----------------------------

  ## Prep TCH data, if needed
  if(!is.null(tch.raster)){
    ## Load TCH calving raster
    tch.calving <- load_spatial(x=tch.raster, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)

    ## Calculate the number of high-quality calving pixels (sensu Johnson et al. 2005).
    tch.hq <- calc_highquality(tch.calving, hq.quant=hq.quant)

    ## If running in parallel below, pack the SpatRaster object so it can be distributed among parallel workers
    if(run.in.parallel) tch.calving <- terra::pack(tch.calving)
  }

  ## Prep WAH data, if needed
  if(!is.null(wah.raster)){
    ## Load WAH calving RSF and weighting rasters
    wah.calving <- load_spatial(x=wah.raster[1], proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)
    wah.hq.weight <- load_spatial(x=wah.raster[2], proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)

    ## Calculate the number of high-quality calving pixels (sensu Johnson et al. 2005).
    wah.hq <- calc_highquality(x=wah.calving, y=wah.hq.weight, wah=TRUE, hq.quant=hq.quant)

    ## If running in parallel below, pack the SpatRaster objects for distribution among parallel workers
    if(run.in.parallel){
      wah.calving <- terra::pack(wah.calving)
      wah.hq.weight <- terra::pack(wah.hq.weight)
    }
  }


  #-------------------------------
  ## For shorebird impact analysis
  #-------------------------------

  ## Prep shorebird data, if needed
  if(!is.null(shorebird.raster)){
    ## Run through each indicated shorebird dataset, loading the rasters, preparing the necessary data,
    ## stacking the results, and preparing a summary data.frame.
    for(ii in 1:length(shorebird.raster)){
      ## Load a shorebird habitat suitability index (HSI) raster and ensure it matches the desired projection.
      sb.hsi <- load_spatial(x=shorebird.raster[ii], proj.info=proj.info, wd.loc=wd.loc, path.in=path.in)

      ## Ensure the correct projection and determine suitabile and high-quality suitable pixels,
      ## returning the number of high-quality suitable pixels.
      sb.hq <- calc_highquality(x=sb.hsi, z=shorebird.threshold[ii], sb=TRUE, hq.quant=hq.quant)

      ## Create a temporary data.frame containing the high-quality habitat information and a blank to
      ## be filled with the impact results. This code assumes that the shorebird raster names each start
      ## with the four-letter species code.
      sb.df.tmp <- data.frame("sp"=substr(shorebird.raster[ii], start=1, stop=4), "quant_hq"=sb.hq$quant_hq,
                              "hq_orig"=sb.hq$highquality_orig, "hq_discount"=NA, "hq_remaining"=NA)

      ## Add the individual species data.frame to a combined shorebird data.frame and the HSI raster to
      ## a combined shorebird raster stack
      if(ii == 1){
        sb.df <- sb.df.tmp
        sb.stack <- sb.hsi
      } else{
        sb.df <- rbind(sb.df, sb.df.tmp)
        sb.stack <- c(sb.stack, sb.hsi)
      }
    }

    ## If running in parallel below, pack the SpatRaster object for distribution among parallel workers
    if(run.in.parallel) sb.stack <- terra::pack(sb.stack)
  }


  #---------------------------
  ## For brant impact analysis
  #---------------------------

  ## Prep brant data, if needed
  if(!is.null(brant.shp)){
    ## Load brant lakes shapefile containing lake boundaries with molting brant data.
    brant.lakes <- load_spatial(x=brant.shp, proj.info=proj.info, wd.loc=wd.loc, path.in=path.in, vec=TRUE)

    ## Buffer by zero to avoid typology issues.
    brant.lakes <- rgeos::gBuffer(sf::as_Spatial(brant.lakes), byid=TRUE, width=0)
    brant.lakes <- sf::st_as_sf(brant.lakes)
  }


  ###########################################################
  ## Run infrastructure simulation and/or impacts analyses ##
  ###########################################################

  ## Loop through each scenario and run the iterations within scenarios in parallel
  for(z in 1:length(scenario)){
    scenario.start <- Sys.time()

    ## Extra preparations if infrastructure simulation is to be run
    if(simulate.inf){
      ## Prepare the scenario-specific inputs
      scen.inputs <- prep_scenario_inputs(scenario = scenario, oil.av = gen.inputs$oil_av,
                                          cost.map.file = cost.map.file, willow.cpf.ras = gen.inputs$willow_cpf_ras, pad.res.file = pad.res.files[z],
                                          road.res.file = road.res.files[z], alt.b.rd.stranded.res.file = alt.b.rd.stranded.res.file,
                                          alt.b.stranded.lease.file = alt.b.stranded.lease.file, alt.c.row.file = alt.c.row.file,
                                          alt.d.north.file = alt.d.north.file, wd.loc = wd.loc, path.in = path.in, proj.info = proj.info, z = z,
                                          run.in.parallel = run.in.parallel)

      ## Reduce memory usage by removing the first two objects from gen.inputs, as these are only used by
      ## prep_scenario_inputs() and are not needed in the impact analysis below
      gen.inputs <- gen.inputs[3:length(gen.inputs)]
    }

    ## Set things up to run in parallel.
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl, cores=n.cores)

    ## Run the infrastructure simulation and/or impact analysis in parallel.
    iter.run.times <- foreach::foreach(i = 1:n.iter, .packages=c("raster", "gdistance", "spatstat", "maptools", "sp",
                                                                 "rgeos", "rgdal", "sf", "terra", "units", "udunits2"), .export=c("footprint_generation", "gen_lcp_rd",
                                                                                                                                  "gen_linkage_rd", "generate_cpf", "generate_sat_rd", "impact_brant", "impact_caribou", "impact_shorebird",
                                                                                                                                  "inf_summary", "infrastructure_exclusion_buffer", "infrastructure_proximity_discounting",
                                                                                                                                  "infrastructure_spacer", "lcp_rds_infield", "lcp_rds_outfield", "load_spatial", "poly_rotate",
                                                                                                                                  "projection_alignment", "pt_to_pad", "raster_to_zero"), .inorder=FALSE) %dorng% {
                                                                                                                                    indiv.start <- Sys.time()

                                                                                                                                    #########################################
                                                                                                                                    ## Simulate infrastructure, if desired ##
                                                                                                                                    #########################################

                                                                                                                                    if(simulate.inf){
                                                                                                                                      ## Generate CPF location(s)
                                                                                                                                      cpf.iter <- generate_cpf(pad.im = scen.inputs$pad_im, n.cpf = n.cpf, d2cpf = d2cpf,
                                                                                                                                                               debug.out = debug.out, wd.loc = wd.loc, path.out = path.out, n.iter = n.iter,
                                                                                                                                                               proj.info = proj.info, scenario = scenario, z = z, i = i)

                                                                                                                                      ## Generate satellite pads and roads for each CPF, including connecting to existing
                                                                                                                                      ## infrastructure, after unpacking SpatRaster objects if necessary
                                                                                                                                      if(run.in.parallel){
                                                                                                                                        scen.inputs$pad_prob <- terra::rast(scen.inputs$pad_prob)
                                                                                                                                        scen.inputs$road_res <- terra::rast(scen.inputs$road_res)
                                                                                                                                        if(!is.null(scen.inputs$alt_B_stranded_leases)) scen.inputs$alt_B_stranded_leases <- terra::rast(scen.inputs$alt_B_stranded_leases)
                                                                                                                                        if(!is.null(scen.inputs$alt_D_north)) scen.inputs$alt_D_north <- terra::rast(scen.inputs$alt_D_north)
                                                                                                                                      }
                                                                                                                                      satrd.iter <- generate_sat_rd(cpf.df = cpf.iter$cpf_df, cpf.sf = cpf.iter$cpf_sf,
                                                                                                                                                                    pad.exist = gen.inputs$pad_exist, exist.coords = gen.inputs$exist_coords,
                                                                                                                                                                    pad.prob = scen.inputs$pad_prob, alt.b.stranded.leases = scen.inputs$alt_B_stranded_leases,
                                                                                                                                                                    tr.cost.alt.c.row = scen.inputs$tr_cost_C_row, alt.d.TLnorth = scen.inputs$alt_D_north,
                                                                                                                                                                    tr.cost.alt = scen.inputs$tr_cost_alt, road.res = scen.inputs$road_res,
                                                                                                                                                                    rd.exist = gen.inputs$rd_exist, tr.cost.alt.b.stranded = scen.inputs$tr_cost_B_stranded,
                                                                                                                                                                    maxd2sat = maxd2sat, mind2sat = mind2sat, n.sat = n.sat, n.iter = n.iter, wd.loc = wd.loc,
                                                                                                                                                                    path.out = path.out, proj.info = proj.info, scenario = scenario, debug.out = debug.out, z = z,
                                                                                                                                                                    i = i)

                                                                                                                                      ## Add the Willow CPF to the cpf.sf object
                                                                                                                                      cpf.sf <- rbind(gen.inputs$willow_cpf[,c("cpf", "geometry")], cpf.iter$cpf_sf)

                                                                                                                                      ## Save off the infrastructure data
                                                                                                                                      utils::write.csv(cpf.iter$cpf_df, file=paste(wd.loc, "/", path.out, "/CPF_data_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), ".csv", sep=""), row.names=FALSE)
                                                                                                                                      utils::write.csv(satrd.iter$sat_df, file=paste(wd.loc, "/", path.out, "/Satellite_pad_data_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), ".csv", sep=""), row.names=FALSE)
                                                                                                                                      utils::write.csv(satrd.iter$rd_df, file=paste(wd.loc, "/", path.out, "/Road_data_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), ".csv", sep=""), row.names=FALSE)
                                                                                                                                      if(shp.out){
                                                                                                                                        sf::st_write(satrd.iter$rd_sf, dsn=paste(wd.loc, path.out, sep="/"), layer=paste("Road_shapefile_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), sep=""), driver="ESRI Shapefile")
                                                                                                                                        sf::st_write(cpf.sf, dsn=paste(wd.loc, path.out, sep="/"), layer=paste("CPF_shapefile_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), sep=""), driver="ESRI Shapefile")
                                                                                                                                        sf::st_write(satrd.iter$sat_sf, dsn=paste(wd.loc, path.out, sep="/"), layer=paste("Satellite_shapefile_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), sep=""), driver="ESRI Shapefile")
                                                                                                                                      }
                                                                                                                                      cpf.df <- cpf.iter$cpf_df
                                                                                                                                      sat.df <- satrd.iter$sat_df
                                                                                                                                      rd.df <- satrd.iter$rd_df
                                                                                                                                      sat.sf <- satrd.iter$sat_sf
                                                                                                                                      rd.sf <- satrd.iter$rd_sf
                                                                                                                                      ## Check whether the output csv folder exists. If not, create it.
                                                                                                                                      dir.create(paste(wd.loc, "/", path.out, "/RData", sep=""), showWarnings=FALSE)
                                                                                                                                      ## Save an .RData file for easier access and reuse
                                                                                                                                      save(cpf.df, sat.df, rd.df, cpf.sf, sat.sf, rd.sf, file=paste(wd.loc, "/", path.out, "/RData/", "DIA_data_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), ".RData", sep=""))
                                                                                                                                    }

                                                                                                                                    ##################################
                                                                                                                                    ## Impacts analyses, if desired ##
                                                                                                                                    ##################################

                                                                                                                                    ## Create an empty output results summary data.frame
                                                                                                                                    out.df <- data.frame(scenario=scenario[z], iteration=i)

                                                                                                                                    ## If any species-specific impact analysis is to be run:
                                                                                                                                    if(!is.null(tch.raster) || !is.null(wah.raster) || !is.null(shorebird.raster) || !is.null(brant.shp)){
                                                                                                                                      ## If infrastructure simulation was not run above, identify and load the
                                                                                                                                      ## previously-generated infrastructure data for the given scenario and iteration.
                                                                                                                                      if(simulate.inf == FALSE){
                                                                                                                                        inf.file <- list.files(path=paste(wd.loc, path.out, "RData", sep="/"), pattern=paste("^DIA_data_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), ".*RData$", sep=""), full.names=TRUE)
                                                                                                                                        load(inf.file)
                                                                                                                                      }

                                                                                                                                      ## Convert the simulated infrastructure locations to a surface disturbance footprint
                                                                                                                                      footprints <- footprint_generation(cpf.sf=cpf.sf, sat.sf=sat.sf, rd.sf=rd.sf, proj.info=proj.info,
                                                                                                                                                                         area.cpf=area.cpf, area.sat=area.sat, road.width=road.width)

                                                                                                                                      ## Calculate infrastructure summary statistics
                                                                                                                                      out.df <- inf_summary(out.df=out.df, cpf.sf=cpf.sf, sat.sf=sat.sf, rd.sf=rd.sf,
                                                                                                                                                            surf.disturb=footprints$surf_disturb, npra.file=npra.file, zoi=zoi, wd.loc=wd.loc, path.in=path.in,
                                                                                                                                                            proj.info=proj.info)

                                                                                                                                      ## Run the desired species-specific impact analyses
                                                                                                                                      if(!is.null(tch.raster) || !is.null(wah.raster)){
                                                                                                                                        ## Unpack SpatRasters, if needed
                                                                                                                                        if(run.in.parallel){
                                                                                                                                          tch.calving <- terra::rast(tch.calving)
                                                                                                                                          wah.calving <- terra::rast(wah.calving)
                                                                                                                                          wah.hq.weight <- terra::rast(wah.hq.weight)
                                                                                                                                        }
                                                                                                                                        ## Run the caribou impact analysis
                                                                                                                                        out.df <- impact_caribou(surf.disturb=footprints$surf_disturb, tch.raster=tch.raster,
                                                                                                                                                                 wah.raster=wah.raster, tch.calving=tch.calving, wah.calving=wah.calving,
                                                                                                                                                                 wah.hq.weight=wah.hq.weight, out.df=out.df, wd.loc=wd.loc, path.out=path.out,
                                                                                                                                                                 proj.info=proj.info, tch.hq=tch.hq, wah.hq=wah.hq, zoi=zoi, caribou.out=caribou.out,
                                                                                                                                                                 scenario=scenario, n.iter=n.iter, z=z, i=i)
                                                                                                                                      }

                                                                                                                                      if(!is.null(shorebird.raster)){
                                                                                                                                        if(run.in.parallel) sb.stack <- terra::rast(sb.stack)
                                                                                                                                        sb.df <- impact_shorebird(sb.stack=sb.stack, rd.sf=rd.sf, sb.df=sb.df)
                                                                                                                                        col.tmp <- ncol(out.df)
                                                                                                                                        out.df <- cbind(out.df, t(sb.df$hq_discount), t(sb.df$hq_remaining))
                                                                                                                                        names(out.df)[(col.tmp+1):ncol(out.df)] <- c(paste(sb.df$sp, "discount", sep="_"), paste(sb.df$sp, "remaining", sep="_"))
                                                                                                                                      }

                                                                                                                                      if(!is.null(brant.shp)){
                                                                                                                                        brant.df <- impact_brant(brant.lakes=brant.lakes, footprints=footprints, land.dist=land.dist,
                                                                                                                                                                 heli.disturb=heli.disturb, proximity.effect=proximity.effect, proj.info=proj.info,
                                                                                                                                                                 brant.out=brant.out, wd.loc=wd.loc, path.out=path.out, scenario=scenario, z=z, i=i,
                                                                                                                                                                 n.iter=n.iter)
                                                                                                                                        out.df <- cbind(out.df, brant.df)
                                                                                                                                      }
                                                                                                                                    }

                                                                                                                                    ## Include the individual run time in the output
                                                                                                                                    indiv.run.time <- Sys.time() - indiv.start
                                                                                                                                    out.df$run_time <- as.numeric(indiv.run.time)
                                                                                                                                    out.df$run_units <- units(indiv.run.time)

                                                                                                                                    ## Write out the results
                                                                                                                                    utils::write.csv(out.df, file=paste(wd.loc, "/", path.out, "/Scenario_impact_data_", scenario[z], "_iter_", formatC(i, width=nchar(n.iter), format="d", flag="0"), "_", Sys.Date(), ".csv", sep=""), row.names=FALSE)

                                                                                                                                    ## Return the run time as the output of the parallelization
                                                                                                                                    return(list(iter=i, time=indiv.run.time))
                                                                                                                                  }

    ## Close the cluster, save summaries, record run timing for the scenario.
    parallel::stopCluster(cl)
    save(iter.run.times, file=paste(wd.loc, "/", path.out, "/RData/Iteration_run_times_Scenario_", scenario[z], "_", Sys.Date(), ".RData", sep=""))
    scenario.run.time <- Sys.time()- scenario.start
    scenario.df <- data.frame(scenario=scenario[z], run_time=as.numeric(scenario.run.time), run_units=units(scenario.run.time))
    utils::write.csv(scenario.df, file=paste(wd.loc, "/", path.out, "/Scenario_run_time_", scenario[z], "_", Sys.Date(), ".csv", sep=""), row.names=FALSE)
    if(!is.null(completion.sound)){
      if(!requireNamespace("beepr", quietly=TRUE)){
        warning("Package \"beepr\" not installed. No sound was made at model completion. To have the model play a sound, please install \"beepr\".", call.=FALSE)
      } else{
        beepr::beep(completion.sound); Sys.sleep(2); beepr::beep(completion.sound); Sys.sleep(2); beepr::beep(completion.sound)
      }
    }
  }
}
