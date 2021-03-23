#' Ensure matching projections of spatial data
#'
#'\code{projection_alignment} is a helper function that ensures the projection
#' information of a spatial object matches exactly that selected for the overall
#' analysis. It checks whether the projection of an object matches the desired
#' projection information, indicated by \code{proj.info} and, if not, projects
#' the object.
#'
#' @param x Raster or vector object of type \code{SpatRaster} from the \href{https://cran.r-project.org/web/packages/terra/index.html}{terra}
#'   package or \code{sf} from the \href{https://cran.r-project.org/web/packages/sf/index.html}{sf} package.
#' @param proj.info Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
#'   common to all spatial objects in the analysis.
#'
#' @details This function is intended to address discrepancies in spatial object
#' projection, not to deal with issues of misaligned extent or cell size for
#' SpatRaster objects. \code{projection_alignment} assumes that the user has already
#' done initial preparation to ensure any raster files have the same extent and
#' resolution. For tools to address such data preparation, see the \href{https://cran.r-project.org/web/packages/terra/index.html}{terra}
#' package.
#'
#' @return Spatial object of the same class as \code{x}, projected to the
#'   desired project-wide projection.
#'
projection_alignment <- function(x, proj.info){
  ## Process for raster objects
  if(inherits(x, "SpatRaster")){
    ## Identify the desired projection
    proj.tmp <- terra::crs(proj.info)
    ## Check whether the object's projection matches the desired projection, if not project it to match
    if(!identical(terra::crs(x), proj.tmp)){
      ## Check if there is already a defined projection and handle accordingly
      if(terra::crs(x) == ""){
        x2 <- x
        terra::crs(x2) <- proj.info
      } else{
        x2 <- terra::project(x=x, y=proj.info)
      }
    } else{
      x2 <- x
    }
  }

  ## Process for Spatial* objects
  if(inherits(x, "sf")){
    ## Identify the desired projection
    proj.tmp <- sf::st_crs(proj.info)$wkt
    ## Check whether the projection matches the desired projection, if not project it to match
    if(!identical(sf::st_crs(x)$wkt, proj.tmp)){
      ## Check if there is already a defined projection and handle accordingly
      if(is.na(sf::st_crs(x))){
        x2 <- x
        sf::st_crs(x2) <- proj.info
      } else{
        x2 <- sf::st_transform(x, proj.info)
      }
    }
    else{
      x2 <- x
    }
  }

  return(x2)
}



#' Load spatial data and ensure projection matching
#'
#'\code{load_spatial} is a helper function that loads spatial data (vector or
#' raster) and ensures the projection information matches exactly that selected
#' for the overall analysis.
#'
#' @param x Character string specifying the file name containing the spatial object.
#' @param proj.info Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
#'   common to all spatial objects in the analysis.
#' @param vec Logical indicator of whether the spatial data to be loaded is in
#'   vector format (\code{TRUE}) or raster data (\code{FALSE}; default).
#' @param wd.loc Character string identifying the base directory containing both
#'   input and output data folders.
#' @param path.in Character string identifying the relative path from the base
#'   directory (\code{wd.loc}) to where input data are stored. Defaults to an
#'   Input_Data folder within the base folder.
#'
#' @return A \code{SpatRaster} or \code{sf} object, projected to the desired project-wide
#'   projection.
#'
load_spatial <- function(x, proj.info, vec = FALSE, wd.loc, path.in){
  if(vec){
    y <- sf::st_read(dsn=paste(wd.loc, path.in, sep="/"), layer=x)
  } else{
    y <- terra::rast(paste(wd.loc, path.in, x, sep="/"))
  }

  ## Ensure this matches the projection of the other spatial data.
  y <- projection_alignment(y, proj.info)

  return(y)
}



#' Set raster values to zero
#'
#'\code{raster_to_zero} is a helper function that sets raster values to zero where
#' they are overlapped by a polygon. Pixels with an \code{NA} value in the input raster
#' will remain \code{NA} in the output.
#'
#' @param x \code{SpatRaster} object.
#' @param y Polgon layer of class \code{SpatVector} representing the area within which
#'   the raster values will be set to zero.
#'
#' @return \code{SpatRaster} object with the same information as \code{x}, but with
#'   areas overlapping \code{y} set to zero.
#'
#' @note \code{raster_to_zero} replaces the \code{infrastructure_overlap_to_zero}
#'   function from \code{dia} version 0.1.0. It also simplifies the code for the
#'   \code{infrastructure_exclusion_buffer} function from that version.
raster_to_zero	<- function(x, y){
  x2 <- terra::mask(x=x, mask=y, updatevalue=0, inverse=TRUE, touches=TRUE)
  x2 <- terra::mask(x2, x)
  return(x2)
}
