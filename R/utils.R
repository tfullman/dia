#' Ensure matching projections of spatial data
#'
#'\code{projection_alignment} is a helper function that ensures the projection
#' information of a spatial object matches exactly that selected for the overall
#' analysis. It checks whether the projection of an object matches the desired
#' projection information, indicated by \code{proj.info} and, if not, projects
#' the object.
#'
#' @param x Raster* or Spatial* object.
#' @param proj.info Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
#'   common to all spatial objects in the analysis.
#'
#' @return Spatial object of the same class as \code{x}, projected to the
#'   desired project-wide projection.
#'
projection_alignment <- function(x, proj.info){
  proj.tmp <- sp::CRS(SRS_string = proj.info)

  ## Check whether the projection matches the desired projection, if not project it to match
  if(!identical(raster::crs(x), proj.tmp)){
    ## Process for raster objects
    if(class(x) == "RasterLayer"){
      ## Check if there is already a defined projection and handle accordingly
      if(is.na(raster::crs(x))){
        x2 <- x
        raster::crs(x2) <- proj.tmp
      } else{
        y <- raster::projectExtent(x, crs=proj.tmp)
        x2 <- raster::projectRaster(x, y)
        x2 <- raster::mask(x2, x)
      }
    } else{
      ## Process for Spatial* objects
      if(is.na(raster::crs(x))){
        x2 <- x
        sp::proj4string(x2) <- proj.tmp
      } else{
        x2 <- sp::spTransform(x, proj.tmp)
      }
    }
  } else{
    x2 <- x
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
#' @return A Raster* or Spatial* object, projected to the desired project-wide
#'   projection.
#'
load_spatial <- function(x, proj.info, vec = FALSE, wd.loc, path.in){
  if(vec){
    y <- rgdal::readOGR(dsn=paste(wd.loc, path.in, sep="/"), layer=x)
  } else{
    y <- raster::raster(paste(wd.loc, path.in, x, sep="/"))
  }

  ## Ensure this matches the projection of the other spatial data.
  y <- projection_alignment(y, proj.info)

  return(y)
}
