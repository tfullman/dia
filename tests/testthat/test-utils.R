context("utils functions")

test_that("projection alignment works for multiple data formats and initial projections", {
  proj.info <- "EPSG:6393"
  vec.geom <- st_sfc(st_polygon(list(rbind(c(0,2150000), c(10000,2150000),
    c(30000,2170000), c(20000,2190000), c(10000,2190000), c(0,2150000)),
    rbind(c(10000,2160000), c(10000,2170000), c(20000,2170000), c(10000,2160000)))))
  vec.nocrs <- st_sf(geometry=vec.geom)
  vec.akalb <- vec.nocrs
  st_crs(vec.akalb) <- proj.info
  vec.wgs84 <- st_transform(vec.akalb, "EPSG:4326")
  ras.nocrs <- terra::rast(vals=runif(n=126266), nrows=311, ncols=406, xmin=-36267.19, xmax=12452.81, ymin=2231350, ymax=2268670)
  ras.akalb <- ras.nocrs
  terra::crs(ras.akalb) <- proj.info
  ras.wgs84 <- terra::project(x=ras.akalb, y="EPSG:4326")

  expect_equal(st_crs(projection_alignment(x=vec.nocrs, proj.info=proj.info)), st_crs(proj.info))
  expect_equal(st_crs(projection_alignment(x=vec.akalb, proj.info=proj.info)), st_crs(proj.info))
  expect_equal(st_crs(projection_alignment(x=vec.wgs84, proj.info=proj.info)), st_crs(proj.info))
  expect_equal(st_crs(projection_alignment(x=ras.nocrs, proj.info=proj.info)), st_crs(proj.info))
  expect_equal(st_crs(projection_alignment(x=ras.akalb, proj.info=proj.info)), st_crs(proj.info))
  expect_equal(st_crs(projection_alignment(x=ras.wgs84, proj.info=proj.info)), st_crs(proj.info))
})

