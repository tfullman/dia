# dia 0.1.0.9000

* Switch `sp` to `sf` code.
* Switch `raster` to `terra` code.
* Add helper functions and consolidate code to reduce repetition.
  * New functions:
    * `raster_to_zero()`
    * `gen_lcp_rd()`
    * `gen_linkage_rd()`
    * `poly_rotate()`
    * `pt_to_pad()`
* Use `units` and `udunits2` packages to automate unit conversion.
* Adds more user options for `dia`:
  * Specify high-quality habitat threshold
  * Make writing out shapefiles optional
  * Add ability to easily change values for `area.cpf`, `area.sat`, and `road.width`
  * Write out shapefile of affected brant molting lakes, if any
* Switch `beepr` from Imports to Suggests.


# dia 0.1.0

* Released initial version of the `dia` package to accompany Fullman et al. (in press) Ecosphere.
* Added a `README.md` file for new users.
* Added a `NEWS.md` file to track changes to the package.
