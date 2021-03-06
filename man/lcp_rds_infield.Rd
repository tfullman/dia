% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_inf.R
\name{lcp_rds_infield}
\alias{lcp_rds_infield}
\title{Connect satellite pads to their CPF via infield roads}
\usage{
lcp_rds_infield(
  x,
  tr,
  cost.map,
  cpf.df,
  min.dist.exist,
  wd.loc,
  path.out,
  scenario,
  n.iter,
  j,
  z = NULL,
  i = NULL,
  debug.out = FALSE
)
}
\arguments{
\item{x}{data.frame of satellite pad locations, created within \code{\link{generate_sat_rd}}.
This may contain all satellite pad locations, or a subset, depending on the
scenario being analyzed.}

\item{tr}{Scenario-specific \code{TransitionLayer} object, used by the \code{gdistance}
package for creating least cost paths.}

\item{cost.map}{RasterLayer indicating land and water areas, used in gravel
pad simulation. Created by \code{\link{prep_general_inputs}}.}

\item{cpf.df}{data.frame of CPF point locations. Created by \code{\link{generate_cpf}}.}

\item{min.dist.exist}{Proximity order-sorted data.frame of distance from CPFs
to existing infrastructure. Created within \code{\link{generate_sat_rd}}.}

\item{wd.loc}{Character string identifying the base directory containing both
input and output data folders.}

\item{path.out}{Character string identifying the relative path from the base
directory (\code{wd.loc}) to where output data will be saved. Defaults to
an "Output_Data" folder within the base folder.}

\item{scenario}{Vector of character strings identifying the scenarios being run.}

\item{n.iter}{Integer indicating the desired number of iterations to be run
for each scenario. Defaults to 100.}

\item{j}{Iterator used within \code{\link{generate_sat_rd}} to indicate the
current CPF being analyzed.}

\item{z}{Iterator value used by \code{\link{dia}} to specify the development
scenario being analyzed.}

\item{i}{Iterator value used by \code{\link{dia}} to run development simulation
and impacts analyses in parallel, indicating the specific iteration being run.}

\item{debug.out}{Logical indicator of whether intermediate infrastructure .csv
files should be written out for the purpose of debugging code. Defaults to
\code{FALSE}.}
}
\value{
List containing two objects, a temporary data.frame of generated road
coordinates (\code{tmp_rd_df}) and a temporary SpatialLines object depicting
generated roads (\code{tmp_rd_sl}). These are fed into \code{\link{lcp_rds_outfield}}.
}
\description{
\code{lcp_rds_infield} is a helper function that uses least cost paths to
connect satellite pads to their parent central processing facility (CPF) via
infield roads. Used by \code{\link{generate_sat_rd}}.
}
