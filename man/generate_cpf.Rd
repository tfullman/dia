% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_inf.R
\name{generate_cpf}
\alias{generate_cpf}
\title{Simulate CPF locations}
\usage{
generate_cpf(
  oil.im,
  n.cpf = c(3, 7),
  d2cpf = 35000,
  debug.out = FALSE,
  wd.loc,
  path.out,
  n.iter = 100,
  proj.info = NULL,
  scenario = NULL,
  z = NULL,
  i = NULL
)
}
\arguments{
\item{oil.im}{Object of type "\code{im}" representing the relative likelihood
of pad development as a function of relative undiscovered oil, development
restrictions, and the presence of waterbodies. Created by \code{\link{prep_scenario_inputs}}.}

\item{n.cpf}{Indicator of the desired number of CPFs to be generated. This
can either be a single integer value, indicating a fixed number of CPFs, or
a range \code{c(min,max)}, from which a random value will be drawn using a
uniform distribution. Defaults to \code{c(3,7)} to reflect a reasonable
range of potential variability, depending on the number and distribution of
undiscovered oil deposits.}

\item{d2cpf}{Numeric value indicating the minimum allowed distance (m) between
two central processing facilities (CPFs). Defaults to 35,000 m.}

\item{debug.out}{Logical indicator of whether intermediate infrastructure .csv
files should be written out for the purpose of debugging code. Defaults to
\code{FALSE}.}

\item{wd.loc}{Character string identifying the base directory containing both
input and output data folders.}

\item{path.out}{Character string identifying the relative path from the base
directory (\code{wd.loc}) to where output data will be saved. Defaults to
an "Output_Data" folder within the base folder.}

\item{n.iter}{Integer indicating the desired number of iterations to be run
for each scenario. Defaults to 100.}

\item{proj.info}{Desired projection string in EPSG code format (\code{"EPSG:XXXX"}),
common to all spatial objects in the analysis.}

\item{scenario}{Vector of character strings identifying the scenarios being run.}

\item{z}{Iterator value used by \code{\link{dia}} to specify the development
scenario being analyzed.}

\item{i}{Iterator value used by \code{\link{dia}} to run development simulation
and impacts analyses in parallel, indicating the specific iteration being run.}
}
\value{
List containing a data.frame of simulated CPF locations (\code{cpf_df}),
SpatialPoints object of CPF locations (\code{cpf_sp}), and integer indicating
the number of simulated CPFs in the current iteration (\code{tmp_n_cpf}).
}
\description{
\code{generate_cpf} simulates the location of central processing facilities
(CPFs) as a function of relative estimated undiscovered oil, development
restrictions, and the presence of waterbodies. CPFs are "the operational
center for long-term production" of an oilfield (BLM 2019a, p. B-6) so
infrastructure simulation starts with CPFs. CPF locations approximate the
locations of newly discovered oil accumulations suitable for development
(Wilson et al. 2013).
}
\references{
BLM [Bureau of Land Management] 2019a. National Petroleum Reserve
in Alaska Draft Integrated Activity Plan and Environmental Impact Statement.
Bureau of Land Management, U.S. Department of the Interior, Anchorage, AK,
USA.

Wilson RR, Liebezeit JR, Loya WM. 2013. Accounting for uncertainty in oil and
gas development impacts to wildlife in Alaska. Conservation Letters
6:350-358.
}
