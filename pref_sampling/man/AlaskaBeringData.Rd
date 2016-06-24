\docType{data}
\name{AlaskaBeringData}
\alias{AlaskaBeringData}
\title{Alaska covariate data for bearded seal analyses}
\description{
  Bering Sea covariate data in spring of 2012 and 2013 for use with bearded seal analyses.  The list \code{Data} includes four objects:
  \itemize{
    \item{\code{Adj}} An adjacency matrix describing the spatial proximity of sampling units (sample units that are queen's move neighbors receive a 1, those that are not neighbors receive a 0)
    \item{\code{Adj2}} An adjacency matrix for a larger spatial extent (not used)
    \item{\code{Grid}} A list vector, containing covariates for 2012 (\code{y2012}) and 2013 (\code{y2013}. Each list element represents a single day - for 2012, day 1 corresponds to April 4. For each day, there is a SpatialPolygonsDataFrame (see \code{sp} package) encoding the spatial position and covariates associated with sample units on each day.  For example, \code{Data$Grid$y2012[[1]]} holds covariate data for day 1 (4 April 2012).  Covariates provided include proportion land cover (\code{land_cover}), standardized distance from sample unit centroid to mainland (\code{dist_mainland}), standardized distance from land (including islands; \code{dist_land}), standardized distance from 1000m shelf break (\code{dist_shelf}), standardized depth (\code{depth}), proportion sea ice coverage (obtained from NSIDC; \code{ice_conc}), standardized distance from 90 percent sea ice contour (\code{dist_contour}), standardized distance from the southern sea ice edge (\code{dist_edge}), and Piatt and Springer ecoregion \code{Ecoregion}.
    \item{\code{Meta}} Provides some basic metadata
  }
}
\author{
  Paul Conn \email{paul.conn@noaa.gov}
}
\keyword{data}

