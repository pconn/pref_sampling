\docType{data}
\name{Bearded_effort}
\alias{Bearded_effort}
\title{Bearded seal counts and effort information for aerial surveys conducted over the Bering Sea, April-May, 2012}
\description{
 The list \code{Bearded_effort} list object includes 5 objects:
  \itemize{
    \item{\code{Mapping}} A matrix giving the cell (first column) and day (second column) that was surveyed.  The cell number corresponds to the cells present in \code{Data$Grid$y2012} objects (see \code{AlaskaBeringData} documentation).  In the day column, the value '1' corresponds to April 10.  
    \item{\code{Count_data}} For each cell surveyed in \code{Bearded_effort$Mapping}, the corresponding row of \code{Bearded_effort$Count_data} gives the grid cell, the fraction of the grid cell that was surveyed, and the bearded seal count.
    \item{\code{Area.hab}} A vector providing the relative size of each grid cell in \code{Data$Grid$y2012[[1]]}.  Cell sizes can be less than 1.0 if the cell includes land.
    \item{\code{Area.trans}} For each surveyed cell, provides the fraction of area that was surveyed
    \item{\code{DayHour}} For each surveyed cell, provides the day and hour when surveying occurred
  }
}
\author{
  Paul Conn \email{paul.conn@noaa.gov}
}
\keyword{data}

