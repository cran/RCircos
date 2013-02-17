\name{RCircos.Base.Plot.Positions}
\alias{RCircos.Base.Plot.Positions}
\title{
Base Positions of Points for Circos Plot
}
\description{
Calculate x- and y-coordinates of points for a circular line which serves as the basis of plot locations for all Crircos tracks and data points.
}
\usage{
RCircos.Base.Plot.Positions(cyto.band, RCircos.Par)
}
\arguments{
  \item{cyto.band}{
A data frame returned from RCircos.Cytoband.Data(cyto.info, chr.exclude, RCircos.Par) and containing plot data for chromosome ideogram .
}
  \item{RCircos.Par}{
A list returned from RCircos.Initialize.Parameters() and containing graphicarameters specific for Circos plots.
}
}
\details{
This function calculates x- and y-coordinates of a set of points that forms a circular line with radius of 1. The coordinates start from (0,1) and move on clockwise. Plot positions for chromosome ideogram, data tracks, data points, and text labels are all derived from these coordinates.
}
\value{
A data frame that contains the x- and y-coordinates of each point on a circular line as well as degrees of text rotation at each point.
}
\author{
Hongen Zhang
}
\keyword{methods}
