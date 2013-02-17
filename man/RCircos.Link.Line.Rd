\name{RCircos.Link.Line}
\alias{RCircos.Link.Line}
\title{
Draw A Circos Link Line
}
\description{
Calculate and draw a quandratic Bezire curve between two points with (0,0) as control points. Graphic device must be initialized before drawing.
}
\usage{
RCircos.Link.Line(point.locations, point.one, point.two, color)
}
\arguments{
  \item{point.locations}{
A data frame returned from RCircos.Base.Plot.Positions(cyto.band, RCircos.Par) containing x- and y-coordinates for a circular line as well as degrees of text rotation at each position.
}
  \item{point.one}{
The first point (start point) of the Bezier curve.
}
  \item{point.two}{
The last point (end point) of the Bezier curve.
}
  \item{color}{
Color for the Bezier curve.
}
}
\author{
Hongen Zhang
}
\keyword{methods}
