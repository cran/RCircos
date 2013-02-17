\name{RCircos.Track.Positions}
\alias{RCircos.Track.Positions}
\title{
Calculate Track Position For A Data Plot Track
}
\description{
Calculate inner and outer plot positions for a data track.
}
\usage{
RCircos.Track.Positions(side = "in", track.num, RCircos.Par)
}
\arguments{
  \item{side}{
Character vector either "in" or "out" representing the position related to chromosome ideogram.
}
  \item{track.num}{
Integer, number of data track to be plotted.
}
  \item{RCircos.Par}{
A list returned from RCircos.Initialize.Parameters() that contains graphics parameters specific for Circos plot.
}
}
\value{
Vector of two float numbers for outer and inner positions of a data track.
}
\author{
Hongen Zhang
}
\keyword{methods}
