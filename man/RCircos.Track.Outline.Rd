\name{RCircos.Track.Outline}
\alias{RCircos.Track.Outline}
\title{
Draw Outline of Data Plot Track
}
\description{
Draw outline with subtracks for one data plot track.
}
\usage{
RCircos.Track.Outline(cyto.band, base.positions, out.pos, in.pos, num.subtrack, RCircos.Par)
}
\arguments{
  \item{cyto.band}{
A data frame containing plot data for chromosome ideogram returned from RCircos.Cytoband.Data(cyto.info, chr.exclude, RCircos.Par).
}
  \item{base.positions}{
A data frame returned from RCircos.Base.Plot.Positions(cyto.band, RCircos.Par) containing x- and y-coordinates for a circular line as well as degrees of text rotation at each position.
}
  \item{out.pos}{
Float number, the outside position of a data track.
}
  \item{in.pos}{
Float number, the inside position of a data track.
}
  \item{num.subtrack}{
Integer, total number of subtrack in a data plot track.
}
  \item{RCircos.Par}{
A list returned from RCircos.Initialize.Parameters() containing graphics paramter specific for Circos plot.
}
}
\author{
Hongen Zhang
}
\keyword{methods}
