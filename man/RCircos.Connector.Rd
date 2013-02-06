\name{RCircos.Connector}
\alias{RCircos.Connector}
\title{
Circos-like Connector Plot
}
\description{
Draw a set of connectors between two tracks (three lines for each connector). Graphics device must be initialized before drawing.
}
\usage{
RCircos.Connector(cyto.band, base.positions, connect.data, track.num, side, RCircos.Par)
}
\arguments{
  \item{cyto.band}{
A data frame returned from RCircos.Cytoband.Data(cyto.info, chr.exclude, RCircos.Par) and containing plot data for chromosome ideogram .
}
  \item{base.positions}{
A data frame returned from RCircos.Base.Plot.Positions(cyto.band, RCircos.Par) and containing x- an dy-coordinates of a circular lines as well as degrees of text rotation at each point.
}
  \item{connect.data}{
A data frame containing index of plot points for outer (e.g., chromosome ideogram) and inner track (e.g., gene names).
}
  \item{track.num}{
Integer, the number of the plot track inside or outside of chromosome ideogram.
}
  \item{side}{
Character vector, either "in" or "out" representing the position related to chromosome ideogram.
}
  \item{RCircos.Par}{
A list returned from RCircos.Initialize.Parameters() containing graphics paramters specific for Circos plot.
}
}
\author{
Hongen Zhang
}
\keyword{methods}
