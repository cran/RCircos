\name{RCircos.Histogram}
\alias{RCircos.Histogram}
\title{
Circos-like Histogram Plot
}
\description{
Draw one track of histogram. Graphic device must be initialized before drawing and the first three columns of input hist.data must be chromosome name, start position, end position.
}
\usage{
RCircos.Histogram(cyto.band, base.positions, hist.data, data.col, track.num, side, RCircos.Par)
}
\arguments{
  \item{cyto.band}{
A data frame containing plot data for chromosoem ideogram returned from RCircos.Cytoband.Data(cyto.info, chr.exclude, RCircos.Par).
}
  \item{base.positions}{
A data frame returned from RCircos.Base.Plot.Positions(cyto.band, RCircos.Par) containing x- and y-coordinates for a circular line as well as degrees of text rotation at each position.
}
  \item{hist.data}{
A data frame, converted from genomic data by RCircos.Get.Plot.Data(genomic.data, cyto.band), containing genomic positions, genomic data values, and Circos plot positions.
}
  \item{data.col}{
Integer, the number of column for the data to be plotted.
}
  \item{track.num}{
Integer, the number of the plot track inside or outside of chromosome ideogram.
}
  \item{side}{
Character vector, either "in" or "out" representing the position related to chromosome ideogram.
}
  \item{RCircos.Par}{
A list returned from RCircos.Initialize.Parameters() containing graphics parameters specific for Circos plot.
}
}
\author{
Hongen Zhang
}
\keyword{methods}
