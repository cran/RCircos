\name{RCircos.ScatterPlot}
\alias{RCircos.ScatterPlot}
\title{
Circos-like Scatter Plot
}
\description{
Draw one track of scatter plot. Graphics device must be initialized before drawing.
}
\usage{
RCircos.ScatterPlot(cyto.band, base.positions, scatter.data, data.col, track.num, side, by.fold = 0, RCircos.Par)
}
\arguments{
  \item{cyto.band}{
  A data frame containing plot data for chromosoem ideogram returned from RCircos.Cytoband.Data(cyto.info, chr.exclude, RCircos.Par).
}
  \item{base.positions}{
  A data frame returned from RCircos.Base.Plot.Positions(cyto.band, RCircos.Par) containing x- and y-coordinates for a circular line as well as degrees of text rotation at each position.
}
  \item{scatter.data}{
  A data frame, converted from genomic data by RCircos.Get.Plot.Data(genomic.data, cyto.band), containing genomic positions, scatterplot data values, and Circos plot positions.
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
  \item{by.fold}{
  Zero or a positive float number, if it is greater than zero, any data point with a value greater or equal to the by.fold will be plot as red color and data points with values less or equal to negative by.fold will be plot as blue colr. Otherwise, black color is used.
}
  \item{RCircos.Par}{
  A list returned from RCircos.Initialize.Parameters() containing graphics paramters specific for Circos plot.
}
}
\author{
Hongen Zhang
}
\keyword{methods}
