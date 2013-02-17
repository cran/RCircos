\name{RCircos.Get.Label.Locations}
\alias{RCircos.Get.Label.Locations}
\title{
Calculate New Position for Gene or Data Point Labels
}
\description{
This function is used to derive a new plot position from genomic position for gene or data point labeling. In case there are too many data points in a genomic interval, text and lines for the lables may become overlapped so that correct numbers and new plot positions are often needed to make the labels readable. This function will remove extra label data points if they cannot fit in the given genomic interval.
}
\usage{
RCircos.Get.Label.Locations(cyto.band, genomic.data, label.type = c("text", "line", "point"), RCircos.Par)
}
\arguments{
  \item{cyto.band}{
A data frame containing plot data for chromosome ideogram returned from RCircos.Cytoband.Data(cyto.info, chr.exclude, RCircos.Par).
}
  \item{genomic.data}{
A data frame containing genomic and Circos plot position data returned from RCircos.Get.Plot.Data(genomic.data, cyto.band).
}
  \item{label.type}{
Character vector, either "text", "line", or "point".
}
  \item{RCircos.Par}{
A list returned from RCircos.Initialize.Parameters() and containing graphics parameters specific for Circos plot.
}
}
\value{
A data frame containing same contents as the genomic.data plus one column with new Circos plot positions.
}
\author{
Hongen Zhang
}
\keyword{methods}
