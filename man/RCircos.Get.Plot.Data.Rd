\name{RCircos.Get.Plot.Data}
\alias{RCircos.Get.Plot.Data}
\title{
Convert Genomic Data to Circos Plot Data
}
\description{
This function calculates the Circos plot position from genomic position for a dataset. The Circos plot position is the row index of matrix that contains x- and y-coordinates of a circular line. A new column for Circos plot positions is appended to the input genomic.data and returned.
}
\usage{
RCircos.Get.Plot.Data(genomic.data, cyto.band)
}
\arguments{
  \item{genomic.data}{
A data frame with the first three columns as chromosome name, start and end position followed by data values for one or more samples.
}
  \item{cyto.band}{
A data frame containing Circos plot data for chromosome ideogram returned from RCircos.Cytoband.Data(cyto.info, chr.exclude, RCircos.Par).
}
}
\value{
A data frame with same contents as the genomic.data plus a new column with Circos plot position.
}
\author{
Hongen Zhang
}
\keyword{methods}
