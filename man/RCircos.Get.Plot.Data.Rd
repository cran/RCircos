\name{RCircos.Get.Plot.Data}
\alias{RCircos.Get.Plot.Data}
\title{
Convert Genomic Data to RCircos Plot Data
}
\description{
This function calculates the RCircos plot positions for genomic positions of the dataset. RCircos plot position is the row index of matrix that contains x- and y-coordinates of a circular line. A new column for RCircos plot positions is appended to the input genomic.data and returned.

This function is for internal use.
}
\usage{
RCircos.Get.Plot.Data(genomic.data, plot.type)
}
\arguments{
  \item{genomic.data}{
A data frame with the first three columns as chromosome name, start and end position followed by data values for one or more samples.
}
  \item{plot.type}{
Character vector, either "plot" or "link".
}
}
\value{
Updated genomic data (with a new column for RCircos plot positions)
}
\author{
Hongen Zhang
}
\keyword{ methods }
