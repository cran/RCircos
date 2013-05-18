\name{RCircos.Get.Gene.Label.Locations}
\alias{RCircos.Get.Gene.Label.Locations}
\title{
Calculate New Position for Gene Labels
}
\description{
This function calculates new plot position from genomic positions for gene laabeling. In case there are too many genes in a genomic interval, the lables may become overlapped so that correct numbers and new plot positions are often needed to make the labels readable. This function will remove extra label data points if they cannot fit in the given genomic interval.

This function is for internal use only.
}
\usage{
RCircos.Get.Gene.Label.Locations(genomic.data)
}
\arguments{
  \item{genomic.data}{
A data frame. The first three columns should be chromosome name, start position, end position.
}
}
\value{
A data frame containing same contents as the input genomic.data plus one column with RCircos plot positions.
}
\author{
Hongen Zhang
}
\keyword{ methods }

