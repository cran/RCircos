\name{RCircos.Data.Point}
\alias{RCircos.Data.Point}
\title{
Convert a Genomic Position to an Index of Circos Plot Position
}
\description{
Calculate the index of Circos plot positions (a set of x- and y-coordinates for a circular line) for a genomic position (chromosome name and start  position). This function is for internal use only.
}
\usage{
RCircos.Data.Point(chromosome, start)
}
\arguments{
  \item{chromosome}{
Character vector, a chromosome name with prefix of "chr", e.g., "chr1".
}
  \item{start}{
Integer, start position of a genomic interval (e.g, a gene or a band) on the chromosome.
}
}
\value{
An integer representing the index of RCircos plot positions.
}
\author{
Hongen Zhang
}
\keyword{ methods }