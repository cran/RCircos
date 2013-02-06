\name{RCircos.Data.Point}
\alias{RCircos.Data.Point}
\title{
Convert Genomic Position to Index of Circos Plot Position
}
\description{
Calculate the index of Circos plot positions (a set of x- and y-coordinates for a circular line) for a genomic position (chromosome and start position).
}
\usage{
RCircos.Data.Point(cyto.band, chromosome, start)
}
\arguments{
  \item{cyto.band}{
A data frame returned from RCircos.Cytoband.Data(cyto.info, chr.exclude, RCircos.Par) and containing plot data for chromosome ideogram.
}
  \item{chromosome}{
Character vector, a chromosome name with  or without prefix of "chr".
}
  \item{start}{
Integer, start position of a genomic interval (e.g, a gene or a band) on the chromosome.
}
}
\value{
An integer representing the index of Circos plot positions.
}
\author{
Hongen Zhang
}
\keyword{methods}
