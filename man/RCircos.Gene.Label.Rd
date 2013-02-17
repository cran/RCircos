\name{RCircos.Gene.Label}
\alias{RCircos.Gene.Label}
\title{
Circos-like Gene Label Plot
}
\description{
Label gene names along chromosome ideogram or a data track.
}
\usage{
RCircos.Gene.Label(base.positions, gene.data, name.col, track.num, side, RCircos.Par)
}
\arguments{
  \item{base.positions}{
A data frame returned from RCircos.Base.Plot.Positions(cyto.band, RCircos.Par) containing x- and y-coordinates for a circular line as well as degrees of text rotation at each position.
}
  \item{gene.data}{
A data frame containing gene annotation data. The first three columns should be chromosome name, start position, and end position. One column is for gene names and the last column is for Circos plot location of each gene.
}
  \item{name.col}{
Integer, the number of column for gene name in the gene.data object.
}
  \item{track.num}{
Integer, the number of plot track inside or outside of chromosome ideogram.
}
\item{side}{
  Character vector, either "in" or "out", representing the position related to chromosome ideogram.
}
  \item{RCircos.Par}{
A list returned from RCircos.Initialize.Parameters() containing graphics parameters specific for Circos plot.
}
}
\author{
Hongen Zhang
}
\keyword{methods}
