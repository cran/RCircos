\name{RCircos.Chromosome.Ideogram}
\alias{RCircos.Chromosome.Ideogram}
\title{
Draw Chromosome Ideogram
}
\description{
Draw a circular chromosome ideogram with all or subset of chromosomes, chromosome names, and chromosome highlights. Graphic device must be initialized before drawing.
}
\usage{
RCircos.Chromosome.Ideogram(cyto.band, base.positions, RCircos.Par)
}
\arguments{
  \item{cyto.band}{
A data frame returned from RCircos.Cytoband.Data(cyto.info, chr.exclude, RCircos.Par) and containing plot data for chromosome ideogram .
}
  \item{base.positions}{
A data frame returned from RCircos.Base.Plot.Positions(cyto.band, RCircos.Par) and containing x- and y-coordinates of a circlular line as well as degrees of text rotation at each point.
}
  \item{RCircos.Par}{
A list returned from RCircos.Initialize.Parameters() containing graphics parameters specific for Circos plot.
}
}
\author{
Hongen Zhang
}
\keyword{methods}
