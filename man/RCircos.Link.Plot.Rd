\name{RCircos.Link.Plot}
\alias{RCircos.Link.Plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ircos-like Link Line Plot
}
\description{
Draw lines(quandratic Bezier curves) between paired genomic positions. Graphic device must be initialize before drawing. Link lines are always in the most inner of data plot area.
}
\usage{
RCircos.Link.Plot(cyto.band, base.positions, link.data, track.num, by.chromosome = FALSE, RCircos.Par)
}
\arguments{
  \item{cyto.band}{
A data frame containing plot data for chromosome ideogram returned from RCircos.Cytoband.Data(cyto.info, chr.exclude, RCircos.Par).
}
  \item{base.positions}{
A data frame returned from RCircos.Base.Plot.Positions(cyto.band, RCircos.Par) containing x- and y-coordinates for a circular line as well as degrees of text rotation at each position.
}
  \item{link.data}{
A data frame that contains data for links between genomic positions. The first three columns of the data frame are chromosome name, start position, and end position for first genomic position followed by the chromosome name, start, and end position for second genomic position.
}
  \item{track.num}{
Integer, the number of the plot track inside or outside of chromosome ideogram.
}
  \item{by.chromosome}{
Boolean, If true, red color will be used for links between positions of same chromosomes and blue color for links between different chromosomes. Otherwise, different colors will be automatically assigned for each line.
}
  \item{RCircos.Par}{
A list returned from RCircos.Initialize.Parameters() containing graphics paramters specific for Circos plot.
}
}
\author{
Hongen Zhang
}
\keyword{methods}
