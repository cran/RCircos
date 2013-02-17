\name{RCircos.Cytoband.Data}
\alias{RCircos.Cytoband.Data}
\title{
Chromosome Ideogram Data
}
\description{
This function sorts chromosome ideogram data by chromosome and cytoband, calculate band length (number of base pairs), number of points needed to cover the band, relavent locations of each band on the circular layout, and assign colors to chromosome highlights and each band.
}
\usage{
RCircos.Cytoband.Data(cytoband, chr.exclude = NULL, RCircos.Par)
}
\arguments{
  \item{cytoband}{
A data frame returned from RCircos.Cytoband.Data(cyto.info, chr.exclude, RCircos.Par) and containing plot data for chromosome ideogram .
}
  \item{chr.exclude}{
Character vector containing of chromosome name(s) that should be excluded from data frame to be returned, or NULL if all chromosomes are included.
}
  \item{RCircos.Par}{
A list returned from RCircos.Initialize.Parameters() and containing graphic parameters specific for Circos plot.
}
}
\value{
 \item{Chromosome}{chromosome name with prefix of "chr"}
 \item{ChromStar}{start position of the band}
 \item{ChromEnd}{end position of the band}
 \item{Band}{band name}
 \item{Stain}{intensity of Giemsa staining}
 \item{Length}{band length in base pairs}
 \item{Unit}{total number of points (chromosome units) to draw for the band}
 \item{Location}{index of last points of the band on the circular line}
 \item{BandColor}{RGB color of the band shown on the Circos image}
 \item{ChrColor}{highlight color of the chromosome}
}
\author{
Hongen Zhang
}
\keyword{methods}
