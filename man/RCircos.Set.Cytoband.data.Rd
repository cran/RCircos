\name{RCircos.Set.Cytoband.data}
\alias{RCircos.Set.Cytoband.data}
\title{
Setup Chromosome Cytoband Data for RCircos Plot
}
\description{
Calculate band length (number of base pairs), number of chromosome units (plot points needed to cover the band), relavent locations of each band on the circular layout, and assign colors to chromosome highlights and each band. This function is for internal use.
}
\usage{
RCircos.Set.Cytoband.data(cyto.band.info)
}
\arguments{
  \item{cyto.band.info}{
A data frame with chromosome ideogram data and columns are ordered by chromosome names, start position, end position, band names, and stain intensity for each chromosome band.
}
}
\author{
Hongen Zhang
}
\keyword{ methods }
