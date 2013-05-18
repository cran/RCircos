\name{RCircos.Validate.Cyto.Info}
\alias{RCircos.Validate.Cyto.Info}
\title{
Validate Chromosome Ideogram Data
}
\description{
Validate chromosome ideogram information including of chromosome names, order of chromosomes, chromosome start and chromosome end positions for each chromosome band. This function is mainly for internal use.
}
\usage{
RCircos.Validate.Cyto.Info(cyto.info, chr.exclude)
}
\arguments{
  \item{cyto.info}{
A data frame with chromosome ideogram data includes chromosome names, start position, end position, band names, and stain intensity for each chromosome band.
}
  \item{chr.exclude}{
Character vector of chromosome names to be excluded from plot, e.g., chr.exclude <- c("chrX", "chrY");
}
}
\value{
Data frame with validated chromosome ideogram information.
}
\author{
Hongen Zhang
}
\keyword{ methods }
