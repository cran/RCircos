\name{RCircos.Validate.Genomic.Data}
\alias{RCircos.Validate.Genomic.Data}
\title{
Validate Genomic Data for RCircos Plot
}
\description{
Validate input dataset for correct chromosome names, chromosome start, and chromosome end positions. Chromosome names will be converted to character vectors if they are factor variables. This function is for internal use.
}
\usage{
RCircos.Validate.Genomic.Data(genomic.data, plot.type = c("plot", "link"))
}
\arguments{
  \item{genomic.data}{
A data frame with the first three columns as chromosome names, start and end positions. Following columns may be another set of chromosome names, start and end positions then zero or more data columns.
}
  \item{plot.type}{
Character vector either "plot" for dataset with one set of genomic position data or "link" for dataset with paired genomic position data . 
}
}
\value{
The validated data frame.
}
\author{
Hongen Zhang
}
\keyword{ methods }

