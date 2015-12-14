\name{RCircos.Reset.Plot.Positions}
\alias{RCircos.Reset.Plot.Positions}
\title{
Reset RCircos Base Position
}
\description{
Reset user defined/modified base plot positions. This is reserved for later use.
}
\usage{
RCircos.Reset.Plot.Positions(plot.positions)
}
\arguments{
  \item{plot.positions}{
A matrix of three columns for x and y coordinates of set of points along a circular line as well as degrees of text should be rotated at each point.
}
}

\value{
None
}
\author{
Hongen Zhang
}
\examples{
\dontrun{plot.positions <- RCircos.Get.Plot.Positions();
plot.positions[1, 3] <- 0;
RCircos.Reset.Plot.Positions(plot.positions);}
}

\keyword{ methods }
