\name{RCircos.Reset.Ideogram.Position}
\alias{RCircos.Reset.Ideogram.Position}
\title{
Reset Graphics Parameters Derived from Radius
}
\description{
Reset parameters derived from radius length such as plot positions for chromosome ideogram, chromosome highlight, chromosome names, inside and outside tracks, and radius of plot area.
}
\usage{
RCircos.Reset.Ideogram.Position(RCircos.Par, radius.len)
}
\arguments{
  \item{RCircos.Par}{
A list returned from RCircos.Initialize.Parameters() containing graphics paramters specific for Circos plot.
}
  \item{radius.len}{
Float number, radius length of the circular line which served as the basis for calculation of Circos plot position.
}
}
\details{
Parameters to be changed and their values: chr.ideog.pos(radius.len + 0.1), highlight.pos(radius.len + 0.2), chr.name.pos(radius.len + 0.3), track.in.start(chr.ideog.pos - 0.05), track.out.start( chr.name.pos + 0.1), and plot.radius (radius.len + 0.5).
}
\value{
A list with updataed values of contents same as the input argument.
}
\author{
Hongen Zhang
}
\keyword{methods}
