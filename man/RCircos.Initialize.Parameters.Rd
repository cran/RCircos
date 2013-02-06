\name{RCircos.Initialize.Parameters}
\alias{RCircos.Initialize.Parameters}
\title{
Initialize Graphics Parameters for Circos Plot
}
\description{
Set up default values for parameters of chromosome ideogram and data plots. All parameters are held in a list and returned for outside use.
}
\usage{
RCircos.Initialize.Parameters()
}
\value{
  \item{base.per.unit }{
Integer, number of base pairs covered by one plot point, default: 3000.
}
  \item{ chrom.paddings }{
Integer, padding length between two chromosomes in chromosome units, default: 3000.
}
  \item{radius.len }{
Float number, radius of circular line to calculate all plot positions, default: 1.0.
}
  \item{chr.ideog.pos }{
Float number, location of chromosome ideogram, default: radius.len + 0.1.
}
  \item{highlight.pos }{
Float number, location of chromosome highlight, default: radius.len + 0.2 .
}
  \item{chr.name.pos }{
Float number, location of chromosome name labels, default: radius.len + 0.3.
}
  \item{plot.radius }{
Float, the radius of maximum plot area, default: radius.len + 0.5.
}
  \item{ track.in.start }{
Float number, location of first track inside of chromosome ideogram, default: chr.ideog.pos - 0.05.
}
  \item{ track.out.start }{
Float number, location of first track outside of chromosome ideigram, default: radius.len + 0.4.
}
  \item{chrom.width }{
Float number, width of chromosomes on ideogram, default: 0.08.
}
  \item{track.padding }{
Float number, padding width beteen two plot tracks, default: 0.02.
}
  \item{ track.height }{
Float number, height of data plot track, default: 0.1.
}
  \item{highlight.width }{
Integer, line type of chromosome highlight, default: 1.
}
  \item{hist.width }{
Integer, width of histogram columns, in chromosome units, default: 1000.
}
  \item{text.size }{
Float number, cex value of gene names for circular label, default: 0.4.
}
  \item{heatmap.width }{
Integer, heatmap cell width, in chromosome units, default: 100.
}
  \item{point.type }{
Character type for point plot, default:  ".".
}
  \item{point.size }{
Integer, character size for point plot, default: 1.
}
  \item{max.layers }{
Integer, maximum number of layers for tile plot,  default: 5.
}
}
\author{
Hongen Zhang
}
\keyword{methods}
