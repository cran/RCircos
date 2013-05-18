\name{RCircos.Heatmap.Plot}
\alias{RCircos.Heatmap.Plot}
\title{
Circos-like Heamap Plot
}
\description{
Draw one track of heatmap. RCircos core components and graphic device must be initialized before drawing.
}
\usage{
RCircos.Heatmap.Plot(heatmap.data, data.col, track.num, side)
}
\arguments{
  \item{heatmap.data}{
A data frame with the first three columns as chromosome names, start position and end position of each data point followed by data values for one or more samples.
}
  \item{data.col}{
Integer, representing the ordinal number of the column in input data set that contains the data to be plotted.
}
  \item{track.num}{
Integer, representing the ordinal number of the plot track where the heatmap will be plotted.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram.
}
}
\author{
Hongen Zhang
}
\examples{
# Load RCircos libaray
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

library(RCircos);


# Load human chromosome ideogram and sample data
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

data(UCSC.HG19.Human.CytoBandIdeogram);
data(RCircos.Heatmap.Data);


# Set RCircso core components
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
chr.exclude <- c("chrX", "chrY");
num.inside <- 5;
num.outside <- 0;
\dontrun{RCircos.Set.Core.Components(cyto.info,  
             chr.exclude, num.inside, num.outside);}


# Initialize graphic device, plot chromosome ideogram
# and one heatmap track
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

\dontrun{RCircos.Set.Plot.Area(); }
\dontrun{RCircos.Chromosome.Ideogram.Plot(); }

data.col <- 5
track.num <- 4;
side <- "in";
\dontrun{RCircos.Heatmap.Plot(RCircos.Heatmap.Data,  
                   data.col, track.num, side); }
}
\keyword{ methods }
