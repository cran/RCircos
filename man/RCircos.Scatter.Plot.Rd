\name{RCircos.Scatter.Plot}
\alias{RCircos.Scatter.Plot}
\title{
Circos-like Scatter Plot
}
\description{
Draw one track of scatter plot. RCircos core components and graphics device must be initialized before drawing.
}
\usage{
RCircos.Scatter.Plot(scatter.data, data.col, track.num,   
                                     side, by.fold = 0)
}
\arguments{
  \item{scatter.data}{
A data frame with the first three columns as chromosome name, start position and end position of each data point followed by data values for one or more samples.
}
  \item{data.col}{
Integer, representing the ordinal number of the column in input dataset that contains the data to be plotted.
}
  \item{track.num}{
Integer, representing the ordinal number of the plot track where the data will be plotted.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram.
}
  \item{by.fold}{
Zero or a positive float number, if it is greater than zero, any data point with a value greater or equal to the by.fold will be plot as red color and data points with values less or equal to negative by.fold will be plot as blue colr. Otherwise, black color is used.
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


# Load human chromosome ideogram
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

data(UCSC.HG19.Human.CytoBandIdeogram);
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


# Set RCircso core components
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

chr.exclude <- c("chrX", "chrY");
tracks.inside <- 5;
tracks.outside <- 0;
\dontrun{RCircos.Set.Core.Components(cyto.info,  
           chr.exclude, num.inside, num.outside); }


# Initialize graphic device, plot chromosome ideogram 
# and scatters
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

\dontrun{RCircos.Set.Plot.Area(); }
\dontrun{RCircos.Chromosome.Ideogram.Plot(); }

data(RCircos.Scatter.Data);
data.col <- 5; 
track.num <- 4;
side <- "in";
by.fold <- 1;
\dontrun{RCircos.Scatter.Plot(RCircos.Scatter.Data,  
              data.col, track.num, side, by.fold);}
}
\keyword{ methods }

