\name{RCircos.Set.Plot.Area}
\alias{RCircos.Set.Plot.Area}
\title{
Setup Plot Area for RCircos Plot
}
\description{
Open an new window for plotting. RCircos core components must be initialized before call this function. All four margins are set to 0.25 and window size are set to plot.radius in RCircos parameters. This function is optional and could be replaced by call par(), plot.new(), and plot.window() if users want control the plot window by themselves.
}
\usage{
RCircos.Set.Plot.Area()
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
      chr.exclude, tracks.inside, tracks.outside); }


# Set plot area, i.e., initialize graphic device 
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

\dontrun{RCircos.Set.Plot.Area(); }
}
\keyword{ methods }
