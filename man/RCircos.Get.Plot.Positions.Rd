\name{RCircos.Get.Plot.Positions}
\alias{RCircos.Get.Plot.Positions}
\title{
Get RCircos Plot Position Object
}
\description{
Get the data frame of base plot positions which are used to derive all other positions of RCircos plot. RCircos core components must be initialized first. This function is mainly for internal use.
}
\usage{
RCircos.Get.Plot.Positions()
}
\value{
A data frame with three columns for x, and y coordinates of a circular line with radius of 1 as well as the degrees used for text rotation at each point of the positions. 
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


# Set RCircso core components
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
chr.exclude <- c("chrX", "chrY");
num.inside <- 5;
num.outside <- 0;
\dontrun{RCircos.Set.Core.Components(cyto.info,  
       chr.exclude, tracks.inside, tracks.outside);}


# Get RCircos base plot positions
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

\dontrun{positions <- RCircos.Get.Plot.Positions();}
\dontrun{head(positions);}
}
\keyword{ methods }
