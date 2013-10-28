\name{RCircos.Clear.Track}
\alias{RCircos.Clear.Track}
\title{
Erase One or More Data Tracks
}
\description{
Erase one data track or all tracks from current track to the center of plot plot area.
}
\usage{
RCircos.Clear.Track(track.num, side, to.center = FALSE)
}
\arguments{
  \item{track.num}{
Integer, representing the ordinal number of the plot track to be erased.
}
  \item{side}{
Character vector, either "in" or "out", representing the position related to chromosome ideogram.
}
  \item{to.center}{
Logic, FALSE for erasing current track only and TRUE for erasing all tracks from current track to the center of plot plot area. 
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
data(RCircos.Link.Data);


# Set RCircso core components
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
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
track.num <- 1;
side <- "in";
by.fold <- 1;
\dontrun{RCircos.Scatter.Plot(RCircos.Scatter.Data,  
              data.col, track.num, side, by.fold);}

track.num <- 3;
\dontrun{RCircos.Link.Plot(RCircos.Link.Data,  
                    track.num, TRUE); }

track.num <- 1;
side <- "in";
to.center <- FALSE;
\dontrun{RCircos.Clear.Track(track.num, side, to.center); }
to.center <- TRUE;
\dontrun{RCircos.Clear.Track(track.num, side, to.center); }

}
\keyword{ methods }

