\name{RCircos.Track.Outline}
\alias{RCircos.Track.Outline}
\title{
Draw Outline for A Data Plot Track
}
\description{
Draw outline with subtracks for one data plot track. RCircos core components and graphic device must be initialized first. This function is mainly for internal use.
}
\usage{
RCircos.Track.Outline(out.pos, in.pos, num.layers)
}
\arguments{
  \item{out.pos}{
Float number, the outside position of a data track.
}
  \item{in.pos}{
Float number, the inside position of a data track.
}
  \item{num.layers}{
Integer, total number of subtrack in a data plot track.
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

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
chr.exclude <- c("chrX", "chrY");
tracks.inside <- 5;
tracks.outside <- 0;
\dontrun{RCircos.Set.Core.Components(cyto.info,  
            chr.exclude, num.inside, num.outside); }



# Initialize graphic device and plot chromosome  
# ideogram
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

\dontrun{RCircos.Set.Plot.Area(); }
\dontrun{RCircos.Chromosome.Ideogram.Plot(); }


# Plot one track outline
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

side <- "in";
track.num <- 1;
locations <- RCircos.Track.Positions(side, track.num);  
out.pos <- locations[1];
in.pos <- locations[2];
subtracks <- 5;
\dontrun{RCircos.Track.Outline(out.pos, in.pos,  
                      subtracks); }
}
\keyword{ methods }
