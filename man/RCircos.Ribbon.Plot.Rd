\name{RCircos.Ribbon.Plot}
\alias{RCircos.Ribbon.Plot}
\title{
Draw Ribbon between Two Genomic Regions
}
\description{
Draw various ribbons instead of lines with uniform thickness. The thickness of ribbons are relative to the size of the genomic regions to be linked. 
}
\usage{
RCircos.Ribbon.Plot(ribbon.data, track.num, by.chromosome = FALSE, twist = FALSE)
}
\arguments{
  \item{ribbon.data}{
A data frame with paired genomic positions in each row. The first three columns of each row are for chromosome name, start and end positions of one genomic position followed by three columns for chromosome name, start and end positions of the second genomic position. 
}
  \item{track.num}{
Integer, representing the ordinal number of the plot track where the ribbons will start.
}
  \item{by.chromosome}{
Logic, If true, red color will be used for ribbons between positions of same chromosomes and blue color for ribbons between different chromosomes. Otherwise, different colors will be automatically assigned for each ribbon.
}
  \item{twist}{
Logic, TRUE for the ribbons with twisted shape and FALSE for no twist.
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
data(RCircos.Link.Data);
data(RCircos.Ribbon.Data);


# Set RCircso core components
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
chr.exclude <- c("chrX", "chrY");
num.inside <- 5;
num.outside <- 0;
\dontrun{RCircos.Set.Core.Components(cyto.info,  
            chr.exclude, num.inside, num.outside); }


# Initialize graphic device, plot chromosome ideogram
# and link lines 
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

\dontrun{RCircos.Set.Plot.Area(); }
\dontrun{RCircos.Chromosome.Ideogram.Plot(); }

track.num <- 4;
\dontrun{RCircos.Link.Plot(RCircos.Link.Data,  
                    track.num, TRUE); }
\dontrun{RCircos.Ribbon.Plot(RCircos.Ribbon.Data,  
                    track.num, FALSE, FALSE); }
}
\keyword{ methods }
