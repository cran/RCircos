\name{RCircos.Gene.Name.Plot}
\alias{RCircos.Gene.Name.Plot}
\title{
Draw Gene Names on a Data Track
}
\description{
Label gene names along chromosome ideogram or a data track. RCircos core components and graphic device must be initialized before drawing.
}
\usage{
RCircos.Gene.Name.Plot(gene.data, name.col, track.num, side)
}
\arguments{
  \item{gene.data}{
A data frame. The first three columns should be chromosome name, start position, end position.
}
  \item{name.col}{
Integer, representing the ordinal number of the column in input data that contains gene names.
}
  \item{track.num}{
Integer, representing the ordinal number of the plot track where the gene names are plotted.
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
data(RCircos.Gene.Label.Data);


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
# and gene names  
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

\dontrun{RCircos.Set.Plot.Area(); }
\dontrun{RCircos.Chromosome.Ideogram.Plot(); }

name.col <- 4; 
track.num <- 2;
side <- "in";
gene.data <- RCircos.Gene.Label.Data;
\dontrun{RCircos.Gene.Name.Plot(gene.data, name.col,  
                                track.num, side);}
}
\keyword{ methods }
