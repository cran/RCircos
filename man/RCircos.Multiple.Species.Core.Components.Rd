\name{RCircos.Multiple.Species.Core.Components}
\alias{RCircos.Multiple.Species.Core.Components}
\title{
Setup RCircos Core Components with Multiple Species
}
\description{
Setup RCircos core components with chromosome ideograms of multiple species for visulization of similarity and difference of genomic data between different species.
}
\usage{
RCircos.Multiple.Species.Core.Components(cyto.info.list,  
                       species.list, chr.exclude = NULL,  
                       tracks.inside, tracks.outside)
}
\arguments{
  \item{cyto.info.list}{
A List contains data frames for chromosome ideogram data of different species. All data frames in the list must have same columns.
}
  \item{species.list}{
Character vector represent multiple species. Each one will be used as prefix of chromosome names for relavant species. The order of species must match the order of the species in cyto.info.list.
}
  \item{chr.exclude}{
Character vector of chromosome names to be excluded from plot, e.g., chr.exclude <- c("chrX", "chrY");
}
  \item{tracks.inside}{
Integer, total number of data tracks inside of chromosome ideogram.
}
  \item{tracks.outside}{
Integer, total number of data tracks outside of chromosome ideogram.
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


# Load mouse and rat chromosome ideograms and make  
# a multiple species cytoband data
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

data(UCSC.Mouse.GRCm38.CytoBandIdeogram);
data(UCSC.Baylor.3.4.Rat.cytoBandIdeogram);
cyto.list <- list(UCSC.Mouse.GRCm38.CytoBandIdeogram,  
               UCSC.Baylor.3.4.Rat.cytoBandIdeogram);
species <- c("M", "R");


# Set RCircso core components with multiple species  
# cytoband data
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

chr.exclude <- c("chrX", "chrY");
tracks.inside <- 5;
tracks.outside <- 0;
\dontrun{RCircos.Multiple.Species.Core.Components(  
       cyto.list, species, chr.exclude,  
       tracks.inside, tracks.outside); }
}
\keyword{ methods }

