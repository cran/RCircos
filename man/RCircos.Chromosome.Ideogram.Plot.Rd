\name{RCircos.Chromosome.Ideogram.Plot}
\alias{RCircos.Chromosome.Ideogram.Plot}
\title{
Draw Chromosome Ideogram
}
\description{
Draw a circular chromosome ideogram with all or subset of chromosomes, chromosome names, and chromosome highlights. RCircos core components and graphic device must be initialized before drawing.
}
\usage{
RCircos.Chromosome.Ideogram.Plot()
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
chr.exclude <- NULL;
num.inside <- 5;
num.outside <- 0;
\dontrun{RCircos.Set.Core.Components(cyto.info,  
	chr.exclude, num.inside, num.outside); }


# Initialize graphic device and plot chromosome  
# ideogram
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

\dontrun{RCircos.Set.Plot.Area();}
\dontrun{RCircos.Chromosome.Ideogram.Plot();}
}
\keyword{methods}

