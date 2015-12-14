\name{RCircos.Add.Ideogram.Tick}
\alias{RCircos.Add.Ideogram.Tick}
\title{
Add Ticks onto Chromosome Ideogram
}
\description{
This function will draw ticks along the outside of chromosome highlight lines. In most cases, ticks are not needed so will be added when necessary.
}
\usage{
RCircos.Add.Ideogram.Tick(tick.interval)
}
\arguments{
  \item{tick.interval}{
Non-negative integer, interval in MB base pairs for the distance between two ticks.
}
}
\value{
None
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
\dontrun{RCircos.Add.Ideogram.Tick(5);}
}
\keyword{ methods }
