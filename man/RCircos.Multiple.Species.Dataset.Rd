\name{RCircos.Multiple.Species.Dataset}
\alias{RCircos.Multiple.Species.Dataset}
\title{
Construct A Dataset from Multiple Genomic Data of Different Species
}
\description{
Combine plot datasets from different species as one data frame for RCircos plot with multiple sets of chromosome ideogram.
}
\usage{
RCircos.Multiple.Species.Dataset(data.list, species)
}
\arguments{
  \item{data.list}{
A List contains data frames for genomic data of different species. All data frames in the list must have same columns.
}
  \item{species}{
Character vector represent multiple species. Each one will be used as prefix of chromosome names for relavant species. The order of species must match the order of the species in data.list.
}
}
\value{
A data frame contains all input plot data.
}
\author{
Hongen Zhang
}
\examples{
# Load RCircos libaray
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

library(RCircos);


# Load mouse and rat gene expression datasets 
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

data(RCircos.Mouse.Expr.Data);
data(RCircos.Rat.Expr.Data);
data.list <- list(RCircos.Mouse.Expr.Data,   
                    RCircos.Rat.Expr.Data);
species <- c("M", "R");


# Generate data frome that contains both mouse and
# rat gene expression data
# __________________________________________________
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

\dontrun{new.data <- RCircos.Multiple.Species.Dataset(  
                          data.list, species);}
}
\keyword{ methods }

