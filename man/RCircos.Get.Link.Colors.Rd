\name{RCircos.Get.Link.Colors}
\alias{RCircos.Get.Link.Colors}
\title{
Set up Colors for Link Lines or Ribbons
}
\description{
Assign colors for link lines or Ribbons. The defult colors are rainbow with lenght of total link lines or ribbons. If by.chromosome is set to TRUE, red color will be used for links between positions of same chromosomes and blue color for links between different chromosomes. Users can also define their own colors by appending a column with color names to link data and give the column name as "PlotColor". 

This function is for internal use only.
}
\usage{
RCircos.Get.Link.Colors(link.data, by.chromosome)
}
\arguments{
  \item{link.data}{
A data frame with paired genomic positions in each row. The first three columns are for chromosome names, start and end position of one genomic position followed by three columns for chromosome names, start and end position of the second genomic position. 
}
  \item{by.chromosome}{
Logic. If true, red color will be used for links between positions of same chromosomes and blue color for links between different chromosomes. Otherwise, user defined or random colors (rainbow) will be automatically assigned for each line.
}
}
\value{
A vector of color names with length same as the rows of input data.
}
\author{
Hongen Zhang
}
\keyword{ methods }
