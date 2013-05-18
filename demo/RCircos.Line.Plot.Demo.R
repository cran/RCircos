#
#	This demo draw chromosome ideogram with padding between
#	chromosomes, highlights, chromosome names, and line plot. 
#
#	Usage:
#
#	library(RCircos);
#	demo("RCircos.Line.Plot.Demo");
#
#	========================================================


RCircos.Line.Plot.Demo<-function()
{
	#	Load RCircos library
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
	library(RCircos);


	#	Load human cytoband data and gene expression data
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	data(RCircos.Line.Data);
	data(UCSC.HG19.Human.CytoBandIdeogram);
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


	#	Setup RCircos core components:
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);



	#	Open the graphic device (here a pdf file)
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
	out.file <- "RCircos.Line.Plot.Demo.pdf";
	pdf(file=out.file, height=8, width=8);

	RCircos.Set.Plot.Area();
	title("RCircos Line Plot Demo");


	#	Draw chromosome ideogram
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Draw chromosome ideogram ...\n");

	RCircos.Chromosome.Ideogram.Plot();
	title("RCircos Line Plot Demo");


	#	Plot lines
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	data.col <- 5; 
	track.num <- 1;
	direction <- "in";
	RCircos.Line.Plot(RCircos.Line.Data, data.col, track.num, "in");


	#	Close the graphic device and clear memory
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	dev.off();
	print("RCircos Line Plot Demo Done!");

	rm(list=ls(all=T));
}
	
RCircos.Line.Plot.Demo();



