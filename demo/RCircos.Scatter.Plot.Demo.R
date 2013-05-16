# ________________________________________________________________________________________
# <><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><>
#
#	This demo draw chromosome ideogram with padding between chromosomes, highlights, 
#	chromosome names, and scatterplot. 
#
#	Usage:
#
#	library(RCircos);
#	demo("RCircos.Scatter.Plot.Demo");
# ________________________________________________________________________________________
# <><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><>


RCircos.Scatter.Plot.Demo<-function()
{
	#	Load RCircos library
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
	library(RCircos);


	#	Load human cytoband data and scatterplot data
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	data(RCircos.Scatter.Data);
	data(UCSC.HG19.Human.CytoBandIdeogram);
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


	#	Setup RCircos core components:
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);


	#	Open the graphic device (here a pdf file)
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
	out.file <- "RCircos.Scatter.Plot.Demo.pdf";
	pdf(file=out.file, height=8, width=8);

	RCircos.Set.Plot.Area();


	#	Draw chromosome ideogram
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Draw chromosome ideogram ...\n");

	RCircos.Chromosome.Ideogram.Plot();
	title("RCircos Scatterplot Demo");



	#	Scatterplot 
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	data.col <- 5;
	track.num <- 6; 
	RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, track.num, "in", 1);



	#	Close the graphic device
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	dev.off();	print("RCircos Scatter Plot Demo Done!");

	rm(list=ls(all=T));
}
	

RCircos.Scatter.Plot.Demo();

