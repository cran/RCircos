# ________________________________________________________________________________________
# <><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><>
#
#	This demo draw chromosome ideogram with padding area between chromosomes, 
#	highlights, and chromosome names, label gene names in the second data track
#	and put connectors between chromosome ideogram and gene names.
#
#	Usage:
#
#	library(RCircos);
#	demo("RCircos.Gene.Connector.Demo");
#
# ________________________________________________________________________________________
# <><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><>



RCircos.Gene.Connector.Demo<-function()
{
	#	Load RCircos package and defined parameters
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
	library(RCircos);


	#	Load gene label data and human cytoband data 
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	data(RCircos.Gene.Label.Data);
	data(UCSC.HG19.Human.CytoBandIdeogram);
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


	#	Setup RCircos core components:
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);


	#	Ready to make Circos plot image. Open the graphic device (pdf file)	
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	out.file <- "RCircos.Gene.Connector.Demo.pdf";
	pdf(file=out.file, height=8, width=8);

	RCircos.Set.Plot.Area();
	title("RCircos Gene Label and Connector Demo");


	#	Draw chromosome ideogram
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Draw chromosome ideogram ...\n");
	RCircos.Chromosome.Ideogram.Plot();


	#	Connectors in first track and gene names in the second track. 
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Add Gene and connector tracks ...\n");
	data(RCircos.Gene.Label.Data);

	track.num <- 1;
	RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, 
			track.num, "in");
	name.col <- 4;
	track.num <- 2;
	RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, 
			track.num, "in");


	#	Close the graphic device and clear memory
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	dev.off();
	cat("R Circos Demo Done ...\n\n");
	rm(list=ls(all=T));
}

RCircos.Gene.Connector.Demo();
