#
#	This demo draw chromosome ideogram with padding between
#	chromosomes, highlights, chromosome names, and tile plot. 
#
#	Usage:
#
#	library(RCircos);
#	demo("RCircos.Tile.Plot.Demo");
#
#	========================================================


RCircos.Tile.Plot.Demo<-function()
{
	#	Load R source files and define parameters
	#	*********************************************
	library(RCircos);
	circpar <- RCircos.Initialize.Parameters();	
	RCircos.List.Parameters(circpar)


	#	Read chromosome cytoband data for the species 
	#	defined in R.Circos.Source.R file
	#	********************************************
	data(UCSC.HG19.Human.CytoBandIdeogram);
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
	cyto.band <- RCircos.Cytoband.Data(cyto.info, 
			chr.exclude=NULL, circpar);


	#	Calculate x and y values for the base circle plot
	#	********************************************
	circle.positions <- RCircos.Base.Plot.Positions(cyto.band, 
			circpar);



	#	Load tile plot data. Tile data are genomic positions
	#	only in the order of chromosome, chromStart, and  
	#	chromEnd. Each row is for one tile.
	#	********************************************
	data(RCircos.Tile.Data);
	tile.data <- RCircos.Get.Plot.Data(RCircos.Tile.Data, cyto.band);
	if(is.null(tile.data)) {  stop("Please check the data!"); }



	#	Open the graphic device (here is png/pdf file)
	#
	#	out.file= "RCircos.Tile.Plot.Demo.png";
	#	png(file=out.file, height=9, width=8, unit="in", 
	#		type="cairo", res=300);
	#
	#	********************************************
	out.file <- "RCircos.Tile.Plot.Demo.pdf";
	pdf(file=out.file, height=9, width=8);

	par(mai=c(0.5, 0.5, 0.5, 0.5));
	plot.new();
	plot.window(c(-1*circpar$plot.radius, circpar$plot.radius), 
		c(-1*circpar$plot.radius, circpar$plot.radius));



	#	Draw chromosome ideogram
	#	********************************************
	RCircos.Chromosome.Ideogram(cyto.band, circle.positions, 
			circpar);
	title("RCircos.Tile.Plot.Demo");


	#	Plot tiles. 
	#	********************************************
	track.num <- 1;
	direction <- "in";
	RCircos.Tile.Plot(cyto.band, circle.positions,  tile.data, 
			track.num, direction, circpar);



	#	Close the graphic device
	#	********************************************
	dev.off();	print("RCircos Tile Plot Demo Done!");

	rm(list=ls(all=T));
}
	

RCircos.Tile.Plot.Demo();

