#
#	This demo draw chromosome ideogram with padding between
#	chromosomes, highlights, chromosome names, and scatterplot. 
#
#	Usage:
#
#	library(RCircos);
#	demo("RCircos.Scatter.Plot.Demo");
#
#	========================================================


RCircos.Scatter.Plot.Demo<-function()
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



	#	Load link plot data. Link data has rows with paired
	#	genomic position data in the order of chromosome a, 
	#	chromStart A, and  chromEnd A, chromosome B,
	#	chromStart B, and chromEnd B. Different from other
	#	plot data, link data loaded should be used without 
	#	extra processing.
	#	********************************************
	data(RCircos.Scatter.Data);
	scatter.data <- RCircos.Get.Plot.Data(RCircos.Scatter.Data, cyto.band);
	if(is.null(scatter.data)) {  stop("Please check the data!"); }



	#	Open the graphic device (here is png/pdf file)
	#	out.file= "RCircos.Scatter.Plot.Demo.png";
	#	png(file=out.file, height=9, width=8, unit="in", 
	#			type="cairo", res=300);
	#
	#	********************************************
	out.file <- "RCircos.Scatter.Plot.Demo.pdf";
	pdf(file=out.file, height=9, width=8);

	par(mai=c(0.5, 0.5, 0.5, 0.5));
	plot.new();
	plot.window(c(-1*circpar$plot.radius, circpar$plot.radius), 
		c(-1*circpar$plot.radius, circpar$plot.radius));



	#	Draw chromosome ideogram
	#	********************************************
	RCircos.Chromosome.Ideogram(cyto.band, circle.positions, 
			circpar);
	title("RCircos.Scatter.Plot.Demo");


	#	Plot link lines. Link lines are always in most inside 
	#	of chromosome ideogram. 
	#	********************************************
	data.col <- 5;
	track.num <- 1;
	direction <- "in";
	RCircos.ScatterPlot(cyto.band, circle.positions,  scatter.data, 
			data.col, track.num, direction, 
			by.fold=0, circpar);



	#	Close the graphic device
	#	********************************************
	dev.off();	print("RCircos Scatter Plot Demo Done!");

	rm(list=ls(all=T));
}
	

RCircos.Scatter.Plot.Demo();

