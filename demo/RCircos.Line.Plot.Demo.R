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



	#	Load line plot data. The first four column of data
	#	file have to be chromosome, chromStart, and  
	#	chromEnd, followed by a column of gene names. 
	#	Note: the hist.data will have one more column 
	#	than the original RCircos.Histogram.Data.
	#	********************************************
	data(RCircos.Line.Data);
	line.data <- RCircos.Get.Plot.Data(RCircos.Line.Data, cyto.band);
	if(is.null(line.data)) {  stop("Please check the data!"); }


	#	Open the graphic device (here is png/pdf file)
	#
	#	out.file= "RCircos.Line.Plot.Demo.png";
	#	png(file=out.file, height=9, width=8, unit="in", 
	#		type="cairo", res=300);
	#
	#	********************************************
	out.file <- "RCircos.Line.Plot.Demo.pdf";
	pdf(file=out.file, height=9, width=8);

	par(mai=c(0.5, 0.5, 0.5, 0.5));
	plot.new();
	plot.window(c(-1*circpar$plot.radius, circpar$plot.radius), 
		c(-1*circpar$plot.radius, circpar$plot.radius));



	#	Draw chromosome ideogram
	#	********************************************
	RCircos.Chromosome.Ideogram(cyto.band, circle.positions, 
			circpar);
	title("RCircos Line Plot Demo");


	#	Plot lines
	#	********************************************
	data.col <- 5; 
	track.num <- 1;
	direction <- "in";
	RCircos.Line.Plot(cyto.band, circle.positions,  line.data, 
			data.col, track.num, direction, circpar);



	#	Close the graphic device
	#	********************************************
	dev.off();	print("RCircos Line Plot Demo Done!");

	rm(list=ls(all=T));
}
	
RCircos.Line.Plot.Demo();



