#
#	This demo draw chromosome ideogram with padding between
#	chromosomes, highlights, chromosome names, and heatmap 
#	with gene expression data.
#
#	Usage:
#
#	library(RCircos);
#	demo("RCircos.Heatmap.Demo");
#
#	========================================================


RCircos.Heatmap.Demo<-function()
{
	#	Load R source files and define parameters
	#	*********************************************
	library(RCircos);
	circpar <- RCircos.Initialize.Parameters();	
	RCircos.List.Parameters(circpar);


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



	#	Load gene expression data. The first four column of
	#	data file have to be chromosome, chromStart, and 
	#	chromEnd, followed by a column of gene names. 
	#	Note: the expr.data will have one more column than
	#	RCircos.Heatmap.Data.
	#	********************************************
	data(RCircos.Heatmap.Data);
	expr.data <- RCircos.Get.Plot.Data(RCircos.Heatmap.Data, 
				cyto.band);
	if(is.null(expr.data)) {  stop("Please check the data!"); }


	#	Open the graphic device (here is png/pdf file)
	#	********************************************
	#	out.file= "RCircos.Heatmap.Demo.png";
	#	png(file=out.file, height=9, width=8, unit="in", 
	#			type="cairo", res=300);
	
	out.file <- "RCircos.Heatmap.Demo.pdf";
	pdf(file=out.file, height=9, width=8);

	par(mai=c(0.5, 0.5, 0.5, 0.5));
	plot.new();
	plot.window(c(-1*circpar$plot.radius, circpar$plot.radius), 
		c(-1*circpar$plot.radius, circpar$plot.radius));



	#	Draw chromosome ideogram
	#	********************************************
	RCircos.Chromosome.Ideogram(cyto.band, circle.positions, 
			circpar);
	title("RCircos Heatmap Demo");



	#	Marking plot areas both inside and outside of
	#	chromosome ideogram ( 5 for each)
	#	********************************************
	total.track <- 5;



	#	plot heatmap Inside of chromosome ideogram
	#	********************************************
	direction <- "in";
	for(a.track in 1:total.track)
	{
		data.col <- a.track + 4;
		RCircos.Heatmap(cyto.band, circle.positions, 
			expr.data, data.col, a.track, 
			direction,  circpar);
	}


	#	Close the graphic device
	#	********************************************
	dev.off();	print("RCircos Heatmap demo done!");

	rm(list=ls(all=T));
}
	

RCircos.Heatmap.Demo();

