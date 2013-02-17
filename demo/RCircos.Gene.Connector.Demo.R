#
#	This demo draw chromosome ideogram with padding between
#	chromosomes, highlights, and chromosome names, label gene
#	names in the second data track and put connectors between 
#	chromosome ideogram and gene names.
#
#	Usage:
#
#	library(RCircos);
#	demo("RCircos.Gene.Connector.Demo");
#
#	=====================================================


RCircos.Gene.Connector.Demo<-function()
{
	#	Load R source files and define parameters
	#	*********************************************
	library(RCircos);
	circpar <- RCircos.Initialize.Parameters();	
	RCircos.List.Parameters(circpar)


	#	Load and processing human cytoband data 
	#	********************************************
	data(UCSC.HG19.Human.CytoBandIdeogram);
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
	cyto.band <- RCircos.Cytoband.Data(cyto.info, 
			chr.exclude=NULL, circpar);


	#	Calculate x and y values for the base circle plot
	#	********************************************
	circle.positions <- RCircos.Base.Plot.Positions(cyto.band, 
			circpar);



	#	Load gene label data. The first four columns of data
	#	are chromosome, chromStart, chromEnd, and gene
	#	names. Note: Plot location of each gene will be 
	#	added to the last column in label.data.
	#	********************************************
	data(RCircos.Gene.Label.Data);
	label.data <- RCircos.Get.Plot.Data(RCircos.Gene.Label.Data, 
			cyto.band);
	if(is.null(label.data)) {  stop("Please check the data!"); }



	#	To avoid overlap of neighbour gene names, reset
	#	the plot position for gene names. 
	#	*********************************************
	gene.data <- RCircos.Get.Label.Locations(cyto.band, label.data, 
			 label.type="text", circpar)
	conn.data <- data.frame(gene.data$Location, 
			gene.data$Label.Position);



	#	Ready to make Circos plot image.
	#	Open the graphic device (here is png/pdf file)
	#	
	#	out.file= "RCircos.Gene.Connector.Demo.png";
	#	png(file=out.file, height=9, width=8, unit="in",
	#		type="windows", res=300);
	#	*********************************************
	out.file <- "RCircos.Gene.Connector.Demo.pdf";
	pdf(file=out.file, height=9, width=8);

	par(mai=c(0.5, 0.5, 0.5, 0.5));
	plot.new();
	plot.window(c(-1*circpar$plot.radius, circpar$plot.radius), 
		c(-1*circpar$plot.radius, circpar$plot.radius));


	#	Draw chromosome ideogram
	#	********************************************
	RCircos.Chromosome.Ideogram(cyto.band, circle.positions, 
			circpar);
	title("RCircos Connector and Gene Label Demo");



	#	Plot connectors in first track and gene names
	#	in the second track. 
	#	********************************************
	name.col <- 4;
	direction <- "in";

	track.num <- 1;
	RCircos.Connector(cyto.band, circle.positions, conn.data, 
		track.num, direction, circpar);

	track.num <- 2;
	RCircos.Gene.Label(circle.positions, gene.data, name.col, 
			track.num, direction,  circpar);



	#	Close the graphic device
	#	********************************************
	dev.off();	
	print("RCircos Connector and Gene Lable Demo Done!");

	rm(list=ls(all=T));
}

RCircos.Gene.Connector.Demo();
