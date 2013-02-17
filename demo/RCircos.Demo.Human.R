#  ========================================================================
#
#	This demo draw human chromosome ideogram and data tracks for:
#
#	1.	Connectors
#	2.	Gene lables
#	3.	Heatmap
#	4.	Scatterplot 
#	5.	Line plot
#	6.	Histogram
#	7.	Tile plot
#	8.	Link lines
#
#  ========================================================================


RCircos.Demo.Human <- function()
{
	#	Load RCircos package and defined parameters
	#  	*********************************************		
	library(RCircos);
	circpar <- RCircos.Initialize.Parameters();


	#	Load and processing human cytoband data 
	#  	*********************************************
	data(UCSC.HG19.Human.CytoBandIdeogram);
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
	cyto.band <- RCircos.Cytoband.Data(cyto.info, 
			chr.exclude=NULL, circpar);


	#	Calculate x and y coordinates for the base position.  
	#	These coordinates will be used to derive plot 
	#	positions for all other plot items
	#  	**********************************************
	cat("Calculate x and y coordinates for base circle ...\n\n");
	circle.positions <- RCircos.Base.Plot.Positions(cyto.band, 
			circpar);


	#	We will plot muttiple tracks. So increse the plot area
	#	by adjust radius.len from 1.0 to 2.5. Other variables 
	#	derive from radius.len also need changing.
	#  	**********************************************
	radius.len <- 2.5;
	circpar <- RCircos.Reset.Ideogram.Position(circpar, 
			radius.len);




	#	Open the graphic device (here a pdf file)
	#
	# 	out.file= "RCircos.Demo.Human.png";
	# 	png(file=out.file, height=10, width=8, unit="in", 
	#		type="cairo", res=300);
	#
 	#	out.file= "RCircos.Demo.Human.tif";
 	#	tiff(file=out.file, height=10, width=8, unit="in", 
	#		type="cairo", res=300);
	#
	#  	**********************************************
	cat("Open graphic device and start plot ...\n");

	out.file= "RCircos.Demo.Human.pdf";
	pdf(file=out.file, height=10, width=8);

	par(mai=c(0.25, 0.25, 0.25, 0.25));
	plot.new();
	plot.window(c(-1*circpar$plot.radius, circpar$plot.radius), 
		c(-1*circpar$plot.radius, circpar$plot.radius+0.25));
	title("RCircos 2D Track Plot with Human Genome");




	#	Draw chromosome ideogram
	#  	**********************************************
	cat("Draw chromosome ideogram ...\n");
	RCircos.Chromosome.Ideogram(cyto.band, circle.positions, 
			circpar);




	#	Plot connectors in first track and gene names
	#	in the second track. 
	#  	**********************************************
	data(RCircos.Gene.Label.Data);
	label.data <- RCircos.Get.Plot.Data(RCircos.Gene.Label.Data, 
			cyto.band);
	if(is.null(label.data)) {  stop("Please check the data!"); }

	#	To avoid overlap of neighbour gene names, reset
	#	the plot position for gene names. 
	#	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	gene.data <- RCircos.Get.Label.Locations(cyto.band, label.data, 
			 label.type="text", circpar)
	conn.data <- data.frame(gene.data$Location, 
			gene.data$Label.Position);
	name.col <- 4;
	direction <- "in";

	track.num <- 1;
	RCircos.Connector(cyto.band, circle.positions, conn.data, 
		track.num, direction, circpar);

	track.num <- 2;
	RCircos.Gene.Label(circle.positions, gene.data, name.col, 
			track.num, direction,  circpar);




	#	Draw Heatmap.  Since some gene names are longer 
	#	than one track height, we skip two tracks 
	#  	**********************************************
	cat("Add heatmap track ...\n");

	data(RCircos.Heatmap.Data);
	expr.data <- RCircos.Get.Plot.Data(RCircos.Heatmap.Data, 
				cyto.band);
	if(is.null(expr.data)) {  stop("Please check the data!"); }

	track.num <- 5;	
	direction <- "in";
	data.col <- 6;
	RCircos.Heatmap(cyto.band, circle.positions, expr.data, 
			data.col, track.num, direction, circpar);




	#	Scatterplot with DNA copy number variation data
	#  	**********************************************
	cat("Add scatterplot track ...\n");

	data(RCircos.Scatter.Data);
	scatter.data <- RCircos.Get.Plot.Data(RCircos.Scatter.Data, cyto.band);
	if(is.null(scatter.data)) {  stop("Please check the data!"); }


	track.num <- 6; 
	direction <- "in";
	data.col <- ncol(scatter.data)-1;

	RCircos.ScatterPlot(cyto.band, circle.positions, scatter.data, 
			data.col, track.num,  direction,  
			by.fold=1, circpar);



	#	Draw Line plot with DNA copy number variation data
	#  	**********************************************
	cat("Add line plot track ...\n");

	data(RCircos.Line.Data);
	line.data <- RCircos.Get.Plot.Data(RCircos.Line.Data, cyto.band);

	track.num <- 7;
	direction <- "in";
	data.col <- ncol(line.data)-1;

	RCircos.Line.Plot(cyto.band, circle.positions, line.data, data.col, 
			track.num, direction,  circpar);




	#	Draw Histogram
	#  	**********************************************
	cat("Add histogram track ...\n");

	data(RCircos.Histogram.Data);
	hist.data <- RCircos.Get.Plot.Data(RCircos.Histogram.Data, 
			cyto.band);

	track.num <- 8; 
	direction <- "in";
	data.col <- 4;
	RCircos.Histogram(cyto.band, circle.positions, hist.data, data.col, 
			track.num, direction, circpar);




	#	Draw Tile plot. Note: tile plot data have chromosome 
	#	locations only and each data file is for one track
	#  	**********************************************
	cat("Add tile track ...\n");

	data(RCircos.Tile.Data);
	tile.data <- RCircos.Get.Plot.Data(RCircos.Tile.Data, cyto.band);

	track.num <- 9;
	direction <- "in";
	RCircos.Tile.Plot(cyto.band, circle.positions, tile.data, 
			track.num, direction,  circpar);




	#	Draw Link lines. Note: Link data has only paired 
	#	chromosome locations in each row an does not
	#	need extra processing. Also, link lines are always 
	#	drawn inside of chromosome ideogram.
	#  	**********************************************
	cat("Add link track ...\n");
	
	data(RCircos.Link.Data);
	track.num <- 11;
	RCircos.Link.Plot(cyto.band, circle.positions, RCircos.Link.Data, 
			track.num,  by.chromosome=FALSE, 
			circpar);




	#	Close the graphic device and clear memory
	#  	**********************************************
	dev.off();

	cat("R Circos Demo Done ...\n\n");
	rm(list=ls(all=T));
}


RCircos.Demo.Human();













