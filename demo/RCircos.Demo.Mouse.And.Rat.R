#
#	This demo draws a heatmap with both mouse and rat gene expresion data 
#	and links between same genes. 
#
#	We will put mouse data on right side and rat data on left side. Since we are 
#	going to draw the image with both mouse and rat chromosome ideogram, 
#	we deduct total points of each ideogram by half while keep chromosome 
#	padding as unchanged. This is done by modify the base.per.unit from 3000 
#	to 6000 .
#
#	We also need to use the last mouse chromosome location as the start of rat 
#	chromosome location. Also, modify chromosome names and their positions 
#	to make sure they will not overlap with other items
#
#  ==========================================================================


RCircos.Demo.Mouse.And.Rat<-function()
{
	#	Load RCircos package and defined parameters. 
	#	********************************************	
	library(RCircos);	
	circpar <- RCircos.Initialize.Parameters();


	#	Read chromosome cytoband data. We generate a
	#	pseudo karyotype by adding the last position of 
	#	mouse chromosome plus chromosome padding 
	#	to each of rat chromosome band.
	#  	*********************************************
	cat("Read in cytoband data ...\n");

	data(UCSC.Mouse.GRCm38.CytoBandIdeogram);
	mouse.cyto.info <- UCSC.Mouse.GRCm38.CytoBandIdeogram;
	mouse.cyto.band <- RCircos.Cytoband.Data(mouse.cyto.info, 
	 	chr.exclude=NULL, circpar);

	data(UCSC.Baylor.3.4.Rat.cytoBandIdeogram);
	rat.cyto.info <- UCSC.Baylor.3.4.Rat.cytoBandIdeogram;
	rat.cyto.band <- RCircos.Cytoband.Data(rat.cyto.info,
		chr.exclude=NULL, circpar);

	mouse.last <- mouse.cyto.band$Location[nrow(mouse.cyto.band)];
	rat.cyto.band$Location <- rat.cyto.band$Location + 
		circpar$chrom.paddings + mouse.last;

	mouse.chrs <- paste0("m", as.character(mouse.cyto.band$Chromosome));
	rat.chrs <- paste0("r", as.character(rat.cyto.band$Chromosome));

	two.cyto.band <- rbind(mouse.cyto.band, rat.cyto.band);
	two.cyto.band$Chromosome <- c(mouse.chrs, rat.chrs);


	#	Calculate x and y coordinates for the two base circles. 
	#	These coordinates will be used to derive plot positions 
	#	for other plot items. Since we are going to draw a circle
	#	with both mouse and rat chromosome ideogram, total 
	#	 points here will be the sum of total points for each 
	#	species plus one more padding area. To improve the 
	#	performance speed, we deduct the totoal points to half.
	#  	*************************************************
	cat("Calculate x and y coordinates for base circle ...\n\n");
	circle.positions <- RCircos.Base.Plot.Positions(two.cyto.band,  
			circpar);


	#	Get the plot data
	#  	*********************************************
	cat("Read in gene expression data ...\n");

	data(RCircos.Mouse.Expr.Data);
	mouse.plot.data <- RCircos.Get.Plot.Data(RCircos.Mouse.Expr.Data, 
			mouse.cyto.band);

	data(RCircos.Rat.Expr.Data);
	rat.plot.data <- RCircos.Get.Plot.Data(RCircos.Rat.Expr.Data, 
			rat.cyto.band);


	#	Initialize graphic device (here a png/pdf/tiff file)
	#
	#	out.file= "RCircos.Demo.Mouse.And.Rat.tif";
	#	tiff(file=out.file, height=10, width=9, unit="in", 
	#		type="cairo", res=300);
	#
	# 	out.file= "RCircos.Demo.Mouse.And.Rat.png";
	# 	png(file=out.file, height=10, width=9, unit="in", 
	#	type="windows", res=300);
	#
	#  	*********************************************
	cat("Initialize graphic device ...\n");

	out.file= "RCircos.Demo.Mouse.And.Rat.pdf";
	pdf(file=out.file, height=10, width=9);

	layout(matrix(data=c(1,2), nrow=1, ncol=2),
		 widths=c(8,1), heights=10);

	par(mai=c(0.5, 0.5, 0.5, 0.5));
	plot.new();
	plot.window(c(-1*circpar$plot.radius, circpar$plot.radius), 
		c(-1*circpar$plot.radius, circpar$plot.radius+0.25));

	title("RCircos Demo: Mouse and Rat Gene Expression");


	#	Plot the mouse chromosome ideogram and 
	#	heatmap at default regions (radius.len <- 1.0
	#	and plot.radius <- 1.4)
	#  	*********************************************
	cat("Draw chromosome ideograms ...\n");
	RCircos.Chromosome.Ideogram(two.cyto.band, 
			circle.positions,  circpar);

	cat("Plot mouse and rat data ...\n");
	track.num <- 1;
	direction <- "in";
	data.col <- 5;

	RCircos.Heatmap(mouse.cyto.band, circle.positions, 
		mouse.plot.data, data.col, 
		track.num, direction,  circpar);
	RCircos.Heatmap(rat.cyto.band, circle.positions, 
		rat.plot.data, data.col, 
		track.num, direction,  circpar);


	#	Generate link data and draw link lines
	#  	*********************************************
	mouse.gene <- as.character(mouse.plot.data$Gene);
	rat.gene   <- as.character(rat.plot.data$Gene);
	mouse.data <- mouse.plot.data[mouse.gene %in% rat.gene,];
	mouse.data <- mouse.data[order(mouse.data$Expr.Mean),];
	last <- nrow(mouse.data);
	first <- last - 50;

	mouse.data <- mouse.data[first:last,];
	mouse.gene <- as.character(mouse.data$Gene);
	rat.data <- rat.plot.data[rat.gene %in% mouse.gene,];

	mouse.data <- mouse.data[order(mouse.data$Gene),];
	mouse.data[,1] <- paste0("m", mouse.data[,1]);
	rat.data <- rat.data[order(rat.data$Gene),];
	rat.data[, 1] <- paste0("r", rat.data[,1]);
	link.data <- data.frame(mouse.data[,1:3], rat.data[,1:3]);

	cat("Draw link lines between same genes\n");
	track.num<-2;
	RCircos.Link.Plot(two.cyto.band, circle.positions, link.data, 
			track.num,  by.chromosome = FALSE, 
			circpar);


	#	Add a legend
	#  	*********************************************
	legend(0.6, 1.5, legend=c("Right: GSE42081 Mouse", 
			"Left: GSE42081 Rat"), cex=0.8);


	#	Add a color key
	#  	*********************************************
	RedRamp <- rgb( seq(1, 1, length=256),  		
                     	           		seq(0, 1, length=256),  		
                                		seq(0, 1, length=256)) ; 		

	BlueRamp <- rgb(seq(0, 1, length=256),  		
                                		seq(0, 1, length=256),  		
                                		seq(1, 1, length=256));  		

	ColorRamp   <- cbind(BlueRamp, rev(RedRamp));
	ColorLevels <- seq(min(min(mouse.data$Expr.Mean), 
			min(rat.data$Expr.Mean)), 
			max(max(mouse.data$Expr.Mean), 
			max(mouse.data$Expr.Mean)), 
			length=length(ColorRamp));

	par(mai=c(3,0.1, 3,0.5));
	image(1, ColorLevels, matrix(data=ColorLevels, 
			ncol=length(ColorLevels), nrow=1),
			col=ColorRamp,  
			xlab="",  ylab="",  xaxt="n");




	#	Done ...
	#	********************************************
	dev.off();
	cat("R Circos Demo with mouse and Rat data done ...\n\n");

	rm(list=ls(all=T));
}


RCircos.Demo.Mouse.And.Rat();