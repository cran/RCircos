#
#	This demo draw chromosome ideogram with padding between
#	chromosomes, highlights, chromosome names, and three empty 
#	tracks inside and outside of chromosome ideogram.
#
#	Usage:
#
#	library(RCircos);
#	demo("RCircos.Layout.Demo");
#
#	========================================================


RCircos.Layout.Demo<-function()
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


	#	Modify circpar$plot.radius here since there will
	#	be plot tracks outside of chromosome ideogram
	#	********************************************
	circpar$plot.radius <- 2;



	#	Open the graphic device (here is png/pdf file)
	#
	#	out.file= "RCircos.Layout.Demo.png";
	#	png(file=out.file, height=9, width=8, unit="in", 
	#		type="cairo", res=300);
	#
	#	********************************************
	out.file <- "RCircos.Layout.Demo.pdf";
	pdf(file=out.file, height=9, width=8);

	par(mai=c(0.5, 0.5, 0.5, 0.5));
	plot.new();
	plot.window(c(-1*circpar$plot.radius, circpar$plot.radius), 
		c(-1*circpar$plot.radius, circpar$plot.radius));



	#	Draw chromosome ideogram
	#	********************************************
	RCircos.Chromosome.Ideogram(cyto.band, circle.positions, 
			circpar);
	title("R Circos Layout Demo");



	#	Marking plot areas both inside and outside of
	#	chromosome ideogram ( 3 for each)
	#	********************************************
	total.track <- 3;
	subtrack <- 5;
	one.track <- circpar$track.height + circpar$track.padding;


	#	Outside of chromosome ideogram
	#	********************************************
	for(a.track in 1:total.track)
	{	
		in.pos    <- circpar$track.out.start  + (a.track-1)*one.track;
		out.pos <-in.pos + circpar$track.height;
		RCircos.Track.Outline(cyto.band, circle.positions, 
			out.pos, in.pos, subtrack, circpar);
	}

	#	Inside of chromosome ideogram
	#	********************************************
	for(a.track in 1:total.track)
	{
		out.pos <- circpar$track.in.start  - (a.track-1)*one.track;
		in.pos    <- out.pos - circpar$track.height;
		RCircos.Track.Outline(cyto.band, circle.positions, 
			out.pos, in.pos, subtrack, circpar);
	}



	#	Close the graphic device
	#	********************************************
	dev.off();	print("RCircos Layout Demo Done!");

	rm(list=ls(all=T));
}
	

RCircos.Layout.Demo();


