# R ____________________________________________________________________________________
# R ************************************************************************************
# R
# R	R Source code for RCircos package
# R	
# R	Date created: January 2, 2013
# R	by Hongen Zhang, Ph.D. (hzhang@mail.nih.gov)
# R
# R	Last revised on March 6, 2013
# R	by Hongen Zhang, Ph.D. (hzhang@mail.nih.gov)
# R
# R	Center for Cancer Research Bioiformatics Program, 
# R	National Cancer Institute
# R	National Institutes of Health
# R	Bethesda, Maryland 20892
# R
# R
# R ____________________________________________________________________________________
# R ************************************************************************************




#  =====================================================================
#	Initialize default global parameters for Circos plot. Must be run first.
#  =====================================================================
#
RCircos.Initialize.Parameters<-function()
{
	RCircos.Par <- list(chrom.paddings=3000, 
			base.per.unit=3000,
			radius.len=1.0,   
			chr.ideog.pos=1.1,
			highlight.pos=1.2,
			chr.name.pos=1.3,
			track.out.start=1.4, 
			plot.radius=1.5,
			track.in.start=1.05, 
			chrom.width=0.08,
			track.padding=0.02, 
			track.height=0.1, 
			highlight.width=1, 
			hist.width=1000, 
			text.size=0.4,
			heatmap.width=100, 
			point.type=".", 
			point.size=1,
			max.layers=5);
	RCircos.List.Parameters(RCircos.Par);

	return (RCircos.Par);
}


#  =====================================================================
#	Print out all paramters
#  =====================================================================
#
RCircos.List.Parameters<-function(RCircos.Par)
{
	cat("Parameters initialized.\n\n");

	cat(paste("radius.len:\t",  RCircos.Par$radius.len, "\n"));
	cat(paste("chr.ideog.pos:\t",  RCircos.Par$chr.ideog.pos, "\n"));
	cat(paste("highlight.pos:\t",  RCircos.Par$highlight.pos, "\n"));
	cat(paste("chr.name.pos:\t",  RCircos.Par$chr.name.pos, "\n"));
	cat(paste("plot.radius:\t",  RCircos.Par$plot.radius, "\n"));

	cat(paste("track.in.start:\t",  RCircos.Par$track.in.start, "\n"));
	cat(paste("track.out.start:",  RCircos.Par$track.out.start, "\n"));

	cat(paste("chrom.width:\t",  RCircos.Par$chrom.width, "\n"));
	cat(paste("track.padding:\t",  RCircos.Par$track.padding, "\n"));
	cat(paste("track.height:\t",  RCircos.Par$track.height, "\n\n"));


	cat(paste("base.per.unit:\t",  RCircos.Par$base.per.unit, "\n\n"));
	cat(paste("chrom.paddings:\t", RCircos.Par$chrom.paddings, "\n"));


	cat(paste("highlight.width:",  RCircos.Par$highlight.width, "\n"));
	cat(paste("hist.width:\t",  RCircos.Par$hist.width, "\n"));
	cat(paste("text.size:\t",  RCircos.Par$text.size, "\n"));
	cat(paste("heatmap.width:\t",  RCircos.Par$heatmap.width, "\n\n"));

	cat(paste("point.type:\t",  RCircos.Par$point.type, "\n"));
	cat(paste("point.size:\t",   RCircos.Par$point.size, "\n"));
	cat(paste("max.layers:\t", RCircos.Par$max.layers, "\n\n"))

	cat("Note: Following parameters are derived form radius.len:\n\n");  
	cat("chr.ideog.pos\nhighlight.pos\nchr.name.pos\ntrack.in.start\n");
	cat("track.out.start\nplot.radius\n\n");

	cat("Reset them with RCircos.Reset.Ideogram.Position(radius.len)\n");
	cat("when radius.len is changed. All other parameters could be \n");
	cat("modified from command line.\n\n");
}



#  =====================================================================
#	Reset plot position for chromosome ideogram and related items when
#	radius.len changed
#  =====================================================================
#
RCircos.Reset.Ideogram.Position<-function(RCircos.Par, radius.len)
{
	 RCircos.Par$radius.len <- radius.len;
	 RCircos.Par$chr.ideog.pos <- radius.len + 0.1;
	 RCircos.Par$highlight.pos <- radius.len + 0.2;
	 RCircos.Par$chr.name.pos <- radius.len + 0.3;

	 RCircos.Par$track.in.start    <- RCircos.Par$chr.ideog.pos - 0.05;
	 RCircos.Par$track.out.start <- RCircos.Par$chr.name.pos + 0.1;
	 RCircos.Par$plot.radius        <- radius.len + 0.5;

	cat("Chromosome ideogram position reset:\n\n");
	cat(paste("radius.len",  RCircos.Par$radius.len, "\n"));
	cat(paste("chr.ideog.pos",  RCircos.Par$chr.ideog.pos, "\n"));
	cat(paste("highlight.pos",  RCircos.Par$highlight.pos, "\n"));
	cat(paste("chr.name.pos", RCircos.Par$ chr.name.pos, "\n\n"));

	cat(paste("track.in.start",  RCircos.Par$track.in.start, "\n"));
	cat(paste("track.out.start",  RCircos.Par$track.out.start, "\n\n"));
	cat(paste("plot.radius",  RCircos.Par$plot.radius, "\n\n"));

	cat("Increase value of plot.radius if there are plot tracks\n");
	cat("outside of chromosome ideogram\n\n");

	return (RCircos.Par);
}



#  =====================================================================
#	Read chromosome ideogram data. The cyto.file could be either the 
#	returned value from RCircos.Cyto.File(species) or a chraracter vector. 
#	The chr.exclude should be a vector of chromosome names, if not null.
#  =====================================================================
#
RCircos.Cytoband.Data<-function(cytoband, chr.exclude=NULL,  
			RCircos.Par)
{
	new.cyto <- cytoband[cytoband[,1]=="chr1",];
	new.cyto <- new.cyto[order(new.cyto$ChromEnd),];

	chroms <- paste("chr", c(2:22, "X", "Y"), sep="");
	for(a.chr in 1:length(chroms))
	{
		new.rows <- cytoband[cytoband$Chromosome==chroms[a.chr],];
		if(nrow(new.rows)>0) 
		{
			new.rows <- new.rows[order(new.rows$ChromEnd),];
			new.cyto <- rbind(new.cyto, new.rows);
		}
	}
	cytoband <- new.cyto;

	#       	Reset colors for chromosome bands 
	#	**********************************
	Color <- cytoband$Stain;
	Color <- gsub("gneg", "white", Color);
	Color <- gsub("acen", "red", Color);
	Color <- gsub("stalk", "steelblue", Color);
	Color <- gsub("gvar", "lightgray", Color);
	Color <- gsub("^gpos$", "black", Color);
	Color <- gsub("gpos100", "black", Color);
	Color <- gsub("gpos75", "gray40", Color);
	Color <- gsub("gpos66", "gray50", Color);
	Color <- gsub("gpos50", "gray60", Color);
	Color <- gsub("gpos33", "gray70", Color);
	Color <- gsub("gpos25", "gray80", Color);
	cytoband["BandColor"] <- Color;

	#	Assign colors to chromosome highlight
	#	************************************
	Color <- cytoband$Chromosome;
	Color <- gsub("chr1$", "red", Color);
	Color <- gsub("chr2$", "seagreen", Color);
	Color <- gsub("chr3", "violetred4", Color);
	Color <- gsub("chr4", "orange", Color);
	Color <- gsub("chr5", "magenta", Color);
	Color <- gsub("chr6", "darkgreen", Color);
	Color <- gsub("chr7", "blue", Color);
	Color <- gsub("chr8", "sienna", Color);
	Color <- gsub("chr9", "palevioletred", Color);
	Color <- gsub("chr10", "mediumseagreen", Color);
	Color <- gsub("chr11", "brown", Color);
	Color <- gsub("chr12", "coral", Color);
	Color <- gsub("chr13", "steelblue", Color);
	Color <- gsub("chr14", "turquoise", Color);
	Color <- gsub("chr15", "purple", Color);
	Color <- gsub("chr16", "green", Color);
	Color <- gsub("chr17", "darkred", Color);
	Color <- gsub("chr18", "cyan", Color);
	Color <- gsub("chr19", "violetred", Color);
	Color <- gsub("chr20", "skyblue", Color);
	Color <- gsub("chr21", "aquamarine", Color);
	Color <- gsub("chr22", "darkorchid", Color);
	Color <- gsub("chrX", "salmon", Color);
	Color <- gsub("chrY", "paleturquoise", Color);
	cytoband["ChrColor"] <- Color;


	#	Delete the unneeded chromosomes
	#	*********************************************
	if(length(chr.exclude)>0) 
	{
		for(a.chr in 1:length(chr.exclude))
		{
			the.chr <- paste(chr.exclude[a.chr], "$", sep="");
			if(length(grep("chr", the.chr))==0)
			{  the.chr <- paste("chr", the.chr, sep=""); }

			the.ind <- grep(the.chr, cytoband$Chromosome);
			if(length(the.ind)==0) {
			        stop("Check chromosomes to be excluded!"); 
			} else { cytoband <- cytoband[-the.ind,]; }
		}
	}


	#	Total base pairs and relative length of each band
	#	************************************************
	band.len <- cytoband$ChromEnd - cytoband$ChromStart;
	cytoband["Length"] <- band.len;
	cytoband["Unit"]   <- round(band.len/RCircos.Par$base.per.unit, digits=0);
	

	#       Relative locations of each band in clockwise
	#      ******************************************************
	Relative.Loc <- cytoband$Unit;
	for(i in 2:length(Relative.Loc))
	{ Relative.Loc[i] <- Relative.Loc[i] + Relative.Loc[i-1]; }
	cytoband["Location"] <- Relative.Loc;

	if( RCircos.Par$chrom.paddings>0) {  
		chroms <- unique(cytoband$Chromosome);
		chroms <- chroms[(chroms=="chr1")==F];
		num.pad <-  RCircos.Par$chrom.paddings;
		for(a.chr in 1:length(chroms))
		{
			index <- grep(paste(chroms[a.chr], "$", sep=""), 
				cytoband$Chromosome);
			cytoband$Location[index] <- num.pad + 
				cytoband$Location[index];
			num.pad <- num.pad +  RCircos.Par$chrom.paddings;
		}
	}
	
	return (cytoband);
}



#  =====================================================================
#	Calculate the x and y coordinates for a circle line. These values will be 
#	used as base values to calculate the locations of plot tracks and 
#	positions of chromosome band. Rotation degrees are also attached for 
#	text labeling at each point
#  =====================================================================
#
RCircos.Base.Plot.Positions<-function(cyto.band,  RCircos.Par)
{
	total.points <- cyto.band$Location[nrow(cyto.band)] +  
			RCircos.Par$chrom.paddings;

	interval <- 2*pi/total.points;
	base.val <- seq(0, 2*pi, interval);

	cor.x <- sin(base.val);
	cor.y <- cos(base.val);

	degree <- rep(0, length(base.val));
	mid <- round((length(base.val)-1)/2, digits=0) + 1;

	for(pt in 1:mid)
	{ degree[pt] <- 90 - (base.val[pt]*180/pi);}
	
	for(pt in (mid+1):length(base.val))
	{ degree[pt] <- 270 - (base.val[pt]*180/pi); }

	return (data.frame(cor.x, cor.y, degree));
}



#  =====================================================================
#	Draw chromosome ideogram. Graphic device must be initialized first
#  =====================================================================
#
RCircos.Chromosome.Ideogram<-function(cyto.band, base.positions,  RCircos.Par)
{
	#	Plot chromosome highlights, chromosome outlines,
	#	and chromosome names
	#	**********************************************
	right.side <- nrow(base.positions)/2;
	outer.location <- RCircos.Par$chr.ideog.pos + RCircos.Par$chrom.width;
	inner.location <- RCircos.Par$chr.ideog.pos

	chroms <- unique(cyto.band$Chromosome);
	for(a.chr in 1:length(chroms))
	{
		the.chr  <- cyto.band[cyto.band$Chromosome==chroms[a.chr],];
		start <- the.chr$Location[1]- the.chr$Unit[1] + 1;
		end   <- the.chr$Location[nrow(the.chr)];
		mid <- round((end-start+1)/2, digits=0)+start;
		
		chr.color <- the.chr$ChrColor[nrow(the.chr)];


		#	draw chromosome outlines
		#	***************************************
		pos.x<- c(base.positions[start:end,1]*outer.location, 
			base.positions[end:start,1]*inner.location);
		pos.y<- c(base.positions[start:end,2]*outer.location, 
			base.positions[end:start,2]*inner.location);
		polygon(pos.x, pos.y);


		#	add chromosome names
		#	***********************************************
		chr.name <- sub(pattern="chr", replacement="", chroms[a.chr]);
		text(base.positions[mid,1]*RCircos.Par$chr.name.pos,
			base.positions[mid,2]*RCircos.Par$chr.name.pos,
			label=chr.name, srt=base.positions$degree[mid]);


		#	add chromosome highlights
		#	*******************************
		lines(base.positions[start:end,]*RCircos.Par$highlight.pos,
			col=chr.color, lwd=RCircos.Par$highlight.width);
	}
	
	#	Add chromosome bands
	#	*********************************************************
	for(a.band in 1:nrow(cyto.band))
	{
		a.color <- cyto.band$BandColor[a.band];
		if(a.color=="white") { next; }

		start <- cyto.band$Location[a.band]-cyto.band$Unit[a.band]+1;
		end   <- cyto.band$Location[a.band];

		pos.x<- c(base.positions[start:end,1]*outer.location, 
			base.positions[end:start,1]*inner.location);
		pos.y<- c(base.positions[start:end,2]*outer.location, 
			base.positions[end:start,2]*inner.location);
		polygon(pos.x, pos.y, col=a.color, border=NA);
		
	}

	return (TRUE);
}



#  =====================================================================
#	Calculate plot coordinates for a chromosome position associated to a 
#	data point (e.g., the chromosome name and start position of a gene).
#  =====================================================================
#
RCircos.Data.Point<-function(cyto.band, chromosome, start)
{
	the.point <- 0;

	if(length(grep("chr", chromosome)) == 0) 
	{chromosome <- paste("chr", chromosome, sep=""); }

	chrom.rows <- grep(paste("^", chromosome, "$", sep=""), 
			cyto.band$Chromosome);
	the.row <- which(cyto.band$ChromStart[chrom.rows] <= start &  
			cyto.band$ChromEnd[chrom.rows] >= start)[1];

	band.length <- cyto.band$Length[chrom.rows[the.row]];
	band.units <- cyto.band$Unit[chrom.rows[the.row]];
	band.location <- cyto.band$Location[chrom.rows[the.row]];

	the.bases <- start - cyto.band$ChromStart[chrom.rows[the.row]] ;
	the.units  <- the.bases/band.length*band.units;

	the.point <- band.location - band.units + the.units;
	
	return (round(the.point, digits=0));
}



#  =====================================================================
#	Validate input dataset for correct chromosome names, chromosome
#	start and chromosome end positions. Chromosome names will be
#	converted to character vectors if they are factor variables.
#  =====================================================================
#
RCircos.Validate.Genomic.Data<-function(genomic.data, cyto.band,
			plot.type=c("plot", "link"))
{
	
	if(plot.type=="plot") { chrom.col <- 1; } else { chrom.col <- c(1,4); }
	
	for (a.col in 1:length(chrom.col))
	{
		the.col <- chrom.col[a.col];

		#	Make sure chromosome names have prefix
		#	**************************************************
		genomic.data[,the.col] <- as.character(genomic.data[,the.col]);
		for(a.row in 1:nrow(genomic.data)) {
			if(length(grep("chr", genomic.data[a.row,the.col]))==0) 
			{ genomic.data[a.row,the.col] <- paste("chr", 
				genomic.data[a.row,the.col], sep=""); }
		}

		#	Make sure chromosomes in input data are all included 
		#	in chromosome ideogram data
		#	*************************************************
		cyto.chroms <- as.character(unique(cyto.band$Chromosome));
		data.chroms <- unique(genomic.data[,the.col]);
		if(sum(data.chroms %in% cyto.chroms) < length(data.chroms)) { 
			stop("Error! Some chromosomes not in cyto.band data."); 
			return (NULL); 
		}
	

		#	Make sure chromosome start and end locations in
		#	input data are greater than 0
		#	************************************************
		if(min(genomic.data[,the.col+1])<0) 
		{ print("Error! chromStart position less than 0."); return (NULL);  }
		if(min(genomic.data[,the.col+2])<0) 
		{ print("Error! chromEnd position less than 0."); return (NULL);  }	


		#	Make sure chromosome start and end locations in
		#	input data are not out of chromosome length
		#	************************************************
		for(a.chr in 1:length(data.chroms))
		{
			the.chr      <- data.chroms[a.chr];
			in.data      <- genomic.data[genomic.data[,the.col]==the.chr,];
			cyto.data <- cyto.band[grep(the.chr, cyto.band[,1]),]

			if(max(in.data[,the.col+1])>max(cyto.data[,3]) | 
				max(in.data[,the.col+2])>max(cyto.data[,3]))
			{  
				stop("Error! Location is outside of chromosome length."); 
				print(paste(the.chr, max(in.data[,2]), max(in.data[,3]))); 
				return (NULL); 
			}
		}

		
		#	Make sure all chromosome start positions are smaller than
		#	their paired chromosome end positions
		#	****************************************************
		for(a.row in 1:nrow(genomic.data))
		{
			if(genomic.data[a.row, the.col+1]>genomic.data[a.row, the.col+2]) 
			{ 
				print("chromStart greater than chromEnd"); 
				print(paste("Row:", a.row, genomic.data[a.row, 2],  
					genomic.data[a.row, 3]));
				return(NULL); 
			}
		}
	}

	cat("Input data validated.\n");
	return (genomic.data);
}


#  =====================================================================
#	Calculate x and y coordinates for each row of a data set.
#  =====================================================================
#	
RCircos.Get.Plot.Data<-function(genomic.data, cyto.band)
{
	#	Check chromosome names,Start, andEnd positions
	#	***********************************************
	genomic.data <- RCircos.Validate.Genomic.Data(genomic.data, 
			cyto.band, plot.type="plot");


	#	Calculate the point index for each chromosome location
	#	**************************************************
	data.points <- rep(0, nrow(genomic.data));
	for(a.row in 1:nrow(genomic.data))
	{
		chromosome <- as.character(genomic.data[a.row, 1]);

		if(genomic.data[a.row, 2]>genomic.data[a.row, 3]) 
		{ 
			print("chromStart greater than chromEnd"); 
			print(paste("Row:", a.row, genomic.data[a.row, 2],  
				genomic.data[a.row, 3]));
			return(NULL); 
		}

		location <- round((genomic.data[a.row, 2]+genomic.data[a.row, 3])/2, 
				digits=0);
		data.points[a.row] <- RCircos.Data.Point(cyto.band, 
				chromosome, location);
	}
	genomic.data["Location"] <- data.points;


	#	Sort the data by chromosome then start position
	#	*******************************************
	genomic.data <- genomic.data[order(genomic.data$Location),];
	cat("Plot data done!\n\n");
	return (genomic.data);
}



#  =====================================================================
#	Calculate inner and outer positions of a track
#  =====================================================================
#
RCircos.Track.Positions<-function(side="in", track.num, RCircos.Par)
{
	#	Plot position for current track. There is padding between 
	#	tracks
	#	**************************************************
	one.track <- RCircos.Par$track.height + RCircos.Par$track.padding;
	side <- tolower(side);
	if(side=="in") {
		out.pos <- RCircos.Par$track.in.start -(track.num-1)*one.track;
		in.pos    <- out.pos - RCircos.Par$track.height;

	} else if(side=="out"){
		in.pos <- RCircos.Par$track.out.start +(track.num-1)*one.track;
		out.pos <-in.pos + RCircos.Par$track.height;

	} else {  stop("Incorrect track location. It must be \"in\" or \"out\" ..."); }

	return (c(out.loc=out.pos, in.loc=in.pos));
}



#  =====================================================================
#	In case of too many gene or genomic positions needed to be labeled
#	or connected, the text or lines may become very crowded. Reset plot 
#	positions to make the image more readable.
#  =====================================================================
RCircos.Get.Label.Locations<-function(cyto.band, genomic.data,  
		label.type=c("text", "line", "point"), RCircos.Par)
{
	#	Check chromosome names, Start, and End positions
	#	********************************************************
	genomic.data <- RCircos.Validate.Genomic.Data(genomic.data, 
			cyto.band, plot.type="plot");


	#	Get maximum number of lables for each chromosome. 
	#	**********************************************************
	chromosomes <- unique(cyto.band$Chromosome);

	label.type <- tolower(label.type);
	if(label.type=="text") {  
		label.width <- RCircos.Par$text.size/0.4*5000; 
	} else if (label.type=="line" || label.type=="point") { 
		label.width <- RCircos.Par$base.per.unit*100; 
	} else { stop("label.type must be \"text\", \"line\", or \"point\" ... ");  }

	max.labels <- rep(0, length(chromosomes));
	start.loc <- rep(0, length(chromosomes));	
	end.loc  <- rep(0, length(chromosomes));	
	for(a.chr in 1:length(chromosomes))
	{
		index <- which(cyto.band$Chromosome==chromosomes[a.chr]);
		total.units <- sum(cyto.band$Unit[index]);
		max.labels[a.chr] <- round(total.units/label.width, digits=0);
		start.loc[a.chr] <- min(cyto.band$Location[index]);
		end.loc[a.chr] <- max(cyto.band$Location[index]);
	}


	#	Attach new column for label locations. 
	#	**********************************************
	label.data <- NULL;
	genomic.data["Label.Position"] <- genomic.data$Location;


	#	Reset label locations
	#	**************************************************
	for(a.chr in 1:length(chromosomes))
	{
		index<- which(genomic.data[,1]==chromosomes[a.chr]);
		if(length(index)==0) { next; }

		#	If there too many gene labels, remove extra 
		#	genes for best visulization
		#	*****************************************
		if(length(index)>=max.labels[a.chr])
		{
			cat(paste("Maximum lables for ", chromosomes[a.chr],
				"is", max.labels[a.chr],". " ));
			cat("Extra ones are ignored.\n");

			the.chr <- genomic.data[index[1:max.labels[a.chr]],];
			the.chr <- the.chr[order(the.chr$Location),];

			for(a.gene in 1:nrow(the.chr))
			{
				new.loc  <- start.loc[a.chr] + (a.gene-1)*
						label.width;  
				the.chr$Label.Position[a.gene] <- new.loc;
			}
		} else {
			#	modify label locations if necessary
			#	***************************************
			the.chr <- genomic.data[index,];
			the.chr <- the.chr[order(the.chr$Location),];
			all.gene <- nrow(the.chr);
			last.pos <- the.chr$Label.Position[1];
			for(a.gene in 1:all.gene)
			{
				curr.loc <-  the.chr$Label.Position[a.gene];

				len.needed <- curr.loc + 
					(all.gene-a.gene)*label.width;
				if(len.needed>end.loc[a.chr]) 
				{  curr.loc  <- curr.loc - (len.needed-end.loc[a.chr]); }

				if(a.gene>1 && (curr.loc-last.pos) < label.width)
				{ curr.loc  <- last.pos + label.width;  }
				
				the.chr$Label.Position[a.gene] <- curr.loc;
				last.pos <- curr.loc;				
			}
		}
		label.data <- rbind(label.data, the.chr); 
	}
	colnames(label.data) <- colnames(genomic.data);

	return (label.data);
}



#  =====================================================================
#	Draw connectors between two tracks. The connect.data argument has 
#	two columns paired point locations on the circle (chromosomes) . 
#	The first column is outer points,  and the second column is inner points. 
#	The points are sorted by relative positions on their chromosomes  and
#	are held in the last two columns of connect.data
#	==============================================================
#
RCircos.Connector<-function(cyto.band, base.positions, connect.data, track.num, 
		side, RCircos.Par)
{
	#	Plot position for current track. 
	#	*****************************************************
	locations <- RCircos.Track.Positions(side, track.num, RCircos.Par);
	out.pos <- locations[1];
	in.pos <- locations[2];

	#	Heights for the two vertical lines of connectors and
	#	the horizontal line range
	#	**************************************************
	v.height <- round((out.pos-in.pos)/10, digits=4);
	h.range <- out.pos - in.pos - 2*v.height;

	genomic.col <- ncol(connect.data) - 1;
	label.col <- ncol(connect.data);

	chroms <- unique(connect.data[,1]);
	for(a.chr in 1:length(chroms))
	{
		chr.row <- which(connect.data[,1]==chroms[a.chr]);
		total <- length(chr.row);

		for(a.point in 1:total)
		{
			top.loc <- out.pos - v.height;
			bot.loc <-  in.pos +  v.height;

			p1 <- connect.data[chr.row[a.point], genomic.col];
			p2 <- connect.data[chr.row[a.point], label.col];


			#	draw top vertical line
			#	************************************
			lines(c(base.positions[p1, 1]*out.pos,
				base.positions[p1,  1]*top.loc),
				c(base.positions[p1,2]*out.pos,
				base.positions[p1, 2]*top.loc));

			#	draw bottom vertical line
			#	************************************
			lines(c(base.positions[p2, 1]*bot.loc, 
				base.positions[p2, 1]*in.pos),
				c(base.positions[p2,2]*bot.loc,
				base.positions[p2, 2]*in.pos));

			#	draw horizontal line
			#	************************************
			lines(c(base.positions[p1,  1]*top.loc, 
				base.positions[p2, 1]*bot.loc),
				c(base.positions[p1, 2]*top.loc,
				base.positions[p2, 2]*bot.loc));
		}
	}
}



#  =====================================================================
#	Draw one track of heatmap with blue and read colors
#  =====================================================================
#
RCircos.Heatmap<-function(cyto.band, base.positions, heatmap.data, data.col, 
		track.num, side,  RCircos.Par)
{
	#	Plot position for current track. 
	#	***********************************************
	locations <- RCircos.Track.Positions(side, track.num, RCircos.Par);
	out.pos <- locations[1];
	in.pos <- locations[2];


	#	outline of chromosomes. No lines inside.
	#	**********************************************
	chroms <- unique(cyto.band$Chromosome);
	for(a.chr in 1:length(chroms))
	{
		the.chr  <- cyto.band[cyto.band[,1]==chroms[a.chr],];
		start <- the.chr$Location[1]- the.chr$Unit[1] + 1;
		end   <- the.chr$Location[nrow(the.chr)];

		polygon.x<- c(base.positions[start:end,1]*out.pos, 
			base.positions[end:start,1]*in.pos);
		polygon.y<- c(base.positions[start:end,2]*out.pos, 
			base.positions[end:start,2]*in.pos);
		polygon(polygon.x, polygon.y, col="white");
	}


	#	Colors for different data values
	#	*********************************************
	RedRamp <- rgb( seq(1, 1, length=256),  seq(0, 1, length=256),  
			seq(0, 1, length=256)) ; 
	BlueRamp <- rgb(seq(0, 1, length=256),  seq(0, 1, length=256),  
			seq(1, 1, length=256));		
	ColorRamp   <- cbind(BlueRamp, rev(RedRamp));

	columns <- 5:(ncol(heatmap.data)-1);
	min.value <- min(as.matrix(heatmap.data[, columns]));
	max.value <- max(as.matrix(heatmap.data[, columns]));
	ColorLevel  <- seq(min.value, max.value, length=length(ColorRamp));


	#	Draw heatmap
	#	***************************************************
	heatmap.value <- as.numeric(heatmap.data[, data.col]);
	for(a.point in 1:length(heatmap.value))
	{
		the.point <- heatmap.data[a.point, ncol(heatmap.data)];
		the.level <- which(ColorLevel>=heatmap.value[a.point]);
		cell.color <- ColorRamp[min(the.level)];

		gene.length <- heatmap.data[a.point, 3] - 
				heatmap.data[a.point, 2];
		cell.width <- gene.length/RCircos.Par$base.per.unit;
		if(cell.width<RCircos.Par$heatmap.width) 
		{ cell.width <- RCircos.Par$heatmap.width; }

		start <- the.point - cell.width/2;
		end <- the.point + cell.width/2

		polygon.x<- c(base.positions[start:end,1]*out.pos, 
			base.positions[end:start,1]*in.pos);
		polygon.y<- c(base.positions[start:end,2]*out.pos, 
			base.positions[end:start,2]*in.pos);
		polygon(polygon.x, polygon.y, col=cell.color, border=NA);
	}
}



#  =====================================================================
#	Draw one track of histogram. 
#  =====================================================================
#
RCircos.Histogram<-function(cyto.band, base.positions, hist.data, data.col, 
			track.num, side, RCircos.Par)
{
	#	Plot position for current track. 
	#	************************************************
	locations <- RCircos.Track.Positions(side, track.num, RCircos.Par);
	out.pos <- locations[1];
	in.pos <- locations[2];


	#	Draw histogram
	#	***********************************************
	num.subtrack <- 5;
	RCircos.Track.Outline(cyto.band, base.positions, out.pos, in.pos, 
			num.subtrack, RCircos.Par);

	for(a.point in 1:nrow(hist.data))
	{
		the.point <- hist.data[a.point, ncol(hist.data)];
		hist.height <- hist.data[a.point, data.col];

		if(the.point>=RCircos.Par$hist.width) { 
			start  <- the.point-RCircos.Par$hist.width; 
		} else { start  <- 1; };

		if((the.point+RCircos.Par$hist.width)<=length(base.positions[,1])) { 
			end   <- the.point+RCircos.Par$hist.width;
		} else { end   <- length(base.positions[,1]); }

		height <- in.pos + RCircos.Par$track.height*hist.height;

		polygon.x<- c(base.positions[start:end,1]*height, 
			base.positions[end:start,1]*in.pos);
		polygon.y<- c(base.positions[start:end,2]*height, 
			base.positions[end:start,2]*in.pos);
		polygon(polygon.x, polygon.y, col="red", border=NA);
	}
}



#  =====================================================================
#	Draw one track of line plot
#  =====================================================================
#
RCircos.Line.Plot<-function(cyto.band, base.positions, line.data, data.col, 
		track.num, side, RCircos.Par)
{
	#	Plot position for current track. 
	#	*****************************************************
	locations <- RCircos.Track.Positions(side, track.num, RCircos.Par);
	out.pos <- locations[1];
	in.pos <- locations[2];


	#	Check if the data has negative values such as copy number
	#	change or log2 of fold change. If yes the zero line will use 
	#	the middle of track height otherwise the inner boundary
	#	****************************************************
	if(min(as.numeric(line.data[,data.col]))>=0) {
		point.bottom <- in.pos; data.ceiling <- 10;
	} else {  
		point.bottom <- in.pos + (RCircos.Par$track.height/2);  
		data.ceiling <- 5;
	}
	sub.height <- out.pos-point.bottom;


	#	Start plotting
	#	***************************************************	
	RCircos.Track.Outline(cyto.band, base.positions, out.pos, in.pos, 
				5, RCircos.Par);

	for(a.point in 1:(nrow(line.data)-1))
	{
		point.one <-  line.data[a.point, ncol(line.data)];
		point.two <- line.data[a.point+1, ncol(line.data)];

		if(line.data[a.point,1]!= line.data[a.point+1,1]) { next;}

		if(line.data[a.point, data.col] >data.ceiling) { 
			value.one<- data.ceiling;
		} else if(line.data[a.point, data.col] <(-1*data.ceiling)) {  
			value.one <- data.ceiling*-1;
		} else { 	value.one <- line.data[a.point, data.col];  }

		height.one <- point.bottom  + value.one/data.ceiling*sub.height;


		if(line.data[a.point+1, data.col] >data.ceiling) { 
			value.two<- data.ceiling;
		} else if(line.data[a.point+1, data.col] <(-1*data.ceiling)) {  
			value.two <- data.ceiling*-1;
		} else { 	value.two <- line.data[a.point+1, data.col];  }

		height.two <- point.bottom  + value.two/data.ceiling*sub.height;

		lines(c(base.positions[point.one ,1]*height.one,
			base.positions[point.two ,1]*height.two),
		       	c(base.positions[point.one , 2]*height.one,
			base.positions[point.two ,2]*height.two));
	}
}



#  =====================================================================
#
#	Draw a quandratic Bezier curve between two chromosome locations 
#	with the equation:  
#
#		B(t) = (1-t) ((1-t)P0 + tP1) + t((1-t)P1 + tP2)
#
#	where P0 is the start point, P2 is the end point, and P1 is the control 
#	point. Since we set P1 to (0,0), the equation become: 
#
#		B(t) =(1-t)^2*P0 + t^2*P2
#
#  ======================================================================
#
RCircos.Link.Line<-function(point.locations, point.one, point.two, color)
{	
	#	Set up the points for Bezure curve
	#	******************************************
	P0 <- as.numeric(point.locations[point.one,]);
	P2 <- as.numeric(point.locations[point.two,]);

	#	Calculate total number of points for Bezuer curve
	#	**********************************************
	bc.point.num <- 1000;
	t <- seq(0, 1, 1/bc.point.num);


	#	Calculate the point values for Bezuer curve using the
	#	equation: (1-t) ((1-t)P0 + tP1) + t((1-t)P1 + tP2). Since
	#	P1 is 0, the equation will become: (1-t)^2*P0 + t^2*P2
	#	***********************************************
	
	link.x <- (1-t)^2*P0[1] + t^2*P2[1];
	link.y <- (1-t)^2*P0[2] + t^2*P2[2];


	#	Draw Bezuer curve
	#	******************************
	lines(link.x, link.y, type="l", col=color);
}



#  ======================================================================
#	Draw link lines between two chromosome locations. Link is always 
#	inside of ideogram.
#  ======================================================================
#
RCircos.Link.Plot<-function(cyto.band, base.positions, link.data, track.num, 
		by.chromosome=FALSE,  RCircos.Par)
{
	#	Check chromosome names, Start, and End positions
	#	********************************************************
	link.data <- RCircos.Validate.Genomic.Data(link.data, 
			cyto.band, plot.type="link");


	#	Plot position for link track.
	#	************************************************
	one.track <- RCircos.Par$track.height + RCircos.Par$track.padding;
	start <- RCircos.Par$track.in.start - (track.num-1)*one.track;
	base.positions <- base.positions*start;

	data.points <- matrix(rep(0, nrow(link.data)*2), ncol=2);
	line.colors <- rep("blue", nrow(link.data));
	for(a.link in 1:nrow(link.data))
	{
		data.points[a.link, 1] <- RCircos.Data.Point(cyto.band,
			link.data[a.link, 1], link.data[a.link, 2]);
		data.points[a.link, 2] <- RCircos.Data.Point(cyto.band,
			link.data[a.link, 4], link.data[a.link, 5]);

		if(data.points[a.link, 1]==0 || data.points[a.link, 2]==0)
		{  print("Error in chromosome locations ...");  break; }

		if(by.chromosome==TRUE) {
			if(link.data[a.link, 1]==link.data[a.link, 4]) {
				line.colors[a.link] <- "red";
			} 
		} else {  line.colors[a.link] <- a.link; }
	}


	#	Draw link lines for each pair of locations
	#	****************************************
	for(a.link in 1:nrow(data.points))
	{  
		point.one <- data.points[a.link, 1];
		point.two <- data.points[a.link, 2];
		if(point.one > point.two)
		{ 
			point.one <- data.points[a.link, 2];
			point.two <- data.points[a.link, 1];
		 }
		RCircos.Link.Line(base.positions, point.one, 
			point.two, line.colors[a.link] ); 
	}
}



#  ======================================================================
#	Draw one track of scatterplot. The scatterplot data should contain
#  ======================================================================
#
RCircos.ScatterPlot<-function(cyto.band, base.positions, scatter.data, data.col, 
		track.num, side, by.fold=0,  RCircos.Par)
{
	#	Plot position for current track. 
	#	*****************************************************
	locations <- RCircos.Track.Positions(side, track.num, RCircos.Par);
	out.pos <- locations[1];
	in.pos <- locations[2];


	#	Check if the data has negative value such as copy number
	#	change or log2 of fold change. If yes the zero line will use 
	#	the middle of track height otherwise the inner boundary
	#	*****************************************************
	if(min(as.numeric(scatter.data[,data.col]))>=0) {
		point.bottom <- in.pos; data.ceiling <- 10;
	} else {  
		point.bottom <- in.pos + (RCircos.Par$track.height/2);  
		data.ceiling <- 5;
	}
	sub.height <- out.pos-point.bottom;

	
	#	Start plotting
	#	***************************************************
	RCircos.Track.Outline(cyto.band, base.positions, out.pos, in.pos, 
			5, RCircos.Par);

	for(a.point in 1:nrow(scatter.data))
	{
		the.point <- scatter.data[a.point, ncol(scatter.data)];

		#	Adjust the data value to avoid overflow
		#	****************************************
		if(scatter.data[a.point, data.col] >data.ceiling) { 
			the.value<- data.ceiling;
		} else if(scatter.data[a.point, data.col] <(-1*data.ceiling)) {  
			the.value <- data.ceiling*-1;
		} else { 	the.value <- scatter.data[a.point, data.col];  }

		color <- "black";
		if(by.fold>0) {
			if(the.value>=by.fold) { color <- "red"; 
			} else if (the.value<=-by.fold) { color <- "blue";
			} else {  color <- "black"; }
		} 

		#	plot a scatter
		#	********************************************
		height <- point.bottom  + the.value/data.ceiling*sub.height;
		points(base.positions[the.point ,1]*height,
		       	base.positions[the.point ,2]*height,
			col=color, pch=RCircos.Par$point.type, 
			cex=RCircos.Par$point.size);
	}
}



#  ======================================================================
#	Label genes beside of track. This is only suitable for small number of 
#	labels. When cex=0.4, each character of the lable will occupy about 5000 
#	units. This is the best visulization for a 8x8 inches image.
#  ======================================================================
#
RCircos.Gene.Label<-function(base.positions, gene.data, name.col, track.num,  
		side,  RCircos.Par)
{
	#	Label positions
	#	***************************************************
	side <- tolower(side);
	locations <- RCircos.Track.Positions(side, track.num, RCircos.Par);

	if(side=="in") {label.pos <- locations[1]; } else {label.pos  <- locations[2];}

	right.side <- nrow(base.positions)/2;


	#	Plot labels
	#	****************************************************
	for(a.text in 1:nrow(gene.data))
	{
		gene.name <- as.character(gene.data[a.text, name.col]);
		the.point <- as.numeric(gene.data[a.text, ncol(gene.data)]);

		rotation <- base.positions$degree[the.point];

		if(side=="in") {
			if(the.point<=right.side) {  text.side<-2; 
			} else {  text.side<-4;  }
		} else {
			if(the.point<=right.side) {  text.side<-4; 
			} else {  text.side<-2;  }
		}

		text(base.positions[the.point,1]*label.pos,
			base.positions[the.point,2]*label.pos,
			label=gene.name, pos=text.side,  
			cex=RCircos.Par$text.size, 
			srt=rotation, offset=0);
	}
}



#  ======================================================================
#	Draws one track of tiles. 
#  ======================================================================
#
RCircos.Tile.Plot<-function(cyto.band, base.positions, tile.data, track.num, 
		side,  RCircos.Par)
{
	#	Plot position for current track. 
	#	*****************************************************
	locations <- RCircos.Track.Positions(side, track.num, RCircos.Par);
	out.pos <- locations[1];
	in.pos <- locations[2];

	
	#	Assign a lay number to each data point and find the maxium
	#	layer number
	#	*****************************************************
	the.layer <- 1;
	the.chr <- tile.data[1,1];
	start <- tile.data[1,2];  
	end  <- tile.data[1,3];
	tile.layers <- rep(1, nrow(tile.data));
	for(a.row in 2:nrow(tile.data))
	{
		#	Meet a new region without overlap with previous or
		#	a different chromosome, reset relevant variables
		#	**********************************************
		if (tile.data[a.row, 2] >= end ) {  
			the.layer <- 1;  
			start <- tile.data[a.row, 2];   
			end <- tile.data[a.row, 3];
		} else if (tile.data[a.row, 1] != the.chr) { 
			the.layer <- 1;  
			the.chr <- tile.data[a.row,1];
			start <- tile.data[a.row, 2]; 
			end <- tile.data[a.row, 3];
		} else {  
			the.layer <- the.layer+1; 
			if(tile.data[a.row, 3]>end) 
			{ end <- tile.data[a.row, 3]; }
		}

		tile.layers[a.row] <- the.layer;
	}
	num.layers <- max(tile.layers);
	if(num.layers>RCircos.Par$max.layers) 
	{ RCircos.Par$track.height <- RCircos.Par$track.height/RCircos.Par$max.layers*num.layers; }
	layer.height <- RCircos.Par$track.height/num.layers;


	#	Start plotting
	#	***********************************************************
	RCircos.Track.Outline(cyto.band, base.positions, out.pos, in.pos, 
			num.layers, RCircos.Par);

	the.loc <- ncol(tile.data);
	for(a.row in 1:nrow(tile.data))
	{
		tile.len <- tile.data[a.row,3] - tile.data[a.row,2];
		tile.range <- round(tile.len/RCircos.Par$base.per.unit/2, digits=0);
		start <- tile.data[a.row,the.loc] - tile.range
		end   <- tile.data[a.row,the.loc] + tile.range;

		layer.bot <- in.pos + layer.height*(tile.layers[a.row]-1);
		layer.top <- layer.bot + layer.height*0.8;

		polygon.x<- c(base.positions[start:end,1]*layer.top, 
			base.positions[end:start,1]*layer.bot);
		polygon.y<- c(base.positions[start:end,2]*layer.top, 
			base.positions[end:start,2]*layer.bot);
		polygon(polygon.x, polygon.y, col="black");
	}
}



#  ======================================================================
#	Draw outline of one plot track. Internal use for histogram plot,
#	scatterplot, line plot, and tile plot.
#  ======================================================================
RCircos.Track.Outline<-function(cyto.band, base.positions, out.pos, in.pos, 
		num.subtrack,  RCircos.Par)
{
	subtrack.height <- RCircos.Par$track.height/num.subtrack;

	chroms <- unique(cyto.band$Chromosome);
	for(a.chr in 1:length(chroms))
	{
		the.chr  <- cyto.band[cyto.band[,1]==chroms[a.chr],];
		start <- the.chr$Location[1]- the.chr$Unit[1] + 1;
		end   <- the.chr$Location[nrow(the.chr)];

		polygon.x<- c(base.positions[start:end,1]*out.pos, 
			base.positions[end:start,1]*in.pos);
		polygon.y<- c(base.positions[start:end,2]*out.pos, 
			base.positions[end:start,2]*in.pos);
		polygon(polygon.x, polygon.y, col="wheat");

		for(a.line in 1:(num.subtrack-1))
		{
			height <- out.pos-a.line*subtrack.height;
			lines(base.positions[start:end,1]*height, 
				base.positions[start:end,2]*height,
				col="gray");
		}
	}
}

