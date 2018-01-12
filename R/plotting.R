#' Plotting function for BreakPointR

#' @param output genomic ranges from the function runBreakpointr or runBreakpointrALL

#' @import ggplot2
#' @import grid
#' @importFrom cowplot plot_grid
#' @author David Porubsky
#' @export

plotBreakpoints <- function(data, plotLibraries = NULL, file=NULL) {

	if (!is.null(plotLibraries)) {
		libs2plot <- plotLibraries
	} else {
		libs2plot <- c(1:length(data$deltas))
	}

	plots <- list()
	for (i in libs2plot) {
		files <- names(data$deltas)
		filename <- basename(files[i])

		message("Preparing plot for ",filename)

		deltas <- data$deltas[[i]]
		breaks <- data$breaks[[i]]
		counts <- data$counts[[i]]
		params <- data$params

		## Calculate used threshold
		#zlim <- params['zlim']

		## Transform coordinates from "chr, start, end" to "genome.start, genome.end
		cum.seqlengths <- cumsum(as.numeric(seqlengths(deltas)))
		cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
		names(cum.seqlengths.0) <- seqlevels(deltas)

		chr.lines <- data.frame( y=cum.seqlengths[-length(cum.seqlengths)] )
		chr.label.pos <- round( cum.seqlengths.0 + (0.5 * seqlengths(deltas) ) )

		transCoord <- function(gr) {
			gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
			gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
			return(gr)
		}

		trans.deltas <- transCoord(deltas)
		trans.breaks <- transCoord(breaks)
		trans.counts <- transCoord(counts)

		dfplot.deltas <- as.data.frame(trans.deltas)
		dfplot.breaks <- as.data.frame(trans.breaks)
		dfplot.counts <- as.data.frame(trans.counts)

		## function to put vector of chromsome lengths into data.frame
		reformat <- function(x) {
			out_list <- list() 

			for ( i in 2:length(x) ) {
				out_list[[i]] <- c(x[i-1], x[i])
			}
			mt <- do.call("rbind",out_list)
			df <- data.frame(mt)
			colnames(df) <- c("start", "end")
			df
		}

		chr.coord <- reformat(c(0,cum.seqlengths))

		## Scaling size of the rectangle to amout of reads in a given region
		scale <- (dfplot.counts[,c('Ws','Cs')] / dfplot.counts$width) * 1000000
	
		## filter regions small regions => hard to see on the plot
		outlier.W <- round(max(stats::runmed(scale$Ws, 3)))
		outlier.C <- round(max(stats::runmed(scale$Cs, 3)))

		set.max.W <- round(max(scale$Ws[scale$Ws < outlier.W]))
		set.max.C <- round(max(scale$Cs[scale$Cs < outlier.C]))	

		scale$Ws[scale$Ws > outlier.W] <- set.max.W
		scale$Cs[scale$Cs > outlier.C] <- set.max.C
		
		## putting scaled and filtered regions into a data frame
		names(scale) <- c('W.scaled', 'C.scaled')
		df.W <- cbind(dfplot.counts, scaled=scale$W.scaled, fill.strand=rep('W', length(scale$W.scaled)) )
		df.C <- cbind(dfplot.counts, scaled=-scale$C.scaled, fill.strand=rep('C', length(scale$W.scaled)) )
		dfplot.counts <- rbind(df.W, df.C)

	
		## make altering colors for every chromosome
		chr.num <- length(levels(dfplot.deltas$seqnames))
		times <- chr.num/2

		if (chr.num %% 2) {
			chr.colors <- c( rep(c('black', 'gray30'), times), 'black')
		} else {
			chr.colors <- rep(c('black', 'gray30'), times)
		}

		#get midpoint values for each breakpoint
		dfplot.breaks$midpoint <- dfplot.breaks$start.genome + ( (dfplot.breaks$end.genome - dfplot.breaks$start.genome) %/% 2)

		my_theme <- theme(
			legend.position="none",
			panel.background=element_blank(),
			panel.border=element_blank(),
			panel.grid.major=element_blank(),
			panel.grid.minor=element_blank(),
			plot.background=element_blank()
		)

		## plotting data	
		ggplt1 <- ggplot(dfplot.deltas) + geom_line(aes_string(y='deltaW', x='start.genome', color='seqnames'), size=0.2) + scale_color_manual(values=chr.colors, guide='none')
		ggplt1 <- ggplt1 + geom_rect(data=dfplot.breaks, aes_string(xmin='start.genome', xmax='end.genome', ymin=0, ymax='Inf'), fill='red', color='red', alpha=0.3, inherit.aes = FALSE)  + xlab(NULL) + ylab("DeltaW") + ggtitle(filename) + scale_x_continuous(expand = c(0,0)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + my_theme

		ggplt2 <- ggplot(dfplot.counts) + geom_rect(aes_string(xmin='start.genome', xmax='end.genome', ymin=0, ymax='scaled', fill='fill.strand')) + scale_fill_manual(values=c("sandybrown","paleturquoise4"))
		ggplt2 <- ggplt2 + geom_linerange(data=chr.lines, aes_string(ymin=-Inf, ymax=Inf, x='y'), col='black')
		ggplt2 <- ggplt2 + geom_point(data=dfplot.breaks, aes(x=midpoint, y=0), size=5, color='red', shape=124, inherit.aes = FALSE) + xlab(NULL) + ylab("Strands") + scale_x_continuous(expand = c(0,0)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + my_theme

		ggplt3 <- ggplot(dfplot.counts) + geom_rect(aes_string(xmin='start.genome', xmax='end.genome', ymin=0, ymax=10, fill='states'))
		ggplt3 <- ggplt3 + geom_linerange(data=chr.lines, aes_string(ymin=0, ymax=10, x='y'), col='black') + xlab("Chromosomes") + ylab("States") + scale_x_continuous(breaks=chr.label.pos, labels=names(chr.label.pos), expand = c(0,0)) + scale_fill_manual(values=c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown",'?'="red")) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + my_theme

		## storing plot for each file in a list
		#theme_set(theme_cowplot(axis.text.x=element_blank(), axis.ticks.x=element_blank()))
		p <- cowplot::plot_grid(ggplt1, ggplt2, ggplt3, ncol=1, align="v", rel_heights = c(3,3,2))
		plots[[length(plots)+1]] <- p
	}

	## printing to a file or returning plot object
	if (!is.null(file)) {
		message("Printing to PDF ",file)

		pdf(file, width=30, height=10)
		bquiet = lapply(plots, print)
		d <- dev.off()
	} else {
		return(plots)
	}

}


################################################################################## some notes TMP
## translating strand values into numbers
#	dfplot.deltas$value <- as.numeric(dfplot.deltas$strand)
#	dfplot.deltas$value[dfplot.deltas$value == 1] <- 10
#	dfplot.deltas$value[dfplot.deltas$value == 2] <- -10

	## subsetting random fraction of original data
#	sub.rows <- nrow(dfplot.deltas) * 0.25
#	sub.dfplot.deltas <- dfplot.deltas[sample(nrow(dfplot.deltas), sub.rows), ]

#	ggplt2 <- ggplot(sub.dfplot.deltas) + geom_linerange(aes_string(ymin=0, ymax='value', x='start.genome', col='strand'), size=0.2, alpha=0.1) + scale_color_manual(values=c("sandybrown","paleturquoise4"), guide='none')
#	ggplt2 <- ggplt2 + geom_linerange(data=chr.lines, aes_string(ymin=-10, ymax=10, x='y'), col='black')
#	ggplt2 <- ggplt2 + geom_point(data=dfplot.breaks, aes(x=midpoint, y=11), size=5, fill='red', shape=25, inherit.aes = FALSE) + xlab("Chromosomes") + ylab("Strand") + scale_x_continuous(breaks=chr.label.pos, labels=names(chr.label.pos), expand = c(0,0)) + my_theme

#grid.newpage()
#pushViewport(viewport(layout = grid.layout(3, 1, heights=c(4,2,1))))
#vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#print(ggplt1, vp = vplayout(1, 1))
#print(ggplt2, vp = vplayout(2, 1)) 
#print(ggplt3, vp = vplayout(3, 1)) 

#grid.arrange(ggplt1,ggplt2,ggplt3, ncol=1)


#+ theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

#minus <- trans.deltas[strand(trans.deltas)=='-']
#plus <- trans.deltas[strand(trans.deltas)=='+']
#cov.minus <- coverage(minus)
#rle.cov.minus <- runValue(cov.minus)
#rle.cov.minus <- unlist(rle.cov.minus, use.names = FALSE)
#ranges.cov.minus <- ranges(cov.minus)
#ranges.cov.minus <- unlist(ranges.cov.minus)
#dfplot.minus <- data.frame(start=start(ranges.cov.minus), end=end(ranges.cov.minus), cov=rle.cov.minus)
#ggplt <- ggplot(dfplot.minus, aes(x=start, y=cov)) + geom_step(col='red')

#cov.plus <- coverage(plus)
