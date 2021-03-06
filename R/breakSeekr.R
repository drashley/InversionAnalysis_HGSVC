#' Find breakpoints from deltaWs
#'
#' Find breakpoints from deltaWs by <here a description of how breakpoints are found>.
#'
#' @param deltaWs A \code{\link{GRanges}} object with metadata column "deltaW" generated by \code{\link{deltaWCalculator}}.
#' @param trim Trim the upper and lower deltaWs in the dataset (e.g. 10 calculates SD of all values btwn 10th and 90th percentile).
#' @param peakTh Determines the fraction of the upper peak to use as the threshold (e.g. 0.33 is 1/3 of max(delta)).
#' @param zlim The degree that the value must be above the threshold to be included in the peaklist.
#' @author Ashley Sanders, David Porubsky, Aaron Taudt
#' @export
breakSeekr<- function(deltaWs, trim=10, peakTh=0.33, zlim=3.291) 
{
 
	if (length(deltaWs)==0) {
		message("No windows here.")
		return()
	}
	deltaW <- deltaWs$deltaW
  th<- max(deltaW)*peakTh # calculates the fraction (e.g. 0.33) of the highest deltaW > this will be the th
  
  ## #group the breakpoints to refine at the peak of the peak:
 
  if (mean(deltaW,na.rm=TRUE) <= 1) { # if mean is so low than there is no peak, skip library
		message('no peaks here!')
		return()
	}
  	
  trimVs<- deltaWs[quantile(deltaW, (trim/100) ) < deltaW & deltaW < quantile(deltaW, ((100-trim)/100) )]$deltaW  # trims the outer 5% limit of the deltaWs to calulate the stdev
  if(length(trimVs) <= 1){trimVs <- c(0,0,0)} #'*'#
  
  # calculate zscores, which determine if deltaW is significantly ABOVE the sd
  if (sd(trimVs) == 0){ 
    zscores <- (deltaW - th) / sd(deltaW) 
  }else{
    zscores <- (deltaW - th) / sd(trimVs)
  }  
  
	## Make peak list
	mask <- zscores > zlim
	rlemask <- rle(mask)
	rlegroup <- rlemask
	rlemaskT <- rlemask$values==TRUE
	rlegroup$values[rlemaskT] <- cumsum(rlemask$values[rlemaskT])
	deltaWs$peakGroup <- inverse.rle(rlegroup)

	### refine breakpoints of the grouped peaks in the peakList
	breaks.unrefined <- deltaWs[deltaWs$peakGroup>0]
	peak.lengths <- list()
	breaks <- GRangesList()
	for (group in unique(breaks.unrefined$peakGroup)) {
		peak <- breaks.unrefined[breaks.unrefined$peakGroup == group]
		if (length(peak) > 1){ # ensures peak has more than 1 window above the threshold
			maxpeakdeltaW <- max(peak$deltaW)
			breaks[[as.character(group)]] <- peak[peak$deltaW == max(peak$deltaW)]
			peak.lengths[[as.character(group)]] <- rep(length(peak),length(breaks[[as.character(group)]]))
		}
	}
	breaks <- unlist(breaks)
	breaks$windows <- unlist(peak.lengths)
	if (length(breaks)>0) {
		return(breaks)
	} else {
		message('no peaks here!')
	}

	# simple plot, if desired
	# plot(deltaWs[,1], deltaWs[,3], type="l")
	# lines(peakList[,1], peakList[,3], col='red')
	
}
    
