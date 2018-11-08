
plotKinematics <- function(input) {
	
	# load csv file from event generator
	data <- read.csv(input, header=TRUE, sep=",")

	 # set bin range and bin width for histogram
	binning <- (-10:10)/10

	# weighted histogram
	weighted.hist(data$costheta, data$weight, breaks=binning, main="Muon Angular Distribution", xlab="costheta", ylab="dsigma/dcostheta")

}
