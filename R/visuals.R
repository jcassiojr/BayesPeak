plot.job <- function(x, raw.out, job, strand = "+", threshold = 0.5, xlim = c(0,1), highlight = TRUE, col.un = "grey", col.enr = "blue", bin = 100L, ...) ##x is from read.bed, out is raw output, job is integer
{
	j.info <- raw.out$QC[job,]

	width <- as.integer(j.info$end - j.info$start)
	start <- as.integer(j.info$start + xlim[1]*width)
	end <- as.integer(j.info$start + xlim[2]*width)
	chr <- as.character(j.info$chr)

	##reduce x to...
	x.sel <- x[chr] ##...correct chromosome
	x.sel <- ranges(x.sel)[[chr]][x.sel$strand == strand]##...correct strand (now just an IRange)
	x.sel <- start(x.sel) ##...only read starts
	x.sel <- x.sel[start < x.sel & x.sel < end] ##...correct region


	peaks.sel <- raw.out$peaks[raw.out$peaks$job == j.info$job & raw.out$peaks$PP > threshold,]

	h <- hist(x.sel, seqcover(start, end, by = bin), plot = FALSE)
	colours <- rep(col.un, length(h$counts))
#	if(highlight)
#	{
		enriched <- (peaks.sel$start + peaks.sel$end)/2
		bins <- cut(enriched, h$breaks)
		colours[bins] <- col.enr
#	}
	plot(h, main = paste("Job " ,j.info$job, ", ", strand, " strand", sep = ""), col = colours, border = colours, xlab = "Base pairs")
	if(highlight)
	{
		h.dup <- h
		h.dup$counts[!(1:length(h$counts) %in% as.numeric(bins))] <- 0
		plot(h.dup, main = paste("Job" ,j.info$job), col = colours, border = colours, add = TRUE)
	}
}

plot.bed <- function(x, chr, start, end, strand = "+", bin = 50L, ...)
{

	##reduce x to...
	x.sel <- x[chr] ##...correct chromosome
	x.sel <- ranges(x.sel)[[chr]][x.sel$strand == strand]##...correct strand (now just an IRange)
	x.sel <- start(x.sel) ##...only read starts
	x.sel <- x.sel[start < x.sel & x.sel < end] ##...correct region

	h <- hist(x.sel, seqcover(start, end, by = bin), ...)

}

plot.PP <- function(x, job, breaks =150L, ...)
{
	bin.width <- x$peaks$end[1] - x$peaks$start[1]

	for(i in job)
	{
		PP <- x$peaks$PP[x$peaks$job == i]

		job.bins <- (x$QC$end[i] - x$QC$start[i])/bin.width

		PP <- c(PP, rep(0, job.bins - length(PP)))

		hist(PP, breaks=breaks, main=paste("PP histogram - Job", i), ...)
	}
}
