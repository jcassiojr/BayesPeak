
##-----------------------------------
##autocovariance - simple estimate of R(1)

autocov <- function(x)
{
	n <- length(x)
	x.shift = x[-1]
	mu = mean(x)
	(sum((x[-n] - mu)*(x.shift - mu)))/(n-1)
}

##-----------------------------------
##seqcover - partition of a region 

seqcover <- function(from, to, by, offset = 0L)
{
	if(offset > by) {stop("Bad seqcover arguments: offset > by.")}
	output <- seq(from = from, to = to, by = by) + offset
	output <- output[from <= output & output <= to]
	if (output[length(output)] < to) {output <- c(output, to)} ##cap at beginning
	if (output[1] > from) {output <- c(from, output)} ##cap at end
	output
}

##-----------------------------------
##bed file manipulation

read.bed <- function(filename, chr)
{
	if(missing(chr)) {chr <- NULL}
	##check format of bed files and prepare to read them in
	temp <- readLines(filename, n = 2)
	skiplines <- as.numeric(substring(temp[1],1,5) == "track") ##should we skip a line?
	temp <- length(strsplit(temp[2], "\\s+")[[1]]) ##count columns
	if(temp < 6) {stop(paste(".bed file has only ", temp, " columns - need 6"))}
	colC <- rep("NULL", temp)
	colC[1:6] <- c("character", 0L, 0L, "NULL", "NULL", "character")

	##read in file using scan
	colC <- as.list(colC)
	for(j in which(colC == "NULL")) {colC[j] <- list(NULL)}
	colC[[2]] <- colC[[3]] <- 0L ##has converted to character? TODO
	bed <- scan(filename, what = colC, skip = skiplines)
	gc()
	bed <- bed[!sapply(bed, is.null)] ##remove nulls
	gc()
	#bed <- do.call(data.frame, bed)

	names(bed) <- c("chr", "start", "end", "strand")

	if(is.null(chr)) ##If chr was unspecified...
	{
		chr <- unique(bed$chr) ##get a list of chromosomes
		#FIXME filter
		message("Chromosomes found:\n")
		print(chr)

	} else {
		sel <- bed$chr %in% chr ##filter by chr
		if(sum(sel) == 0) ##specified non-existant chr!
		{
			chr <- unique(bed$chr)
			stop("Specified chromosome(s) not in .bed file.")
		}
		bed <- lapply(bed, function(x){x <- x[sel]})
	}
	gc()
	IR <- IRanges::IRanges(start = bed$start, end = bed$end)
	gc()
	IRanges::RangedData(IR, space = bed$chr, strand = bed$strand)
}

strand.split <- function(bed) ##data frame, headings "chr", "start", "end", "strand"
{
	if(is.null(bed)) {return(NULL)}
	sel <- bed$strand == "+"

	if(class(bed) == "data.frame")
	{
		strand <- list("+" <- bed[sel, c("chr","start")], 
		               "-" <- bed[!sel, c("chr","end")])

		colnames(strand$"+") <- c("chr", "x")
		colnames(strand$"-") <- c("chr", "x")
	} else if (class(bed) == "RangedData")
	{
		strand <- list("+" = data.frame(chr = space(bed)[sel], x = start(bed)[sel]),
		               "-" = data.frame(chr = space(bed)[!sel], x = end(bed)[!sel]))
	}

	strand ##list "+" "-"
}


bin.strand <- function(strand, chr, region = NULL, bin.size = 100L) ## bed is: $chr, $x
{
	sel <- strand$chr == chr
	x.sel <- strand[sel,]

	if((is.null(region)))
	{
		region = c(0L, max(x.sel$x))
	} else {
		sel <- (x.sel$x <= region[2])
		if(region[1] > 0) {sel <- (sel & x.sel$x >= region[1])}
	}

	breaks <- seqcover(from = region[1], to = region[2], by = bin.size, offset = 0L)
	breaks.offset <- seqcover(from = region[1], to = region[2], by = bin.size, offset = floor(bin.size/2))

	##binning
	bin <- list(norm = hist(unique(x.sel$x[sel]), breaks, plot = FALSE)$counts,
	            off = hist(unique(x.sel$x[sel]), breaks.offset, plot = FALSE)$counts)

	bin ## output is list of 2 vectors
}

##-----------------------------------
##bayespeak - main function

bayespeak <- function(treatment, control, chr = NULL, start, end, bin.size = 100L, iterations = 10000L, repeat.offset = TRUE, into.jobs = TRUE, job.size = 6E6L, job.overlap = 20L, use.multicore = FALSE, mc.cores = getOption("cores"), prior = c(5, 5, 10, 5, 25, 4, 0.5, 5), report.p.samples = TRUE)
{
	if(missing(start)) {start <- NA}
	if(missing(end)) {end <- NA}

	##A little bit of parameter checking:
	if(!any(is.na(start)) && !any(is.na(end))) ##don't check start and end if they need to be estimated
	{
		if(!length(start) == length(end)) {stop("Lengths of start and end do not match.")}
		if(any(end < start + bin.size)) {stop("end too small compared to start.")}
		if(!(is.null(chr))) ##must match start, end with chr (if specified
		{
			if(!(length(start) %in% c(1, length(chr)))) {stop("start length must be 1 or the same as chr")}
			if(!(length(end) %in% c(1, length(chr)))) {stop("stop length must be 1 or the same as chr")}
		}
	}
	if(use.multicore && !("multicore" %in% names(sessionInfo()$otherPkgs)))
	{
		message("\nPackage 'multicore' is not loaded - parallel processing disabled. Please load multicore with library(multicore). (See ?bayespeak for more information.)\n")
		use.multicore = FALSE
	}
	
	##check prior (if it's very likely that a_0, b_0, a_1, b_1 < 0.0001 then bayespeak will hang)
	##we are enforcing E(a_i) > 0.0001, E(b_i) > 0.0001
	sel <- matrix(prior, ncol = 2, byrow = TRUE)
	sel <- apply(sel, 1, prod) < 0.0001
	if(any(sel))
	{
		
		errormsg <- paste(c("a0", "b0", "a1", "b1")[sel], collapse = ", ")
		errormsg <- paste("Bad prior, since", errormsg, "have mean < 0.0001. Did you supply rate instead of scale?") 
		stop(errormsg)
	}	

	if(missing(control))
	{
		control <- NULL
	}

	##TODO convert function to S4 method, then can remove this
	##check inputs and fetch data if necessary
	for(input in c("treatment", "control"))
	{
		inputdata <- get(input) ##fetch named object - what class is it?
		if(class(inputdata) == "character")
		{
			##we got a file location string - replace with the file
			eval(parse(text = paste(input, "<- read.bed(inputdata, chr)")))
		}
		if(class(inputdata) == "data.frame")
		{
			##check we have the required columns
			mynames = names(inputdata)
			if(all(c("chromosome","start","end","strand") %in% mynames))
			{
				stop("Bad input - ", input, " does not have column names c(\"chr\", \"start\", \"end\", \"strand\") - instead it has column names c(\"", paste(names(inputdata), collapse = "\", \""), "\").")
			}
		}
		if(class(inputdata) == "RangedData")
		{
			##check that strand info is present
			if(!("strand" %in% colnames(inputdata)))
			{
				stop("Bad input - ", input, "does not have a strand column.")
			}
		}
	}


	##treatment, control: location OR data frame of case, control
	##chr: chromosome string(s)

	##list of parameters in C code output
	para.names <- c("p","theta","lambda0","lambda1","gamma", "loglhood")
	para.list <- rep(0, length(para.names))
	names(para.list) <- para.names

	##get chromosomes
	if(class(treatment) == "data.frame")
	{
		chr <- levels(treatment$chr)
		if(is.null(chr)) {chr <- unique(treatment$chr)} ##for non-factors
	} else if(class(treatment) == "RangedData")
	{
		chr <- unique(space(treatment))
	}

	##check start/end parameters were appropriate
	if(!(length(start) %in% c(1, length(chr)))) {stop("start has length", length(start),"- must be either 1 or equal to no. of chromosomes (", length(chr) ,")")}
	if(!(length(end) %in% c(1, length(chr)))) {stop("end has length", length(end),"- must be either 1 or equal to no. of chromosomes (", length(chr) ,")")}

	##begin
	##split into strands
	treatment <- strand.split(treatment) #$"+", $"-"
	gc()
	control <- strand.split(control)
	gc()

	##prepare for data and QC info
	QC <- data.frame()
	peaks <- list()
	p.samples <- list()

	##collect start and end co-ordinates (TODO Take centromere into account)
	jobs <- rep(0L, length(chr))
	region.full <- cbind(start, end, jobs)
	colnames(region.full) = c("start", "end", "jobs") ##information about each chromosome to be used later

	if(any(is.na(start) | is.na(end)))
	{
		for(i in 1:length(chr))
		{
			ch = chr[i]
			if(any(is.na(start))) ##autofit start location
			{
#				temp <- c(sapply(treatment, function(y){min(y$x[y$chr == ch])}), sapply(control, function(y){min(y$x[y$chr == ch])}))
#				region.full[i, 1] <- min(unlist(temp))
				temp <- sapply(treatment, function(y){min(y$x[y$chr == ch])})
				region.full[i, 1] <- min(unlist(temp))
			}
			if(any(is.na(end))) ##autofit end location
			{
#				temp <- c(sapply(treatment, function(y){max(y$x[y$chr == ch])}), sapply(control, function(y){max(y$x[y$chr == ch])}))
#				region.full[i, 2] <- max(unlist(temp))
				temp <- sapply(treatment, function(y){max(y$x[y$chr == ch])})
				region.full[i, 2] <- max(unlist(temp))
			}
		}
	}
	##check for non-absurd autofit co-ordinates
	sel <- (region.full[,1] + bin.size > region.full[,2])
	if(any(sel))
	{
		warning(paste("No data to analyse on chromosomes:", paste(chr[sel], collapse = ", "), "(perhaps start is too large or end is too small?)"))
	}
	chr <- chr[!sel]

	##calculate how many jobs are expected on each chromosome
	region.full[,"jobs"] <- ceiling((region.full[,"end"] - region.full[,"start"]) / job.size)
	if(repeat.offset == TRUE) {region.full[,"jobs"] <- region.full[,"jobs"] * 2}
	

	for(i in 1:length(chr))
	{
		output <- list()
		ch <- chr[i]
		##get region
		region <- region.full[i,1:2]

		##percentage done...
		if(i == 1)
		{
			percentdone <- 0
		} else {
			percentdone <- 100*(sum(region.full[1:(i-1),"jobs"])/sum(region.full[,"jobs"]))
			percentdone <- round(percentdone, 1)
		}
		message("\n[", percentdone, "% done] Starting ", ch, ":", sep = "", appendLF = FALSE)

		message(region[1], "-", region[2], ":", sep = "", appendLF = FALSE)

		##Put bin counts in bin[[]]
		##NB [1:4] = default (ChIP +, ChIP -, ctrl +, ctrl -)
		##   [5:8] = offset by half of bin.size

		##Bin treated sample
		tr <- lapply(treatment, bin.strand, chr = ch, region = region, bin.size = bin.size)
		gc()

		##Bin control sample, but only if we have any!
		if(is.null(control))
		{
			w <- list(norm = rep(0L, length(tr$"+"$norm)),
			          off = rep(0L, length(tr$"+"$off)))

		} else {
			co <- lapply(control, bin.strand, chr = ch, region = region, bin.size = bin.size)

			##calculate w from the control data - summation over "tetris shape"
			w <- list(norm = co$"+"$norm + co$"-"$norm
				+ c(0L, co$"+"$norm[-length(co$"+"$norm)]) + c(co$"-"$norm[-1], 0L),
			          off = co$"+"$off + co$"-"$off
				+ c(0L, co$"+"$off[-length(co$"+"$off)]) + c(co$"-"$off[-1], 0L))
		} 
		gc()


		if(into.jobs)
		{
			jobs <- NULL ##this will be a list of .Call job expressions, will lapply eval() to it.

			##Break up region into e.g. 6Mb jobs. (seqcover takes into account if region needs to be shorter than 6MB)
			n <- length(tr$"+"$norm)
			boundaries <- floor(seqcover(1, n, job.size/bin.size))
			boundaries[length(boundaries)] = boundaries[length(boundaries)] + 1
			
			##assign bins to each job
			bin.start <- boundaries[-length(boundaries)]
			bin.end <- boundaries[-1] - 1
			##jobs need to overlap, so dilate the jobs a bit
			bin.start[-1] <- bin.start[-1] - job.overlap
			bin.end[-length(bin.start)] <- bin.end[-length(bin.start)] + job.overlap
			##but stop them dilating too far...
			bin.end <- pmin(bin.end, rep(n, length(bin.end)))

			##convert bins -> chromosome loci
			job.bins <- paste(bin.start, ":", bin.end, sep = "")
			job.start <- region[1] + bin.size * (bin.start - 1)
			job.end <- region[1] + bin.size * (bin.end)
			n.jobs <- length(job.start)

			##Get var, autocov (we calculate autocorr later)
			##NB: Variance and autocovariance are bad statistics to use to compare jobs, because of difficulty choosing degrees of freedom
			##(highly dependent on which job - no. of peaks, job may contain centromere, etc...)
			##Assuming the d.f. is the same for variance and autocovariance, it (i.e. (n-k)) cancels out when we divide one by the other to get autocorr.
			job.var <- apply(cbind(bin.start, bin.end), 1, function(x) {0.5*(var(tr$"+"$norm[x[1]:x[2]]) + var(tr$"-"$norm[x[1]:x[2]]))} )
			job.autocov <- apply(cbind(bin.start, bin.end), 1, function(x) {0.5*(autocov(tr$"+"$norm[x[1]:x[2]]) + autocov(tr$"-"$norm[x[1]:x[2]]))} )

			##form jobs
			jobs <- paste(
					'.Call("bayespeak",	as.integer(tr$"+"$norm[', job.bins, ']),
					as.integer(tr$"-"$norm[', job.bins ,']),
					as.integer(w$norm[', job.bins, ']),
					as.double(w$norm[', job.bins, ']),
					as.integer(max(w$norm[', job.bins, '])),
					as.integer(', boundaries[-1] - boundaries[-length(boundaries)], '),
					as.character(ch),
					start = as.integer(', job.start, '),
					end = as.integer(', job.end, '),
					as.integer(bin.size),
					as.integer(iterations),
					para = as.double(para.list),
					as.double(prior))',
					sep = ""
				)

			if(repeat.offset)
			{
				##offset peaks
				n <- length(tr$"+"$off)
				boundaries <- floor(seqcover(1, n, job.size/bin.size))
				boundaries[length(boundaries)] = boundaries[length(boundaries)] + 1

				bin.start <- boundaries[-length(boundaries)]
				bin.end <- boundaries[-1] - 1
				##jobs need to overlap, so dilate the jobs a bit
				bin.start[-1] <- bin.start[-1] - job.overlap
				bin.end[-length(bin.start)] <- bin.end[-length(bin.start)] + job.overlap

				job.bins <- paste(bin.start, ":", bin.end, sep = "")
				job.start <- region[1] + bin.size * (bin.start - 1) - floor(bin.size/2)
				job.end <- region[1] + bin.size * (bin.end) - floor(bin.size/2)

				n.jobs <- length(job.start)

				##*append* var, autocov to previous. autocorr is quotient of this
				job.var <- c(job.var, apply(cbind(bin.start, bin.end), 1, function(x) {0.5*(var(tr$"+"$off[x[1]:x[2]]) + var(tr$"-"$off[x[1]:x[2]]))} ))
				job.autocov <- c(job.autocov, apply(cbind(bin.start, bin.end), 1, function(x) {0.5*(autocov(tr$"+"$off[x[1]:x[2]]) + autocov(tr$"-"$off[x[1]:x[2]]))} ))
				job.autocorr <- job.autocov/job.var

				jobs <- c(jobs, paste(
						'.Call("bayespeak",	as.integer(tr$"+"$off[', job.bins, ']),
						as.integer(tr$"-"$off[', job.bins ,']),
						as.integer(w$off[', job.bins, ']),
						as.double(w$off[', job.bins, ']),
						as.integer(max(w$off[', job.bins, '])),
						as.integer(', boundaries[-1] - boundaries[-length(boundaries)], '),
						as.character(ch),
						start = as.integer(', job.start, '),
						end = as.integer(', job.end, '),
						as.integer(bin.size),
						as.integer(iterations),
						para = as.double(para.list),
						as.double(prior))',
						sep = ""
					))
			}

			jobs <- as.list(jobs)

			##go for it
			if(use.multicore)
			{
				output <- multicore::mclapply(jobs, function(x){eval(parse(text = x))}, mc.cores = mc.cores)
			} else {
				output <- lapply(jobs, function(x){eval(parse(text = x))})
			}

		} else {
			##(into.jobs == FALSE)

			##Collect information on variance and autocorrelation
			job.var <- 0.5 * c( var(tr$"+"$norm) + var(tr$"-"$norm), var(tr$"+"$off) + var(tr$"-"$off))
			job.autocov <- 0.5 * c( autocov(tr$"+"$norm) + autocov(tr$"-"$norm), autocov(tr$"+"$off) + autocov(tr$"-"$off))
			job.autocorr <- job.autocov/job.var

			##generate peaks file
			output[[1]] <- .Call("bayespeak",
				as.integer(tr$"+"$norm),	##score2pos
				as.integer(tr$"-"$norm),	##score2neg
				as.integer(w$norm),	##w
				as.double(w$norm),	##w_dbl
				as.integer(max(w$norm)),	##w_max
				as.integer(length(tr$"+"$norm)),	##passtoN
				as.character(ch),	##chr
				start = as.integer(region[1]),	##start
				end = as.integer(region[2]),	##end
				as.integer(bin.size),	##win
				as.integer(iterations),	##iterations
				para = as.double(para.list),
				as.double(prior))

			if(repeat.offset)
			{
				##generate offset peaks file
				output[[2]] <- .Call("bayespeak",
					as.integer(tr$"+"$off),	##score2pos
					as.integer(tr$"-"$off),	##score2neg
					as.integer(w$off),	##w
					as.double(w$off),	##w_dbl
					as.integer(max(w$off)),	##w_max
					as.integer(length(tr$"+"$off)),	##passtoN
					as.character(ch),	##chr
					start = as.integer(region[1] - floor(bin.size/2)),	##start
					end = as.integer(region[2] + ceiling(bin.size/2)),	##end
					as.integer(bin.size),	##win
					as.integer(iterations),	##iterations
					para = as.double(para.list),
					as.double(prior))
			}
		}
		message("Done.")

		##QC: collect parameters
		temp <- cbind(t(sapply(output, function(x){c(x$jobstart, x$jobend, x$para)})))

		if(repeat.offset)
		{
			temp <- data.frame(ch, temp, job.var, job.autocorr, status = rep(c("normal", "offset"), each = nrow(temp)/2))
		} else {
			temp <- data.frame(ch, temp, job.var, job.autocorr, status = "normal")
		}
		QC <- rbind(QC, temp)

		##peaks: collect hits
		temp <- lapply(output, 
		function(x)
			{
				if(length(x$start) > 0) {data.frame(x$chr, x$start, x$end, x$PP)} else {data.frame(NULL)}
			}
		)
		peaks <- c(peaks, as.list(temp))

		##p.samples: collect parameters passed across
		if(report.p.samples)
		{
			p.samples <- c(p.samples, lapply(output, function(x) {cbind(p = x$p, theta = x$theta, a0 = x$a0, b0 = x$b0, lambda0 = x$lambda0, a1 = x$a1, b1 = x$b1, lambda1 = x$lambda1, loglhood = x$loglhood)}))
		}
	}

	##peaks: add job label and names (X)
	for(i in 1:length(peaks))
		{
			if(nrow(peaks[[i]]) > 0)
			{
				peaks[[i]] <- cbind(peaks[[i]], job = i)
				colnames(peaks[[i]]) = c("chr", "start", "end", "PP", "job")
			}
		}
		
	##QC: rename cols, add label by job
	colnames(QC) <- c("chr", "start", "end", para.names, "var", "autocorr", "status")
	QC <- cbind(score = sapply(peaks, function(x) {sum(x$PP > 0.5)/length(x$PP)}), QC)
	QC <- cbind(job = 1:nrow(QC), calls = sapply(peaks, nrow) , QC) ##FIXME START, END info
	if(is.null(control)) {QC[,"gamma"] = 1} ##we had no control - gamma is irrelevant

	##peaks: list => data.frame, and sort (else offset causes problems)
	peaks <- do.call("rbind", peaks)
	for(ch in chr)
	{
		sel <- peaks$chr == ch
		peaks[sel,] <- peaks[sel,][order(peaks$start[sel]),]
	}

	message("\n")
	list(peaks = peaks, QC = QC, p.samples = p.samples, call = match.call())
}

##-----------------------------------
##summarise.peaks

summarise.peaks <- function(x, threshold = 0.5, method = c("lowerbound", "max"), exclude.jobs = NULL)
{
	method <- match.arg(method)

	##check threshold arg is OK
	if(!(length(threshold) %in% c(1, nrow(x$peaks))))
	{
		stop("Bad 'threshold' argument - length of vector given was ", length(threshold),". This should be 1 or equal to the number of jobs in 'x' argument (", nrow(x),")")
	} 

	##force exclude.jobs arg to integers
	if(is.logical(exclude.jobs)){exclude.jobs <- which(exclude.jobs)}

	output <- NULL
	chr <- unique(x$QC$chr)

	##remove excluded jobs from peak list
	x$peaks <- x$peaks[!(x$peaks$job %in% exclude.jobs),]

	##vector for unenriched chromosomes
	unenriched.chr <- character(0)

	for(j in 1:length(chr))
	{
		##get peaks on this chr
		ch = chr[j]
		sel <- x$peaks[,1] == ch
		x.sel <- x$peaks[sel,]

		##threshold
		if(length(threshold) == 1)
		{
			thr <- threshold
		} else {
			thr <- threshold[x.sel$job]
		}
		x.sel <- x.sel[x.sel$PP > thr, ]

		if(is.unsorted(x.sel$start)) {x.sel <- x.sel[order(x.sel),]}##ensure output is sorted

		##if there are duplicates, then only keep bin with largest PP
		dup <- which(duplicated(x.sel$start))
		dup.all <- sort(c(dup - 1, dup))
		dup.start <- x.sel$start[dup.all]
		for(i in unique(dup.start))
		{
			sel <- dup.all[dup.start == i]
			##overwrite all PP values with max, then delete later
			temp <- max(x.sel$PP[sel])
			x.sel$PP[sel] <- temp
		}
		if(length(dup)>0) {x.sel <- x.sel[-dup,]}

		##exit if chromosome has been reduced to no reads
		if(nrow(x.sel) == 0)
		{
			#warning("No enrichment found on ", ch, " at PP > ", threshold)
			unenriched.chr <- c(unenriched.chr, as.character(ch))
			next
		}

		boundaries <- take.union(x.sel[,2:3])

		start <- which(x.sel$start %in% boundaries$start) ##OUTPUT
		end <- which(x.sel$end %in% boundaries$end) ##OUTPUT

		sel <- apply(rbind(cbind(start, end), 0), 1, function(y){y[1]:y[2]}) ##build a list of regions
		sel <- sel[-length(sel)]

		##combine probabilities - choice of technique
		if(method == "max")
		{
			p <- unlist(lapply(sel, function(y){max(x.sel$PP[y])})) ##build PP vals

		} else if (method == "lowerbound") {

			##a non overlapping set has probability 1 - product (1-P). Find the set with the highest lower bound via dynamic programming.
			p <- rep(0, length(sel))
			halfwidth <- (x.sel$end[1] - x.sel$start[1])/2

			for(i in 1:length(sel)) ##in each region...
			{
				temp <- x.sel[sel[[i]],c("start","PP")]
				temp$start = floor((temp$start - temp$start[1])/halfwidth + 1)
				Q = rep(1, max(temp$start)) ##can fill gaps with Q = 1-P = 1
				Q[temp$start] = 1 - temp$PP

				best <- c(1,1)
				for(j in 1:length(Q))
				{
					slot <- (j %% 2) + 1
					best[slot] <- min(best[3-slot], Q[j]*best[slot])
				}
				p[i] = 1 - min(best)
			}

		} else {
			stop(paste("method = '", method, "' not supported", sep = ""))
		}

		temp <- data.frame(chr = ch, start = boundaries$start, end = boundaries$end, PP = p)
		output <- rbind(output, temp) ##append
	}

	if(length(unenriched.chr) > 0) {warning("No enrichment found on at PP > ", threshold , " on the following chromosomes: ", paste(unenriched.chr, collapse = ", "))}

	##output
	if(is.null(output)) {return(NULL)}
	IRanges::RangedData(IRanges(start = output$start, end = output$end), PP = output$PP, space = output$chr)
}

##-----------------------------------
##take.union - take the union of a matrix of intervals (ie on a single chromosome)

take.union <- function(input)
{
	##if(nrow(input) == 0) {return(list())}
	##put all starts and ends in order
	starts <- cbind(input[,1], 1)
	ends <- cbind(input[,2], -1)
	full <- rbind(starts, ends) ##NB when a start and an end have same value, the start comes first (important)
	full <- full[order(full[,1]),]

	##taking union: we are therefore interested in any region where an interval begins, until all intervals close again.
	endpoints <- cumsum(full[,2]) == 0
	startpoints <- c(TRUE, endpoints[-length(endpoints)])

	##output start & endpoints
	list(start = full[startpoints,1], end = full[endpoints,1])
}

###-----------------------------------
###show results

#show.bayespeak <- function(x, job = , )
#{
#	
#}

