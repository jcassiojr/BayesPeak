\name{plot.job}
\alias{plot.job}
\title{Plot the called peaks in a job.}

\description{
  Plot the distribution of reads in a .bed file, with BayesPeak's calls highlighted.
}

\usage{
plot.job(x, raw.out, job, strand = "+", threshold = 0.5,
	xlim = c(0,1), highlight = TRUE, col.un = "grey",
	col.enr = "blue", bin = 50L, ...)
}


\arguments{

	\item{x}{
		RangedData (from the IRanges package) with a value column entitled \code{"strand"}. The .bed file to plot. This could have been read in with \code{\link{read.bed}}, or alternatively by the \code{import} function in the rtracklayer library.
	}

	\item{raw.out}{
		Raw output from the \code{\link{bayespeak}} function.
	}

	\item{job}{
		Integer. The number of the job to plot.
	}

	\item{strand}{
		Character. Strand to plot - usually either "+" or "-". If the \code{bed} argument uses a different convention in its strand column, this can be used instead.
	}

	\item{threshold}{
		Numeric. Bins with a PP higher than this value will be classed as enriched.
	}

	\item{xlim}{
		Numeric vector. This controls which part of the job is plotted. For example, \code{c(0.5,1)} would plot the last half of the job, whereas \code{c(0.25,0.75)} would plot the central half.
	}

	\item{highlight}{
		Logical. FIXME
	}

	\item{col.un}{
		The colour used to plot counts in unenriched bins.
	}

	\item{col.enr}{
		The colour used to plot counts in enriched bins.
	}

	\item{bin}{
		What sized bin should be used? Currently, this value should be the same as the value used in \code{\link{bayespeak}}. The function will behave strangely if this parameter is altered, particularly to pathological values.
	}

	\item{...}{
		Additional arguments to be passed through to \code{\link{hist}}.
	}
}


\value{
	Plots a histogram on the active graphical device.
}


\details{
	Similar to \code{\link{plot.bed}}, \code{plot.job} takes the reads in a bed file, and plots a histogram of their locations - i.e. plots the bin counts. It then goes on to highlight the histogram bins that have been made in the raw.output from \code{\link{bayespeak}}.

	It is worth bearing in mind that BayesPeak takes the information on both strands into account when calling peaks, and therefore judgements based on a one-stranded view of the data should be treated with caution. For a better picture of what is going on, both strands should be viewed simultaneously, as is done in the examples below.
}

\author{
	Jonathan Cairns
}


\references{
	Spyrou C, Stark R, Lynch AG, Tavare S
	BayesPeak: Bayesian analysis of ChIP-seq data, BMC Bioinformatics 2009, 10:299 doi:10.1186/1471-2105-10-299
}


\seealso{
	\code{\link{bayespeak}}, \code{\link{read.bed}}, \code{\link{plot.bed}}.
}


\examples{
##get the ChIP .bed file
dir <- system.file("extdata", package="BayesPeak")
treatment <- file.path(dir, "H3K4me3reduced.bed")
bed <- read.bed(treatment)

##get the corresponding raw.output object
data(raw.output.H3K4me3)

##plot job 1, + and - strand
par(mfrow = c(2,1))
plot.job(bed, raw.output.H3K4me3, 1)
plot.job(bed, raw.output.H3K4me3, 1, "-")

##zoom in for a closer look...
plot.job(bed, raw.output.H3K4me3, 1, xlim = c(0.58,0.6))
plot.job(bed, raw.output.H3K4me3, 1, "-", xlim = c(0.58,0.6))

}
%\keyword{}
