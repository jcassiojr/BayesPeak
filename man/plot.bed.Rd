\name{plot.bed}
\alias{plot.bed}
\title{Plot bed file.}

\description{
Plot the distribution of reads in a .bed file.
}

\usage{
plot.bed(x, chr, start, end, strand = "+", bin = 50L, ...)
}


\arguments{
	\item{x}{
	RangedData (from the IRanges package) with a value column entitled \code{"strand"}. The .bed file to plot. This could have been read in with \code{\link{read.bed}}, or alternatively by the \code{import} function in the rtracklayer library.
	}

	\item{chr}{
	Character. Chromosome to plot - should be exactly the same as a space in \code{bed}.
	}

	\item{start, end}{
	Integer. Start and end locations on chromosome.
	}

	\item{strand}{
	Character. Strand to plot - usually either "+" or "-". If the \code{bed} argument uses a different convention in its strand column, this can be used instead.
	}

	\item{bin}{
	Integer. The width of bin to plot.
	}

	\item{...}{
	Additional arguments, passed to \code{\link{hist}}.
	}
}


\value{
	Plots a histogram on the active graphical device.
}


\details{
	This function takes the reads in a bed file, and plots a histogram of their locations. This allows us to see the peaks present in the data.
}

\author{
	Jonathan Cairns
}


\references{
	Spyrou C, Stark R, Lynch AG, Tavare S
	BayesPeak: Bayesian analysis of ChIP-seq data, BMC Bioinformatics 2009, 10:299 doi:10.1186/1471-2105-10-299
}


\seealso{
	\code{\link{read.bed}}
}


\examples{

dir <- system.file("extdata", package="BayesPeak")
treatment <- file.path(dir, "H3K4me3reduced.bed")
bed <- read.bed(treatment)

##look at region chr16:91000000-94000000
plot.bed(bed, "chr16", 9.1E7, 9.4E7)
}
%\keyword{}
