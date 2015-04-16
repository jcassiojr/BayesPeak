\name{read.bed}
\alias{read.bed}
\title{BayesPeak - Bayesian analysis of ChIP-seq data}

\description{
  Read a .bed file into a data frame, but only the chr, start, end and strand columns.
}

\usage{
read.bed(filename, chr)
}


\arguments{
	\item{filename}{
		Character - The path to the .bed file in question.
	}

	\item{chr}{
		Character vector, specifying which chromosomes to read in. Chromosome names must be specified exactly as they appear in the .bed files.

		If \code{chr} is missing, then \code{read.bed} will read in the entire data set.
	}
}


\value{
  A \link{RangedData} object, split into spaces by chromosome. This object has a "strand" data track.

	See the IRanges package vignette for more information.
}


\details{
	The purpose of this function is to extract 4 columns from a bed file: chromosome, start, end and strand. These are assumed to be in columns 1, 2, 3 and 6 respectively.

	If the first line begins with "track" then it will be skipped.

%For example:
%
%\samp{track name="080211-TC1-H3K4me3" visibility=4 itemRgb="On"
%chr16 3016155 3016190 080229-TC1-H3K4me3s_4_sorted.txt1600869 0 - - - 255,0,0
%chr16 3016438 3016473 080229-TC1-H3K4me3s_4_sorted.txt1600870 0 + - - 0,0,255
%...}

	The strand sense is expected to be given as "+"/"-".
}

\author{
	Jonathan Cairns
}


\references{
	UCSC BED format FAQ - http://genome.ucsc.edu/FAQ/FAQformat.html#format1
}


\seealso{
	\code{\link{bayespeak}}.
}


\examples{
dir <- system.file("extdata", package="BayesPeak")
file <- file.path(dir, "H3K4me3reduced.bed")

treatment <- read.bed(file)
treatment
}
%\keyword{}
