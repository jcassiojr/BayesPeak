\name{summarize.peaks}
\alias{summarise.peaks}
\alias{summarize.peaks}
\title{Summarize Peaks}

\description{
  Summarize Peaks - Combine the potentially enriched bins found by BayesPeak into contiguous peaks, and associate each with a posterior probability.
  \code{summarise.peaks} is an alias for \code{summarize.peaks}
}

\usage{
summarize.peaks(x, threshold = 0.5, method = c("lowerbound", "max"), exclude.jobs = NULL)
summarise.peaks(x, threshold = 0.5, method = c("lowerbound", "max"), exclude.jobs = NULL)
}


\arguments{
	\item{x}{
		Raw output from the function \code{\link{bayespeak}}.
	}

	\item{threshold}{
		Numeric vector. \code{threshold} must have length equal to either 1 or \code{nrow(x$QC)} (i.e. the number of jobs).			
\itemize{
			\item If \code{threshold} is of length 1, then for each job, all bins with a posterior probability (PP) lower than \code{threshold} will be discarded before summarising.
			\item If \code{threshold} is of length \code{nrow(x$QC)}, then jobs are taken to have separate thresholds - in other words, bins in job i will be discarded if they have a PP less than threshold[i]. Note that this behaviour is irrespective of how many jobs are excluded (see the exclude.jobs argument below) - excluded jobs are still assigned a PP threshold, which is essentially ignored.
		}
	}

	\item{method}{
	The method used to combine the posterior probabilities of multiple peaks. Current methods are:

		\itemize{
			\item \code{lowerbound}:
			A lower bound is found for the posterior probability of the region containing a peak. In contiguous regions of moderately high probability, this method should report a fairer value than \code{method = max}.

			Suppose we have a set of \eqn{n} non-intersecting calls within our region, with posterior probabilities \eqn{p_1} to \eqn{p_n} respectively of containing peaks. Then the probability of there being a peak in this region is at least \eqn{1 - \prod_{i=1}^n(1-p_i)}{1 - (1-p_1)...(1-p_n)}. We maximise this over all possible sets of non-intersecting calls within the region. (Usually, this will simply be a choice between exclusively using the offset or the non-offset analyses.)
			\item \code{max}:
			Combined region has probability equal to the maximum posterior probability over all the peaks it contains.

		}
	}

	\item{exclude.jobs}{
		A vector of integers, denoting jobs to be excluded from later analysis. Alternatively, a logical vector (to be passed through \code{which()}).
	}
}


\value{
  A RangedData object corresponding to the peaks called - each range has an associated PP (Posterior Probability) value.
}


%\details{
%	
%}

\author{
	Jonathan Cairns
}


%\references{
%	Spyrou C, Stark R, Lynch AG, Tavare S
%	BayesPeak: Bayesian analysis of ChIP-seq data, BMC Bioinformatics 2009, 10:299 doi:10.1186/1471-2105-10-299
%}


\seealso{
  \code{\link{bayespeak}}.
}


\examples{
dir <- system.file("extdata", package="BayesPeak")
treatment <- file.path(dir, "H3K4me3reduced.bed")
input <- file.path(dir, "Inputreduced.bed")

##look at specific region 92-95Mb on chromosome 16
##(we've used half the number of iterations here to reduce the time this example takes)
raw.output <- bayespeak(treatment, input, chr = "chr16", start = 9.2E7, end = 9.5E7, iterations = 5000L, use.multicore = TRUE)
output <- summarize.peaks(raw.output)
output

##higher threshold
output.ht <- summarize.peaks(raw.output, threshold = 0.9)
output.ht

}

%\keyword{}
