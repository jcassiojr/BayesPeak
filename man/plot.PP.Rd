\name{plot.PP}
\alias{plot.PP}
\title{Job PP plot.}

\description{
  Plot the PP profile of a particular job.
}

\usage{
plot.PP(x, job, breaks = 150L, ...)
}


\arguments{
	\item{x}{
		Raw output from the \code{\link{bayespeak}} function.
	}

	\item{job}{
		Integer. Which of the jobs in the output should have its PP values plotted?
	}

	\item{breaks}{
		Integer. Analogous to the breaks argument in \code{\link{hist}}.
	}

	\item{...}{
		Additional arguments passed to the \code{\link{hist}} function.
	}
}


\value{
	Plots a histogram on the active graphical device.
}


\details{
	\code{plot.PP} plots a histogram of the PP values returned in a particular BayesPeak job. This can be used to identify overfitting in a particular job. It is not suitable for identifying the prevalence of overfitting in all of the jobs in a genome-wide analysis - for that, please see \code{\link{plot.overfitdiag}}.
}

\author{
	Jonathan Cairns
}


\references{
	Spyrou C, Stark R, Lynch AG, Tavare S
	BayesPeak: Bayesian analysis of ChIP-seq data, BMC Bioinformatics 2009, 10:299 doi:10.1186/1471-2105-10-299
}


\seealso{
	\code{\link{bayespeak}}.
}


\examples{

##recreation of the plots in the vignette

data(raw.output) ##output from bayespeak()
plot.PP(raw.output, 324, ylim = c(0,50))
plot.PP(raw.output, 325, ylim = c(0,50))

}
%\keyword{}
