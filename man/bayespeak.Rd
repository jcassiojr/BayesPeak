\name{bayespeak}
\alias{bayespeak}
\title{BayesPeak - Bayesian analysis of ChIP-seq data}

\description{
  BayesPeak - Bayesian analysis of ChIP-seq data. This function divides the genome into jobs, and performs the BayesPeak algorithm on each using a C backend. The jobs can be performed in parallel, using the package \code{multicore}. Results are returned in R.
}

\usage{
bayespeak(treatment, control, chr = NULL, start,
	end, bin.size = 100L, iterations = 10000L,
	repeat.offset = TRUE, into.jobs = TRUE, job.size = 6E6L,
	job.overlap = 20L, use.multicore = FALSE,
	mc.cores = getOption("cores"),
	prior = c(5, 5, 10, 5, 25, 4, 0.5, 5),
	report.p.samples = TRUE)
}


\arguments{
	\item{treatment, control}{
		These arguments should contain the treated ChIP-seq data and the control data, respectively.

		Each of these arguments can be:

			\itemize{
			\item a path to a .bed file (this file will be read in as per \code{\link{read.bed}}).
			\item OR a \code{data.frame}, which should have columns \code{"chr", "start", "end", "strand"}.
			\item OR a \code{\link{RangedData}} object. This object is expected to be split into spaces by chromosome, and should have a data track labelled "strand".
			}

		The \code{control} argument is entirely optional. (Mathematically, leaving this argument out is equivalent to setting gamma = 1 in the model.)

		Strand information is expected to be given as "+" or "-".
	}

	\item{chr}{
		Character vector, specifying which chromosomes to restrict analysis to. Chromosome names must be specified exactly as they appear in the treatment and control arguments.

		If left as the default value \code{chr = NULL}, then BayesPeak will find all chromosomes present in the \code{treatment} file.
	}

	\item{start, end}{
		Numeric. Locations on the chromosome to start and end at, respectively. If unspecified, then the algorithm will start and end at the minimum and maximum reads found in the data, respectively.
	}

	\item{bin.size}{
		Numeric. Reads are collected into bins. This parameter controls the width of each bin. The bin size is related to the mean fragment length in the library being sequenced, and thus a smaller mean fragment may merit a smaller bin size - please see Spyrou et al. (2009) for more information.
	}

	\item{iterations}{
		Numeric. Number of iterations to run the Monte Carlo analysis for.
	}

	\item{repeat.offset}{
		Logical. If \code{TRUE}, the algorithm is run a second time, this time with the bins offset by \code{floor(window/2)}.
	}

	\item{into.jobs}{
		Logical. By default, BayesPeak will divide a large region into smaller jobs and analyse each one separately. To prevent this behaviour, set \code{into.chunks = FALSE}. This may put BayesPeak at increased risk of overflow and underflow issues, and will additionally prevent usage of the parallel processing options.
	}

	\item{job.size}{
		Numeric. The size of the jobs in base pairs, as described above.
	}

	\item{job.overlap}{
		Numeric. Jobs are expanded to overlap each other. This is prevent peaks on the boundary between two jobs being missed. \code{job.overlap} corresponds to the number of bins by which each job is expanded.
	}

	\item{use.multicore}{
		Logical. If \code{use.multicore = TRUE}, then the individual chunks will be processed in parallel, using the \code{multicore} package. The \code{multicore} package must be installed for this feature to be enabled. At time of writing, it can be downloaded from \url{http://www.rforge.net/multicore/index.html}.
	}

	\item{mc.cores}{
		Numeric. The number of cores to be used for parallel processing. This argument is passed directly to the \code{\link[multicore]{mclapply}} function.
	}

	\item{prior}{
		Numeric. A vector, specifying the prior on the hyperparameters as follows. We have lambda_0 ~ gamma(alpha_0, beta_0) and lambda_1 ~ gamma(alpha_1, beta_1). Additionally, we have that alpha_0, alpha_1, beta_0, beta_1 all have gamma priors. This argument should be c(alpha_0 shape, alpha_0 scale, beta_0 shape, beta_1 scale, alpha_1 shape, alpha_1 scale, beta_1 shape, beta_1 scale).
	}

	\item{report.p.samples}{
		Logical. If FALSE, do not collect information required for the parameter samples reported in the output. Thus, output\$p.samples will be an empty list. If this information is not required, setting this parameter to FALSE will reduce memory usage.
	}
}


\value{

	A list of 4 objects:

	\itemize{
		\item peaks: A \code{data.frame} corresponding to the bins that BayesPeak has identified as potentially being enriched. \code{chr, start, end} give the genomic co-ordinates of the bin. \code{PP} refers to the posterior probability of the bin being enriched. \code{job} is the number of the job within which the bin was called, which corresponds to a row in the QC data.frame (see below).

		\item QC: details of each individual job, listed in columns as follows:
			\itemize{
				\item \code{calls} is the number of potentially enriched bins identified in a job (i.e. bins with \code{PP} > 0.01).
				\item \code{score} is simply the proportion of potentially enriched bins with a PP value above 0.5. Intuitively, a larger score is "better", as it indicates that more of the PP values have tended to 0 or 1.
				\item \code{chr, start, end} are the genomic co-ordinates of the job.
				\item We report the average value, across iterations of the algorithm, of the important parameters \code{p, theta, lambda0, lambda1, gamma} and the average log likelihood \code{loglhood}.
				\item \code{var} is the variance of the bin counts.
				\item \code{autocorr} is an estimate of the first order autocorrelation of bin counts.
				\item \code{status} indicates whether the job was normal, or offset by half a bin width.
			}

		\item call: the line of code used to run BayesPeak.

		\item p.samples: A \code{list} of \code{matrix} objects, each containing parameter samples from the MCMC runs. \code{p.samples[[i]]} corresponds to the samples taken in job \code{i}. Samples are taken every 10 iterations, with the first half of the run being discarded, and to avoid using too much memory, not all parameters are given. This output can be used to assess convergence e.g. using the CRAN packages coda or boa. (see the vignette - \code{vignette("BayesPeak")})
	}

		Note that the raw output of this function is not intended to be used directly as results - the output should be summarized using the \code{\link{summarize.peaks}} function before using it in later analysis.

}


\details{
	BayesPeak uses a fully Bayesian hidden Markov model to detect enriched locations in the genome. The structure accommodates the natural features of the Solexa/Illumina sequencing data and allows for overdispersion in the abundance of reads in different regions. Markov chain Monte Carlo algorithms are applied to estimate the posterior distributions of the model parameters, and posterior probabilities are provided for the sites of interest.
}

\author{
	Christiana Spyrou and Jonathan Cairns
}


\references{
	Spyrou C, Stark R, Lynch AG, Tavare S
	BayesPeak: Bayesian analysis of ChIP-seq data, BMC Bioinformatics 2009, 10:299 doi:10.1186/1471-2105-10-299
}


\seealso{
	\code{\link{read.bed}}, \code{\link{summarize.peaks}}.
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

\dontrun{
##analyse all data in file
raw.output.wg <- bayespeak(treatment, input, use.multicore = TRUE)
output <- summarize.peaks(raw.output.wg)
	}
}
%\keyword{}
