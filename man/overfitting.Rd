\name{overfittingDiagnostics}
\alias{plot.overfitdiag}
\alias{identify.overfitdiag}
\alias{region.overfitdiag}


\title{
Overfitting diagnostic functions
}
\description{
Three functions that provide diagnostic plots and tools to mitigate the effects of overfitting.
}
\usage{
plot.overfitdiag(x, whatX = "lambda1", whatY = "score",
	logX = TRUE, logY = FALSE,
	main = "Overfitting diagnostic", ...)
identify.overfitdiag(x, whatX = "lambda1", whatY = "score",
	logX = TRUE, logY = FALSE,
	main = "Overfitting diagnostic", ...)
region.overfitdiag(x, ...)
}

\arguments{
	\item{x}{
		Raw output from the \code{\link{bayespeak}} function.
	}

	\item{whatX, whatY}{
		Character. The quantities to plot on the X and Y axes. Common choices would be \code{"lambda1"}, \code{"score"}, \code{"calls"}. Any choice in \code{names(raw.output$QC)} is, in theory, acceptable (except for \code{"chr"} and \code{"status"}, which do not correspond to numeric quantities).
	}

	\item{logX, logY}{
		Logical. If TRUE, the quantity on the corresponding axis undergoes a log transformation before being plotted.
	}

	\item{main}{
		Title of plot (corresponds to \code{main} argument in \code{\link{plot}} function).
	}

	\item{...}{
		Further arguments.
		\itemize{
			\item\code{plot.overfitdiag} passes these through to \code{\link{plot}}.
			\item\code{identify.overfitdiag} passes these through to \code{identify}.
			\item\code{region.overfitdiag} passes these through to \code{plot.overfitdiag}.
		}
	}
}


\details{
	These three functions are used to investigate the prevalence of overfitting in a data set, and to aid selection of sensible criteria for performing overfitting corrections.

	\code{plot.overfitdiag} provides a scatterplot of the key parameters associated with jobs. Please see section 9 of the vignette for an description of how to interpret this information.

	\code{identify.overfitdiag} is used after a call \code{plot.overfitdiag}, with the same arguments, to find out which job was plotted at a particular location. The interface is operated in the same manner as \code{\link{identify}} - left-click on the plot to label the job closest to that point, and right-click on the plot to end this process.

	\code{region.overfitdiag} is used to define an overfit region on the plot, and return the jobs in the region. The function is used in the same manner as \code{\link{locator}} - left-click on the plot to define the vertices of a polygon, and then right click anywhere to close the polygon (there is no need to left-click on the first vertex again). The area selected will be filled in with red hatching. The function then returns the IDs of the jobs in the hatched area. Typically, this output will be used as an \code{exclude.jobs} argument in \code{\link{summarise.peaks}}
}

\value{
	All three functions output to the active graphical device. In addition, \code{identify.overfitdiag} and \code{region.overfitdiag} return integer vectors corresponding to the jobs selected on the plot.
}

\author{
	Jonathan Cairns
}

\examples{

data(raw.output)

plot.overfitdiag(raw.output)

##recreate figures in vignette
plot.overfitdiag(raw.output, whatX="calls", logX = TRUE, whatY = "lambda1", logY = TRUE)
plot.overfitdiag(raw.output, whatX="calls", logX = TRUE, whatY = "score", logY = TRUE)

\dontrun{

##identify particular jobs in the plot
plot.overfitdiag(raw.output, whatX="calls", logX = TRUE, whatY = "score", logY = TRUE)
identify.overfitdiag(raw.output, whatX="calls", logX = TRUE, whatY = "score", logY = TRUE)

##define an overfit region
##left-click to define the polygon vertices, right-click to close the polygon
sel <- region.overfitdiag(raw.output, whatX="calls", logX = TRUE, whatY = "score", logY = TRUE)
output <- summarise.peaks(raw.output, exclude.jobs = sel)

}

}


