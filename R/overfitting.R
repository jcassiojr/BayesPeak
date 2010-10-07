##function used to aid region.overfitdiag
IsInPolygon <- function(x, y, px, py)
{
	#px, py are vectors containing the x and y co-ords of the polygon

	n <- length(px)
	stopifnot(n == length(py))
	stopifnot(n > 2)
	stopifnot(length(x) == length(y))
	stopifnot(is.numeric(x))
	stopifnot(is.numeric(y))

	nodes <- rep(FALSE, length(x))

	for(i in 1:n)
	{
		j <- i%%n + 1

		##which points are on the LHS of this line, and have y values in the same range? Flip truth of each one.
		nodes <- xor(nodes, (px[i] + (px[j] - px[i])*abs((y-py[i])/(py[j]-py[i])) > x & max(py[i], py[j]) > y & min(py[i], py[j]) <= y))
	}
	return(nodes)
}


##Plot the overfitting diagnostic (e.g. plot each job's lambda1 and score, or calls)

plot.overfitdiag <- function(x, whatX = "lambda1", whatY = "score", logX = TRUE, logY = FALSE, main = "Overfitting diagnostic", ...)
{
	raw <- x

	x <- raw$QC[[whatX]]
	y <- raw$QC[[whatY]]

	if(logX == TRUE) x = log(x)
	if(logY == TRUE) y = log(y)

	plot(x,y, main=main, xlab = ifelse(logX, paste("log(",whatX,")", sep = ""), whatX),
						ylab = ifelse(logY, paste("log(",whatY,")", sep = ""), whatY), ...)
}

identify.overfitdiag <- function(x, whatX = "lambda1", whatY = "score", logX = TRUE, logY = FALSE, main = "Overfitting diagnostic", ...)
{
	raw <- x

	x <- raw$QC[[whatX]]
	y <- raw$QC[[whatY]]

	if(logX == TRUE) x = log(x)
	if(logY == TRUE) y = log(y)

	identify(x,y, ...)
}


region.overfitdiag <- function(x, whatX = "lambda1", whatY = "score", logX = TRUE, logY = FALSE, main = "Overfitting diagnostic", ...)
{
	##How to use this function:
	##Left click to define the vertices of a polygon. To close polygon off, right click.
	##Red hatched area shows the selection
	##Output can be fed into summarise.peaks as exclude.jobs parameter

	raw <- x

	vertices = NULL
	options(locatorBell = FALSE)

	plot.overfitdiag(raw, whatX = whatX, whatY = whatY, logX = logX, logY = logY, main = main, ...)
	temp <- NULL
	temp <- locator(512, type = "o", pch = 20, col = "red", lty = 1)
	if(length(temp$x) == 0)
		{	
			temp$x <- 0
			temp$y <- 0
		}
	vertices <- temp

	if(length(vertices$x) > 2)
	{
		#rescale from imageplot to (x,y) co-ords
		vertices$x <- vertices$x
		vertices$y <- vertices$y

		x <- raw$QC[[whatX]]
		y <- raw$QC[[whatY]]

		if(logX == TRUE) x = log(x)
		if(logY == TRUE) y = log(y)

		##find selected beads using IsInPolygon
		sel <- IsInPolygon(x, y, vertices$x, vertices$y)
		sel[is.na(sel)] <- FALSE

		#display
		polygon(vertices$x, vertices$y, density = c(10,20), col = "Red")
		return(which(sel))
	}

	else{message("Not a polygon - NULL returned.")}
}

