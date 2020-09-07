plotMaxent <- function(coords, color, palette, pres, pres.pts=FALSE, method=c("points", "raster"), res, pch=15, cex=1, pch.pr=pch, cex.pr=cex, col.pr=2, bg.pr="transparent", lwd.pr=1, outline, xlim, ylim, axes=FALSE, axs="i", mai, device=c("quartz","x11","pdf"), width=7, height=7, file, scaleparams=NULL, ...) {
	method <- method[1]
	if (!is.null(device)) {
		device <- device[1]
	}
	if (missing(palette)) {
		palette <- attr(color, "palette")
	}
	if (missing(mai)) {
		mai <- rep(0.02, 4)
	}
	if (missing(xlim)) {
		if (!missing(outline)) {
			xlim <- extent(outline)[1:2]
		} else {
			xlim <- range(coords[,1])
		}
	}
	if (missing(ylim)) {
		if (!missing(outline)) {
			ylim <- extent(outline)[3:4]
		} else {
			ylim <- range(coords[,2])
		}
	}
	if (missing(pres)) {
		pres.pts <- FALSE
	} else if (isTRUE(pres.pts) & pch.pr == pch & cex.pr == cex) {
		pres.pts <- FALSE
		color[pres] <- col.pr
	}

	if (method == "raster") {
		if (!missing(outline)) {
			crs <- crs(outline)
		} else {
			crs <- NA
		}
		if (missing(res)) {
			res <- c(min(diff(sort(unique(coords[,1])))), min(diff(sort(unique(coords[,2])))))
			res <- ifelse(res < 1, 1 / round(1 / res), res)
		} else {
			res <- rep_len(res, 2)
		}
		ext <- extent(t(apply(coords, 2, range) + rbind(-res, res) / 2))
		ras <- raster(ext=ext, resolution=res, vals=NULL, crs=crs)
		val <- as.numeric(RANN::nn2(data=coordinates(ras), query=coords, k=1, searchtype="standard", treetype="kd", radius=0, eps=0)$nn.idx)
		values(ras)[val] <- as.integer(match(color, palette))
		
		ras <- list(x=sort(unique(coordinates(ras)[,1])), y=sort(unique(coordinates(ras)[,2])), z=as.matrix(t(ras))[,nrow(ras):1])
	}
	
	if (!is.null(device)) {
		if (device == "pdf") {
			match.fun(device)(file=file, width=width, height=height)
		} else if (tolower(device) %in% c("quartz", "x11")) {
			match.fun(device)(width=width, height=height)
		}	
	}
	axt <- c("n", "s")[axes + 1]
	par(mai=mai)
	if (method == "raster") {
		graphics::image(ras, col=palette, xlim=xlim, ylim=ylim, asp=1, bty="n", xaxs=axs, yaxs=axs, xaxt=axt, yaxt=axt)
	} else if (!missing(outline)) {
		plot(outline, xlim=xlim, ylim=ylim, asp=1, ...)
		points(coords, pch=pch, col=color, cex=cex, ...)
	} else {
		plot(coords, xlim=xlim, ylim=ylim, asp=1, pch=pch, col=color, cex=cex, bty="n", xaxs=axs, yaxs=axs, xaxt=axt, yaxt=axt, ...)
	}
	if (!missing(outline)) {
		plot(outline, xlim=xlim, ylim=ylim, add=TRUE, ...)
	}
	if (isTRUE(pres.pts)) {
		points(coords[pres,], pch=pch.pr, lwd=lwd.pr, col=col.pr, bg=bg.pr, cex=cex.pr)
	}
	if (!is.null(scaleparams)) {
		if (!"digits" %in% names(scaleparams)) {
			scaleparams[["digits"]] <- 2
		}
		if (!"cex" %in% names(scaleparams)) {
			scaleparams[["cex"]] <- 1
		}
		if (!"adj" %in% names(scaleparams)) {
			scaleparams[["adj"]] <- 0.5
		}
		if (!"offset" %in% names(scaleparams)) {
			scaleparams[["offset"]] <- 0
		}
		if (!"boundary" %in% names(scaleparams)) {
			scaleparams[["boundary"]] <- NA
		}
		putColorScale(breaks=attr(color, "breaks"), palette=attr(color, "palette"), threshold=attr(color, "threshold"),
			xlim=scaleparams[["xlim"]], ylim=scaleparams[["ylim"]],
			digits=scaleparams[["digits"]], cex=scaleparams[["cex"]], adj=scaleparams[["adj"]], offset=scaleparams[["offset"]], boundary=scaleparams[["boundary"]])
	}
	if (!is.null(device)) {
		if (device == "pdf") {
			dev.off()
		}
	}

}



adjustColorScale <- function(ror, extremes=c("#FFFFFF", "#00A600"), breakpoint="#E6E600", nbg=length(ror), qqplot=TRUE) {
	ror <- sort(ror / sum(ror))

	x <- seq(nbg) / nbg
	y <- cumsum(ror)
	thr <- min(which((nbg * ror) > 1))

	add <- sum(ror[seq(thr - 1)]) / (nbg - thr + 1)
	ref <- ror * (seq(nbg) >= thr) + rep(c(0, add), c(thr - 1, nbg - thr + 1)) / (nbg - thr + 1)
	reratio <- sum((ror * log(nbg * ror))[ror > 0]) / sum((ref * log(nbg * ref))[ref > 0])
	bpratio <- (sum(ror[ror > 1 / nbg]) - 0.5) / 0.5

	mincolor <- grDevices::col2rgb(extremes[1])
	maxcolor <- grDevices::col2rgb(extremes[2])
	bpcolor <- grDevices::col2rgb(breakpoint)	
	midcolor <- 0.5 * (mincolor + maxcolor)
	maxcolor <- midcolor + reratio * (maxcolor - midcolor)
	bpcolor <- midcolor + bpratio * (bpcolor - midcolor)
	mincolor <- grDevices::rgb(mincolor[1], mincolor[2], mincolor[3], max=255)
	maxcolor <- grDevices::rgb(maxcolor[1], maxcolor[2], maxcolor[3], max=255)
	midcolor <- grDevices::rgb(midcolor[1], midcolor[2], midcolor[3], max=255)
	bpcolor <- grDevices::rgb(bpcolor[1], bpcolor[2], bpcolor[3], max=255)

	if (isTRUE(qqplot)) {
		par(mar=c(5.1, 5.1, 2.1, 2.1))
		plot(x, y, xlim=c(0,1), ylim=c(0,1), xlab="Rank", ylab="ROR", cex.lab=1.5)
		points(cbind(c(0, 0.5, 1, 1), c(0, 0.5, 1, 0)), pch=21, cex=3, bg=c(extremes[1], midcolor, extremes[2], breakpoint))
		points(x[thr], y[thr], pch=21, cex=3, bg=bpcolor)
		rr <- 0.5 + reratio * 0.5
		points(rr, rr, pch=21, cex=3, bg=maxcolor)
		arrows(0.975, 0.025, x[thr] + 0.025, y[thr] - 0.025, length=0.15)
		arrows(0.975, 0.975, rr + 0.025, rr + 0.025, length=0.15)
	}

	return(list(extremes=c(mincolor, maxcolor), breakpoint=bpcolor, H=reratio, C=bpratio))
}


makeColorScale <- function(pred, ncat=20, extremes, breakpoint, type=c("ror", "prob"), normalized=TRUE, nbg=length(pred)) {
	if (missing(extremes)) {
		extremes <- c("#F2F2F2", "#00A600")
	}
	if (missing(extremes)) {
		breakpoint <- "#E6E600"
	}
	colors <- c(extremes[1], breakpoint, extremes[2])
	
	type <- type[1]
	if (type == "prob") {
		thr <- 0.5
		upp <- 1
	} else if (normalized == TRUE) {
		thr <- 1
		upp <- nbg
	} else {
		thr <- 1 / nbg
		upp <- 1
	}

	qprior <- sum(pred <= thr) / length(pred)
	below <- round(ncat * qprior)
	above <- ncat - below
	below <- seq(quantile(log(pred), 0.05), quantile(log(pred), qprior), length=below + 1)[-(below + 1)]
	above <- seq(quantile(log(pred), qprior), quantile(log(pred), 0.95), length=above + 1)[-1]
	breaks <- c(below, log(thr), above)
	breaks[1] <- min(log(pred))
	breaks[ncat + 1] <- log(upp) 
	breaks <- exp(breaks)
	breaks[1] <- 0
	less <- grDevices::colorRampPalette(colors[1:2])
	more <- grDevices::colorRampPalette(colors[2:3])
	palette <- c(less(length(below)), more(length(above)+1)[-1])
	col <- palette[.bincode(pred, breaks, include.lowest=TRUE)]
	attr(col, "breaks") <- breaks
	attr(col, "palette") <- palette
	attr(col, "threshold") <- thr

	return(col)
}


applyColorScale <- function(pred, color) {
	col <- attr(color, "palette")[.bincode(pred, attr(color, "breaks"), include.lowest=TRUE)]
	attr(col, "breaks") <- attr(color, "breaks")
	attr(col, "palette") <- attr(color, "palette")
	attr(col, "threshold") <- attr(color, "threshold")
	return(col)
}



putColorScale <- function(breaks, palette, threshold, xlim, ylim, digits=2, cex=1, adj=0.5, offset=0, boundary=NA) {
	ncat <- length(breaks) - 1
	vals <- which(breaks == 1)
	vals <- ncat - c(rev(seq(vals, ncat, by=2)), seq(vals, 2, by=-2)[-1]) + 2
	breakvals <- rev(formatC(breaks, digits=digits, format="f"))
	step <- diff(ylim) / ncat
	bar <- cbind(xlim[c(1,1,2,2,1)], ylim[rep(2,5)] - c(0,step,step,0,0))
	polygon(bar, border=NA, col=rev(palette)[1])
	for (p in 2:ncat) {
		bar[,2] <- bar[,2] - step
		polygon(bar, border=NA, col=rev(palette)[p])
		if (p %in% vals) {
			text(xlim[1] - offset, max(bar[,2]), labels=breakvals[p], cex=cex, adj=c(adj,0.5))
		}
	}
	if (!is.na(boundary)) {
		polygon(cbind(xlim[c(1,1,2,2,1)], ylim[c(1,2,2,1,1)]), border=boundary, col="transparent")
	}
}



ror_weighted_curves <- function(scurves, ror, pres=NULL, minprop=0.05, curve=1, rows=TRUE, main=NULL, xaxis="months", xlab=NULL, ylab=NULL, ylim, col.ror="#00A600", col.pr="gold") {
	curves <- scurves$curves[[curve]][rows,]
	sorder <- scurves$reverse[rows,]
	for (i in seq(nrow(curves))) {
		curves[i,] <- curves[i,][sorder[i,]]
	}
	regimes <- sorder[,1]
	weight <- tapply(ror, regimes, sum)
	included <- weight > minprop
	sorder <- sorder[match(as.numeric(names(weight)[included]), regimes),,drop=FALSE]
	weight <- weight[included]
	ror_curves <- sorder
	for (i in seq(nrow(sorder))) {
		ii <- regimes == sorder[i,1]
		ror_curves[i,] <- colSums(diag(ror[ii] / sum(ror[ii])) %*% curves[ii,])
	}
	if (missing(ylim)) {
		 ylim <- range(curves)
	}
	if (xaxis == "months") {
		 xaxis <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
	}
	# orig. colors: col.pr="darkolivegreen2", col.ror="firebrick2"
	plot(1, xlim=c(0.75, ncol(sorder) + 0.25), ylim=ylim, type="n", xaxt="n", xlab=xlab, ylab=ylab, main=main)
	axis(1, at=1:ncol(sorder), labels=xaxis, las=2)
	for (i in 1:ncol(sorder)) {
		vioplot::vioplot(curves[,i], at=i, col="NA", border="grey", lwd=3, pchMed=NA, drawRect=FALSE, add=TRUE)
	}
	if (!is.null(pres)) {
		for (i in which(pres == 1)) {
			lines(1:ncol(sorder), curves[i,], col=col.pr, lwd=0.5)
		}
	}
	for (i in seq(nrow(ror_curves))) {
		lines(1:ncol(sorder), ror_curves[i,], col=col.ror, lwd=weight[i]*12)
	}
}


