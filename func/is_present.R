# BACKGROUND POINTS HAVING PRESENCE RECORD
# bg - background points
# pr - presence points

is_present <- function(bg, pr, convert=FALSE, check=TRUE, tol=0) {
	if (isTRUE(convert)) {
		zero <- colMeans(rbind(as.matrix(pr), as.matrix(bg)))
		pr <- geodXY(lat=pr[,2], lon=pr[,1], lat0=zero[2], lon0=zero[1], R=6371)
		bg <- geodXY(lat=bg[,2], lon=bg[,1], lat0=zero[2], lon0=zero[1], R=6371)
	}
	closest <- RANN::nn2(data=bg, query=pr, k=1, searchtype="standard", treetype="kd", radius=0, eps=0)
	if (isTRUE(check)) {
		missing <- closest$nn.dists > tol
	} else {
		missing <- rep(FALSE, length(closest$nn.dists))
	}	
	record <- seq(nrow(bg)) %in% closest$nn.idx[!missing]
	attr(record, "missing") <- which(missing)
	return(record)
}


# SPATIAL NEAREST NEIGHBORS
# based on: https://www.gis-blog.com/nearest-neighbour-search-for-spatial-points-in-r/
# dependency: RANN

spatialnn <- function(query, reference, k, searchtype=c("standard", "priority", "radius"), treetype=c("kd", "bd"), radius=0, eps=0, convert=FALSE, R=6371) {
	if (isTRUE(convert)) {
		zero <- colMeans(rbind(as.matrix(query), as.matrix(reference)))
		query <- geodXY(lat=query[,2], lon=query[,1], lat0=zero[2], lon0=zero[1], R=R)
		reference <- geodXY(lat=reference[,2], lon=reference[,1], lat0=zero[2], lon0=zero[1], R=R)
	}
	closest <- RANN::nn2(data=reference, query=query, k=k, searchtype=searchtype[1], treetype=treetype[1], radius=radius, eps=eps)
	return(closest$nn.idx)
}
