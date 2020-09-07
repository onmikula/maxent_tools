### TOOLS
library(raster)
library(maxnet)
source("func/is_present.R")
source("func/plotMaxent.R")
standard <- function(x) x / sum(x)


### DATA
# read in climatic predictors
bioclim <- readRDS("data/helio_ex_bioclim.rds")
# get background
vals <- !is.na(values(bioclim[[1]]))
bg <- round(coordinates(bioclim[[1]])[vals,], 5)
bioclim <- as.data.frame(do.call(cbind, lapply(bioclim, values))[vals,])
# obtain presence records
records <- read.delim("data/heliophobius.txt")
pr <- is_present(bg, records[,3:2])


### MODEL FITTING
# fit MaxEnt
enm <- maxnet::maxnet(p=pr, data=bioclim, f=maxnet.formula(p=pr, data=bioclim, classes="l"))
# predict RORs
ror <- standard(predict(enm, newdata=bioclim, type="exponential", clamp=FALSE))


### PLOTS
# original palette
palette0 <- list(extremes=c("#FFFFFF", "#00A600"), breakpoint="#E6E600")

# adjusted palette
palette1 <- adjustColorScale(ror)

# corresponding colors scales
color0 <- makeColorScale(pred=nrow(bg)*ror, ncat=20, extremes=palette0$extremes, breakpoint=palette0$breakpoint, type="ror", normalized=TRUE)
color1 <- makeColorScale(pred=nrow(bg)*ror, ncat=20, extremes=palette1$extremes, breakpoint=palette1$breakpoint, type="ror", normalized=TRUE)


### FIGURE 1
adjust <- adjustColorScale(ror, extremes=c("#FFFFFF", "#00A600"), breakpoint="#E6E600", qqplot=FALSE)
x <- seq_along(ror) / length(ror)
y <- cumsum(sort(ror))
thr <- min(which((length(ror) * sort(ror)) > 1))
midcolor <- 0.5 * (grDevices::col2rgb("#FFFFFF") + grDevices::col2rgb("#00A600"))
midcolor <- grDevices::rgb(midcolor[1], midcolor[2], midcolor[3], max=255)

mat <- matrix(1, 3, 5)
mat[2,2] <- 2
mat[2,4] <- 3
layout(mat, widths=c(0.1, 1/3, 0.01, 1/3, 0.2233333), heights=c(0.05, 0.3560607, 0.5939394))

par(mar=c(5.1, 5.1, 2.1, 2.1))
plot(x, y, xlim=c(0,1), ylim=c(0,1), xlab="Rank", ylab="ROR", cex.lab=2.5, cex=5/3, cex.axis=5/3)
points(cbind(c(0, 0.5, 1, 1), c(0, 0.5, 1, 0)), pch=21, cex=5, bg=c("#FFFFFF", midcolor, "#00A600", "#E6E600"))
points(x[thr], y[thr], pch=21, cex=5, bg=adjust$breakpoint)
rr <- 0.5 + adjust$H * 0.5
points(rr, rr, pch=21, cex=5, bg=adjust$extremes[2])
arrows(0.975, 0.025, x[thr] + 0.025, y[thr] - 0.025, length=0.15)
arrows(0.975, 0.975, rr + 0.025, rr + 0.025, length=0.15)
plotMaxent(coords=bg, color=color0, pres=pr, pres.pts=FALSE, method="points", res=1/6, cex=0.9, mai=c(0,0,0,0), device=NULL)
plotMaxent(coords=bg, color=color1, pres=pr, pres.pts=FALSE, method="points", res=1/6, cex=0.9, mai=c(0,0,0,0), device=NULL)

