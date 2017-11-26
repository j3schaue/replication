pirateplot_custom <- function (formula, data, line.fun = mean, pal = "appletv", back.col = gray(1), 
    point.cex = 1, point.pch = 1, point.lwd = 1, cut.min = NULL, 
    cut.max = NULL, width.min = 0.3, width.max = NA, bean.o = NULL, 
    point.o = NULL, bar.o = NULL, hdi.o = NULL, line.o = NULL, 
    theme.o = 1, hdi.iter = 1000, jitter.val = 0.03, line.lwd = 4, 
    gl.col = NULL, ylim = NULL, xlim = NULL, xlab = NULL, ylab = NULL, 
    main = NULL, yaxt = NULL, at = NULL, bw = "nrd0", adjust = 1, 
    add = F, y.levels = NULL, bar.border.col = NULL, ...) 
{
    data.2 <- model.frame(formula = formula, data = data)
    dv.name <- names(data.2)[1]
    dv.v <- data.2[, 1]
    iv.levels <- lapply(2:ncol(data.2), FUN = function(x) {
        unique(data.2[, x])
    })
    iv.lengths <- sapply(1:length(iv.levels), FUN = function(x) {
        length(iv.levels[[x]])
    })
    iv.names <- names(data.2)[2:ncol(data.2)]
    n.iv <- length(iv.levels)
    if (is.na(width.max)) {
        if (n.iv == 1) {
            width.max <- 0.45
        }
        if (n.iv == 2) {
            width.max <- 0.5
        }
    }
    bean.mtx <- expand.grid(iv.levels)
    names(bean.mtx) <- names(data.2)[2:ncol(data.2)]
    n.beans <- nrow(bean.mtx)
    bean.mtx$bean.num <- 1:nrow(bean.mtx)
    if (is.null(at)) {
        bean.loc <- 1:n.beans
        group.spacing <- 1
        if (n.iv == 2) {
            bean.loc <- bean.loc + rep(group.spacing * (0:(iv.lengths[2] - 
                1)), each = iv.lengths[1])
        }
    }
    if (!is.null(at)) {
        bean.loc <- rep(at, length.out = n.beans)
    }
    bean.mtx$x.loc <- bean.loc
    data.2 <- merge(data.2, bean.mtx)
    n.cols <- iv.lengths[1]
    if (theme.o == 1) {
        point.o <- ifelse(is.null(point.o), 0.3, point.o)
        bean.o <- ifelse(is.null(bean.o), 0.1, bean.o)
        hdi.o <- ifelse(is.null(hdi.o), 0, hdi.o)
        line.o <- ifelse(is.null(line.o), 0.5, line.o)
        bar.o <- ifelse(is.null(bar.o), 0.5, bar.o)
    }
    if (theme.o == 2) {
        point.o <- ifelse(is.null(point.o), 0.8, point.o)
        bean.o <- ifelse(is.null(bean.o), 0.5, bean.o)
        hdi.o <- ifelse(is.null(hdi.o), 0, hdi.o)
        line.o <- ifelse(is.null(line.o), 0.1, line.o)
        bar.o <- ifelse(is.null(bar.o), 0.1, bar.o)
    }
    if (theme.o == 3) {
        point.o <- ifelse(is.null(point.o), 0.2, point.o)
        bean.o <- ifelse(is.null(bean.o), 0.2, bean.o)
        hdi.o <- ifelse(is.null(hdi.o), 0.8, hdi.o)
        line.o <- ifelse(is.null(line.o), 1, line.o)
        bar.o <- ifelse(is.null(bar.o), 0.1, bar.o)
    }
    if (theme.o == 0) {
        point.o <- ifelse(is.null(point.o), 0, point.o)
        bean.o <- ifelse(is.null(bean.o), 0, bean.o)
        hdi.o <- ifelse(is.null(hdi.o), 0, hdi.o)
        line.o <- ifelse(is.null(line.o), 0, line.o)
        bar.o <- ifelse(is.null(bar.o), 0, bar.o)
    }
    if (mean(pal %in% piratepal(action = "p")) == 1) {
        col.vec <- rep(piratepal(palette = pal, length.out = n.cols))
        point.col <- piratepal(palette = pal, length.out = n.cols, 
            trans = 1 - point.o)
        bean.border.col <- piratepal(palette = pal, length.out = n.cols, 
            trans = 1 - bean.o)
        hdi.band.col <- piratepal(palette = pal, length.out = n.cols, 
            trans = 1 - hdi.o)
        average.line.col <- piratepal(palette = pal, length.out = n.cols, 
            trans = 1 - line.o)
        bar.col <- piratepal(palette = pal, length.out = n.cols, 
            trans = 1 - bar.o)
    }
    if (mean(pal %in% piratepal(action = "p")) != 1) {
        if (length(pal) < n.cols) {
            pal <- rep(pal, n.cols)
        }
        col.vec <- pal
        point.col <- sapply(1:length(pal), function(x) {
            transparent(pal[x], trans.val = 1 - point.o)
        })
        bean.border.col <- sapply(1:length(pal), function(x) {
            transparent(pal[x], trans.val = 1 - bean.o)
        })
        hdi.band.col <- sapply(1:length(pal), function(x) {
            transparent(pal[x], trans.val = 1 - hdi.o)
        })
        average.line.col <- sapply(1:length(pal), function(x) {
            transparent(pal[x], trans.val = 1 - line.o)
        })
        bar.col <- sapply(1:length(pal), function(x) {
            transparent(pal[x], trans.val = 1 - bar.o)
        })
    }
    if (n.iv == 2) {
        col.vec <- rep(col.vec, times = iv.lengths[2])
        point.col <- rep(point.col, times = iv.lengths[2])
        bean.border.col <- rep(bean.border.col, times = iv.lengths[2])
        hdi.band.col <- rep(hdi.band.col, times = iv.lengths[2])
        average.line.col <- rep(average.line.col, times = iv.lengths[2])
        bar.col <- rep(bar.col, times = iv.lengths[2])
    }
    if (is.null(bar.border.col)) {
        bar.border.col <- bar.col
    }
    if (is.null(bar.border.col) == F) {
        bar.border.col <- rep(bar.border.col, times = length(bar.col))
    }
    if (is.null(ylim) == TRUE & is.null(y.levels) == TRUE) {
        steps.p <- c(1, 2, 5, 10, 25, 50, 100, seq(100, 1000, 
            100), seq(1000, 10000, 1000), seq(10000, 1e+05, 10000), 
            seq(1e+05, 1e+06, 1e+05), seq(1e+06, 1e+07, 1e+06), 
            seq(1e+07, 1e+08, 1e+07), seq(1e+08, 1e+09, 1e+08))
        range <- max(dv.v) - min(dv.v)
        steps.p.m <- range/steps.p
        best.step.size <- steps.p[which(abs(steps.p.m - 10) == 
            min(abs(steps.p.m - 10)))]
        plot.min <- floor(min(dv.v)/best.step.size) * best.step.size
        plot.max <- plot.min + 10 * best.step.size
        ylim <- c(plot.min, plot.max)
        y.levels <- seq(plot.min, plot.max, by = best.step.size)
    }
    if (is.null(ylim) == FALSE & is.null(y.levels) == TRUE) {
        steps.p <- c(1, 2, 5, 10, 25, 50, 100, seq(100, 1000, 
            100), seq(1000, 10000, 1000), seq(10000, 1e+05, 10000), 
            seq(1e+05, 1e+06, 1e+05), seq(1e+06, 1e+07, 1e+06), 
            seq(1e+07, 1e+08, 1e+07), seq(1e+08, 1e+09, 1e+08))
        range <- ylim[2] - ylim[1]
        steps.p.m <- range/steps.p
        best.step.size <- steps.p[which(abs(steps.p.m - 10) == 
            min(abs(steps.p.m - 10)))]
        plot.min <- floor(ylim[1]/best.step.size) * best.step.size
        plot.max <- ylim[1] + 10 * best.step.size
        y.levels <- seq(ylim[1], ylim[2], by = best.step.size)
    }
    if (is.null(ylim) == TRUE & is.null(y.levels) == FALSE) {
        ylim <- c(min(y.levels), max(y.levels))
    }
    if (is.null(xlim)) {
        xlim <- c(min(bean.loc) - 0.5, max(bean.loc) + 0.5)
    }
    if (is.null(xlab)) {
        xlab <- iv.names[1]
    }
    if (is.null(ylab)) {
        ylab <- dv.name
    }
    if (add == F) {
        #par(mar = c(5, 4, 4, 1) + 0.1)
        plot(1, xlim = xlim, ylim = ylim, type = "n", xaxt = "n", 
            yaxt = "n", xlab = xlab, ylab = ylab, main = main, 
            yaxt = yaxt, ...)
        if (is.null(yaxt)) {
            axis(side = 2, at = y.levels, las = 1, lwd = 1, lwd.ticks = 1)
        }
        rect(-1000, -1000, 10000, 10000, col = back.col)
        if (is.null(gl.col) == F) {
            abline(h = seq(min(y.levels), max(y.levels), length.out = 21), 
                lwd = c(0.75, 0.25), col = gl.col)
        }
    }
    for (bean.i in 1:n.beans) {
        dv.i <- data.2[data.2$bean.num == bean.i, dv.name]
        fun.val <- line.fun(dv.i)
        x.loc.i <- bean.mtx$x.loc[bean.i]
        rect(x.loc.i - width.max, 0, x.loc.i + width.max, fun.val, 
            col = bar.col[bean.i], border = bar.border.col[bean.i], 
            lwd = 1)
        {
            if (length(dv.i) > 5) {
                dens.i <- density(dv.i, bw, adjust)
                dens.y.i <- dens.i$y
                dens.x.i <- dens.i$x
                if (max(dens.y.i) < width.min) {
                  dens.y.i <- dens.y.i/max(dens.y.i) * width.min
                }
                if (max(dens.y.i) > width.max) {
                  dens.y.i <- dens.y.i/max(dens.y.i) * width.max
                }
                dens.x.plot.i <- dens.x.i
                dens.y.plot.i <- dens.y.i
                if (is.null(cut.min) == F) {
                  dens.x.plot.i <- dens.x.i[dens.x.i > cut.min]
                  dens.y.plot.i <- dens.y.i[dens.x.i > cut.min]
                }
                if (is.null(cut.max) == F) {
                  dens.x.plot.i <- dens.x.i[dens.x.i < cut.max]
                  dens.y.plot.i <- dens.y.i[dens.x.i < cut.max]
                }
                polygon(c(x.loc.i - dens.y.plot.i[1:(length(dens.x.plot.i))], 
                  x.loc.i + rev(dens.y.plot.i[1:(length(dens.x.plot.i))])), 
                  c(dens.x.plot.i[1:(length(dens.x.plot.i))], 
                    rev(dens.x.plot.i[1:(length(dens.x.plot.i))])), 
                  col = gray(1, alpha = bean.o), border = bean.border.col[bean.i], 
                  lwd = 2)
            }
        }
        segments(x.loc.i - width.max, fun.val, x.loc.i + width.max, 
            fun.val, col = average.line.col[bean.i], lwd = line.lwd, 
            lend = 3)
        if (hdi.o > 0) {
            hdi.i <- BEST::hdi(BEST::BESTmcmc(dv.i, numSavedSteps = hdi.iter, 
                verbose = F))
            hdi.lb <- hdi.i[1, 1]
            hdi.ub <- hdi.i[2, 1]
            dens.hdi.x <- dens.x.i[dens.x.i >= hdi.lb & dens.x.i <= 
                hdi.ub]
            dens.hdi.y <- dens.y.i[dens.x.i >= hdi.lb & dens.x.i <= 
                hdi.ub]
            band.type <- "wide"
            if (band.type == "constrained") {
                polygon(c(x.loc.i - dens.hdi.y, x.loc.i + rev(dens.hdi.y)), 
                  c(dens.hdi.x, rev(dens.hdi.x)), col = hdi.band.col[bean.i], 
                  border = NA, lwd = 2)
            }
            if (band.type == "wide") {
                rect(x.loc.i - width.max, hdi.lb, x.loc.i + width.max, 
                  hdi.ub, col = hdi.band.col[bean.i], border = NA)
            }
        }
        points(rep(x.loc.i, length(dv.i)) + rnorm(length(dv.i), 
            mean = 0, sd = jitter.val), dv.i, pch = point.pch, 
            col = point.col[bean.i], cex = point.cex, lwd = point.lwd)
    }
    mtext(bean.mtx[, 1], side = 1, at = bean.mtx$x.loc, line = 0.5)
    if (n.iv == 2) {
        text.loc <- (iv.lengths[1] + 1)/2 * (2 * (1:iv.lengths[2]) - 
            1)
        mtext(text = paste(names(bean.mtx)[2], "=", unique(bean.mtx[, 
            2])), side = 1, line = 2, at = text.loc)
    }
	
	return(xlim)
}