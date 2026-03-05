# color palettes
.pal_kid <- .pal <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

.pal_sid <- unname(pals::trubetskoy())
.pal_div <- c("turquoise4", "turquoise", "grey95", "magenta", "magenta4")

.pal_sub <- c(
    epi="lightgrey",
    mye="turquoise", 
    nkt="gold", 
    str="deeppink", 
    tum="darkslateblue")
.pal_typ <- c(LN="#56B4E9", SO="#009E73", NP="#F0E442", EY="#D55E00", TO="#0072B2")
.pal_gid <- c(m="lightsteelblue", n="khaki", d="indianred", t="darkseagreen")
.pal_pid <- c(pals::tol(12), "grey50")

.pal_ctx <- pals::brewer.set3(10)

# thresholded z-normalization
.z <- \(x, th=2.5) {
    if (is.null(dim(x))) {
        x[x < 0] <- 0
        sd <- sd(x, na.rm=TRUE)
        if (is.na(sd)) sd <- 1
        x <- x-mean(x, na.rm=TRUE)
        if (sd != 0) x <- x/sd
    } else {
        mus <- colMeans(x)
        sds <- colSds(x)
        x <- sweep(x, 2, mus, `-`)
        x <- sweep(x, 2, sds, `/`)
    }
    x[x > +th] <- +th
    x[x < -th] <- -th
    return(x)
}

# upper/lower quantile scaling
.q <- \(x, margin=1, q=0.01) {
    if (length(q) == 1)
        q <- c(q, 1-q)
    if (!is.matrix(x)) {
        qs <- quantile(x, q, na.rm=TRUE)
        x <- (x-qs[1])/diff(qs)
    } else {
        qs <- c(rowQuantiles, colQuantiles)[[margin]]
        qs <- matrix(qs(x, probs=q), ncol=2)
        x <- switch(margin, 
            `1`=(x-qs[, 1])/(qs[, 2]-qs[, 1]), 
            `2`=t((t(x)-qs[, 1])/(qs[, 2]-qs[, 1])))
    }
    x[x < 0] <- 0
    x[x > 1] <- 1
    return(x)
}

# hierarchical clustering
.xo <- \(.) rownames(.)[hclust(dist(.))$order]
.yo <- \(.) colnames(.)[hclust(dist(t(.)))$order]

# px/mm conversion
.px2mm <- \(.) .*0.00012028
.mm2px <- \(.) ./0.00012028

# compute each cell's distance (by default, in mm) 
# to t(op), b(ottom), r(ight), l(eft) FOV borders
.d2b <- \(dat) {
    require(sce, quietly=TRUE)
    se <- is(dat, "SingleCellExperiment")
    df <- if (se) data.frame(colData(dat)) else dat
    xy <- "Center(X|Y)_local_mm"
    xy <- grep(xy, names(df))
    x <- df[[xy[1]]]
    y <- df[[xy[2]]]
    ds <- data.frame(
        l=x-min(x), r=max(x)-x,
        b=y-min(y), t=max(y)-y)
    if (!se) return(as.matrix(ds))
    colData(dat)[names(ds)] <- ds
    return(dat)
}

# principal component regression
.pcr <- \(sce, x, n=NULL) {
    require(SingleCellExperiment)
    y <- reducedDim(sce, "PCA")
    if (!is.null(n)) y <- y[, seq(n)]
    lapply(x, \(.) {
        z <- summary(lm(y ~ sce[[.]]))
        r2 <- sapply(z, \(.) .$adj.r.squared)
        data.frame(x=., pc=seq_along(r2), r2)
    }) |> do.call(what=rbind)
}

# feature selection
.sel <- \(sce, lab="lab", top=100) {
    df <- scran::findMarkers(sce, groups=sce[[lab]], direction="up")
    lapply(df, \(df) rownames(df)[df$Top <= top])
}

# get reference profiles
.pbs <- \(sce, ids, bp) {
    # dependencies
    library(scuttle)
    library(BiocParallel)
    library(SingleCellExperiment)
    # aggregate counts by both
    ids <- colData(sce)[ids]
    pbs <- aggregateAcrossCells(sce, ids, BPPARAM=bp)
    # library size normalization
    sizeFactors(pbs) <- NULL
    pbs <- logNormCounts(pbs, log=FALSE)
    # average across second
    aggregateAcrossCells(pbs, 
        ids=pbs[[names(ids)[1]]], BPPARAM=bp,
        statistics="mean", use.assay.type="normcounts") 
}

# run 'InSituType' (gs = features to use, nk = number of clusters)
.ist <- \(sce, nk, gs=TRUE, pbs=NULL, bkg=TRUE, ns=c(1e4, 2e4, 1e5)) {
    # dependencies
    library(InSituType)
    library(SingleCellExperiment)
    # load counts
    mtx <- counts(sce[gs, ])
    mtx <- as(t(mtx), "dgCMatrix")
    # cohorting based on IF data
    j <- names(cd <- colData(sce))
    i <- grep("^Mean", j, value=TRUE)
    i <- setdiff(i, "Mean.G")
    i <- c("Area", "AspectRatio", i)
    coh <- fastCohorting(as.matrix(cd[i]))
    # background estimation
    neg <- grep("^neg", altExpNames(sce), ignore.case=TRUE, value=TRUE)
    neg <- sce$nCount_negprobes/nrow(altExp(sce, neg))
    # update reference profiles
    pbs <- if (!is.null(pbs)) {
        bkg <- if (bkg) {
            rna <- sce$nCount_RNA
            rna*mean(neg)/mean(rna) 
        }
        updateReferenceProfiles(
            reference_profiles=pbs, reference_sds=NULL,
            counts=mtx, neg=neg, bg=bkg)$updated_profiles
    }
    # clustering
    insitutype(mtx, 
        reference_profiles=pbs,
        update_reference_profiles=FALSE,
        neg=neg, cohort=coh, n_clusts=nk,
        n_chooseclusternumber=ns[1],
        n_benchmark_cells=ns[1],
        n_phase1=ns[1],
        n_phase2=ns[2],
        n_phase3=ns[3])
}

# relabel 'InSituType' clustering results
.jst <- \(x, y, na="xxx") {
    df <- if (is.list(y)) {
        old <- unlist(y)
        int <- sapply(y, length)
        new <- rep.int(names(y), int)
        data.frame(old, new)
    } else y
    i <- x$clust
    j <- df[[2]][match(i, df[[1]])]
    j[.] <- i[. <- is.na(j)]
    j[j == na] <- NA
    names(j) <- names(i)
    x$clust <- j
    
    i <- colnames(x$logliks)
    j <- df[[2]][match(i, df[[1]])]
    j[.] <- i[. <- is.na(j)]
    colnames(x$logliks) <- j
    colnames(x$profiles) <- j
    
    mx <- x$profiles
    ex <- colnames(mx) == na
    x$profiles <- mx[, !ex]
    
    # i <- split(seq_along(j), j)
    # x$logliks <- sapply(i, \(.) rowMeans(x$logliks[, ., drop=FALSE]))
    # x$profiles <- sapply(i, \(.) rowMeans(x$profiles[, ., drop=FALSE]))
    return(x)
}

# spatial smoothing
.smooth_xy <- \(x, a=1, r=0.01, k=101) {
    cd <- SingleCellExperiment::colData(x)
    xy <- as.matrix(cd[grep("global_mm", names(cd))])
    nn <- RANN::nn2(xy, k=k, searchtype="radius", r=r)
    is <- nn$nn.idx; is[is == 0] <- NA
    apply(assay(x, a), 1, \(y) {
        z <- matrix(y[c(is)], nrow(is), ncol(is))
        rowMeans(z, na.rm=TRUE) 
    }) |> t() |> `colnames<-`(colnames(x))
}

# plt ----

# dependencies
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(SingleCellExperiment)
})

# plot dimensionality reduction
.plt_dr <- \(x, c, id=NULL, dr="UMAP", s=0.1, 
    t=c("n", "z", "q"), th=2.5, qs=0.01, a=NULL) {
    library(ggplot2)
    library(ggrastr)
    library(SingleCellExperiment)
    if (is.data.frame(x)) {
        df <- x
    } else {
        f <- \(.) inherits(tryCatch(as.vector(.), error=\(.) .), "error")
        i <- !sapply(as.list(colData(x)), f)
        xy <- colnames(df <- reducedDim(x, dr)[, c(1, 2)])
        df <- data.frame(colData(x)[i], df, check.names=FALSE)
        if (length(c) == 1 && c %in% rownames(x)) {
            if (is.null(a)) a <- "logcounts"
            df[[c]] <- assay(x, a)[c, ]
        }
    }
    df <- df[sample(nrow(df)), ]
    if (is.numeric(df[[c]])) {
        aes <- switch(match.arg(t), 
            n={ # no transformation
                scale_color_gradientn(colors=pals::jet())
            },
            z={ # thresholded z-normalization
                df[[c]] <- .z(df[[c]], th)
                scale_color_gradient2(low="turquoise", mid="lavender", high="deeppink")
            },
            q={ # upper/lower quantile scaling
                df[[c]] <- .q(df[[c]], qs)
                scale_color_gradientn(
                    limits=c(0, 1), 
                    breaks=c(0, 1),
                    colors=pals::jet())
            })
        aes <- list(aes, .thm_xy_c(s))
        df <- df[order(abs(df[[c]]), na.last=FALSE), ]
    } else {
        aes <- list(.thm_xy_d(s), scale_color_manual(values=.pal))
    }
    suppressMessages(
    ggplot(df, aes(.data[[xy[1]]], .data[[xy[2]]], col=.data[[c]])) +
        (if (!is.null(id)) ggtitle(.lab(id, nrow(df)))) + aes + 
        scale_x_continuous(expand=expansion(0.05)) +
        scale_y_continuous(expand=expansion(0.05)) +
        coord_cartesian(expand=TRUE) + theme(
            aspect.ratio=1,
            plot.margin=margin(2,2,2,2),
            panel.border=element_rect(linewidth=0.2, fill=NA, color="grey"))
    )
}

# dotplot of selected genes
.plt_dp <- \(sce, id, gs, hc=FALSE) {
    gs <- lapply(gs, intersect, rownames(sce))
    ns <- sapply(gs, length)
    gt <- unlist(gs)
    mu <- scran::summaryMarkerStats(sce[gt, ], sce[[id]])
    mu <- lapply(mu, \(.) data.frame(g=rownames(.), .))
    ys <- "self.average"
    mu <- mu |>
        bind_rows(.id="k") |>
        group_by(g) |>
        mutate_at(ys, .z)
    if (hc) {
        mx <- mu |>
            select(k, g, ys) |>
            pivot_wider(
                names_from="g",
                values_from=ys)
        my <- as.matrix(mx[, -1])
        rownames(my) <- mx[[1]]
        xo <- .yo(my)
        yo <- .xo(my)
        box <- NULL
    } else {
        xo <- gt
        yo <- rev(names(gs))
        s <- cumsum(ns)
        n <- length(ns)
        bb <- 0.5+data.frame(
            xmin=rev(c(0,s[-n])), xmax=rev(s),
            ymin=seq(n), ymax=c(0,seq(n)[-n]))
        box <- geom_rect(
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), bb,
            inherit.aes=FALSE, fill=NA, col="black", linewidth=0.4)
    }
    mu <- mu |>
        mutate_at("g", factor, xo) |>
        mutate_at("k", factor, yo)
    ggplot(mu, aes(g, k, col=self.average, size=self.detected)) + 
        scale_color_gradient2(low="blue3", mid="grey95", high="red3") +
        scale_size_continuous(breaks=seq(0, 1, .2), range=c(1, 5)) +
        guides(col=guide_colorbar(order=1)) +
        geom_point() + coord_equal() + box +
        .thm_fig_c("minimal") + theme(
            panel.grid=element_blank(),
            axis.title=element_blank(),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
}

# gene x cluster heatmaps including look-up, joint & split markers
.plt_de <- \(x, id="", z=TRUE, gs=NULL, n1=30, n2=80) {
    library(dplyr)
    library(tidyr)
    library(scran)
    library(scuttle)
    library(ggplot2)
    # DGE analysis
    ids <- factor(. <- x$clust, sort(unique(.)))
    nk <- length(names(ks) <- ks <- names(ns <- table(ids)))
    # selection 
    es <- x$profiles
    es <- normalizeCounts(es)
    # duplicates
    nms <- colnames(es)
    dup <- duplicated(nms)
    dup <- unique(nms[dup])
    for (d in dup) {
        i <- colnames(es) == d
        mu <- rowMeans(es[, i])
        es <- cbind(es[, !i], mu)
        colnames(es)[ncol(es)] <- d
    }
    # selection
    fcs <- lapply(ks, \(i) {
        j <- setdiff(ks, i)
        a <- es[, i]
        b <- rowMeans(es[, j])
        c <- (c <- a/b)[!is.na(c) & is.finite(c)]
        c*rowDiffs(rowRanges(es))[names(c), ]
    })
    gs <- intersect(gs, rownames(es))
    sel <- \(n) lapply(ks, \(k) names(tail(sort(fcs[[k]]), n)))
    top <- unique(c(gs, unlist(sel(n1)))); all <- sel(n2)
    # aesthetics
    pal <- if (z) {
        es <- t(apply(es, 1, .z))
        scale_fill_gradientn(
            "z-scaled\nmean expr.", 
            colors=.pal_div,
            limits=c(-2.5, 2.5), breaks=seq(-2, 2, 2)) 
    } else {
        scale_fill_gradientn(
            "mean\nexpr.", limits=c(0, NA), 
            colors=c("navy", "red", "gold"))
    }
    aes <- list(
        coord_equal(4/3, expand=FALSE),
        geom_tile(), theme_bw(6), theme(
            legend.position="bottom",
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            plot.background=element_blank(),
            legend.key=element_blank(),
            legend.background=element_blank(),
            legend.title=element_text(vjust=1.5),
            legend.key.width=unit(0.8, "lines"),
            legend.key.height=unit(0.4, "lines"),
            axis.text.y=element_text(size=8),
            axis.text.x=element_text(size=6, 
                angle=90, hjust=1, vjust=0.5)))
    # hierarchical clustering
    .x <- \(.) rownames(.)[hclust(dist(.))$order]
    .y <- \(.) colnames(.)[hclust(dist(t(.)))$order]
    df <- data.frame(gene=rownames(es), es, check.names=FALSE)
    df <- pivot_longer(df, -gene, names_to="k", values_to="y") 
    # plotting
    p1 <- ggplot(
        df[df$gene %in% top, ],
        aes(gene, k, fill=y)) + 
        scale_x_discrete(limits=.x(es[top, ])) + 
        scale_y_discrete(limits=.y(es[top, ])) +
        ggtitle(.lab(id, length(x$clust))) +
        aes + pal 
    ps <- lapply(ks, \(k) ggplot(
        df[df$gene %in% (gs <- all[[k]]), ], 
        aes(gene, k, fill=y)) + 
            scale_x_discrete(limits=.x(es[gs, ])) + 
            scale_y_discrete(limits=.y(es[gs, ])) +
            ggtitle(.lab(k, ns[[k]])) +
            aes + pal)
    c(list(p1), ps)
}

# df = `arrow::Table` containing cell boundaries
# sce = corresponding 'SingleCellExperiment'
# c = character string; feature name or 'colData' to color points by
# t = "n"(o transformation), "z"(-normalization), or "q"(uantile) scaling
# th = scalar numeric; threshold to use when 't == "z"'
# qs = scalar or length-2 numeric; quantiles to use when 't == "q"'
# hl = logical/character vector; cells to highlight (others are 'blacked out')
.plt_ps <- \(df, sce=NULL, c="white", a=1,
    t=c("n", "z", "q"), th=2.5, qs=0.01, 
    hl=NULL, lw=0.1, lc="lightgrey", id=NULL) {
    library(dplyr)
    library(ggplot2)
    library(SingleCellExperiment)
    # filter for cells present in object
    cs <- pull(df, "cell", as_vector=TRUE)
    cs <- which(cs %in% sce$cell)
    df <- df[cs, ] |>
        mutate(x=.px2mm(x_global_px)) |>
        mutate(y=.px2mm(y_global_px)) 
    # join cell metadata & polygons data
    i <- match(pull(df, "cell", as_vector=TRUE), sce$cell)
    j <- setdiff(names(colData(sce)), names(df))
    df <- cbind(as.data.frame(df), colData(sce)[i, j])
    if (c %in% rownames(sce)) {
        df[[c]] <- logcounts(sce)[c, i]
        # continuous coloring
        pal <- switch(match.arg(t), 
            n={ # no transformation
                scale_fill_gradientn(colors=pals::jet())
            },
            z={ # thresholded z-normalization
                df[[c]] <- .z(df[[c]], th)
                scale_fill_gradientn(
                    colors=pals::coolwarm(),
                    limits=c(-th, th), n.breaks=5)
            },
            q={ # lower/upper quantile scaling
                df[[c]] <- .q(df[[c]], qs)
                scale_fill_gradientn(
                    colors=pals::jet(),
                    limits=c(0, 1), n.breaks=5)
            })
        thm <- list(pal, .thm_fig_c("void"))
    } else if (c %in% names(df)) {
        if (!is.numeric(df[[c]])) {
            # discrete coloring
            n <- nlevels(df[[c]])
            pal <- if (n == 5) .pal_sub else 
                if (length(.pal) >= n) .pal else
                    grDevices::colorRampPalette(.pal)(n) 
            if (is.null(names(pal))) {
                names(pal) <- levels(df[[c]])
                pal <- pal[!is.na(names(pal))]
            }
            pal <- scale_fill_manual(NULL, values=pal, 
                na.value="lightgrey", breaks=names(pal))
            thm <- list(pal, .thm_fig_d("void", "f"))
        } else {
            pal <- scale_fill_gradientn(colors=pals::jet())
            thm <- list(pal, .thm_fig_c("void"))
        }
    } else {
        df[[c <- "foo"]] <- c
        thm <- list(scale_fill_identity(NULL))
    }
    # highlighting
    if (!is.null(hl)) {
        if (is.logical(hl)) 
            hl <- sce$cell[hl]
        df[[c]][!df$cell %in% hl] <- NA
    }
    # plotting
    ggplot(df, aes(x, y, fill=.data[[c]], group=cell)) + 
        ggrastr::rasterize(dpi=300,
        geom_polygon(col=lc, alpha=a, linewidth=lw, 
            show.legend=TRUE, key_glyph="point")) +
        coord_equal(expand=FALSE) + thm +
        if (!is.null(id)) ggtitle(.lab(id, length(unique(df$cell)))) 
}

# spatial plot
.plt_xy <- \(x, k, id="", s=NULL, split=FALSE, na=FALSE) {
    # dependencies
    library(ggplot2)
    library(ggrastr)
    library(SingleCellExperiment)
    # wrangling
    if (is.logical(k)) split <- FALSE
    cd <- if (is.data.frame(x)) x else colData(x)
    xy <- "Center(X|Y)_global_mm"
    xy <- grep(xy, names(cd))
    names(cd)[xy] <- c("x", "y")
    if (length(names(k))) {
        cs <- match(colnames(x), names(k))
        ko <- order(tolower(ks <- unique(k)))
        nk <- length(ks <- levels(k <- factor(k[cs], ks[ko])))
    }
    # aesthetics
    df <- data.frame(cd, k)
    dx <- diff(range(df$x))
    dy <- diff(range(df$y))
    pt <- if (is.null(s)) min(dx, dy)/100/2 else s
    # plotting
    if (!is.numeric(df$k)) {
        fd <- if (na) df else df[!is.na(df$k), ]
        p0 <- ggplot(fd, aes(x, y, col=k)) + .thm_xy_d(pt) +
            scale_color_manual(NULL, drop=FALSE, values=.pal, 
                na.value="lightgrey", breaks=\(.) setdiff(., NA)) +
            ggtitle(.lab(id, sum(!is.na(fd$k))))
        if (!split) return(p0)
        ps <- if (split) lapply(c(setNames(c(ks, NA), c(ks, "nan"))), \(k) {
            df$. <- if (is.na(k)) is.na(df$k) else grepl(sprintf("^%s$", k), df$k)
            ggplot(df[order(df$.), ], aes(x, y, col=.)) + 
                .thm_xy_d(pt) + theme(legend.position="none") +
                scale_color_manual(NULL, values=c(c("lavender", "darkslateblue"))) +
                ggtitle(.lab(k, sum(df$.)))
        })
        c(list(all=p0), ps)
    } else {
        df <- df[order(df$k, na.last=FALSE), ]
        ggplot(df, aes(x, y, col=k)) + .thm_xy_c(pt) +
            scale_color_gradientn(NULL, colors=pals::jet(), na.value="lightgrey") +
            ggtitle(.lab(id, nrow(df)))
    }
}

.plt_rgb <- \(x, id, s=NULL) {
    # dependencies
    library(ggplot2)
    library(ggrastr)
    library(SingleCellExperiment)
    # wrangling
    xy <- "Center(X|Y)_global_mm"
    xy <- grep(xy, names(colData(x)))
    names(colData(x))[xy] <- c("x", "y")   
    y <- reducedDim(x, "PCA")[, seq_len(3)]
    z <- sweep(y, 1, rowMins(y), `-`)
    z <- sweep(z, 1, rowMaxs(z), `/`)
    z <- apply(z, 1, \(.) rgb(.[1], .[2], .[3]))
    df <- data.frame(colData(x), y, z)
    # aesthetics
    dx <- diff(range(df$x))
    dy <- diff(range(df$y))
    pt <- if (is.null(s)) min(dx, dy)/100/2 else s
    # plotting
    p0 <- ggplot(df, aes(x, y, col=z)) + 
        ggtitle(.lab(id, nrow(df))) +
        scale_color_identity() + 
        .thm_xy_d(pt)
    ps <- lapply(colnames(y), \(.) {
        ggplot(df, aes(x, y, col=.q(.data[[.]]))) + 
            scale_color_gradientn(
                colors=pals::jet(), n.breaks=6,
                paste0("q-scaled\n", ., " value")) +
            ggtitle(.lab(id, nrow(df))) + .thm_xy_c(pt)
    })
    c(list(p0), ps)
}

# save series of spatial plots,
# adjusting for dimensions
.pdf <- \(ps, nm, sf=1) {
    tf <- replicate(length(ps), tempfile(fileext=".pdf"), FALSE)
    for (. in seq_along(ps)) {
        df <- ps[[.]]$data
        dx <- sf*diff(range(df$x))
        dy <- sf*diff(range(df$y))
        pdf(tf[[.]], 
            width=(2+dx)/2.54, 
            height=(0.5+dy)/2.54)
        print(ps[[.]]); dev.off()
    }
    qpdf::pdf_combine(unlist(tf), output=nm)
}

# 'flightpath plot', i.e., embedding of 
# 'InSituType' assignment likelihoods
.plt_fp <- \(x, id=NULL, s=4) {
    library(dplyr)
    library(ggplot2)
    library(InSituType)
    ks <- x$clust
    ll <- x$logliks
    ps <- x$profiles
    set.seed(194849)
    f <- \(.) sample(., min(1e4, length(.)))
    i <- unlist(lapply(split(seq(nrow(ll)), ks), f))
    y <- flightpath_layout(logliks=ll[i, ], profiles=ps)
    df <- data.frame(y$cellpos, k=ks[i])[sample(i), ]
    fd <- summarize(group_by(df, k), across(c("x", "y"), median))
    ggplot(df, aes(x, y, col=k)) + 
        scale_color_manual(values=.pal) + 
        .thm_xy_d(0.1) + theme(aspect.ratio=1, legend.position="none") +
        geom_text(data=fd, aes(label=k), size=s, col="black") +
        if (!is.null(id)) ggtitle(.lab(id, sum(!is.na(ks))))
}

# compositional barplot
.plt_fq <- \(z, x, y, id=NULL, by=NULL, hc=TRUE, h=FALSE, a=1, ...) {
    library(ggplot2)
    library(SingleCellExperiment)
    # tabulate cell counts
    if (is(z, "SingleCellExperiment"))
        z <- data.frame(colData(z))
    ns <- table(z[[x]], z[[y]])
    df <- as.data.frame(ns)
    i <- match(df[[1]], z[[x]])
    j <- setdiff(names(z), c(x, y))
    df <- cbind(df, z[i, j])
    if (!is.null(by)) df[[by]] <- z[[by]][i]
    # hierarchical clustering
    xo <- if (hc) {
        ds <- dist(prop.table(ns, 1))
        ds[is.na(ds)] <- 0
        hc <- hclust(ds)
        hc$labels[hc$order]
    } else rownames(ns)
    # plotting
    aes <- if (h) {
        list(
            coord_flip(expand=FALSE), theme(
                axis.text.x=element_blank(), 
                axis.title.x=element_blank()))
    } else {
        list(coord_cartesian(expand=FALSE), theme(
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)))
    }
    ggplot(df, aes(Var1, Freq, fill=Var2)) +
        (if (!is.null(by)) facet_wrap(by, ...)) +
        geom_col(
            position="fill", col="white", alpha=a,
            linewidth=0.1, width=1, key_glyph="point") +
        labs(x=x, y="frequency", fill=y) +
        scale_fill_manual(values=.pal) +
        scale_x_discrete(limits=xo) +
        .thm_fig_d("minimal", "f") + 
        aes + theme(
            panel.grid=element_blank(),
            axis.title=element_blank(),
            axis.ticks=element_blank()) +
        if (!is.null(id)) ggtitle(.lab(id, sum(ns)))
}

# aes ----

# prettified plot title in the style of
# 'title (N = count)' with bold 'title'
.lab <- \(x, n=NULL, m=FALSE) {
    if (is.null(n)) {
        if (is.null(x)) "" else # blanc
            bquote(bold(.(x)))  # 'x' only
    } else {
        n <- format(n, big.mark=",")
        if (is.null(x)) {
            # 'n' only
            bquote("N ="~.(n))
        } else {
            # both
            if (!m) {
                # single line
                bquote(bold(.(x))~"("*.(n)*")")
            } else {
                # w/ linebreak
                bquote(atop(bold(.(x)),"("*.(n)*")"))
            }
        }
    }
}

# base figure theme
.thm_fig <- \(.="minimal") {
    thm <- get(paste0("theme_", .))(9)
    list(thm, theme(
        legend.key=element_blank(),
        plot.background=element_blank(),
        panel.background=element_blank(),
        strip.background=element_blank(),
        panel.grid.minor=element_blank(),
        legend.background=element_blank(),
        plot.title=element_text(hjust=0.5)))
}

# discrete coloring
.thm_fig_d <- \(., l=c("c", "f")) {
    aes <- switch(match.arg(l),
        c=list(alpha=1, shape=19, size=1),
        f=list(alpha=1, shape=21, stroke=0, col=NA, size=2))
    thm <- list(theme(
        legend.key.size=unit(0, "lines")),
        guides(col=guide_legend(ncol=1, override.aes=aes)),
        guides(fill=guide_legend(ncol=1, override.aes=aes)))
    c(.thm_fig(.), list(thm))
}

# continuous coloring
.thm_fig_c <- \(.) {
    thm <- theme(
        legend.key.width=unit(0.2, "lines"),
        legend.key.height=unit(0.4, "lines"))
    c(.thm_fig(.), list(thm))
}

# theme for spatial plots
.thm_xy <- \(s=0.1) list(
    ggrastr::geom_point_rast(show.legend=TRUE,
        shape=16, stroke=0, size=s, raster.dpi=100),
    scale_x_continuous(expand=expansion(0, 0.1)),
    scale_y_continuous(expand=expansion(0, 0.1)),
    coord_equal(), theme(
        plot.margin=margin(),
        plot.title=element_text(hjust=0.5),
        panel.background=element_rect(fill=NA)))
.thm_xy_d <- \(s=0.1) c(.thm_fig_d("void"), .thm_xy(s))
.thm_xy_c <- \(s=0.1) c(.thm_fig_c("void"), .thm_xy(s))

# mgs ----

.mgs <- list(
    nkt=list(
        Thn=c(
            "CD2", "CD3D", "CD3E", "CD3G", "CD4", "CD5", "CD6", "CD28", "CD44", # canonical
            "LEF1", "TCF7", "LTB"), # naivety
        Tcm=c("IL7R", "CCR7", "KLF2", "KLF3", "BACH2", "FOXP1", "IKZF1", "SATB1"),
        Tfh=c("ICOS", "BCL6", "CXCL13", "CD40LG"),
        Treg=c("TNFRSF4", "IL2RA", "FOXP3", "CTLA4"),
        Tcn=c("CD8A", "CD8B"),
        Tem=c("STAT4", "CX3CR1", "JUN", "JUNB", "CCL5", "IL10RA"),
        Tex=c("HAVCR2", "TIGIT", "LAG3", "TOX", "PDCD1", "IFNG", 
            "KLRG1", "PRF1", "CXCR6", "GLNY", "DUSP4", "EOMES",
            "GZMA", "GZMB", "GZMH", "GZMK"),
        NK=c("ENTPD1", "NKG7", "TBX21", "KLRB1", "KLRC3", "IL2RB", "CD7", "FOS", "FOSB")
    ),
    mye=list(
        macro=c(
            "CD14", "CD68", "CD163",
            "MS4A4A", "MS4A6A", "MRC1", "SLC40A1", 
            "IGF1", "SELENOP", "APOE"),
        moDC=c(
            "CSF1R", "TLR2", "TLR4", "ITGA9",
            "MMP9", "MMP12", "MMP14", "LST1",
            "ITGAX", "CTSD", "TREM2", "FCER1G"),
        DC=c(
            "CD1C", "XCR1", "LAMP3", "IDO1", "CLEC9A", "CLEC10A", "BATF3", 
            "HLA-DPA1", "HLA-DPB1", "HLA-DQB1", "CIITA", "MARCH1"),
        pDC=c("IRF4", "IRF7", "IRF8", "CR2", "FCRL2", "FCRL5", "CLEC4C"),
        mono.nc=c("MARCO", "C1QA", "C1QB", "C1QC", 
            "B2M", "AIF1", "BLVRB", "FCGR3A", "CX3CR1", "LILRB2"),
        mono.c=c("FBP1", "PCK2", "LYZ", "S100A8","S100A9")
    ),
    str=list(
        fib=c(
            "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL5A1", 
            "COL5A2", "COL6A3", "COL12A1", "COL14A1", "COL18A1"),
        mCAF=c("FN1", "LAMA4", "DCN", "MMP1", "MMP2", "IGF2", "POSTN", "IGFBP5"),
        FRCtcz=c("DNM1", "COL27A1", "FDCSP", "CLU", "CR2", "CXCL13"),
        FRCcts=c("CCN1", "JUNB", "JUN", "FOS", "FOSB", "CCL2", "EGR1", "EGR3", "GBP1"),
        FRCpv=c(
            "PDGFRA", "PDGFRB", "VCAM1", "VEGFA", "NDRG2", 
            "MYL9", "CCL19", "CCL21", "CXCL12", "NOTCH3"),
        BEC=c(
            "NOTCH4", "FLT1", "CLEC14A", "CDH5", "CD34", "CD93", 
            "VWF", "HEY1", "DLL4", "CAV1"),
        LEC=c("LYVE1", "PROX1", "CXCL1", "CXCL2"),
        epi=c("KRT5", "KRT7", "KRT14", "KRT16", "KRT17", "KRT19")
    ),
    tum=list(
        tum=c("CCND1", 
            "CD79A", "CD79B", "BCL2", "CD5", 
            "BCL6", "SOX10", "SOX11", "MYC", "TCL1A", 
            "CIITA", "CD19", "CD83", "MS4A1", "BCL7A",
            "CDK1", "MKI67", "TOP2A", "TUBB4A", "TUBB4B"),
        PC=c("MZB1", "XBP1", "JCHAIN", "IGHA1", "IGHA2")
    )
)
