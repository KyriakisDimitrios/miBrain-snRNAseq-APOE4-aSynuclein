# =============================================================================
# 06_gsea_clusterProfiler.R
# Targeted gene set enrichment analysis (GSEA) using curated GO terms and
# custom pathways for Neurons and Astrocytes in A53T miBRAIN snRNA-seq.
#
# Author: Dimitrios Kyriakis
#
# Usage:
#   Rscript 06_gsea_clusterProfiler.R \
#       --dreamlet-rds /path/to/miBRAIN_CLASS_cova_Diff.rds \
#       --outdir       /path/to/output
#
# Arguments:
#   --dreamlet-rds  dreamlet model fit RDS (output of step 05)
#   --outdir        Output directory for GSEA results and figures
#                   [default: current working directory]
# =============================================================================

library(optparse)

option_list <- list(
  make_option("--dreamlet-rds", type="character", default=NULL,
              help="Path to dreamlet RDS from step 05"),
  make_option("--outdir",       type="character", default=".",
              help="Output directory [default: .]")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt[["dreamlet-rds"]])) stop("--dreamlet-rds is required")

DREAMLET_RDS <- opt[["dreamlet-rds"]]
OUTDIR       <- opt$outdir
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

file.path(R.home("bin"), "R")

getwd()

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(RColorBrewer)
library(GOSemSim)

sessionInfo()

# 1. Proteostasis (autophagy / ubiquitin / ER stress / UPR)
# Autophagy: ATG5, ATG7, MAP1LC3B, SQSTM1
# Ubiquitin–proteasome: UBB, UBC, PSMA1, PSMB5
# ER stress / UPR: HSPA5 (BiP), DDIT3 (CHOP), ATF4, XBP1
# 2. Synapse (presynaptic + postsynaptic)
# Presynaptic: SNAP25, STX1A, VAMP2, SYN1, RAB3A
# Postsynaptic: DLG4 (PSD95), SHANK3, HOMER1, GRIN1
# 3. Protein phosphorylation / kinase activity
# Kinases: PLK2, GRK5, CSNK2A1, LRRK2, GSK3B
# Phosphatases: PPP2CA, PPP1CA
#
# Astrocytes (E3 vs E4)
# 1. Lipid metabolism
# Core: APOE, ABCA1, ABCG1
# Regulation: SREBF1, SREBF2
# Lipid handling: FASN, SCD
# 2. Endosomal–lysosomal system
# Endosomal trafficking: RAB5A, RAB7A, RAB11A
# Lysosome function: LAMP1, LAMP2, CTSD, CTSB
# MVB / ESCRT: TSG101, VPS28, CHMP4B, PDCD6IP (ALIX)
# 3. Exosomes / extracellular vesicles
# Markers: CD63, CD81, CD9
# Secretion regulators: RAB27A, RAB27B, RAB35
# Link to lipids: NPC1, NPC2
# 4. Reactive astrocytes (A1 vs A2)
# A1: C3, SERPING1, GBP2
# A2: S100A10, PTX3, EMP1

# -------------------------------------------------------------------
# Neurons (CTR vs E3)
# -------------------------------------------------------------------
neuron_go_targets <- c(
  # 1. Proteostasis
  "Proteostasis_Autophagy"       = "GO:0006914", # autophagy
  "Proteostasis_Ubiquitin"       = "GO:0006511", # ubiquitin-dependent protein catabolic process
  "Proteostasis_ER_Stress"       = "GO:0034976", # response to endoplasmic reticulum stress
  
  # 2. Synapse
  "Synapse_Presynaptic"          = "GO:0099504", # synaptic vesicle cycle (captures SNAP25, VAMP2)
  "Synapse_Postsynaptic"         = "GO:0099084", # postsynaptic specialization organization (captures DLG4, SHANK3)
  
  # 3. Protein phosphorylation / kinase activity
  # "Phosphorylation_Kinases"      = "GO:0006468", # protein phosphorylation
  # "Phosphorylation_Phosphatases" = "GO:0006470",  # protein dephosphorylation
  "Protein_dephosphorylation"    = "GO:0006470"
)

# -------------------------------------------------------------------
# Astrocytes (E3 vs E4)
# -------------------------------------------------------------------
astrocyte_go_targets <- c(
  # 1. Lipid metabolism
  "Lipid_Core"                   = "GO:0006629", # lipid metabolic process (captures APOE, ABCA1, ABCG1)
  "Lipid_Regulation"             = "GO:0019216", # regulation of lipid metabolic process (captures SREBFs)
  "Lipid_Handling"               = "GO:0006633", # fatty acid biosynthetic process (captures FASN, SCD)
  
  # 2. Endosomal-lysosomal system
  "Endolysosomal_Trafficking"    = "GO:0016197", # endosomal transport (captures RAB5/7/11)
  "Endolysosomal_Lysosome"       = "GO:0007040", # lysosome organization (captures LAMPs, CTSD/B)
  "Endolysosomal_MVB"            = "GO:0071985", # multivesicular body sorting pathway (captures TSG101, ESCRT)
  
  # 3. Exosomes / extracellular vesicles
  "Exosome_Markers"              = "GO:0140112", # extracellular vesicle assembly (captures CD63/81/9)
  "Exosome_Secretion"            = "GO:1990182", # exosomal secretion (captures RAB27A/B)
  "Exosome_Lipid_Link"           = "GO:0030301",  # cholesterol transport (captures NPC1/2)
  "protein_phosphorylation"      = "GO:0006468",
  "protein_dephosphorylation"    = "GO:0006470"
    
)

library(dreamlet)

#' @rdname gseaplot
#' @exportMethod gseaplot
setMethod(
    "gseaplot",
    signature(x = "gseaResult"),
    function(
        x,
        geneSetID,
        by = "all",
        title = "",
        color = 'black',
        color.line = "green",
        color.vline = "#FA5860",
        ...
    ) {
        gseaplot.gseaResult(
            x,
            geneSetID = geneSetID,
            by = by,
            title = title,
            color = color,
            color.line = color.line,
            color.vline = color.vline,
            ...
        )
    }
)

#' @rdname gseaplot
#' @param color color of line segments
#' @param color.line color of running enrichment score line
#' @param color.vline color of vertical line indicating the
#' maximum/minimal running enrichment score
#' @return ggplot2 object
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_linerange
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplotGrob
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 ggplot_gtable
#' @importFrom ggplot2 ggplot_build
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 rel
#' @importFrom aplot plot_list
#' @author Guangchuang Yu
gseaplot.gseaResult <- function(
    x,
    geneSetID,
    by = "all",
    title = "",
    color = 'black',
    color.line = "green",
    color.vline = "#FA5860",
    ...
) {
    by <- match.arg(by, c("runningScore", "preranked", "all"))
    gsdata <- gsInfo(x, geneSetID)
    p <- ggplot(gsdata, aes(x = .data$x)) +
        theme_dose() +
        xlab("Position in the Ranked List of Genes")
    if (by == "runningScore" || by == "all") {
        p.res <- p +
            geom_linerange(
                aes(ymin = .data$ymin, ymax = .data$ymax),
                color = color
            )
        p.res <- p.res +
            geom_line(
                aes(y = .data$runningScore),
                color = color.line,
                linewidth = 1
            )
        enrichmentScore <- x@result[geneSetID, "enrichmentScore"]
        es.df <- data.frame(
            es = which.min(abs(p$data$runningScore - enrichmentScore))
        )
        p.res <- p.res +
            geom_vline(
                data = es.df,
                aes(xintercept = .data$es),
                colour = color.vline,
                linetype = "dashed"
            )
        p.res <- p.res + ylab("Running Enrichment Score")
        p.res <- p.res + geom_hline(yintercept = 0)
    }
    if (by == "preranked" || by == "all") {
        df2 <- data.frame(x = which(p$data$position == 1))
        df2$y <- p$data$geneList[df2$x]
        p.pos <- p +
            geom_segment(
                data = df2,
                aes(x = .data$x, xend = .data$x, y = .data$y, yend = 0),
                color = color
            )
        p.pos <- p.pos +
            ylab("Ranked List Metric") +
            xlim(0, length(p$data$geneList))
    }
    if (by == "runningScore") {
        return(p.res + ggtitle(title))
    }
    if (by == "preranked") {
        return(p.pos + ggtitle(title))
    }

    p.pos <- p.pos +
        xlab(NULL) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p.pos <- p.pos +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5, size = rel(2)))
    #plot_list(gglist =  list(p.pos, p.res), ncol=1)

    aplot::gglist(gglist = list(p.pos, p.res), ncol = 1)
}


#' extract gsea result of selected geneSet
#'
#'
#' @title gsInfo
#' @param object gseaResult object
#' @param geneSetID gene set ID
#' @return data.frame
#' @author Guangchuang Yu
## @export
gsInfo <- function(object, geneSetID) {
    geneList <- object@geneList

    if (is.numeric(geneSetID)) {
        geneSetID <- object@result[geneSetID, "ID"]
    }

    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- enrichit::gseaScores(geneList, geneSet, exponent, fortify = TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore)) / 20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList
    if (length(object@gene2Symbol) == 0) {
        df$gene <- names(geneList)
    } else {
        df$gene <- object@gene2Symbol[names(geneList)]
    }

    df$Description <- object@result[geneSetID, "Description"]
    return(df)
}


get_gsdata <- function(x, geneSetID) {
    if (length(geneSetID) == 1) {
        gsdata <- gsInfo(x, geneSetID)
        return(gsdata)
    }

    lapply(geneSetID, gsInfo, object = x) |>
        yulab.utils::rbindlist()
}

#' Horizontal plot for GSEA result
#'
#'
#' @title hplot
#' @param x gseaResult object
#' @param geneSetID gene set ID
#' @return horizontal plot
#' @export
#' @author Guangchuang Yu
hplot <- function(x, geneSetID) {
    if (!inherits(x, "gseaResult")) {
        stop("hplot only work for GSEA result")
    }

    gsdata <- get_gsdata(x, geneSetID)

    ggplot(gsdata, aes(.data$x, .data$runningScore)) +
        ggHoriPlot::geom_horizon(origin = 'min', horizonscale = 4) +
        facet_grid(Description ~ .) +
        #ggHoriPlot::scale_fill_hcl(palette = 'Peach', reverse = TRUE) +
        ggHoriPlot::scale_fill_hcl(palette = 'BluGrn', reverse = TRUE) +
        theme_minimal() +
        ggfun::theme_noyaxis() +
        theme(
            panel.spacing.y = unit(0, "lines"),
            strip.text.y = element_text(angle = 0),
            legend.position = 'none',
            panel.border = element_blank(),
            panel.grid = element_blank(),
        ) +
        xlab(NULL) +
        ylab(NULL)
}

#' GSEA plot that mimic the plot generated by broad institute's GSEA software
#'
#'
#' @title gseaplot2
#' @param x gseaResult object
#' @param geneSetID gene set ID
#' @param title plot title
#' @param color color of running enrichment score line
#' @param base_size base font size
#' @param rel_heights relative heights of subplots
#' @param subplots which subplots to be displayed
#' @param pvalue_table whether add pvalue table
#' @param pvalue_table_columns selected columns to be plotted in the `pvalue_table`
#' @param pvalue_table_rownames selected column as the rownames of the `pvalue_table`. If set to NULL, no rownames will be displayed.
#' @param ES_geom geom for plotting running enrichment score,
#' one of 'line' or 'dot'
#' @return plot
#' @export
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 theme_void
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 margin
#' @importFrom ggplot2 annotation_custom
#' @importFrom stats quantile
#' @importFrom RColorBrewer brewer.pal
#' @author Guangchuang Yu
gseaplot2 <- function(
    x,
    geneSetID,
    title = "",
    color = "green",
    base_size = 11,
    rel_heights = c(1.5, .5, 1),
    subplots = 1:3,
    pvalue_table = FALSE,
    pvalue_table_columns = c("pvalue", "p.adjust"),
    pvalue_table_rownames = "Description",
    ES_geom = "line"
) {
    ES_geom <- match.arg(ES_geom, c("line", "dot"))

    geneList <- position <- NULL ## to satisfy codetool

    gsdata <- get_gsdata(x, geneSetID)

    p <- ggplot(gsdata, aes(x = .data$x)) +
        xlab(NULL) +
        theme_classic(base_size) +
        theme(
            panel.grid.major = element_line(colour = "grey92"),
            panel.grid.minor = element_line(colour = "grey92"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
        ) +
        scale_x_continuous(expand = c(0, 0))

    if (ES_geom == "line") {
        es_layer <- geom_line(
            aes(y = .data$runningScore, color = .data$Description),
            linewidth = 1
        )
    } else {
        es_layer <- geom_point(
            aes(y = .data$runningScore, color = .data$Description),
            size = 1,
            data = subset(gsdata, position == 1)
        )
    }

    p.res <- p +
        es_layer +
        theme(
            legend.position = "inside",
            legend.position.inside = c(.8, .8),
            legend.title = element_blank(),
            legend.background = element_rect(fill = "transparent")
        )

    p.res <- p.res +
        ylab("Running Enrichment Score") +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            plot.margin = margin(t = .2, r = .2, b = 0, l = .2, unit = "cm")
        )

    # Vectorized ymin/ymax assignment
    terms <- unique(gsdata$Description)
    term_indices <- match(gsdata$Description, terms) - 1
    idx <- which(gsdata$ymin != 0)
    gsdata[idx, "ymin"] <- term_indices[idx]
    gsdata[idx, "ymax"] <- term_indices[idx] + 1
    p2 <- ggplot(gsdata, aes(x = .data$x)) +
        geom_linerange(aes(
            ymin = .data$ymin,
            ymax = .data$ymax,
            color = .data$Description
        )) +
        xlab(NULL) +
        ylab(NULL) +
        theme_classic(base_size) +
        theme(
            legend.position = "none",
            plot.margin = margin(t = -.1, b = 0, unit = "cm"),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line.x = element_blank()
        ) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))

    if (length(geneSetID) == 1) {
        ## geneList <- gsdata$geneList
        ## j <- which.min(abs(geneList))
        ## v1 <- quantile(geneList[1:j], seq(0,1, length.out=6))[1:5]
        ## v2 <- quantile(geneList[j:length(geneList)], seq(0,1, length.out=6))[1:5]

        ## v <- sort(c(v1, v2))
        ## inv <- findInterval(geneList, v)

        v <- seq(1, sum(gsdata$position), length.out = 9)
        inv <- findInterval(rev(cumsum(gsdata$position)), v)
        if (min(inv) == 0) {
            inv <- inv + 1
        }

        col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))

        ymin <- min(p2$data$ymin)
        yy <- max(p2$data$ymax - p2$data$ymin) * .3
        xmin <- which(!duplicated(inv))
        xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
        d <- data.frame(
            ymin = ymin,
            ymax = yy,
            xmin = xmin,
            xmax = xmax,
            col = col[unique(inv)]
        )
        p2 <- p2 +
            geom_rect(
                aes(
                    xmin = .data$xmin,
                    xmax = .data$xmax,
                    ymin = .data$ymin,
                    ymax = .data$ymax,
                    fill = I(col)
                ),
                data = d,
                alpha = .9,
                inherit.aes = FALSE
            )
    }

    ## p2 <- p2 +
    ## geom_rect(aes(xmin=x-.5, xmax=x+.5, fill=geneList),
    ##           ymin=ymin, ymax = ymin + yy, alpha=.5) +
    ## theme(legend.position="none") +
    ## scale_fill_gradientn(colors=color_palette(c("blue", "red")))

    df2 <- p$data #data.frame(x = which(p$data$position == 1))
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p +
        geom_segment(
            data = df2,
            aes(x = .data$x, xend = .data$x, y = .data$y, yend = 0),
            color = "grey"
        )
    p.pos <- p.pos +
        ylab("Ranked List Metric") +
        xlab("Rank in Ordered Dataset") +
        theme(
            plot.margin = margin(t = -.1, r = .2, b = .2, l = .2, unit = "cm")
        )

    if (!is.null(title) && !is.na(title) && title != "") {
        p.res <- p.res + ggtitle(title)
    }

    if (length(color) == length(geneSetID)) {
        p.res <- p.res + scale_color_manual(values = color)
        if (length(color) == 1) {
            p.res <- p.res + theme(legend.position = "none")
            p2 <- p2 + scale_color_manual(values = "black")
        } else {
            p2 <- p2 + scale_color_manual(values = color)
        }
    }

    if (pvalue_table) {
        pd <- x[geneSetID, pvalue_table_columns]
        # pd <- pd[order(pd[,1], decreasing=FALSE),]
        if (is.null(pvalue_table_rownames)) {
            rows <- NULL
        } else {
            # rownames(pd) <- pd$Description
            if (length(pvalue_table_rownames) != 1) {
                stop(
                    "the length of `pvalue_table_rownames` should be equal to 1"
                )
            }

            rows <- x[geneSetID, pvalue_table_rownames]
        }

        # pd <- round(pd, 4)
        for (i in seq_len(ncol(pd))) {
            pd[, i] <- format(pd[, i], digits = 4)
        }
        tp <- tableGrob2(d = pd, p = p.res, rows = rows)

        p.res <- p.res +
            theme(legend.position = "none") +
            annotation_custom(
                tp,
                xmin = quantile(p.res$data$x, .5),
                xmax = quantile(p.res$data$x, .95),
                ymin = quantile(p.res$data$runningScore, .75),
                ymax = quantile(p.res$data$runningScore, .9)
            )
    }

    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] +
        theme(
            axis.line.x = element_line(),
            axis.ticks.x = element_line(),
            axis.text.x = element_text()
        )

    if (length(subplots) == 1) {
        return(
            plotlist[[1]] +
                theme(
                    plot.margin = margin(
                        t = .2,
                        r = .2,
                        b = .2,
                        l = .2,
                        unit = "cm"
                    )
                )
        )
    }

    if (length(rel_heights) > length(subplots)) {
        rel_heights <- rel_heights[subplots]
    }

    # aplot::plot_list(gglist = plotlist, ncol=1, heights=rel_heights)
    aplot::gglist(gglist = plotlist, ncol = 1, heights = rel_heights)
}


#' plot ranked list of genes with running enrichment score as bar height
#'
#'
#' @title gsearank
#' @param x gseaResult object
#' @param geneSetID gene set ID
#' @param title plot title
#' @param output one of 'plot' or 'table' (for exporting data)
#' @return ggplot object
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 theme_minimal
#' @export
#' @author Guangchuang Yu
gsearank <- function(x, geneSetID, title = "", output = "plot") {
    output <- match.arg(output, c("plot", "table"))

    position <- NULL
    gsdata <- gsInfo(x, geneSetID)
    gsdata <- subset(gsdata, position == 1)

    if (output == "table") {
        res <- gsdata[, c("gene", "x", "runningScore")]
        if (x[geneSetID, "NES"] > 0) {
            res$core <- "NO"
            res$core[1:which.max(gsdata$runningScore)] <- "YES"
        } else {
            res$core <- "NO"
            res$core[which.min(gsdata$runningScore):nrow(res)] <- "YES"
        }
        names(res) <- c(
            "gene",
            "rank in geneList",
            "running ES",
            "core enrichment"
        )
        rownames(res) <- NULL
        return(res)
    }

    p <- ggplot(gsdata, aes(x = .data$x, y = .data$runningScore)) +
        geom_segment(aes(xend = .data$x, yend = 0)) +
        ggtitle(title) +
        xlab("Position in the Ranked List of Genes") +
        ylab("Running Enrichment Score") +
        theme_minimal()
    return(p)
}


#' label genes in running score plot
#'
#'
#' @title geom_gsea_gene
#' @param genes selected genes to be labeled
#' @param mapping aesthetic mapping, default is NULL
#' @param geom geometric layer to plot the gene labels, default is geom_text
#' @param ... additional parameters passed to the 'geom'
#' @param geneSet choose which gene set(s) to be label if the plot contains multiple gene sets
#' @return ggplot object
#' @importFrom rlang .data
#' @export
#' @author Guangchuang Yu
geom_gsea_gene <- function(
    genes,
    mapping = NULL,
    geom = ggplot2::geom_text,
    ...,
    geneSet = NULL
) {
    default_mapping <- aes(
        x = .data$x,
        y = .data$runningScore,
        label = .data$gene
    )
    if (is.null(mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, mapping)
    }
    if (is.null(geneSet)) {
        data <- ggtree::td_filter(.data$gene %in% genes)
    } else {
        data <- ggtree::td_filter(
            .data$gene %in% genes & .data$Description %in% geneSet
        )
    }

    geom(mapping = mapping, data = data, ...)
}

#' plot table
#'
#'
#' @title ggtable
#' @param d data frame
#' @param p ggplot object to extract color to color rownames(d), optional
#' @importFrom rlang check_installed
#' @return ggplot object
#' @export
#' @author guangchuang yu
ggtable <- function(d, p = NULL) {
    # has_package("ggplotify")
    rlang::check_installed('ggplotify', 'for `ggtable()`.')
    ggplotify::as.ggplot(tableGrob2(d, p))
}

#' @importFrom grid gpar
#' @importFrom ggplot2 ggplot_build
#' @importFrom rlang check_installed
tableGrob2 <- function(d, p = NULL, rows=NULL) {
    # has_package("gridExtra")
    order_index <- order(rownames(d))
    d <- d[order_index,]
    if (!is.null(rows)) {
        rows <- rows[order_index]
    }
    rlang::check_installed('gridExtra', 'for `tableGrob2()`.')
    tp <- gridExtra::tableGrob(d, rows=rows)
    if (is.null(p) || is.null(rows)) {
        return(tp)
    }

    # Fix bug: The 'group' order of lines and dots/path is different
    p_data <- ggplot_build(p)$data[[1]]
    # pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
    p_data <- p_data[order(p_data[["group"]]), ]
    pcol <- unique(p_data[["colour"]])
    ## This is fine too
    ## pcol <- unique(p_data[["colour"]])[unique(p_data[["group"]])]  
    j <- which(tp$layout$name == "rowhead-fg")

    for (i in seq_along(pcol)) {
        tp$grobs[j][[i+1]][["gp"]] <- gpar(col = pcol[i])
    }
    return(tp)
}

run_standard_gsea_unfiltered <- function(ranked_genes, ont = "BP") {
  log_msg("INFO", "Initiating gseGO (Ontology: BP) with unfiltered p-values...")
  
  set.seed(42)
  res <- gseGO(
    geneList     = ranked_genes,
    ont          = ont,
    OrgDb        = org.Hs.eg.db,
    keyType      = "SYMBOL", 
    minGSSize    = 5,
    maxGSSize    = 1000,
      eps          = 0, 
      nPermSimple  = 10000,
    # CRITICAL CHANGE: Keep all calculated pathways in the object so gseaplot2 can find them
    pvalueCutoff = 1,   
    verbose      = FALSE
  )
  
  log_msg("INFO", "Unfiltered gseGO completed.")
  return(res)
}

plot_automatic_3panel <- function(gsea_result_obj, pathway_id) {
  #' Uses enrichplot's native function to generate a 3-panel GSEA plot.
  #' 
  #' @param gsea_result_obj The complete output object from gseGO() or GSEA()
  #' @param pathway_id The specific ID to plot (e.g., "GO:0006914")
  
  log_msg("INFO", sprintf("Generating 3-panel gseaplot2 for %s...", pathway_id))
  
  # gseaplot2 automatically handles the curve, barcode, and rank metric.
  # pvalue_table = TRUE overlays the actual NES and p-values on the plot.
  p <- gseaplot2(
    gsea_result_obj, 
    geneSetID = pathway_id, 
    color = "firebrick",      # Color of the enrichment curve
    pvalue_table = TRUE,      # Adds a nice stats table inside the plot
    subplots = 1:3            # Ensures all 3 panels (Curve, Barcode, Rank) are drawn
  )
  
  return(p)
}



suppressPackageStartupMessages({
  library(enrichplot)
  library(ggplot2)
})

# Standard pipeline logging
log_msg <- function(level, message) {
  cat(sprintf("%s | %s | %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, message), file=stderr())
}

plot_labeled_3panel_safe <- function(gsea_result_obj, pathway_id, condition_positive = "A53T", condition_negative = "Ctrl") {
  
  log_msg("INFO", sprintf("Generating labeled 3-panel plot for %s...", pathway_id))
  
  # 1. Base plot
  p <- gseaplot2(
    gsea_result_obj,
    title=pathway_id, 
    geneSetID = pathway_id, 
    subplots = 1:2, # Show enrichment score, hits, and rank
    color = "forestgreen",
    pvalue_table = TRUE,
    pvalue_table_rownames = NULL,
    pvalue_table_columns = c("ID", "NES", "p.adjust"))        
    
  p[[2]] <- p[[2]] + labs(x = axis_label) +
    theme(
      axis.title.x = element_text(face = "bold", size = 12),
      legend.position = "none" 
    )

  
  # # 2. Extract rankings directly from the GSEA object (Version-proof)
  # ranks_vec <- gsea_result_obj@geneList
  # safe_metric_data <- data.frame(
  #   x_rank = 1:length(ranks_vec),
  #   y_metric = as.numeric(ranks_vec)
  # )
  
  # # 3. Construct label
  # axis_label <- sprintf("← Upregulated in %s (Red)                   Upregulated in %s (Blue) →", 
  #                       condition_positive, condition_negative)
  
  # log_msg("INFO", "Injecting color gradient cleanly...")
  
  # # 4. Inject using inherit.aes = FALSE to sever ties with the locked enrichplot environment
  # p[[2]] <- p[[2]] + 
  #   geom_col(
  #     data = safe_metric_data, 
  #     aes(x = x_rank, y = y_metric, fill = y_metric), 
  #     inherit.aes = FALSE  # <- This is the magic bullet that stops the crash
  #   ) + 
  #   scale_fill_gradient2(
  #     low = "dodgerblue",   
  #     mid = "white",        
  #     high = "firebrick",   
  #     midpoint = 0
  #   ) +
  #   labs(x = axis_label) +
  #   theme(
  #     axis.title.x = element_text(face = "bold", size = 12),
  #     legend.position = "none" 
  #   )
  
  return(p)
}

dreamlet_fit <- readRDS(DREAMLET_RDS)

# # NEURONS

set.seed(20210224)

cluster = "Neurons"
COEF <- 'E3vsCntrl'
condition_positive='E3'
condition_negative ='Ctrl'
file_cl = 'Neurons'
# GO BP
axis_label <- sprintf("<- Upregulated in %s (Red)                   Upregulated in %s (Blue) ->", 
                        condition_positive, condition_negative)

ranked_ExN <- topTable(dreamlet_fit, coef = COEF, number = Inf) %>%
  as.data.frame() %>%
  dplyr::filter(assay == cluster) %>%
  dplyr::arrange(desc(t)) %>%
  { setNames(.$t, .$ID) }


library(clusterProfiler)
library(enrichplot)
library(dplyr)

# Standard logging
message(sprintf("[%s] | INFO | Preparing custom gene sets for GSEA", Sys.time()))

# Define custom pathways from your image
custom_pathways_list <- list(
  Synuclein_Pathology = c("ABCA1", "SNCA", "CTSB", "SNCAIP", "ATG7", "PRKN"),
  Synaptic_Vesicle_Maturation = c("UNC13C", "UNC13A", "UNC13B", "RIMS1", "RIMS2", 
                                  "RIMBP2", "STX1A", "SNAP25", "VAMP2", "SYT1", 
                                  "SYT2", "SYN1", "SYN2", "RAB3A", "RAB3B")
    
)

# Convert list to a long-format data frame (TERM2GENE)
# Column 1: Pathway Name, Column 2: Gene Symbol
term2gene <- stack(custom_pathways_list) %>%
  dplyr::select(term = ind, gene = values)

# Run GSEA
# Using pvalueCutoff = 1 to ensure we can plot even if not significant
gsea_custom <- GSEA(ranked_ExN, 
                    TERM2GENE = term2gene, 
                    minGSSize = 1,      # Important for small custom sets
                    maxGSSize = 500,
                    eps          = 0,
                    nPermSimple  = 10000,
                    pvalueCutoff = 1, 
                    verbose = FALSE)

set.seed(20210224)
# Generate gseaplot2
# You can plot one or both together
message(sprintf("[%s] | INFO | Generating gseaplot2 for %s", Sys.time(), cluster))
options(repr.plot.width = 10, repr.plot.height = 5, repr.plot.res = 300)

for (custo_term in names(custom_pathways_list)){
    p <- gseaplot2(gsea_custom, 
                    geneSetID = custo_term, 
                    title = paste(cluster, ":",custo_term),
                    subplots = 1:2, # Show enrichment score, hits, and rank
                    color = "forestgreen",
                    pvalue_table = TRUE,
                      pvalue_table_rownames = NULL,
                    pvalue_table_columns = c("ID", "NES", "p.adjust"))
    
    # 3. Construct label
    axis_label <- sprintf("<- Upregulated in %s (Red)                   Upregulated in %s (Blue) ->", 
                        condition_positive, condition_negative)
    
    
    # 4. Inject using inherit.aes = FALSE to sever ties with the locked enrichplot environment
    p[[2]] <- p[[2]] + 
        labs(x = axis_label) +
        theme(
          axis.title.x = element_text(face = "bold", size = 12),
          legend.position = "none" 
        )
    
    file_name_pdf <- sprintf("figures/GSEA/GSEA_%s_%s_%s.pdf",file_cl,COEF, custo_term)
    pdf(file = file_name_pdf, width = 10, height = 5)
    print(p)
    dev.off()

}

# p <- gseaplot2(gsea_custom, 
#                 geneSetID = "Synuclein_Pathology", 
#                 title = paste(cluster, ": Synuclein Pathology"),
#                 subplots = 1:2, # Show enrichment score, hits, and rank
#                 color = "forestgreen",
#                 pvalue_table = TRUE,
#                   pvalue_table_rownames = NULL,
#                 pvalue_table_columns = c("ID", "NES", "p.adjust"))

# # 3. Construct label
# axis_label <- sprintf("← Upregulated in %s (Red)                   Upregulated in %s (Blue) →", 
#                     condition_positive, condition_negative)


# # 4. Inject using inherit.aes = FALSE to sever ties with the locked enrichplot environment
# p[[2]] <- p[[2]] + 
#     labs(x = axis_label) +
#     theme(
#       axis.title.x = element_text(face = "bold", size = 12),
#       legend.position = "none" 
#     )

# p

# file_name_pdf <- sprintf("figures/GSEA/GSEA_%s_%s_%s.pdf",file_cl,COEF, target_name)
# pdf(file = file_name_pdf, width = 10, height = 5)
# print(final_plot)
# dev.off()

file_name_pdf

set.seed(20210224)
# 1. Rerun the GSEA to keep all results
gsea_neurons_all <- run_standard_gsea_unfiltered(ranked_ExN)

set.seed(20210224)

# Ensure output directory exists before saving
dir.create("figures/GSEA", recursive = TRUE, showWarnings = FALSE)

options(repr.plot.width = 10, repr.plot.height = 5, repr.plot.res = 300)
for (target_name in names(neuron_go_targets)) {
target_id <- neuron_go_targets[[target_name]]

# The Safety Gate: Check if the ID survived the min/max size filters
if (!target_id %in% rownames(gsea_neurons_all@result)) {
    cat(sprintf("Skipping %s (%s): Filtered out (likely didn't meet minGSSize=15)\n", target_name, target_id))
    next # Skip to the next target in the list
    }
    
    # If it exists, generate the plot
    cat(sprintf("Plotting %s (%s)...\n", target_name, target_id))
    final_plot <- plot_labeled_3panel_safe(
        gsea_neurons_all, 
        pathway_id = target_id, 
        condition_positive = condition_positive, 
        condition_negative = condition_negative
    )
    
    
    # ---------------------------------------------------------
    # THE FIX: Base R graphic device instead of ggsave
    # ---------------------------------------------------------
    file_name <- sprintf("figures/GSEA/GSEA_%s_%s_%s.png",file_cl,COEF, target_name)
    file_name_pdf <- sprintf("figures/GSEA/GSEA_%s_%s_%s.pdf",file_cl,COEF, target_name)
    
    # 1. Open a blank PNG file with your desired dimensions
    # png(filename = file_name, width = 10, height = 5, units = "in", res = 300)
    
    # 2. Draw the plot onto the canvas
    pdf(file = file_name_pdf, width = 10, height = 5)
    print(final_plot)
    dev.off()
    


cat(sprintf("Saved: %s\n\n", file_name))

}

# # Astrocytes

cluster = "Astrocytes"
COEF <- 'E4vsE3'
condition_positive='E4'
condition_negative ='E3'
file_cl = 'Astro'

set.seed(20210224)
# GO BP
ranked_Astro <- topTable(dreamlet_fit, coef = COEF, number = Inf) %>%
  as.data.frame() %>%
  dplyr::filter(assay == cluster) %>%
  dplyr::arrange(desc(t)) %>%
  { setNames(.$t, .$ID) }

# 1. Rerun the GSEA to keep all results
gsea_astro_all <- run_standard_gsea_unfiltered(ranked_Astro)

library(clusterProfiler)
library(enrichplot)
library(dplyr)
set.seed(20210224)

# Standard logging
message(sprintf("[%s] | INFO | Preparing custom gene sets for GSEA", Sys.time()))

# Define custom pathways from your image
custom_pathways_list <- list(
    Lipid_droplet_biogenesis = c("PLIN1","PLIN2","PLIN3","PNPLA2","DGAT1","DGAT2","CIDEC","CIDEA","ATGL","SCD","FABP5","G0S2","ACAT1","ACAT2"),   
    Matrisome = c("LAMA2","LAMA4","LAMB1","LAMC1","FN1","SPARC","TGFBI","LOX","THBS1","THBS2","SERPINE1","MMP2","MMP14","TIMP1","TIMP3"),
    Synuclein_Pathology = c("ABCA1","SNCA","CTSB","SNCAIP","ATG7", "PRKN"),
    ASTRO_INFLAMMATORY_IL1 = c("IL1A","IL1B","IL1R1","IL1RAP","IL1RAPL1","IL1RN","MYD88","IRAK1","IRAK4","TRAF6", "NFKB1","RELA","MAPK1","MAPK14","GFAP","VIM","SERPINA3","IL6","CXCL10","CCL2","C3","C1S","C1R","STAT3","SOCS3","TNFAIP3")
    
)

# Convert list to a long-format data frame (TERM2GENE)
# Column 1: Pathway Name, Column 2: Gene Symbol
term2gene <- stack(custom_pathways_list) %>%
  dplyr::select(term = ind, gene = values)

# Run GSEA
# Using pvalueCutoff = 1 to ensure we can plot even if not significant
gsea_custom <- GSEA(ranked_Astro, 
                    TERM2GENE = term2gene, 
                    minGSSize = 1,      # Important for small custom sets
                    maxGSSize = 500,
                     eps          = 0,
                    nPermSimple  = 10000,
                    pvalueCutoff = 1, 
                    verbose = FALSE)

set.seed(20210224)
# Generate gseaplot2
# You can plot one or both together
message(sprintf("[%s] | INFO | Generating gseaplot2 for %s", Sys.time(), cluster))
options(repr.plot.width = 10, repr.plot.height = 5, repr.plot.res = 300)

for (custo_term in names(custom_pathways_list)){
    p <- gseaplot2(gsea_custom, 
                    geneSetID = custo_term, 
                    title = paste(cluster, ":",custo_term),
                    subplots = 1:2, # Show enrichment score, hits, and rank
                    color = "forestgreen",
                    pvalue_table = TRUE,
                      pvalue_table_rownames = NULL,
                    pvalue_table_columns = c("ID", "NES", "p.adjust"))
    
    # 3. Construct label
    axis_label <- sprintf("<- Upregulated in %s (Red)                   Upregulated in %s (Blue) ->", 
                        condition_positive, condition_negative)
    
    
    # 4. Inject using inherit.aes = FALSE to sever ties with the locked enrichplot environment
    p[[2]] <- p[[2]] + 
        labs(x = axis_label) +
        theme(
          axis.title.x = element_text(face = "bold", size = 12),
          legend.position = "none" 
        )
    
    file_name_pdf <- sprintf("figures/GSEA/GSEA_%s_%s_%s.pdf",file_cl,COEF, custo_term)
    pdf(file = file_name_pdf, width = 10, height = 5)
    print(p)
    dev.off()
    

}

astrocyte_go_targets

gsea_astro_all@result %>% filter(ID=='GO:0006470')

gsea_astro_all@result %>% filter(ID=='GO:0006468')


# Ensure output directory exists before saving
dir.create("figures/GSEA", recursive = TRUE, showWarnings = FALSE)

options(repr.plot.width = 10, repr.plot.height = 5, repr.plot.res = 300)
for (target_name in names(astrocyte_go_targets)) {
target_id <- astrocyte_go_targets[[target_name]]

# The Safety Gate: Check if the ID survived the min/max size filters
if (!target_id %in% rownames(gsea_astro_all@result)) {
    cat(sprintf("Skipping %s (%s): Filtered out (likely didn't meet minGSSize=15)\n", target_name, target_id))
    next # Skip to the next target in the list
    }
    
    # If it exists, generate the plot
    cat(sprintf("Plotting %s (%s)...\n", target_name, target_id))
    final_plot <- plot_labeled_3panel_safe(
        gsea_astro_all, 
        pathway_id = target_id, 
        condition_positive = condition_positive, 
        condition_negative = condition_negative
    )
    
    # print(final_plot)
    
    # Save the plot with a dynamically generated filename
    # file_name <- sprintf("figures/GSEA/GSEA_Astrocytes_E4_vs_E3_%s.png", target_name)
    # ggsave(filename = file_name, plot = final_plot, width = 10, height = 5, dpi = 300)
    # cat(sprintf("Saved: %s\n\n", file_name))
    
    # ---------------------------------------------------------
    # THE FIX: Base R graphic device instead of ggsave
    # ---------------------------------------------------------
    file_name <- sprintf("figures/GSEA/GSEA_%s_%s_%s.png",cluster,COEF, target_name)
    file_name_pdf <- sprintf("figures/GSEA/GSEA_%s_%s_%s.pdf",cluster,COEF, target_name)

    # 1. Open a blank PNG file with your desired dimensions
    # png(filename = file_name, width = 10, height = 5, units = "in", res = 300)
    
    # 2. Draw the plot onto the canvas
    # print(final_plot)
    pdf(file = file_name_pdf, width = 10, height = 5)
    print(final_plot)
    dev.off()

# 3. Close the canvas and save the file to disk
# invisible(dev.off())
# ---------------------------------------------------------

cat(sprintf("Saved: %s\n\n", file_name_pdf))

}

# # Generate gseaplot2
# # You can plot one or both together
# message(sprintf("[%s] | INFO | Generating gseaplot2 for %s", Sys.time(), cluster))
# options(repr.plot.width = 10, repr.plot.height = 5, repr.plot.res = 300)

# p <- gseaplot2(gsea_custom, 
#                 geneSetID = "Lipid_droplet_biogenesis", 
#                 title = paste(cluster, ": Lipid_droplet_biogenesis"),
#                 subplots = 1:2, # Show enrichment score, hits, and rank
#                 color = "forestgreen",
#                 pvalue_table = TRUE,
#                   pvalue_table_rownames = NULL,
#                 pvalue_table_columns = c("ID", "NES", "qvalue"))

# # 3. Construct label
# axis_label <- sprintf("← Upregulated in %s (Red)                   Upregulated in %s (Blue) →", 
#                     condition_positive, condition_negative)


# # 4. Inject using inherit.aes = FALSE to sever ties with the locked enrichplot environment
# p[[2]] <- p[[2]] + 
#     labs(x = axis_label) +
#     theme(
#       axis.title.x = element_text(face = "bold", size = 12),
#       legend.position = "none" 
#     )

# p

# p <- gseaplot2(gsea_custom, 
#                 geneSetID = "Matrisome", 
#                 title = paste(cluster, ": Matrisome"),
#                 subplots = 1:2, # Show enrichment score, hits, and rank
#                 color = "forestgreen",
#                 pvalue_table = TRUE,
#                   pvalue_table_rownames = NULL,
#                 pvalue_table_columns = c("ID", "NES", "qvalue"))

# # 3. Construct label
# axis_label <- sprintf("← Upregulated in %s (Red)                   Upregulated in %s (Blue) →", 
#                     condition_positive, condition_negative)


# # 4. Inject using inherit.aes = FALSE to sever ties with the locked enrichplot environment
# p[[2]] <- p[[2]] + 
#     labs(x = axis_label) +
#     theme(
#       axis.title.x = element_text(face = "bold", size = 12),
#       legend.position = "none" 
#     )

# p


# p <- gseaplot2(gsea_custom, 
#                 geneSetID = "Matrisome", 
#                 title = paste(cluster, ": Matrisome"),
#                 subplots = 1:2, # Show enrichment score, hits, and rank
#                 color = "forestgreen",
#                 pvalue_table = TRUE,
#                   pvalue_table_rownames = NULL,
#                 pvalue_table_columns = c("ID", "NES", "qvalue"))

# # 3. Construct label
# axis_label <- sprintf("← Upregulated in %s (Red)                   Upregulated in %s (Blue) →", 
#                     condition_positive, condition_negative)


# # 4. Inject using inherit.aes = FALSE to sever ties with the locked enrichplot environment
# p[[2]] <- p[[2]] + 
#     labs(x = axis_label) +
#     theme(
#       axis.title.x = element_text(face = "bold", size = 12),
#       legend.position = "none" 
#     )

# p


# p <- gseaplot2(gsea_custom, 
#                 geneSetID = "Matrisome", 
#                 title = paste(cluster, ": Matrisome"),
#                 subplots = 1:2, # Show enrichment score, hits, and rank
#                 color = "forestgreen",
#                 pvalue_table = TRUE,
#                   pvalue_table_rownames = NULL,
#                 pvalue_table_columns = c("ID", "NES", "qvalue"))

# # 3. Construct label
# axis_label <- sprintf("← Upregulated in %s (Red)                   Upregulated in %s (Blue) →", 
#                     condition_positive, condition_negative)


# # 4. Inject using inherit.aes = FALSE to sever ties with the locked enrichplot environment
# p[[2]] <- p[[2]] + 
#     labs(x = axis_label) +
#     theme(
#       axis.title.x = element_text(face = "bold", size = 12),
#       legend.position = "none" 
#     )

# p
