## The following list of commands will generate the plots created in iSEE
## Copy them into a script or an R session containing your SingleCellExperiment.
## All commands below refer to your SingleCellExperiment object as `se`.

se <- sc_data
colormap <- ExperimentColorMap()
colormap <- synchronizeAssays(colormap, se)
all_contents <- list()

################################################################################
# Defining brushes
################################################################################

all_active <- list()
all_active[['ReducedDimensionPlot1']] <- list()
all_active[['FeatureAssayPlot1']] <- list()
all_active[['ColumnDataPlot1']] <- list()
all_active[['RowDataPlot1']] <- list()
all_active[['SampleAssayPlot1']] <- list()

################################################################################
## Reduced dimension plot 1
################################################################################


red.dim <- reducedDim(se, "PCA");
plot.data <- data.frame(X=red.dim[, 1], Y=red.dim[, 2], row.names=colnames(se));

# Avoid visual biases from default ordering by shuffling the points
set.seed(1221);
plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];

dot.plot <- ggplot() +
    geom_point(aes(x=X, y=Y), alpha=1, plot.data, color='#000000', size=1) +
    labs(x="Dimension 1", y="Dimension 2", title="PCA") +
    coord_cartesian(xlim=range(plot.data$X, na.rm=TRUE),
                    ylim=range(plot.data$Y, na.rm=TRUE), expand=TRUE) +
    theme_bw() +
    theme(legend.position='bottom', legend.box='vertical', legend.text=element_text(size=9), legend.title=element_text(size=11),
          axis.text=element_text(size=10), axis.title=element_text(size=12), title=element_text(size=12))

# Saving data for transmission
all_contents[['ReducedDimensionPlot1']] <- plot.data

################################################################################
## Row data table 1
################################################################################


tab <- as.data.frame(rowData(se));

# Saving data for transmission
all_contents[['RowDataTable1']] <- tab

################################################################################
## Feature assay plot 1
################################################################################


plot.data <- data.frame(Y=assay(se, "logcounts")["S100A9", ], row.names=colnames(se))
plot.data$X <- colData(se)[, "cluster_louvain"];

plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y,
                                                width=0.4, varwidth=FALSE, adjust=1,
                                                method='quasirandom', nbins=NULL);

# Avoid visual biases from default ordering by shuffling the points
set.seed(1221);
plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];

dot.plot <- ggplot() +
    geom_violin(aes(x=X, y=Y, group=GroupBy), alpha=0.2, data=plot.data, scale='width', width=0.8) +
    geom_point(aes(y=Y, x=jitteredX), alpha=1, plot.data, color='#000000', size=1) +
    labs(x="cluster_louvain", y="S100A9 (logcounts)", title="S100A9 vs cluster_louvain") +
    coord_cartesian(ylim=range(plot.data$Y, na.rm=TRUE), expand=TRUE) +
    scale_x_discrete(drop=FALSE) +
    theme_bw() +
    theme(legend.position='bottom', legend.text=element_text(size=9),
          legend.title=element_text(size=11), legend.box='vertical',
          axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=12), title=element_text(size=12))

# Saving data for transmission
all_contents[['FeatureAssayPlot1']] <- plot.data

################################################################################
## Column data plot 1
################################################################################


plot.data <- data.frame(Y=colData(se)[, "Sample"], row.names=colnames(se));
plot.data$X <- factor(character(ncol(se)))

plot.data[["Y"]] <- factor(plot.data[["Y"]]);

set.seed(100);
j.out <- iSEE:::jitterSquarePoints(plot.data$X, plot.data$Y);
summary.data <- j.out$summary;
plot.data$jitteredX <- j.out$X;
plot.data$jitteredY <- j.out$Y;

# Avoid visual biases from default ordering by shuffling the points
set.seed(1221);
plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];

dot.plot <- ggplot(plot.data) +
    geom_tile(aes(x=X, y=Y, height=2*YWidth, width=2*XWidth, group=interaction(X, Y)),
              summary.data, color='black', alpha=0, size=0.5) +
    geom_point(aes(x=jitteredX, y=jitteredY), alpha=1, plot.data, color='#000000', size=1) +
    labs(x="", y="Sample", title="Sample ") +
    scale_x_discrete(drop=FALSE) +
    scale_y_discrete(drop=FALSE) +
    theme_bw() +
    theme(legend.position='bottom', legend.text=element_text(size=9),
          legend.title=element_text(size=11), legend.box='vertical',
          axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=12), title=element_text(size=12))

# Saving data for transmission
all_contents[['ColumnDataPlot1']] <- plot.data

################################################################################
## Row data plot 1
################################################################################


plot.data <- data.frame(Y=rowData(se)[, "ID"], row.names=rownames(se));
plot.data$X <- factor(character(nrow(se)))

plot.data$Y <- as.numeric(as.factor(plot.data$Y));

plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y,
                                                width=0.4, varwidth=FALSE, adjust=1,
                                                method='quasirandom', nbins=NULL);

# Avoid visual biases from default ordering by shuffling the points
set.seed(33538);
plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];

dot.plot <- ggplot() +
    geom_violin(aes(x=X, y=Y, group=GroupBy), alpha=0.2, data=plot.data, scale='width', width=0.8) +
    geom_point(aes(y=Y, x=jitteredX), alpha=1, plot.data, color='#000000', size=1) +
    labs(x="", y="ID", title="ID ") +
    coord_cartesian(ylim=range(plot.data$Y, na.rm=TRUE), expand=TRUE) +
    scale_x_discrete(drop=FALSE) +
    theme_bw() +
    theme(legend.position='bottom', legend.text=element_text(size=9),
          legend.title=element_text(size=11), legend.box='vertical',
          axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=12), title=element_text(size=12))

# Saving data for transmission
all_contents[['RowDataPlot1']] <- plot.data

################################################################################
## Sample assay plot 1
################################################################################


plot.data <- data.frame(Y=assay(se, "logcounts")[,"AAACCCAAGGAGAGTA-1"], row.names=rownames(se));
plot.data$X <- factor(character(nrow(se)));

plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y,
                                                width=0.4, varwidth=FALSE, adjust=1,
                                                method='quasirandom', nbins=NULL);

# Avoid visual biases from default ordering by shuffling the points
set.seed(33538);
plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];

dot.plot <- ggplot() +
    geom_violin(aes(x=X, y=Y, group=GroupBy), alpha=0.2, data=plot.data, scale='width', width=0.8) +
    geom_point(aes(y=Y, x=jitteredX), alpha=1, plot.data, color='#000000', size=1) +
    labs(x="", y="AAACCCAAGGAGAGTA-1 (logcounts)", title="AAACCCAAGGAGAGTA-1") +
    coord_cartesian(ylim=range(plot.data$Y, na.rm=TRUE), expand=TRUE) +
    scale_x_discrete(drop=FALSE) +
    theme_bw() +
    theme(legend.position='bottom', legend.text=element_text(size=9),
          legend.title=element_text(size=11), legend.box='vertical',
          axis.text.x=element_text(angle=90, size=10, hjust=1, vjust=0.5),
          axis.text.y=element_text(size=10),
          axis.title=element_text(size=12), title=element_text(size=12))

# Saving data for transmission
all_contents[['SampleAssayPlot1']] <- plot.data

################################################################################
## Column data table 1
################################################################################


tab <- as.data.frame(colData(se));

# Saving data for transmission
all_contents[['ColumnDataTable1']] <- tab

################################################################################
## Complex heatmap 1
################################################################################


.heatmap.rows <- "MIR1302-2HG";
.heatmap.columns <- colnames(se);
plot.data <- assay(se, "logcounts")[.heatmap.rows, .heatmap.columns, drop=FALSE]
plot.data <- as.matrix(plot.data);

.assay_colors <- assayColorMap(colormap, "logcounts", discrete=FALSE)(21L)
.assay_colors <- circlize::colorRamp2(breaks = seq(0, 1, length.out = 21L), colors = .assay_colors)

hm <- ComplexHeatmap::Heatmap(matrix=plot.data, col=.assay_colors,
                              cluster_rows=FALSE, cluster_columns=FALSE, name="logcounts",
                              show_row_names=TRUE, show_column_names=FALSE,
                              heatmap_legend_param=list(direction="horizontal"))

ComplexHeatmap::draw(hm, heatmap_legend_side="bottom", annotation_legend_side="bottom")

# Saving data for transmission
all_contents[['ComplexHeatmapPlot1']] <- plot.data

################################################################################
## To guarantee the reproducibility of your code, you should also
## record the output of sessionInfo()
sessionInfo()
