library(dplyr)
library(ggplot2)
library(Maaslin2)
library(SIAMCAT)
library(caret)
library(Boruta)
library(tidyr)
library(vegan)
library(zCompositions)
library(rlist)
library(ggrepel)


# TODO: Load domain/signalp data
# TODO: try deleting unused columns
# signalp_res <- read.table("~/Documents/lab/pgn_hydrolases/annotations/data/db_annotations/signalp/pooled_prediction_results.tsv.csv",
#                           sep = "\t", fill = TRUE)

# TODO: apply HUMANn3/DB to healthy people, maybe from different parts of the world

relab_abun_df <- function(abun_df) {

  rownames(abun_df) <- abun_df$Gene.Family
  abun_df$Gene.Family <- NULL

  colnames(abun_df) <- unlist(lapply(colnames(abun_df), function(col_name) {
    strsplit(col_name, "_")[[1]][1]
  }))

  data.frame(t(abun_df))

}

unstrat_abun_df <- function(abun_df) {
  abun_df[, grep(".", colnames(abun_df), fixed = TRUE, invert = TRUE)]
}

strat_abun_df <- function(abun_df) {
  abun_df[, c(1, grep(".", colnames(abun_df), fixed = TRUE))]
}

pooled_abun_df <- function(abun_df) {

  classes <- unlist(lapply(colnames(abun_df), function(col_name) {
    strsplit(col_name, "_")[[1]][1]
  }))

  pooled_abun_df <- aggregate(t(abun_df), by = list(classes), FUN = sum)

  rownames(pooled_abun_df) <- pooled_abun_df$Group.1
  pooled_abun_df$Group.1 <- NULL

  pooled_abun_df <- data.frame(t(pooled_abun_df))

  colnames(pooled_abun_df) <- unlist(lapply(colnames(pooled_abun_df), function(col_name) {
    sub(".", "-", col_name, fixed = TRUE)
  }))

  pooled_abun_df

}

# domain_pooled_abun_df <- function(abun_df) {
#
# }

# Apply CZM and CLR to a count df
# - count_df: the count df that is to be transformed, with samples by row
czm_clr <- function(count_df) {

  set.seed(101)
  czm <- cmultRepl(count_df, method = "CZM",
                   suppress.print = TRUE, z.warning = 1)
  clr <- t(apply(czm, 1, function(x){log(x) - mean(log(x))}))

  clr

}

study_wise_permanova <- function(count_df, meta_df) {
  count_df <- count_df[rownames(meta_df), ]

  set.seed(101)
  adonis2(czm_clr(count_df[, colSums(count_df) > 0]) ~ meta_df[, "status"],
          method = "euclidean")
}



status_maaslin2 <- function(count_df, meta_df, maaslin2_reference) {

  set.seed(101)
  maaslin2_res <- Maaslin2(input_data = count_df,
                           input_metadata = meta_df,
                           output = "tmp",
                           normalization = "CLR", transform = "none",
                           fixed_effects = "status",
                           cores = 6, #plot_heatmap = FALSE, plot_scatter = FALSE,
                           reference = maaslin2_reference)[[1]]

  maaslin2_res <- maaslin2_res[, c("feature", "coef", "stderr", "pval")]
  colnames(maaslin2_res) <- c("feature", "estimate", "std_error", "p_val")

  maaslin2_res$adj_p_val <- p.adjust(maaslin2_res$p_val, method = "BH")

  maaslin2_res

}

rf <- function(count_df, meta_df) {

  labels <- meta_df[rownames(count_df), "status"]
  labels[labels == "No GvHD"] <- "No_GvHD"

  const_mtry <- floor(sqrt(ncol(count_df)))

  caret::train(x = count_df,
               y = labels,
               method = "rf", ntree = 501,
               metric = "ROC", tuneGrid = data.frame(mtry = const_mtry),
               trControl = trainControl(method = "repeatedcv", number = 2, repeats = 5,
                                        summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE))

}

boruta <- function(count_df, meta_df) {

  set.seed(101)
  Boruta(czm_clr(count_df),
         as.factor(meta_df[rownames(count_df), "status"]),
         maxRuns = 1000)

}

delta_time_plot <- function(count_df, meta_df) {

  plotting_df <- merge(count_df[, colnames(count_df) != "UNMAPPED"],
                       meta_df[, c("patient", "status", "day")], by = 0)

  plotting_df$Row.names <- NULL

  plotting_df <- pivot_longer(plotting_df, !c("patient", "status", "day"),
                              names_to = "pg_hydrolase",
                              values_to = "abun")

  ggplot(plotting_df, aes(x = day, y = abun, color = status)) +
    geom_line(aes(group = patient)) +
    geom_point() +
    labs(x = "Day Post-HSCT", y = "CLR Relative Abundance") +
    facet_wrap(~pg_hydrolase, scales = "free") +
    theme_linedraw() +
    theme(panel.grid = element_blank(),
          # legend.key.size = unit("0.35", "cm"),
          legend.position = c(0.775, 0.14),
          strip.text = element_text(color = "black",
                                    margin = margin(t = 0, r = 0, b = 2, l = 0, unit = "pt")),
          strip.background = element_rect(fill = "white", color = "white"))

}

compute_deltas <- function(count_df, meta_df, base_status, alt_status) {

  patients <- unique(meta_df$patient)

  delta_count_df <- list.rbind(lapply(patients, function(patient) {

    print(patient)
    initial_t <- min(meta_df[meta_df$patient == patient, "day"])
    final_t <- max(meta_df[meta_df$patient == patient, "day"])

    count_df[rownames(meta_df[meta_df$patient == patient & meta_df$day == final_t, ]), ] -
      count_df[rownames(meta_df[meta_df$patient == patient & meta_df$day == initial_t, ]), ]

  }))
  rownames(delta_count_df) <- patients

  plotting_df <- merge(delta_count_df[, colnames(delta_count_df) != "UNMAPPED"],
                       unique(meta_df[c("patient", "status")]), by.x = 0,
                       by.y = "patient")

  plotting_df$Row.names <- NULL

  plotting_df <- pivot_longer(plotting_df, !"status",
                              names_to = "pg_hydrolase",
                              values_to = "delta")

  label_df <- data.frame(list.rbind(lapply(unique(plotting_df$pg_hydrolase), function(pg_hydrolase) {
    list(pg_hydrolase,
         wilcox.test(plotting_df$delta[plotting_df$pg_hydrolase == pg_hydrolase & plotting_df$status == alt_status],
                     plotting_df$delta[plotting_df$pg_hydrolase == pg_hydrolase & plotting_df$status == base_status])$p.value)
  })))

  colnames(label_df) <- c("pg_hydrolase", "delta_p_val")
  label_df$pg_hydrolase <- as.character(label_df$pg_hydrolase)
  label_df$delta_p_val <- as.numeric(label_df$delta_p_val)

  list(plotting_df, label_df)

}


pca_biplot <- function(count_df, meta_df, scaling_factor, labels) {

  count_df <- count_df[rownames(meta_df), ]

  clr <- czm_clr(count_df[, colSums(count_df) > 0])
  pca <- prcomp(clr)

  sample_positions <- data.frame(pca[["x"]])
  sample_positions <- merge(sample_positions,
                            meta_df[, "status", drop = FALSE],
                            by = 0)

  feature_positions <- data.frame(pca[["rotation"]])

  ggplot() +
    geom_point(data = sample_positions, aes(x = PC1, y = PC2, fill = status),
               pch = 21, alpha = 0.8, size = 2.5, stroke = 0.3) +
    geom_segment(data = feature_positions, aes(x = 0, y = 0, xend = scaling_factor * PC1, yend = scaling_factor * PC2),
                 arrow = arrow(length = unit(0.2, 'picas')),
                 color = "grey74",
                 size = 0.3) +
    geom_text_repel(data = feature_positions[labels, ], aes(x = scaling_factor * PC1, y = scaling_factor * PC2),
                    label = rownames(feature_positions[labels, ]), force_pull = 100, force = 10, point.size = NA,
                    color = "black",
                    size = 2.5) +
    labs(x = paste0("PC1: ", round(pca$sdev[1]^2/sum(pca$sdev^2),3) * 100, "%"),
         y = paste0("PC2: ", round(pca$sdev[2]^2/sum(pca$sdev^2),3) * 100, "%")) +
    theme_linedraw() +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          legend.key.size = unit("0.35", "cm"),
          legend.position = "bottom",
          legend.margin = margin(-0.25, 0, 0, -0.2, unit = "cm"),
          legend.justification = "left")

}

maaslin2_plot <- function(l_clr_dr, maaslin2_df, label_x_coord) {

  ggplot(l_clr_dr[l_clr_dr$pg_hydrolase != "UNMAPPED", ], aes(x = abun, y = pg_hydrolase)) +
    geom_point(aes(color = status, group = status), alpha = 0.8,  size = 2, stroke = 0.3,
               position = position_jitterdodge(jitter.height = 0, jitter.width = 0.15)) +
    stat_summary(aes(group = status),fun.min = function(x) {quantile(x, 0.25)}, fun.max = function(x) {quantile(x, 0.75)},
                 size = 0.2, width = 0.3, geom = "errorbar", position = position_dodge(width = 0.8)) +
    stat_summary(aes(group = status), fun.min = median, fun.max = median,
                 size = 0.2, width = 0.6, geom = "errorbar", position = position_dodge(width = 0.8)) +
    geom_label(data = maaslin2_df[maaslin2_df$feature != "UNMAPPED", ],
               aes(x = label_x_coord, y = feature, label = paste0("italic('P')==", round(adj_p_val, 3))), size = 2.5, parse = TRUE,
               hjust = 0, label.padding = unit(0.25, "lines"), label.size = NA, fill = NA) +
    scale_color_manual(values = c("#cf5154", "#4e88aa"), name = NULL) +
    labs(x = "CLR Relative Abundance", y = NULL) +
    theme_linedraw() +
    theme(panel.grid = element_blank(),
          legend.key.size = unit("0.35", "cm"),
          legend.position = "bottom",
          legend.margin = margin(-0.25, 0, 0, -0.2, unit = "cm"),
          legend.justification = "left")

}



##### PRJNA283642 #####
PRJNA283642_meta_df <- read.csv("/home/ccarr/Documents/lab/pgn_hydrolases/shotgun/data/gvhd_meta/PRJNA283642_meta.csv")

PRJNA283642_abun_df <- read.table("/home/ccarr/Documents/lab/pgn_hydrolases/shotgun/data/gvhd_humann/PRJNA283642_humann_genefamilies_relab.tsv",
                                  sep = "\t", header = TRUE)

PRJNA283642_abun_df <- relab_abun_df(PRJNA283642_abun_df)

tmp_PRJNA283642_abun_df <- aggregate(PRJNA283642_abun_df[PRJNA283642_meta_df$sample, ], by = list(PRJNA283642_meta_df$patient), FUN = mean)

rownames(tmp_PRJNA283642_abun_df) <- tmp_PRJNA283642_abun_df$Group.1
tmp_PRJNA283642_abun_df$Group.1 <- NULL

PRJNA283642_abun_df <- tmp_PRJNA283642_abun_df
rm(tmp_PRJNA283642_abun_df)

PRJNA283642_meta_df <- unique(PRJNA283642_meta_df[, c("patient", "status")])
rownames(PRJNA283642_meta_df) <- PRJNA283642_meta_df$patient
PRJNA283642_meta_df <- PRJNA283642_meta_df[rownames(PRJNA283642_abun_df), ]

# PRJNA283642_unstrat_abun_df <- unstrat_abun_df(PRJNA283642_abun_df)
# PRJNA283642_strat_abun_df <- strat_abun_df(PRJNA283642_abun_df)
PRJNA283642_pooled_abun_df <- pooled_abun_df(PRJNA283642_abun_df)
# PRJNA283642_domain_pooled_abun_df <- domain_pooled_abun_df(PRJNA283642_abun_df)

study_wise_permanova(PRJNA283642_pooled_abun_df, PRJNA283642_meta_df)

tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/gvhd_late_pca_plot.tiff",
     units = "in", width = 3, height = 3, res = 600)
b <- pca_biplot(PRJNA283642_pooled_abun_df, PRJNA283642_meta_df, 4,
           c("DL-endopeptidase", "LD-endopeptidase", "Muramidase", "Glucosaminidase", "LD-carboxypeptidase")) +
  scale_fill_manual(values = c("#cf5154", "#4e88aa"), name = NULL) +
  geom_label(aes(x = -Inf, y = Inf, label = "italic('R')^2 == 0.174~(italic('P')==0.065)"), size = 2.5, parse = TRUE,
             hjust = 0, vjust = 1, label.padding = unit(0.25, "lines"), label.size = NA, fill = NA)
dev.off()

# MaAsLin2
# TODO: try other parameters
# TODO: check all transforms are OK on rel. abun. data
# + = higher in "No GvHD"
PRJNA283642_pooled_maaslin2_df <- status_maaslin2(PRJNA283642_pooled_abun_df, PRJNA283642_meta_df, "status,No GvHD")
PRJNA283642_pooled_maaslin2_df$feature <- unlist(lapply(PRJNA283642_pooled_maaslin2_df$feature, function(col_name) {
  sub(".", "-", col_name, fixed = TRUE)
}))

PRJNA283642_pooled_clr_df <- data.frame(czm_clr(PRJNA283642_pooled_abun_df))
colnames(PRJNA283642_pooled_clr_df) <- unlist(lapply(colnames(PRJNA283642_pooled_clr_df), function(col_name) {
  sub(".", "-", col_name, fixed = TRUE)
}))

tmp_PRJNA283642_pooled_clr_df <- PRJNA283642_pooled_clr_df
tmp_PRJNA283642_pooled_clr_df$patient <- rownames(tmp_PRJNA283642_pooled_clr_df)
l_PRJNA283642_pooled_clr_df <- pivot_longer(tmp_PRJNA283642_pooled_clr_df, !"patient",
                                            names_to = "pg_hydrolase",
                                            values_to = "abun")
rm(tmp_PRJNA283642_pooled_clr_df)

l_PRJNA283642_pooled_clr_df <- merge(l_PRJNA283642_pooled_clr_df, PRJNA283642_meta_df,
                                     by = "patient")

l_PRJNA283642_pooled_clr_df$pg_hydrolase <- factor(l_PRJNA283642_pooled_clr_df$pg_hydrolase,
                                                   levels = rev(unique(l_PRJNA283642_pooled_clr_df$pg_hydrolase)))

tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/gvhd_late_maaslin2_plot.tiff",
     units = "in", width = 4, height = 3, res = 600)
c <- maaslin2_plot(l_PRJNA283642_pooled_clr_df, PRJNA283642_pooled_maaslin2_df, 3.5) +
  xlim(NA, 5.75)
dev.off()

#

PRJNA283642_pooled_rf_res <- rf(PRJNA283642_pooled_abun_df, PRJNA283642_meta_df)

library(plotROC)
tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/gvhd_late_auroc_plot.tiff",
     units = "in", width = 3.5, height = 3, res = 600)
a <- ggplot() +
  geom_abline(slope = 1, intercept = 0, size = 0.75, color = "grey74", linetype = "dashed") +
  geom_roc(data = PRJNA283642_pooled_rf_res$pred, aes(m = No_GvHD, d = obs),
           n.cuts = 0, size = 0.75, color = "#4e88aa") +
  labs(x = "True Positive Rate", y = "False Positive Rate",
       caption = "AUROC = 0.887") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.0015))) +
  scale_y_continuous(expand = expansion(mult = c(0.0015, 0.01))) +
  theme_linedraw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank())
dev.off()

PRJNA283642_pooled_boruta_res <- boruta(PRJNA283642_pooled_abun_df, PRJNA283642_meta_df)

PRJNA283642_pooled_boruta_imp_df <- data.frame(PRJNA283642_pooled_boruta_res[["ImpHistory"]])

PRJNA283642_pooled_boruta_imp_df <- PRJNA283642_pooled_boruta_imp_df[, !colnames(PRJNA283642_pooled_boruta_imp_df) %in%
                                                                       c("UNMAPPED", "shadowMin", "shadowMean", "shadowMax")]
colnames(PRJNA283642_pooled_boruta_imp_df) <- unlist(lapply(colnames(PRJNA283642_pooled_boruta_imp_df), function(col_name) {
  sub(".", "-", col_name, fixed = TRUE)
}))

PRJNA283642_pooled_boruta_imp_df$run <- rownames(PRJNA283642_pooled_boruta_imp_df)

l_PRJNA283642_pooled_boruta_imp_df <- pivot_longer(PRJNA283642_pooled_boruta_imp_df,
                                                   !"run",
                                                   names_to = "pg_hydrolase",
                                                   values_to = "importance")

l_PRJNA283642_pooled_boruta_imp_df <- l_PRJNA283642_pooled_boruta_imp_df[l_PRJNA283642_pooled_boruta_imp_df$importance != "-Inf", ]
l_PRJNA283642_pooled_boruta_imp_df$importance <- as.numeric(l_PRJNA283642_pooled_boruta_imp_df$importance)

mean_PRJNA283642_pooled_boruta_imp_df <- l_PRJNA283642_pooled_boruta_imp_df %>%
  group_by(pg_hydrolase) %>%
  summarise(mean = mean(importance))

l_PRJNA283642_pooled_boruta_imp_df$pg_hydrolase <- factor(l_PRJNA283642_pooled_boruta_imp_df$pg_hydrolase,
                                                          levels = mean_PRJNA283642_pooled_boruta_imp_df$pg_hydrolase[order(mean_PRJNA283642_pooled_boruta_imp_df$mean)])



tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/gvhd_late_boruta_plot.tiff",
     units = "in", width = 3, height = 3, res = 600)
ggplot(l_PRJNA283642_pooled_boruta_imp_df, aes(x = importance, y = pg_hydrolase)) +
  geom_vline(xintercept = 0, size = 0.2, color = "grey74", linetype = "dashed") +
  stat_boxplot(geom = "errorbar", width = 0.25, size = 0.3) +
  geom_boxplot(color = "black", size = 0.3, width = 0.55,
               outlier.shape = NA) +
  labs(x = "Importance", y = NULL) +
  theme_linedraw() +
  theme(panel.grid = element_blank())
dev.off()

library(cowplot)
tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/fig_2.tiff",
     units = "in", width = 10.5, height = 3, res = 600)
plot_grid(a, b, c, ncol = 3, rel_widths = c(3.5, 3, 4),
          labels = c("a", "b", "c"), label_size = 12)
dev.off()

##### PRJNA525982 #####
PRJNA525982_meta_df <- read.csv("/home/ccarr/Documents/lab/pgn_hydrolases/shotgun/data/gvhd_meta/PRJNA525982_meta.csv")

rownames(PRJNA525982_meta_df) <- PRJNA525982_meta_df$sample
PRJNA525982_meta_df$sample <- NULL

PRJNA525982_meta_df$day <- unlist(lapply(PRJNA525982_meta_df$timepoint, function(timepoint) {
  if (timepoint == "Pre hematopoietic stem cell transplantation") {
    -1
  } else {
    as.numeric(strsplit(timepoint, " ")[[1]][1])
  }
}))

PRJNA525982_abun_df <- read.table("/home/ccarr/Documents/lab/pgn_hydrolases/shotgun/data/gvhd_humann/PRJNA525982_humann_genefamilies_relab.tsv",
                                  sep = "\t", header = TRUE)

PRJNA525982_abun_df <- relab_abun_df(PRJNA525982_abun_df)

# PRJNA525982_unstrat_abun_df <- unstrat_abun_df(PRJNA525982_abun_df)
# PRJNA525982_strat_abun_df <- strat_abun_df(PRJNA525982_abun_df)
PRJNA525982_pooled_abun_df <- pooled_abun_df(PRJNA525982_abun_df)
# PRJNA525982_domain_pooled_abun_df <- domain_pooled_abun_df(PRJNA525982_abun_df)

study_wise_permanova(PRJNA525982_pooled_abun_df,
                     PRJNA525982_meta_df[PRJNA525982_meta_df$timepoint == "Pre hematopoietic stem cell transplantation", ])

tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/gvhd_baseline_pca_plot.tiff",
     units = "in", width = 3, height = 3, res = 600)
pca_biplot(PRJNA525982_pooled_abun_df,
         PRJNA525982_meta_df[PRJNA525982_meta_df$timepoint ==
                               "Pre hematopoietic stem cell transplantation", ],
         3, colnames(PRJNA525982_pooled_abun_df)[colnames(PRJNA525982_pooled_abun_df) != "UNMAPPED"]) +
  scale_fill_manual(values = c("#cf5154", "#4e88aa"), name = NULL) +
  geom_label(aes(x = Inf, y = Inf, label = "italic('R')^2 == 0.463~(italic('P')==0.029)"), size = 2.5, parse = TRUE,
             hjust = 1, vjust = 1, label.padding = unit(0.25, "lines"), label.size = NA, fill = NA)
dev.off()

PRJNA525982_pooled_maaslin2_df <- status_maaslin2(PRJNA525982_pooled_abun_df, PRJNA525982_meta_df, "status,No GvHD")
PRJNA525982_pooled_maaslin2_df$feature <- unlist(lapply(PRJNA525982_pooled_maaslin2_df$feature, function(col_name) {
  sub(".", "-", col_name, fixed = TRUE)
}))

PRJNA525982_pooled_clr_df <- data.frame(czm_clr(PRJNA525982_pooled_abun_df))
colnames(PRJNA525982_pooled_clr_df) <- unlist(lapply(colnames(PRJNA525982_pooled_clr_df), function(col_name) {
  sub(".", "-", col_name, fixed = TRUE)
}))

tmp_PRJNA525982_pooled_clr_df <- PRJNA525982_pooled_clr_df
tmp_PRJNA525982_pooled_clr_df$sample <- rownames(tmp_PRJNA525982_pooled_clr_df)
l_PRJNA525982_pooled_clr_df <- pivot_longer(tmp_PRJNA525982_pooled_clr_df, !"sample",
                                            names_to = "pg_hydrolase",
                                            values_to = "abun")
rm(tmp_PRJNA525982_pooled_clr_df)

l_PRJNA525982_pooled_clr_df <- merge(l_PRJNA525982_pooled_clr_df, PRJNA525982_meta_df,
                                     by.x = "sample", by.y = 0)

l_PRJNA525982_pooled_clr_df$pg_hydrolase <- factor(l_PRJNA525982_pooled_clr_df$pg_hydrolase,
                                                   levels = rev(unique(l_PRJNA525982_pooled_clr_df$pg_hydrolase)))

tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/gvhd_baseline_maaslin2_plot.tiff",
     units = "in", width = 4, height = 3, res = 600)
maaslin2_plot(l_PRJNA525982_pooled_clr_df[l_PRJNA525982_pooled_clr_df$timepoint == "Pre hematopoietic stem cell transplantation", ],
              PRJNA525982_pooled_maaslin2_df, 1.3) +
  xlim(NA, 2.65)
dev.off()

#

tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/gvhd_delta_time_plot.tiff",
     units = "in", width = 6, height = 4, res = 600)
d <- delta_time_plot(PRJNA525982_pooled_clr_df, PRJNA525982_meta_df) +
  scale_color_manual(values = c("#cf5154", "#4e88aa"), name = NULL)
dev.off()

PRJNA525982_pooled_delta_dfs <- compute_deltas(PRJNA525982_pooled_clr_df, PRJNA525982_meta_df,
                                               "No GvHD", "GvHD")

tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/gvhd_delta_boxplot.tiff",
     units = "in", width = 1.75, height = 2.25, res = 600)
e <- ggplot(PRJNA525982_pooled_delta_dfs[[1]][PRJNA525982_pooled_delta_dfs[[1]]$pg_hydrolase == "DL-endopeptidase", ],
       aes(x = status, y = delta)) +
  stat_boxplot(geom = "errorbar", width = 0.25, size = 0.3) +
  geom_boxplot(color = "black", size = 0.3, width = 0.55,
               outlier.shape = NA) +
  geom_jitter(aes(fill = status), shape = 21, stroke = 0.3, size = 2.5,
              height = 0, width = 0.1) +
  geom_segment(aes(x = 1, xend = 2, y = 2.75, yend = 2.75), size = 0.2) +
  geom_text(data = PRJNA525982_pooled_delta_dfs[[2]][PRJNA525982_pooled_delta_dfs[[2]]$pg_hydrolase == "DL-endopeptidase", ],
             aes(x = 1.5, y = 3, label = paste0("italic('P') == ", round(delta_p_val, 3))), parse = TRUE,
             size = 2.5) +
  scale_fill_manual(values = c("#cf5154", "#4e88aa")) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  labs(x = NULL, y = expression(Delta~"CLR Relative Abundance")) +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/fig_3.tiff",
     units = "in", width = 7.75, height = 3.75, res = 600)
plot_grid(d, e, ncol = 2, rel_widths = c(6, 1.75),
          labels = c("a", "b"), label_size = 12)
dev.off()

##### PRJEB40960 #####
PRJEB40960_meta_df <- read.csv("/home/ccarr/Documents/lab/pgn_hydrolases/shotgun/data/gvhd_meta/PRJEB40960_meta.csv")

PRJEB40960_meta_df$patient.sample <- unlist(lapply(PRJEB40960_meta_df$patient.sample, function(patient_sample) {
  sub("a|b", "", patient_sample)
}))

PRJEB40960_patient_meta_df <- read.csv("/home/ccarr/Documents/lab/pgn_hydrolases/shotgun/data/gvhd_meta/PRJEB40960_patient_meta.csv")
PRJEB40960_time_meta_df <- read.csv("/home/ccarr/Documents/lab/pgn_hydrolases/shotgun/data/gvhd_meta/PRJEB40960_time_meta.csv")

PRJEB40960_meta_df <- merge(PRJEB40960_meta_df, PRJEB40960_time_meta_df,
                            by = "patient.sample", all.x = TRUE)

PRJEB40960_meta_df$patient <- unlist(lapply(PRJEB40960_meta_df$patient.sample, function(patient_sample) {
  sub("p0", "", strsplit(patient_sample, ".", fixed = TRUE)[[1]][1])
}))
PRJEB40960_meta_df$sample <- unlist(lapply(PRJEB40960_meta_df$patient.sample, function(patient_sample) {
  strsplit(patient_sample, ".", fixed = TRUE)[[1]][2]
}))

PRJEB40960_meta_df <- merge(PRJEB40960_meta_df, PRJEB40960_patient_meta_df,
                            by = "patient", all.x = TRUE)

PRJEB40960_meta_df$sample <- NULL
PRJEB40960_meta_df$patient_day <- paste0(PRJEB40960_meta_df$patient, "_", PRJEB40960_meta_df$day)

colnames(PRJEB40960_meta_df)[3] <- "sample"
colnames(PRJEB40960_meta_df)[6] <- "status"

PRJEB40960_abun_df <- read.table("/home/ccarr/Documents/lab/pgn_hydrolases/shotgun/data/gvhd_humann/PRJEB40960_humann_genefamilies_relab.tsv",
                                  sep = "\t", header = TRUE)

PRJEB40960_abun_df <- relab_abun_df(PRJEB40960_abun_df)

tmp_PRJEB40960_abun_df <- aggregate(PRJEB40960_abun_df[PRJEB40960_meta_df$sample, ], by = list(PRJEB40960_meta_df$patient_day), FUN = mean)

rownames(tmp_PRJEB40960_abun_df) <- tmp_PRJEB40960_abun_df$Group.1
tmp_PRJEB40960_abun_df$Group.1 <- NULL

PRJEB40960_abun_df <- tmp_PRJEB40960_abun_df
rm(tmp_PRJEB40960_abun_df)

PRJEB40960_meta_df <- unique(PRJEB40960_meta_df[, !colnames(PRJEB40960_meta_df) %in% c("sample", "patient.sample")])

colnames(PRJEB40960_meta_df)[14] <- "sample"
rownames(PRJEB40960_meta_df) <- PRJEB40960_meta_df$sample

# Pool UC118 with DL-endos here and elsewhere

# permanova (maybe last pre tx, first post tx, and last?)
# PCA plots to explore (patient and instrument esp), and for response
# delta plots to see if loss of dls is even before gvhd (and look at other classes too)
# maaslin2 for response
# rf for response

# PRJEB40960_unstrat_abun_df <- unstrat_abun_df(PRJEB40960_abun_df)
# PRJEB40960_strat_abun_df <- strat_abun_df(PRJEB40960_abun_df)
PRJEB40960_pooled_abun_df <- pooled_abun_df(PRJEB40960_abun_df)
# PRJEB40960_domain_pooled_abun_df <- domain_pooled_abun_df(PRJEB40960_abun_df)

PRJEB40960_pooled_abun_df$tmp_dl_endopeptidase <- PRJEB40960_pooled_abun_df$`DL-endopeptidase` + PRJEB40960_pooled_abun_df$UC118
PRJEB40960_pooled_abun_df$`DL-endopeptidase` <- NULL
PRJEB40960_pooled_abun_df$UC118 <- NULL
colnames(PRJEB40960_pooled_abun_df)[9] <- "DL-endopeptidase"

PRJEB40960_pooled_abun_df <- PRJEB40960_pooled_abun_df[, c(1:3, 9, 4:8)]

# # study_wise_permanova(PRJEB40960_pooled_abun_df,
# #                      PRJEB40960_meta_df[PRJEB40960_meta_df$timepoint == "Pre hematopoietic stem cell transplantation", ])
#
# # tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/gvhd_baseline_pca_plot.tiff",
# #      units = "in", width = 3, height = 3, res = 600)
# # pca_biplot(PRJEB40960_pooled_abun_df,
# #            PRJEB40960_meta_df[PRJEB40960_meta_df$timepoint ==
# #                                  "Pre hematopoietic stem cell transplantation", ],
# #            3, colnames(PRJEB40960_pooled_abun_df)[colnames(PRJEB40960_pooled_abun_df) != "UNMAPPED"]) +
# #   scale_fill_manual(values = c("#cf5154", "#4e88aa"), name = NULL) +
# #   geom_label(aes(x = Inf, y = Inf, label = "italic('R')^2 == 0.463~(italic('P')==0.029)"), size = 2.5, parse = TRUE,
# #              hjust = 1, vjust = 1, label.padding = unit(0.25, "lines"), label.size = NA, fill = NA)
# # dev.off()
#
# # PRJEB40960_pooled_maaslin2_df <- status_maaslin2(PRJEB40960_pooled_abun_df, PRJEB40960_meta_df, "status,No GvHD")
# # PRJEB40960_pooled_maaslin2_df$feature <- unlist(lapply(PRJEB40960_pooled_maaslin2_df$feature, function(col_name) {
# #   sub(".", "-", col_name, fixed = TRUE)
# # }))

PRJEB40960_pooled_clr_df <- data.frame(czm_clr(PRJEB40960_pooled_abun_df))
colnames(PRJEB40960_pooled_clr_df) <- unlist(lapply(colnames(PRJEB40960_pooled_clr_df), function(col_name) {
  sub(".", "-", col_name, fixed = TRUE)
}))

tmp_PRJEB40960_pooled_clr_df <- PRJEB40960_pooled_clr_df
tmp_PRJEB40960_pooled_clr_df$sample <- rownames(tmp_PRJEB40960_pooled_clr_df)
l_PRJEB40960_pooled_clr_df <- pivot_longer(tmp_PRJEB40960_pooled_clr_df, !"sample",
                                           names_to = "pg_hydrolase",
                                           values_to = "abun")
rm(tmp_PRJEB40960_pooled_clr_df)

l_PRJEB40960_pooled_clr_df <- merge(l_PRJEB40960_pooled_clr_df, PRJEB40960_meta_df,
                                    by.x = "sample", by.y = 0)

l_PRJEB40960_pooled_clr_df$pg_hydrolase <- factor(l_PRJEB40960_pooled_clr_df$pg_hydrolase,
                                                  levels = rev(unique(l_PRJEB40960_pooled_clr_df$pg_hydrolase)))

plotting_df <- merge(PRJEB40960_pooled_clr_df[, colnames(PRJEB40960_pooled_clr_df) != "UNMAPPED"],
                     PRJEB40960_meta_df[, c("patient", "status", "day")], by = 0)

plotting_df$Row.names <- NULL

plotting_df <- pivot_longer(plotting_df, !c("patient", "status", "day"),
                            names_to = "pg_hydrolase",
                            values_to = "abun")
plotting_df <- plotting_df[plotting_df$pg_hydrolase == "DL-endopeptidase", ]
plotting_df <- plotting_df[!plotting_df$patient %in% c(21, 42, 61, 70), ]

gvhd_onset_df <- data.frame(patient = c("16", "19", "27", "40", "44"),
                            day = c(35, 45, 78, 109, 21))

# TODO: make by defining patient and day cols, then patient_day, then merge to get abun
delta_points_df <- data.frame(patient = rep(c("16", "19", "27", "40", "44"), each = 2),
                               abun = c(-0.2752777, -1.9542991, -1.3241578, -2.4201503,
                                        -1.5129521, -2.5271854, -1.7428919, -1.8104829,
                                        -1.4155316, -1.4424114),
                               day = c(-1, 28, -27, 48, -17, 28, -3, 39, -3, 1))

tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/test.tiff",
     units = "in", width = 6, height = 1.75, res = 2000)
g <- ggplot() +
  geom_vline(data = gvhd_onset_df, aes(xintercept = day),
             linetype = "dashed", size = 0.3, alpha = 0.75) +
  geom_line(data = plotting_df, color = "#cf5154", aes(x = day, y = abun, group = patient)) +
  geom_point(data = plotting_df, color = "#cf5154", aes(x = day, y = abun)) +
  geom_point(data = delta_points_df, aes(x = day, y = abun), color = "black", shape = 8, size = 0.75, stroke = 0.5) +
  labs(x = "Day Post-HSCT", y = "CLR Relative\nAbundance") +
  facet_wrap(~patient, scales = "free_x", nrow = 1) +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        strip.text = element_text(color = "black",
                                  margin = margin(t = 0, r = 0, b = 2, l = 0, unit = "pt")),
        strip.background = element_rect(fill = "white", color = "white"))
dev.off()

delta_df <- as.data.frame(list.rbind(lapply(unique(delta_points_df$patient), function(patient) {
  tmp_df <- delta_points_df[delta_points_df$patient == patient, ]
  c(patient, as.numeric(tmp_df$abun[2] - tmp_df$abun[1]))
})))
colnames(delta_df) <- c("patient", "delta")
delta_df$delta <- as.numeric(delta_df$delta)

tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/test2.tiff",
     units = "in", width = 1.75, height = 2, res = 600)
f <- ggplot(data = delta_df, aes(x = patient, y = delta)) +
  geom_col(color = "#cf5154", fill = "#cf5154", width = 0.6) +
  theme_linedraw() +
  labs(x = "Patient", y = expression(Delta~"CLR Rel. Abun.")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0)),
                     breaks = c(-0.5, -1, -1.5)) +
  theme(panel.grid = element_blank())
dev.off()

tiff("~/Documents/lab/pgn_hydrolases/shotgun/figures/fig_4.tiff",
     units = "in", width = 7.75, height = 2, res = 2000)
plot_grid(f, g, ncol = 2, rel_widths = c(1.75, 6),
          labels = c("a", "b"), label_size = 12)
dev.off()

# TODO: are patients with shortest follow-up the ones with the least delta?
# TODO: can analyze which specific dl-endos (taxonomically and in terms of secretion/domain)

##### PRJEB10878 #####




