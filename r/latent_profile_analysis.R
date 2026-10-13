# Latent profile analysis (Gaussian mixture modelling) of the behavioral scores of the
# General Stroke Cohort and the Pure Aphasia Cohort.
#
# Steps (Supplementary Section S5):
#   1. read the 30 behavioral variables (OCS and NIHSS subtests) of both cohorts
#   2. invert the variables in `invert_cols` so that higher = better performance
#   3. z-score both cohorts with the mean and SD of the General Stroke Cohort
#   4. fit Gaussian mixture models with 1 to `max_clusters` clusters (mclust), choose by BIC
#   5. Sankey diagrams: cohort -> cluster, and clinical syndrome -> cluster (Pure Aphasia Cohort)
#   6. radar plots of the behavioral profile of each cluster
# Optional (not part of the original analysis): adjusted Rand index between clinical
# syndromes and clusters, with a permutation test.
#
# No patient data are included in this repository. Edit the settings below.

# ---------------------------------------------------------------- settings
cfg <- list(
  data_dir         = ".",                              # folder with the two Excel files
  general_file     = "general_stroke_behavior.xlsx",   # General Stroke Cohort
  aphasia_file     = "pure_aphasia_behavior.xlsx",     # Pure Aphasia Cohort
  general_cols     = 4:33,   # columns of the 30 behavioral variables (General Stroke Cohort file)
  aphasia_cols     = 9:39,   # columns of the Pure Aphasia Cohort file that contain the variables
  aphasia_drop_col = 19,     # column (within aphasia_cols) not used  -> DOCUMENT WHY
  aphasia_drop_row = 5,      # patient removed from the Pure Aphasia Cohort file -> DOCUMENT WHY
  label_column     = "label",# clinical syndrome of each Pure Aphasia Cohort patient
  invert_cols      = 16:30,  # variables inverted so that higher = better performance
  missing          = "zero", # "zero": missing values set to 0 (original analysis)
                             # "drop": patients with missing values are excluded
  max_clusters     = 5,      # number of clusters tried (1 to max_clusters); try 10 as a check
  seed             = 1,
  output_dir       = "lpa_output"
)

# ---------------------------------------------------------------- packages
need <- function(p) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(c("readxl", "mclust", "networkD3", "htmlwidgets", "jsonlite", "fmsb", "scales"), need))
suppressPackageStartupMessages({
  library(readxl); library(mclust); library(networkD3)
  library(htmlwidgets); library(jsonlite); library(fmsb)
})
set.seed(cfg$seed)
dir.create(cfg$output_dir, showWarnings = FALSE)

# ---------------------------------------------------------------- 1) data
gen  <- read_excel(file.path(cfg$data_dir, cfg$general_file))
diss <- read_excel(file.path(cfg$data_dir, cfg$aphasia_file))

if (!is.null(cfg$aphasia_drop_row)) diss <- diss[-cfg$aphasia_drop_row, ]
syndrome <- gsub("_", " ", diss[[cfg$label_column]])
diss <- diss[, cfg$aphasia_cols]
if (!is.null(cfg$aphasia_drop_col)) diss <- diss[, -cfg$aphasia_drop_col]
gen  <- gen[, cfg$general_cols]
stopifnot(ncol(gen) == ncol(diss), ncol(gen) == 30)

to_matrix <- function(x) apply(as.data.frame(x), 2, as.numeric)
gen_mat  <- to_matrix(gen)
diss_mat <- to_matrix(diss)
var_names <- names(gen)

message("Missing values: General Stroke Cohort = ", sum(is.na(gen_mat)),
        ", Pure Aphasia Cohort = ", sum(is.na(diss_mat)))
if (cfg$missing == "zero") {
  gen_mat[is.na(gen_mat)]   <- 0
  diss_mat[is.na(diss_mat)] <- 0
} else if (cfg$missing == "drop") {
  keep_g <- complete.cases(gen_mat); keep_d <- complete.cases(diss_mat)
  gen_mat <- gen_mat[keep_g, , drop = FALSE]
  diss_mat <- diss_mat[keep_d, , drop = FALSE]
  syndrome <- syndrome[keep_d]
} else stop("cfg$missing must be 'zero' or 'drop'")

# ---------------------------------------------------------------- 2) invert
for (i in cfg$invert_cols) {
  m <- max(gen_mat[, i], na.rm = TRUE)       # maximum of the General Stroke Cohort
  gen_mat[, i]  <- m - gen_mat[, i]
  diss_mat[, i] <- m - diss_mat[, i]
}

# ---------------------------------------------------------------- 3) standardize
mu    <- colMeans(gen_mat)
sigma <- apply(gen_mat, 2, sd)
gen_z  <- scale(gen_mat,  center = mu, scale = sigma)
diss_z <- scale(diss_mat, center = mu, scale = sigma)
gen_z[is.na(gen_z)]   <- 0                    # variables without variance
diss_z[is.na(diss_z)] <- 0

X <- as.data.frame(rbind(gen_z, diss_z))
names(X) <- var_names
cohort <- c(rep("General Stroke Cohort", nrow(gen_z)), rep("Pure Aphasia Cohort", nrow(diss_z)))
message("Patients: General Stroke Cohort = ", nrow(gen_z), ", Pure Aphasia Cohort = ", nrow(diss_z))

# ---------------------------------------------------------------- 4) mixture models
fit <- Mclust(X, G = 1:cfg$max_clusters)
print(summary(fit))                           # model type, number of clusters, BIC
if (fit$G == cfg$max_clusters)
  warning("The best BIC is at the upper limit of the search (", cfg$max_clusters,
          " clusters): increase max_clusters to check that the solution is stable.")

pdf(file.path(cfg$output_dir, "BIC.pdf")); plot(fit, what = "BIC"); dev.off()
cluster <- fit$classification
write.csv(aggregate(X, by = list(Cluster = cluster), FUN = mean),
          file.path(cfg$output_dir, "cluster_means.csv"), row.names = FALSE)
write.csv(data.frame(cohort = cohort, cluster = cluster),
          file.path(cfg$output_dir, "cluster_assignments.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 5) Sankey diagrams
# 5a) cohort -> cluster
flows <- as.data.frame(table(source = cohort, target = paste("Cluster", cluster)))
names(flows) <- c("source", "target", "value")
flows$source <- as.character(flows$source); flows$target <- as.character(flows$target)
nodes <- data.frame(name = unique(c(flows$source, flows$target)))
flows$IDsource <- match(flows$source, nodes$name) - 1
flows$IDtarget <- match(flows$target, nodes$name) - 1
sankey_cohorts <- sankeyNetwork(Links = flows, Nodes = nodes, Source = "IDsource",
                                Target = "IDtarget", Value = "value", NodeID = "name",
                                fontSize = 14, nodeWidth = 30)
saveWidget(sankey_cohorts, file.path(cfg$output_dir, "sankey_cohorts.html"), selfcontained = TRUE)

# 5b) clinical syndrome -> cluster (Pure Aphasia Cohort only)
idx <- which(cohort == "Pure Aphasia Cohort")
stopifnot(length(idx) == length(syndrome))
df <- data.frame(stage1 = "Pure Aphasia Cohort",
                 stage2 = paste0("Aphasia ", syndrome),
                 stage3 = paste0("Cluster ", cluster[idx]))
flow1 <- aggregate(list(value = rep(1, nrow(df))), by = list(from = df$stage1, to = df$stage2), FUN = sum)
flow2 <- aggregate(list(value = rep(1, nrow(df))), by = list(from = df$stage2, to = df$stage3), FUN = sum)
all_flows <- rbind(flow1, flow2)
nodes2 <- data.frame(name = unique(c(all_flows$from, all_flows$to)))
all_flows$source <- match(all_flows$from, nodes2$name) - 1
all_flows$target <- match(all_flows$to,   nodes2$name) - 1
node_colors <- sapply(nodes2$name, function(x) {
  if (x == "Pure Aphasia Cohort") "red" else if (grepl("^Aphasia ", x)) "green" else if (grepl("^Cluster ", x)) "orange" else "grey"
})
sankey_syndromes <- sankeyNetwork(
  Links = all_flows, Nodes = nodes2, Source = "source", Target = "target",
  Value = "value", NodeID = "name", sinksRight = FALSE, fontSize = 12, nodeWidth = 30,
  colourScale = JS(sprintf("d3.scaleOrdinal().range(%s)", jsonlite::toJSON(node_colors))))
saveWidget(sankey_syndromes, file.path(cfg$output_dir, "sankey_syndromes.html"), selfcontained = TRUE)

# ---------------------------------------------------------------- optional: agreement between syndromes and clusters
ari <- adjustedRandIndex(syndrome, cluster[idx])
perm <- replicate(10000, adjustedRandIndex(sample(syndrome), cluster[idx]))
p_ari <- (1 + sum(perm >= ari)) / (1 + length(perm))
message(sprintf("Adjusted Rand index (syndromes vs clusters) = %.3f, permutation p = %.4f", ari, p_ari))

# ---------------------------------------------------------------- 6) radar plots
profiles <- aggregate(X, by = list(Cluster = cluster), FUN = mean)
prof <- profiles[, -1]
normalize <- function(x) (x - min(x)) / (max(x) - min(x) + 1e-8)   # per variable, across clusters
prof_norm <- as.data.frame(lapply(prof, normalize))
pdf(file.path(cfg$output_dir, "cluster_profiles.pdf"), width = 7, height = 7)
for (k in seq_len(nrow(prof_norm))) {
  df_plot <- rbind(rep(1, ncol(prof_norm)), rep(0, ncol(prof_norm)), prof_norm[k, ])
  rownames(df_plot) <- c("max", "min", paste("Cluster", profiles$Cluster[k]))
  par(mar = c(2, 2, 4, 2))
  radarchart(df_plot, axistype = 1, pcol = "darkred", pfcol = scales::alpha("darkred", 0.4),
             plwd = 3, title = paste("Cluster", profiles$Cluster[k], "behavioral profile"),
             cglcol = "grey80", cglty = 1, cglwd = 0.8, axislabcol = "grey40",
             vlcex = 0.8, seg = 5, caxislabels = seq(0, 1, 0.2))
}
dev.off()
