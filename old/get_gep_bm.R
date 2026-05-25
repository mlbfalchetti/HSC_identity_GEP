library(data.table); library(igraph); library(dplyr)

#--
my_files <- list.files(
  path = "/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/cnmf_outputs_1.5/sk_2022_adult_bm/", 
  pattern = ".gene_spectra_score\\.k_5\\.dt_0_5\\.txt$", 
  recursive = TRUE, 
  full.names = TRUE
)
my_samples <- c("sk_2022_adult_bm_01_hpc_old", 
                "sk_2022_adult_bm_01_hpc_young", 
                "sk_2022_adult_bm_01_hsc_old", 
                "sk_2022_adult_bm_01_hsc_young", 
                "sk_2022_adult_bm_02_hpc_old", 
                "sk_2022_adult_bm_02_hpc_young", 
                "sk_2022_adult_bm_02_hsc_young")
k <- 10
my_list_02 <- list()
for (j in 1:length(my_samples)) {
  my_paths <- paste("/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/cnmf_outputs_1.5/sk_2022_adult_bm/", my_samples[j], ".gene_spectra_score.k_", k, ".dt_0_5.txt", sep = "")
  W <- fread(file = my_paths, header = TRUE, data.table = FALSE)
  rownames(W) <- W[[1]]
  W <- W[, -1]
  my_list_01 <- list()
  for (i in 1:k) {
    tmp <- data.frame(t(W[i, ]))
    colnames(tmp) <- "values"
    tmp$gene_symbol <- rownames(tmp) 
    tmp$values <- scale(tmp$values)
    my_list_01[[i]] <- sort(tmp[which(tmp$values > 2), ]$gene_symbol)
  }
  names(my_list_01) <- 1:10
  my_list_02[[j]] <- my_list_01
}
names(my_list_02) <- my_samples
my_list <- my_list_02

#--
all_programs <- unlist(my_list, recursive = FALSE)
patterns <- c(
  "sk_2022_adult_bm_01_hpc_young")
matches <- Reduce(`|`, lapply(patterns, function(p) grepl(p, names(all_programs), fixed = TRUE)))
all_programs <- all_programs[names(all_programs) %in% c(names(all_programs)[matches])]

load("C:/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/seurat_obj_1.3/sk_2022_adult_bm_01_hpc_young.RData")
my_list_01 <- list()
my_list_02 <- list()
for (i in 1:length(all_programs)) {
  seurat <- AddModuleScore(seurat, features = list(all_programs[[i]]), name = "target")
  test <- t.test(seurat@meta.data[which(seurat@meta.data$hsc == "hsc"), ]$target1, 
                 seurat@meta.data[which(seurat@meta.data$hsc != "hsc"), ]$target1)
  my_list_01[[i]] <- as.numeric(test$statistic)
  my_list_02[[i]] <- as.numeric(test$p.value)
}
tmp <- data.frame(name = names(all_programs), 
                  statistic = unlist(my_list_01), 
                  p = unlist(my_list_02))
tmp[which(tmp$statistic > 0 & tmp$p < 0.05),]

seurat <- AddModuleScore(seurat, features = list(all_programs$sk_2022_adult_bm_01_hpc_young.1), name = "target")
data <- FetchData(seurat, vars = c("hsc", "target1"))
p_value <- 0.05
my_n_01 <- list()
my_n_02 <- list()
for (n in c(10, 25, 50, 100)) {
  my_j_01 <- logical(100)
  my_j_02 <- logical(100)
  for (j in 1:100) {
    set.seed(j)
    subset_data <- data %>%
      group_by(hsc) %>%
      slice_sample(n = n, replace = TRUE) %>%
      ungroup()
    test <- t.test(
      subset_data[subset_data$hsc == "hsc", ]$`target1`,
      subset_data[subset_data$hsc != "hsc", ]$`target1`
    )
    my_j_01[j] <- as.numeric(test$statistic) > 0 & test$p.value < p_value
    my_j_02[j] <- as.numeric(test$statistic) < 0 & test$p.value < p_value
  }
  my_n_01[[as.character(n)]] <- sum(my_j_01, na.rm = TRUE)
  my_n_02[[as.character(n)]] <- sum(my_j_02, na.rm = TRUE)
}
my_n_01 <- my_n_01[c("10", "25", "50", "100")]
my_n_02 <- my_n_02[c("10", "25", "50", "100")]
unlist(my_n_01)

# Cluster 1
# all_programs$sk_2022_adult_bm_01_hpc_young.1

#--
all_programs <- unlist(my_list, recursive = FALSE)
patterns <- c(
  "sk_2022_adult_bm_02_hpc_young")
matches <- Reduce(`|`, lapply(patterns, function(p) grepl(p, names(all_programs), fixed = TRUE)))
all_programs <- all_programs[names(all_programs) %in% c(names(all_programs)[matches])]

load("C:/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/seurat_obj_1.3/sk_2022_adult_bm_01_hpc_young.RData")
my_list_01 <- list()
my_list_02 <- list()
for (i in 1:length(all_programs)) {
  seurat <- AddModuleScore(seurat, features = list(all_programs[[i]]), name = "target")
  test <- t.test(seurat@meta.data[which(seurat@meta.data$hsc == "hsc"), ]$target1, 
                 seurat@meta.data[which(seurat@meta.data$hsc != "hsc"), ]$target1)
  my_list_01[[i]] <- as.numeric(test$statistic)
  my_list_02[[i]] <- as.numeric(test$p.value)
}
tmp <- data.frame(name = names(all_programs), 
                  statistic = unlist(my_list_01), 
                  p = unlist(my_list_02))
tmp[which(tmp$statistic > 0 & tmp$p < 0.05),]

seurat <- AddModuleScore(seurat, features = list(all_programs$sk_2022_adult_bm_02_hpc_young.10), name = "target")
data <- FetchData(seurat, vars = c("hsc", "target1"))
p_value <- 0.05
my_n_01 <- list()
my_n_02 <- list()
for (n in c(10, 25, 50, 100)) {
  my_j_01 <- logical(100)
  my_j_02 <- logical(100)
  for (j in 1:100) {
    set.seed(j)
    subset_data <- data %>%
      group_by(hsc) %>%
      slice_sample(n = n, replace = TRUE) %>%
      ungroup()
    test <- t.test(
      subset_data[subset_data$hsc == "hsc", ]$`target1`,
      subset_data[subset_data$hsc != "hsc", ]$`target1`
    )
    my_j_01[j] <- as.numeric(test$statistic) > 0 & test$p.value < p_value
    my_j_02[j] <- as.numeric(test$statistic) < 0 & test$p.value < p_value
  }
  my_n_01[[as.character(n)]] <- sum(my_j_01, na.rm = TRUE)
  my_n_02[[as.character(n)]] <- sum(my_j_02, na.rm = TRUE)
}
my_n_01 <- my_n_01[c("10", "25", "50", "100")]
my_n_02 <- my_n_02[c("10", "25", "50", "100")]
unlist(my_n_01)

#--
all_programs <- unlist(my_list, recursive = FALSE)
patterns <- c(
  "sk_2022_adult_bm_01_hsc_young")
matches <- Reduce(`|`, lapply(patterns, function(p) grepl(p, names(all_programs), fixed = TRUE)))
all_programs <- all_programs[names(all_programs) %in% c(names(all_programs)[matches])]

load("C:/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/seurat_obj_1.3/sk_2022_adult_bm_01_hpc_young.RData")
my_list_01 <- list()
my_list_02 <- list()
for (i in 1:length(all_programs)) {
  seurat <- AddModuleScore(seurat, features = list(all_programs[[i]]), name = "target")
  test <- t.test(seurat@meta.data[which(seurat@meta.data$hsc == "hsc"), ]$target1, 
                 seurat@meta.data[which(seurat@meta.data$hsc != "hsc"), ]$target1)
  my_list_01[[i]] <- as.numeric(test$statistic)
  my_list_02[[i]] <- as.numeric(test$p.value)
}
tmp <- data.frame(name = names(all_programs), 
                  statistic = unlist(my_list_01), 
                  p = unlist(my_list_02))
tmp[which(tmp$statistic > 0 & tmp$p < 0.05),]

seurat <- AddModuleScore(seurat, features = list(all_programs$sk_2022_adult_bm_01_hsc_young.1), name = "target")
data <- FetchData(seurat, vars = c("hsc", "target1"))
p_value <- 0.05
my_n_01 <- list()
my_n_02 <- list()
for (n in c(10, 25, 50, 100)) {
  my_j_01 <- logical(100)
  my_j_02 <- logical(100)
  for (j in 1:100) {
    set.seed(j)
    subset_data <- data %>%
      group_by(hsc) %>%
      slice_sample(n = n, replace = TRUE) %>%
      ungroup()
    test <- t.test(
      subset_data[subset_data$hsc == "hsc", ]$`target1`,
      subset_data[subset_data$hsc != "hsc", ]$`target1`
    )
    my_j_01[j] <- as.numeric(test$statistic) > 0 & test$p.value < p_value
    my_j_02[j] <- as.numeric(test$statistic) < 0 & test$p.value < p_value
  }
  my_n_01[[as.character(n)]] <- sum(my_j_01, na.rm = TRUE)
  my_n_02[[as.character(n)]] <- sum(my_j_02, na.rm = TRUE)
}
my_n_01 <- my_n_01[c("10", "25", "50", "100")]
my_n_02 <- my_n_02[c("10", "25", "50", "100")]
unlist(my_n_01)

# Cluster 1
# all_programs$sk_2022_adult_bm_01_hsc_young.1

#--
all_programs <- unlist(my_list, recursive = FALSE)
patterns <- c(
  "sk_2022_adult_bm_02_hsc_young")
matches <- Reduce(`|`, lapply(patterns, function(p) grepl(p, names(all_programs), fixed = TRUE)))
all_programs <- all_programs[names(all_programs) %in% c(names(all_programs)[matches])]

load("C:/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/seurat_obj_1.3/sk_2022_adult_bm_01_hpc_young.RData")
my_list_01 <- list()
my_list_02 <- list()
for (i in 1:length(all_programs)) {
  seurat <- AddModuleScore(seurat, features = list(all_programs[[i]]), name = "target")
  test <- t.test(seurat@meta.data[which(seurat@meta.data$hsc == "hsc"), ]$target1, 
                 seurat@meta.data[which(seurat@meta.data$hsc != "hsc"), ]$target1)
  my_list_01[[i]] <- as.numeric(test$statistic)
  my_list_02[[i]] <- as.numeric(test$p.value)
}
tmp <- data.frame(name = names(all_programs), 
                  statistic = unlist(my_list_01), 
                  p = unlist(my_list_02))
tmp[which(tmp$statistic > 0 & tmp$p < 0.05),]

seurat <- AddModuleScore(seurat, features = list(all_programs$sk_2022_adult_bm_02_hsc_young.1), name = "target")
data <- FetchData(seurat, vars = c("hsc", "target1"))
p_value <- 0.05
my_n_01 <- list()
my_n_02 <- list()
for (n in c(10, 25, 50, 100)) {
  my_j_01 <- logical(100)
  my_j_02 <- logical(100)
  for (j in 1:100) {
    set.seed(j)
    subset_data <- data %>%
      group_by(hsc) %>%
      slice_sample(n = n, replace = TRUE) %>%
      ungroup()
    test <- t.test(
      subset_data[subset_data$hsc == "hsc", ]$`target1`,
      subset_data[subset_data$hsc != "hsc", ]$`target1`
    )
    my_j_01[j] <- as.numeric(test$statistic) > 0 & test$p.value < p_value
    my_j_02[j] <- as.numeric(test$statistic) < 0 & test$p.value < p_value
  }
  my_n_01[[as.character(n)]] <- sum(my_j_01, na.rm = TRUE)
  my_n_02[[as.character(n)]] <- sum(my_j_02, na.rm = TRUE)
}
my_n_01 <- my_n_01[c("10", "25", "50", "100")]
my_n_02 <- my_n_02[c("10", "25", "50", "100")]
unlist(my_n_01)

# Cluster 1
# all_programs$sk_2022_adult_bm_02_hsc_young.1

all_programs <- unlist(my_list, recursive = FALSE)
patterns <- c(
  "sk_2022_adult_bm_01_hpc_young",
  "sk_2022_adult_bm_01_hsc_young",
  "sk_2022_adult_bm_02_hsc_young")
matches <- Reduce(`|`, lapply(patterns, function(p) grepl(p, names(all_programs), fixed = TRUE)))
all_programs <- all_programs[names(all_programs) %in% c(names(all_programs)[matches])]
df <- data.frame(table(c(all_programs$sk_2022_adult_bm_01_hpc_young.1,
                         all_programs$sk_2022_adult_bm_01_hsc_young.1,
                         all_programs$sk_2022_adult_bm_02_hsc_young.1)))
features <- sort(df$Var1)
save(features, file = "C:/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/useful_data_1.5/identity_gep_bms.RData")

# hk_genes_01 <- sort(unique(trimws(sub("\\t.*", "", readLines("http://www.tau.ac.il/~elieis/HKG/HK_genes.txt")))))
# features <- setdiff(features, hk_genes_01)
# length(features)
load("/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/useful_data/Housekeeping_Genes_Human_HTA.RData")
hk_genes_02 <- Housekeeping_Genes
features <- setdiff(features, unique(sort(as.character(hk_genes_02$Gene.name))))
length(features)
save(features, file = "/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/useful_data_1.5/identity_gep_bms_without_hks.RData")
