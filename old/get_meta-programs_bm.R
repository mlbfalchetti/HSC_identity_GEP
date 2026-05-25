library(Seurat)

#--
load(file = "/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/results_programs_1.5/cnmf_top_100_adult_bm.RData")

load("C:/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/seurat_obj_1.3/sk_2022_adult_bm_01_hpc_young.RData")
my_list_01 <- list()
my_list_02 <- list()
for (i in 1:10) {
  seurat <- AddModuleScore(seurat, features = list(my_list$sk_2022_adult_bm_01_hpc_young$n_components_10[, i]), name = "target")
  test <- t.test(seurat@meta.data[which(seurat@meta.data$hsc == "hsc"), ]$target1, 
                 seurat@meta.data[which(seurat@meta.data$hsc != "hsc"), ]$target1)
  my_list_01[[i]] <- as.numeric(test$statistic)
  my_list_02[[i]] <- as.numeric(test$p.value)
}
data.frame(statistic = unlist(my_list_01), 
           p = unlist(my_list_02))
# 1

#--
load("C:/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/seurat_obj_1.3/sk_2022_adult_bm_01_hpc_young.RData")
my_list_01 <- list()
my_list_02 <- list()
for (i in 1:length(mp_adult_bm)) {
  seurat <- AddModuleScore(seurat, features = list(mp_adult_bm[[i]]$gene), name = "target")
  test <- t.test(seurat@meta.data[which(seurat@meta.data$hsc == "hsc"), ]$target1, 
                 seurat@meta.data[which(seurat@meta.data$hsc != "hsc"), ]$target1)
  my_list_01[[i]] <- as.numeric(test$statistic)
  my_list_02[[i]] <- as.numeric(test$p.value)
}
data.frame(statistic = unlist(my_list_01), 
           p = unlist(my_list_02))
# 1

seurat <- AddModuleScore(seurat, features = list(mp_adult_bm[[1]]$gene), name = "target")
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
load("C:/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/seurat_obj_1.3/sk_2022_adult_bm_02_hpc_young.RData")
my_list_01 <- list()
my_list_02 <- list()
for (i in 1:length(mp_adult_bm)) {
  seurat <- AddModuleScore(seurat, features = list(mp_adult_bm[[i]]$gene), name = "target")
  test <- t.test(seurat@meta.data[which(seurat@meta.data$hsc == "hsc"), ]$target1, 
                 seurat@meta.data[which(seurat@meta.data$hsc != "hsc"), ]$target1)
  my_list_01[[i]] <- as.numeric(test$statistic)
  my_list_02[[i]] <- as.numeric(test$p.value)
}
data.frame(statistic = unlist(my_list_01), 
           p = unlist(my_list_02))
# 1

seurat <- AddModuleScore(seurat, features = list(mp_adult_bm[[1]]$gene), name = "target")
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
load(file = "/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/results_meta-programs_1.5/mp_adult_bm_old.Rdata")

load("C:/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/seurat_obj_1.3/sk_2022_adult_bm_01_hpc_old.RData")
DefaultAssay(seurat) <- "SCT"
my_list_01 <- list()
my_list_02 <- list()
for (i in 1:length(mp_adult_bm)) {
  seurat <- AddModuleScore(seurat, features = list(mp_adult_bm[[i]]$gene), name = "target")
  test <- t.test(seurat@meta.data[which(seurat@meta.data$hsc == "hsc"), ]$target1, 
                 seurat@meta.data[which(seurat@meta.data$hsc != "hsc"), ]$target1)
  my_list_01[[i]] <- as.numeric(test$statistic)
  my_list_02[[i]] <- as.numeric(test$p.value)
}
data.frame(statistic = unlist(my_list_01), 
           p = unlist(my_list_02))
# 1

seurat <- AddModuleScore(seurat, features = list(mp_adult_bm[[1]]$gene), name = "target")
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

#-
load("C:/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/seurat_obj_1.3/sk_2022_adult_bm_02_hpc_old.RData")
DefaultAssay(seurat) <- "SCT"
my_list_01 <- list()
my_list_02 <- list()
for (i in 1:length(mp_adult_bm)) {
  seurat <- AddModuleScore(seurat, features = list(mp_adult_bm[[i]]$gene), name = "target")
  test <- t.test(seurat@meta.data[which(seurat@meta.data$hsc == "hsc"), ]$target1, 
                 seurat@meta.data[which(seurat@meta.data$hsc != "hsc"), ]$target1)
  my_list_01[[i]] <- as.numeric(test$statistic)
  my_list_02[[i]] <- as.numeric(test$p.value)
}
data.frame(statistic = unlist(my_list_01), 
           p = unlist(my_list_02))
# 1

seurat <- AddModuleScore(seurat, features = list(mp_adult_bm[[1]]$gene), name = "target")
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
load(file = "/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/results_meta-programs_1.5/mp_adult_bm.Rdata")
set_i <- mp_adult_bm$meta_1$gene
load(file = "/Users/marce/OneDrive/Área de Trabalho/bch/m_falchetti/results_meta-programs_1.5/mp_adult_bm_old.Rdata")
set_j <- mp_adult_bm$meta_1$gene
intersection <- length(intersect(set_i, set_j))
union <- length(union(set_i, set_j))
intersection / union # 0.30
