load("/Users/james/Downloads/Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata")
sce_new <- sce[, sce$sample_name == 151673]
library(scater)
library(assertthat)

set.seed(101)
dlpfc <- scater::logNormCounts(sce_new)  # 生成 logcounts assay
dec <- scran::modelGeneVar(dlpfc)  # 此时默认使用 logcounts
top <- scran::getTopHVGs(dec, n = 2000)

set.seed(102)
dlpfc <- scater::runPCA(dlpfc, subset_row=top)

q <- 7  # Number of clusters
d <- 15  # Number of PCs
## Run BayesSpace clustering
colnames(colData(dlpfc))[4] <- "array_row"
colnames(colData(dlpfc))[5] <- "array_col"

source("spatialCluster.R")

system.time({
  cluster <- spatialCluster(dlpfc, q=q, "PCA", d=d, platform = "Visium",
                          nrep = 2000, burn.in = 1000, thin = 100,gamma = 3)
})
source("draw.R")
labels <- dplyr::recode(cluster$spatial.cluster, 3, 4, 5, 6, 2, 7, 1)
cdata <- data.frame(colData(cluster))
coord.multiplier <- list(
  x = 1,
  y = -1
)
vertices <- .make_hex_spots(cdata, labels, coord.multiplier)
splot <- ggplot(vertices, aes(x = x.vertex, y = y.vertex, group = spot, fill = factor(fill))) +
  geom_polygon() +           
  labs(fill = "Cluster") +                     # 图例标题
  coord_equal() +                              # 保持x/y轴比例一致
  theme_void() +                               # 移除背景和坐标轴
  scale_fill_viridis_d(option = "D")           # 使用 viridis 颜色方案（可选）
splot + scale_fill_viridis_d(option = "A", labels = 1:7) +
  labs(title="BayesSpace")

