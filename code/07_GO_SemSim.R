library(GOSemSim)
library(ape)
library(ggtree)
library(clusterProfiler)
library(reshape2)
library(ggplot2)

human_genes.df = read.csv(file = '../data/08_GO_gene_subset/ages_lists/human_genes_id.csv',
                       col.names = c('SYMBOL', 'ENTREZID'),
                       colClasses = "character")
eg <- human_genes.df['ENTREZID'][, 1]

# MF
hsGO_MF <- godata('org.Hs.eg.db', ont="MF", computeIC = FALSE)


sim_MF <- mgeneSim(eg, semData = hsGO_MF, measure = "Wang", drop = "IEA", combine = "BMA")
rownames(sim_MF) <- human_genes[rownames(sim_MF)]
colnames(sim_MF) <- human_genes[colnames(sim_MF)]

write.table(sim_MF, file="../data/07_GO/human_sim_MF.txt")

sim_MF <- as.matrix(read.table("../data/07_GO/human_sim_MF.txt", header=TRUE, row.names=1))
melted_sim_MF <- melt(sim_MF)

png("../data/07_GO_plots/human_MF_01.png", width = 5200, height = 5000, units = "px")
heatmap(sim_MF, Rowv = NA)
dev.off()
gc()
#DOSE::simplot(sim_MF) # not enough space

png("../data/07_GO_plots/human_MF_02.png")
ggplot(melted_sim_MF, aes(x=Var1, y=Var2, fill=value)) +  
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "red", limit = c(0,1))
ggsave("../data/07_GO_plots/human_MF_02.png", width = 5200, height = 5000, units = "px")
gc()

# delete
rm(hsGO_MF)
rm(melted_sim_MF)
rm(sim_MF)

#############################################################################
# HUMAN BP

hsGO_BP <- godata('org.Hs.eg.db', ont="BP", computeIC = FALSE)
sim_BP <- mgeneSim(eg, semData = hsGO_BP, measure = "Wang", drop = "IEA", combine = "BMA")
rownames(sim_BP) <- human_genes[rownames(sim_BP)]
colnames(sim_BP) <- human_genes[colnames(sim_BP)]

write.table(sim_BP, file="../data/07_GO/human_sim_BP.txt")

## read and plot
sim_BP <- as.matrix(read.table("../data/07_GO/human_sim_BP.txt", header=TRUE, row.names=1))
melted_sim_BP <- melt(sim_BP)

# Hierarchical Cluster Analysis
dist_BP <- dist(sim_BP)
# https://uc-r.github.io/hc_clustering#algorithms
hc.complete <- hclust(dist_BP, method = 'complete') 
hc.ward.D2 <- hclust(dist_BP, method = "ward.D2")
## save
dend.complete <- as.phylo(hc.complete)
dend.ward.D2 <- as.phylo(hc.ward.D2)
write.tree(dend.complete, file = "../data/07_GO/dendrograms/humanBP_dend_complete.tree")
write.tree(dend.ward.D2, file = "../data/07_GO/dendrograms/humanBP_dend_ward.tree")

# plots
png("../data/07_GO_plots/human_BP_01.png", width = 5200, height = 5000, units = "px")
heatmap(sim_BP, Rowv = NA)
dev.off()
gc()
#DOSE::simplot(sim_BP) # not enough space

png("../data/07_GO_plots/human_BP_02.png")
ggplot(melted_sim_BP, aes(x=Var1, y=Var2, fill=value)) +  
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "red", limit = c(0,1))
ggsave("../data/07_GO_plots/human_BP_02.png", width = 5200, height = 5000, units = "px")
gc()

rm(hsGO_BP)
rm(melted_sim_BP)
rm(sim_BP)

#############################################################################
# HUMAN CC

hsGO_CC <- godata('org.Hs.eg.db', ont="CC", computeIC = FALSE)
sim_CC <- mgeneSim(eg, semData = hsGO_CC, measure = "Wang", drop = "IEA", combine = "BMA")
rownames(sim_CC) <- human_genes[rownames(sim_CC)]
colnames(sim_CC) <- human_genes[colnames(sim_CC)]

write.table(sim_CC, file="../data/07_GO/human_sim_CC.txt")

sim_CC <- as.matrix(read.table("../data/07_GO/human_sim_CC.txt", header=TRUE, row.names=1))
melted_sim_CC <- melt(sim_CC)

png("../data/07_GO_plots/human_CC_01.png", width = 5200, height = 5000, units = "px")
heatmap(sim_CC, Rowv = NA)
dev.off()
gc()
#DOSE::simplot(sim_CC) # not enough space

png("../data/07_GO_plots/human_CC_02.png")
ggplot(melted_sim_CC, aes(x=Var1, y=Var2, fill=value)) +  
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "red", limit = c(0,1))
ggsave("../data/07_GO_plots/human_CC_02.png", width = 5200, height = 5000, units = "px")
gc()

rm(hsGO_CC)
rm(melted_sim_CC)
rm(sim_CC)

#############################################################################
#############################################################################
#############################################################################
#############################################################################

# MOUSE MF
mouse_genes.df = read.csv(file = '../data/08_GO_gene_subset/ages_lists/mouse_genes_id.csv',
                          col.names = c('SYMBOL', 'ENTREZID'),
                          colClasses = "character")
eg <- mouse_genes.df['ENTREZID'][, 1]

# MF
hsGO_MF <- godata('org.Mm.eg.db', ont="MF", computeIC = FALSE)


sim_MF <- mgeneSim(eg, semData = hsGO_MF, measure = "Wang", drop = "IEA", combine = "BMA")
rownames(sim_MF) <- mouse_genes[rownames(sim_MF)]
colnames(sim_MF) <- mouse_genes[colnames(sim_MF)]

write.table(sim_MF, file="../data/07_GO/mouse_sim_MF.txt")

sim_MF <- as.matrix(read.table("../data/07_GO/mouse_sim_MF.txt", header=TRUE, row.names=1))
melted_sim_MF <- melt(sim_MF)

png("../data/07_GO_plots/mouse_MF_01.png", width = 5200, height = 5000, units = "px")
heatmap(sim_MF, Rowv = NA)
dev.off()
gc()
#DOSE::simplot(sim_MF) # not enough space

png("../data/07_GO_plots/mouse_MF_02.png")
ggplot(melted_sim_MF, aes(x=Var1, y=Var2, fill=value)) +  
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "red", limit = c(0,1))
ggsave("../data/07_GO_plots/mouse_MF_02.png", width = 5200, height = 5000, units = "px")
gc()

# delete
rm(hsGO_MF)
rm(melted_sim_MF)
rm(sim_MF)

#############################################################################
# mouse BP

hsGO_BP <- godata('org.Mm.eg.db', ont="BP", computeIC = FALSE)
sim_BP <- mgeneSim(eg, semData = hsGO_BP, measure = "Wang", drop = "IEA", combine = "BMA")
rownames(sim_BP) <- mouse_genes[rownames(sim_BP)]
colnames(sim_BP) <- mouse_genes[colnames(sim_BP)]

write.table(sim_BP, file="../data/07_GO/mouse_sim_BP.txt")

sim_BP_mouse <- as.matrix(read.table("../data/07_GO/mouse_sim_BP.txt", header=TRUE, row.names=1))
melted_sim_BP <- melt(sim_BP_mouse)

# Hierarchical Cluster Analysis
dist_BP <- dist(sim_BP)
# https://uc-r.github.io/hc_clustering#algorithms
hc.complete <- hclust(dist_BP, method = 'complete') 
hc.ward.D2 <- hclust(dist_BP, method = "ward.D2")
## save
dend.complete <- as.phylo(hc.complete)
dend.ward.D2 <- as.phylo(hc.ward.D2)
write.tree(dend.complete, file = "../data/07_GO/dendrograms/mouseBP_dend_complete.tree")
write.tree(dend.ward.D2, file = "../data/07_GO/dendrograms/mouseBP_dend_ward.tree")

# plots
png("../data/07_GO_plots/mouse_BP_01.png", width = 5200, height = 5000, units = "px")
heatmap(sim_BP, Rowv = NA)
dev.off()
gc()
#DOSE::simplot(sim_BP) # not enough space

png("../data/07_GO_plots/mouse_BP_02.png")
ggplot(melted_sim_BP, aes(x=Var1, y=Var2, fill=value)) +  
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "red", limit = c(0,1))
ggsave("../data/07_GO_plots/mouse_BP_02.png", width = 5200, height = 5000, units = "px")
gc()

rm(hsGO_BP)
rm(melted_sim_BP)
rm(sim_BP)

#############################################################################
# mouse CC

hsGO_CC <- godata('org.Mm.eg.db', ont="CC", computeIC = FALSE)
sim_CC <- mgeneSim(eg, semData = hsGO_CC, measure = "Wang", drop = "IEA", combine = "BMA")
rownames(sim_CC) <- mouse_genes[rownames(sim_CC)]
colnames(sim_CC) <- mouse_genes[colnames(sim_CC)]

write.table(sim_CC, file="../data/07_GO/mouse_sim_CC.txt")

sim_CC <- as.matrix(read.table("../data/07_GO/mouse_sim_CC.txt", header=TRUE, row.names=1))
melted_sim_CC <- melt(sim_CC)

png("../data/07_GO_plots/mouse_CC_01.png", width = 5200, height = 5000, units = "px")
heatmap(sim_CC, Rowv = NA)
dev.off()
gc()
#DOSE::simplot(sim_CC) # not enough space

png("../data/07_GO_plots/mouse_CC_02.png")
ggplot(melted_sim_CC, aes(x=Var1, y=Var2, fill=value)) +  
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "red", limit = c(0,1))
ggsave("../data/07_GO_plots/mouse_CC_02.png", width = 5200, height = 5000, units = "px")
gc()

rm(hsGO_CC)
rm(melted_sim_CC)
rm(sim_CC)




