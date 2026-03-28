# file: phylo_pipeline_fixed.R

# ---- Install if needed ----
install.packages(c("ape", "seqinr", "ggplot2"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("Biostrings", "DECIPHER", "ggtree", "msa"))
install.packages("pheatmap")
library(pheatmap)

# ---- Load libraries ----
library(Biostrings)
library(DECIPHER)
library(seqinr)
library(ape)
library(ggtree)
library(ggplot2)
library(msa)

# ---- Read sequences ----
seqs <- readDNAStringSet("C:/Users/ANVITHA/Downloads/sequence.fasta")

# ---- Orientation (FIX: correct function) ----
seqs <- OrientNucleotides(seqs)

# ---- Alignment ----
aligned <- AlignSeqs(seqs)

# ---- Save alignment ----
writeXStringSet(aligned, file="Amanita_aligned.fasta")

# ---- Read alignment (seqinr format) ----
dna <- read.alignment("Amanita_aligned.fasta", format = "fasta")

# ---- Distance matrix ----
D <- dist.alignment(dna, matrix = "identity")  # FIX: "similarity" → "identity"

# ---- Clean matrix ----
temp <- as.matrix(D)
temp[is.na(temp)] <- 0
temp[is.infinite(temp)] <- 0
any(is.na(D))
D_mat <- as.matrix(D)


# small correction
D_mat[is.na(D_mat)] <- 1 
# ---- SAFE heatmap (FIX: remove table.paint) ----
pheatmap(temp)

# ---- Tree construction ----
tre <- nj(as.dist(D_mat))   # FIX: ensure correct type
tre <- ladderize(tre)

# ---- Base plots ----
plot.phylo(
  tre,
  cex = 0.7,
  no.margin = TRUE
)

plot(tre, cex = 0.6)
title("similarity in Amanita (ITS)")

# ---- Clustering ----
D_mat <- as.matrix(D)

D_mat[is.na(D_mat)] <- max(D_mat, na.rm = TRUE)
D_mat[is.infinite(D_mat)] <- max(D_mat, na.rm = TRUE)

D_clean <- as.dist(D_mat)

h_cluster <- hclust(D_clean, method = "average")
plot(h_cluster)
h_cluster <- hclust(D, method = "average", members = NULL)
# method = average is used for UPGMA, members can be equal to NULL or a vector with a length of size D

# ---- ggtree plots ----
p <- ggtree(tre) +
  geom_tiplab(hjust = -0.3, size = 3) +
  xlim(0, max(node.depth.edgelength(tre)) * 1.2)

print(p)

# ---- Highlight clusters (safe example nodes) ----
p2 <- ggtree(tre) +
  geom_tiplab(hjust = -0.3, size = 3) +
  geom_hilight(node = Nnode(tre) + 1, fill = "purple", alpha = 0.2)

print(p2)

# ---- Match labels safely ----
tre$tip.label <- names(seqs)

# ---- Alignment + tree ----
msaplot(
  p = ggtree(tre),
  fasta = "Amanita_aligned.fasta",
  window = c(1, 50)
)