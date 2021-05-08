chooseCRANmirror()
install.packages("BiocManager")
BiocManager::install("microbiome")
library(microbiome)
library(ggplot2)
library(dplyr)
library(IRanges)
peerj32 = "https://peerj.com/articles/32/"
data("peerj32")
pseq <- peerj32$phyloseq
p <- boxplot_abundance(pseq, x = "time", y = "Akkermansia", line = "subject") + scale_y_log10()
print(p)
dfs <- meta(pseq)
dfs$signal <- abundances(pseq)["Akkermansia", rownames(dfs)]
BiocManager::install("lme4")
library(lme4)
out <- lmer(signal ~ group + (1|subject), data = dfs)
out0 <- lmer(signal ~ (1|subject), data = dfs)
comp <- anova(out0, out)
pv <- comp[["Pr(>Chisq)"]][[2]]
print(pv)
pseq.rel <- microbiome::transform(pseq, "compositional")
data(peerj32)
pseq <- peerj32$phyloseq
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)
p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "group", size = 3)
print(p)
library(vegan)
permanova <- adonis(t(otu) ~ group, data = meta, permutations=99, method = "bray")
print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])
dist <- vegdist(t(otu))
anova(betadisper(dist, meta$group))
coef <- coefficients(permanova)["group1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")