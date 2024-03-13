library(corrplot)

# For mRNA and Protein  individual tabs:
# mat1 = # tumor samples (rows)
# mat2 = # normal samples (columns)

# For mRNA/Protein mix tab:
# mat1 = # mRNA (rows)
# mat2 = # Protein (columns)


cor.mat <- cor(mat1, mat2, use = "pairwise.complete.obs", method = "spearman")
p.mat <- cor.mtest(mat = cor.mat, method = "spearman")


# Bubble correlogram with significant astherisk annotation
par(bg = "transparent") # makes background transparent
corrplot(cor.mat,
  col = c("blue", "red"),
  tl.col = "black",
  method = "circle",
  # bg = "transparent",
  # addCoef.col = "darkgray",
  # number.cex = 0.6,
  type = "upper",
  tl.pos = "lt",
  p.mat = p.mat,
  sig.level = 0.01,
  insig = "label_sig",
  # pch = 4,
  pch.col = "white",
  pch.cex = 1,
  tl.srt = 45
)


# Numbers correlogram only showing significant rho values
par(bg = "transparent")
corrplot(cor.mat,
  method = "number",
  number.cex = 1,
  col = c("blue", "red"),
  tl.col = "black",
  type = "lower",
  tl.pos = "lt",
  p.mat = p.mat,
  sig.level = 0.01,
  insig = "blank",
  tl.srt = 45
)
