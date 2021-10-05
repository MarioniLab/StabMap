# install StabMap package
# where tokenstring is your own personal access token from GitHub
# (group access to MarioniLab required)

# devtools::install_github("MarioniLab/StabMap",
#                          ref="main",
#                          auth_token = "tokenstring")

# vignette example code
############################
library(StabMap)

# simulate some data
full_expr = matrix(rnorm(300*150), nrow = 300, ncol = 150,
                   dimnames = list(paste0("gene_", 1:300),
                                   paste0("cell_", 1:150)))

# build list of discrete assays with non-overlapping features
assay_list = list(
  D_R = full_expr[1:150, 1:50],
  D_j = full_expr[76:215, 51:100],
  D_i = full_expr[151:300, 101:150]
)

# assign labels to one group of cells
labels_list = list(
  D_R = rep(letters[1:5], each = 10)
)

# whether the data should be treated as reference
reference_list = list(
  D_R = TRUE,
  D_j = FALSE,
  D_i = TRUE
)

# first examine the feature relationships:
plotFeatureOverlaps(assay_list)
plot(featureNetwork(assay_list))
out = stabMap(assay_list,
              labels_list = labels_list,
              reference_list = reference_list,
              ncomponentsReference = 20, ncomponentsSubset = 20)

# inspect the output object, a mixture of PCs and (supervised) LDs
head(out)

# look at the scale of each component and discriminant
boxplot(out, las = 2, outline = FALSE)

# re-weight embedding for equal contribution from LDs
out_reweighted = reWeightEmbedding(out)

# look at the new scale of each component and discriminant
boxplot(out_reweighted, las = 2, outline = FALSE)

# testing imputeEmbedding
imp = imputeEmbedding(assay_list,
                      embedding = out_reweighted,
                      reference = colnames(assay_list[[1]]),
                      query = colnames(assay_list[[2]]))

# inspect the imputed values
lapply(imp, dim)
head(imp[[1]])
