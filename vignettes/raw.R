# install StabMap package
# where tokenstring is your own personal access token from GitHub
# (group access to MarioniLab required)

# devtools::install_github("MarioniLab/StabMap",
#                          ref="main",
#                          auth_token = "tokenstring")

# vignette example code
############################
library(StabMap)

set.seed(2021)

# simulate some data as an assay list
assay_list = mockMosaicData()

# assign labels to one group of cells
labels_list = list(
  D1 = rep(letters[1:5], length.out = ncol(assay_list[["D1"]]))
)

# whether the data should be treated as reference
reference_list = list(
  D1 = TRUE,
  D2 = FALSE,
  D3 = TRUE
)

# first examine the feature relationships:
mosaicDataUpSet(assay_list)
plot(mosaicDataTopology(assay_list))
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
