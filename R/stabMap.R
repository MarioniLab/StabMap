#' Stabilised mosaic single cell data integration using unshared features
#'
#' stabMap performs mosaic data integration by first building a mosaic data
#' topology, and for each reference dataset, traverses the topology to
#' project and predict data onto a common principal component (PC) or linear
#' discriminant (LD) embedding.
#'
#' @param assay_list A list of data matrices with rownames (features) specified.
#' @param labels_list (optional) named list containing cell labels
#' @param reference_list Named list containing logical values whether the data
#' matrix should be considered as a reference dataset, alternatively a
#' character vector containing the names of the reference data matrices.
#' @param reference_features_list List of features to consider as reference data
#' (default is all available features).
#' @param reference_scores_list Named list of reference scores (default NULL). If
#' provided, matrix of cells (rows with rownames given) and dimensions (columns
#' with colnames given) are used as the reference low-dimensional embedding to
#' target, as opposed to performing PCA or LDA on the input reference data.
#' @param ncomponentsReference Number of principal components for embedding
#' reference data, given either as an integer or a named list for each
#' reference dataset.
#' @param ncomponentsSubset Number of principal components for embedding query
#' data prior to projecting to the reference, given either as an integer or a
#' named list for each reference dataset.
#' @param suppressMessages Logical whether to suppress messages (default TRUE).
#' @param projectAll Logical whether to re-project reference data along with
#' query (default FALSE).
#' @param restrictFeatures logical whether to restrict to features used in
#' dimensionality reduction of reference data (default FALSE). Overall it's
#' recommended that this be FALSE for single-hop integrations and TRUE for
#' multi-hop integrations.
#' @param maxFeatures Maximum number of features to consider for predicting
#' principal component scores (default 1000).
#' @param plot Logical whether to plot mosaic data UpSet plot and mosaic data
#' topology networks (default TRUE).
#' @param scale.center Logical whether to re-center data to a mean of 0 (default
#' FALSE).
#' @param scale.scale Logical whether to re-scale data to standard deviation of
#' 1 (default FALSE).
#'
#' @return matrix containing common embedding with rows corresponding to cells,
#' and columns corresponding to PCs or LDs for reference dataset(s).
#'
#' @examples
#' set.seed(2021)
#' assay_list = mockMosaicData()
#' lapply(assay_list, dim)
#'
#' # specify which datasets to use as reference coordinates
#' reference_list = c("D1", "D3")
#'
#' # specify some sample labels to distinguish using linear discriminant
#' # analysis (LDA)
#' labels_list = list(
#' D1 = rep(letters[1:5], length.out = ncol(assay_list[["D1"]]))
#' )
#'
#' # examine the topology of this mosaic data integration
#' mosaicDataUpSet(assay_list)
#' plot(mosaicDataTopology(assay_list))
#'
#' # stabMap
#' out = stabMap(assay_list,
#'               reference_list = reference_list,
#'               labels_list = labels_list,
#'               ncomponentsReference = 20,
#'               ncomponentsSubset = 20)
#'
#' head(out)
#'
#' @export
stabMap = function(assay_list,
                   labels_list = NULL,
                   reference_list = sapply(names(assay_list), function(x) TRUE, simplify = FALSE),
                   reference_features_list = lapply(assay_list, rownames),
                   reference_scores_list = NULL,
                   ncomponentsReference = 50,
                   ncomponentsSubset = 50,
                   suppressMessages = TRUE,
                   projectAll = FALSE,
                   restrictFeatures = FALSE,
                   maxFeatures = 1000,
                   plot = TRUE,
                   scale.center = TRUE,
                   scale.scale = TRUE) {

  require(igraph)
  require(scater)

  # check various things and error if not:

  # the columns of each assay_list (cells) should all have different names
  stopifnot("columns of each assay_list (cells) should all have different names"=
              !any(duplicated(unlist(lapply(assay_list, colnames)))))

  # each of the assays should have rownames and colnames defined
  stopifnot("each assay in assay_list must have colnames specified"=
              all(unlist(lapply(assay_list, function(x) !is.null(colnames(x))))))

  stopifnot("each assay in assay_list must have rownames specified"=
              all(unlist(lapply(assay_list, function(x) !is.null(rownames(x))))))

  # assay_list should have names
  stopifnot("each assay in assay_list must be named"=
              !is.null(names(assay_list)) & !any(names(assay_list) == ""))

  # and the assay_list names should be unique
  stopifnot("each assay in assay_list must have a unique name"=
              !any(duplicated(names(assay_list))))

  # remove features with zero variance
  assay_list <- lapply(assay_list, function(x){
    x[rowVars(x) > 0,]
  })

  # if labels_list given the entries should match the ncol(assay_list)
  # to-do

  ## remove messages from calculatePCA
  if (suppressMessages) {
    sm = function(expr) {
      suppressWarnings(suppressMessages(eval(expr)))
    }
  } else {
    sm = function(expr) {
      eval(expr)
    }
  }

  ## reference_list should have names, if it's not a list, then it should be
  ## a character vector of the names of assay_list which should be included
  if (is.character(reference_list)) {
    reference_list <- sapply(names(assay_list), function(x) x %in% reference_list, simplify = FALSE)
  }

  # if ncomponentsReference given as integer convert to a list
  if (is.numeric(ncomponentsReference)) {
    ncomponentsReference <- as.list(rep(ncomponentsReference, length(reference_list)))
    names(ncomponentsReference) <- names(reference_list)
  }

  # if ncomponentsSubset given as integer convert to a list
  if (is.numeric(ncomponentsSubset)) {
    ncomponentsSubset <- as.list(rep(ncomponentsSubset, length(reference_list)))
    names(ncomponentsSubset) <- names(reference_list)
  }

  if (plot) mosaicDataUpSet(assay_list)

  assay_network = mosaicDataTopology(assay_list)
  if (plot) plot(assay_network)

  ## check whether the network is a connected component
  ## the number of components should be 1:
  if (components(assay_network)$no != 1) {
    stop("feature network is not connected, features must overlap in some way via rownames")
  }

  ## if needed, scale the data
  if (any(c(scale.center, scale.scale))) {
    assay_list <- lapply(
      assay_list,
      function(x) t(scale(t(x), center = scale.center, scale = scale.scale))
    )
  }

  all_embeddings_list = list()

  for (reference_dataset in names(assay_list)) {

    ## if not a reference go next
    if (!reference_list[[reference_dataset]]) next

    message(paste0("treating \"", reference_dataset, "\" as reference"))

    ## when the graph has a weight, then by default it will use them
    # shortest path is weighted by the number of shared features
    to_nodes = names(sort(distances(assay_network, to = reference_dataset)[,reference_dataset]))
    all_paths = lapply(all_shortest_paths(assay_network,
                                          from = reference_dataset)$res, names)
    names(all_paths) <- unlist(lapply(all_paths, function(x) rev(x)[1]))
    all_paths <- all_paths[to_nodes]

    ## the PC space of the reference dataset
    nPC = min(ncomponentsReference[[reference_dataset]], nrow(assay_list[[reference_dataset]]))

    for (projectionType in c("PC", "LD")) {

      if (projectionType %in% "PC") {

        if (!is.null(reference_scores_list[[reference_dataset]])) {

          message("reference scores given, using these for mosaic integration")
          reference_scores = reference_scores_list[[reference_dataset]]
          restrictFeatures = FALSE

        } else {

        reference_scores_raw = sm(calculatePCA(assay_list[[reference_dataset]][reference_features_list[[reference_dataset]],],
                                               ncomponents = nPC,
                                               scale = FALSE))

        attr(reference_scores_raw, "rotation") <- list(attr(reference_scores_raw, "rotation"),
                                                   setNames(rep(0,nrow(attr(reference_scores_raw, "rotation"))),
                                                            rownames(attr(reference_scores_raw, "rotation"))))

        loadings_reference = attr(reference_scores_raw, "rotation")

        reference_scores = as.matrix(t(assay_list[[reference_dataset]])) %*1% loadings_reference
        }

        d_nPC = diag(nPC)
        colnames(d_nPC) <- paste0(reference_dataset, "_", colnames(reference_scores))

        P_0 = d_nPC
      }

      if (projectionType %in% "LD") {
        # linear discriminants

        if (is.null(labels_list[[reference_dataset]])) next

        require(MASS)

        message(paste0("labels provided for \"", reference_dataset, "\", adding LD components"))

        features = Reduce(intersect, lapply(assay_list, rownames))
        if (length(features) == 0) {
          features = reference_features_list[[reference_dataset]]
        }

        if (length(features) > maxFeatures) {
          message("more input features than maxFeatures, subsetting features using variance ranking")
          require(scran)
          genevars = modelGeneVar(assay_list[[reference_dataset]][features,])
          genevars_sorted = genevars[order(genevars$bio, decreasing = TRUE),]
          features <- rownames(genevars_sorted)[seq_len(maxFeatures)]
        }

        labels_train = labels_list[[reference_dataset]]

        require(Matrix)
        ## remove features with zero variance for LDA
        vars = rowMaxs(apply(fac2sparse(labels_train), 1, function(x)
          rowWeightedVars(assay_list[[reference_dataset]][features,],x)), na.rm = TRUE)
        if (any(vars == 0)) message("removing features with zero intra-class variance")
        features <- features[vars > 0]

        data_train = t(assay_list[[reference_dataset]][features,])

        lda.fit = sm(lda(data_train[!is.na(labels_train),],
                         grouping = labels_train[!is.na(labels_train)]))
        colnames(lda.fit$scaling) <- paste0(reference_dataset, "_", colnames(lda.fit$scaling))

        reference_scores = t(assay_list[[reference_dataset]][features,]) %projpred% lda.fit

        d_nLD = diag(ncol(reference_scores))
        colnames(d_nLD) <- colnames(reference_scores)

        P_0 = d_nLD

      }


      embedding_list = list()

      for (path in all_paths) {

        message(paste0("generating embedding for path with reference \"",
                       reference_dataset,
                       "\": ",
                       paste0(rev(paste0("\"", path, "\"")), collapse = " -> ")))

        if (identical(as.character(path), reference_dataset)) {
          embedding_list[[reference_dataset]] <- reference_scores %projpred% P_0
          next
        }

        path_current = path

        P = P_0

        obj = list(P)
        ops = list("%projpred%")

        while(length(path_current) > 1) {

          features_current = Reduce(intersect, lapply(assay_list[path_current[1:2]], rownames))
          if (path_current[1] == reference_dataset) {
            current_scores = as.matrix(reference_scores)
            # edit by shila to replace by intersecting among the loadings features when nearest the reference
            if (projectionType == "PC" & restrictFeatures) {
              features_current = intersect(rownames(loadings_reference[[1]]), rownames(assay_list[[path_current[2]]]))
              if (length(features_current) == 0) {
                message("No common features when using restrictFeatures, switching to intersection")
                features_current = Reduce(intersect, lapply(assay_list[path_current[1:2]], rownames))
              }
            }
          } else {
            current_obj = obj[[length(obj)]]
            if (is.list(current_obj)) current_obj <- current_obj[[1]]
            scores_features = setdiff(rownames(current_obj), "intercept")
            current_scores = as.matrix(t(assay_list[[path_current[1]]][scores_features,]))
          }

          if (length(path_current) > 2) {
            nPC_sub = min(ncomponentsSubset[[reference_dataset]], length(features_current))

            dimred_current = sm(calculatePCA(assay_list[[path_current[1]]][features_current,],
                                             ncomponents = nPC_sub,
                                             scale = FALSE))
            attr(dimred_current, "rotation") <- list(attr(dimred_current, "rotation"),
                                                     rowMeans(assay_list[[path_current[1]]][features_current,]))
            loadings_current = attr(dimred_current, "rotation")

            coef = lm.fit(cbind(intercept = 1, dimred_current), current_scores)$coefficients
            coef <- na.omit(coef)

            obj[[length(obj) + 1]] <- coef

            obj[[length(obj) + 1]] <- loadings_current
            ops <- c("%*1%", ops)
            ops <- c("%**%", ops)
          } else {

            ## if there are more than maxFeatures in features_current,
            ## then restrict to HVGs
            if (length(features_current) > maxFeatures) {
              message("more input features than maxFeatures, subsetting features using variance ranking")
              require(scran)
              genevars = modelGeneVar(assay_list[[path_current[1]]][features_current,])
              genevars_sorted = genevars[order(genevars$bio, decreasing = TRUE),]
              features_current <- rownames(genevars_sorted)[seq_len(maxFeatures)]
            }

            coef = lm.fit(
              cbind(
                intercept = 1,
                as.matrix(t(assay_list[[path_current[1]]][features_current,]))
              ),
              current_scores
            )$coefficients
            coef <- na.omit(coef)
            obj[[length(obj) + 1]] <- coef
            ops <- c("%*1%", ops)

          }
          ## if length(path_current) == 2 then this is the last step,
          ## i.e. can use all genes in the remaining assay for coefficient estimation
          path_previous <- path_current[1]
          path_current <- path_current[-1]
        }

        ## now that the path is just length 1, perform the projection
        obj[[length(obj) + 1]] <- t(assay_list[[path_current]])
        if (length(obj) - length(ops) != 1) {
          ops <- c("%**%",ops)
        }
        embedding_list[[path_current]] <- .runOps(rev(obj), ops, leftToRight = FALSE)

        ## also project the prior data too
        if (projectAll) {
          obj[[length(obj)]] <- t(assay_list[[path_previous]])

          embedding_list[[path_previous]] <- .runOps(rev(obj), ops, leftToRight = FALSE)

          for (path_previous_previous in path) {

            if (path_previous_previous == path_previous) break

            features_previous_previous = rownames(assay_list[[path_previous_previous]])

            previous_previous_ind = max(which(unlist(lapply(obj, function(x){
              if (is.list(x)) {
                x <- x[[1]]
              }
              any(features_previous_previous %in% rownames(x))
            }))))

            obj_previous_previous = obj[1:previous_previous_ind]
            ops_previous_previous = rev(rev(ops)[1:(length(obj_previous_previous) - 1)])

            obj_previous_previous[[length(obj_previous_previous) + 1]] <- t(assay_list[[path_previous_previous]])
            ops_previous_previous <- c("%**%",ops_previous_previous)

            embedding_list[[path_previous_previous]] <- .runOps(rev(obj_previous_previous), ops_previous_previous, leftToRight = FALSE)
          }

        }
      }

      embedding = as.matrix(do.call(rbind, embedding_list))

      ## re-centre the embedding
      embedding <- t(t(embedding) - colMeans(embedding))

      all_embeddings_list[[paste0(reference_dataset, "_", projectionType)]] <- embedding
    }

  }

  all_cells = rownames(all_embeddings_list[[1]])

  all_embeddings = do.call(cbind, lapply(all_embeddings_list, "[", all_cells, ))

  return(all_embeddings)
}
