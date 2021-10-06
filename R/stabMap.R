#' stabMap
#'
#' stabMap
#'
#' @param assay_list assay_list
#' @param labels_list labels_list
#' @param reference_list reference_list
#' @param reference_features_list reference_features_list
#' @param ncomponentsReference ncomponentsReference
#' @param ncomponentsSubset ncomponentsSubset
#' @param suppressMessages suppressMessages
#' @param projectAll projectAll
#' @param maxFeatures maxFeatures
#' @param plot plot
#' @param scale.center scale.center
#' @param scale.scale scale.scale
#'
#' @return matrix
#'
#' @examples
#'
#' @export
stabMap = function(assay_list,
                              labels_list = NULL,
                              reference_list = sapply(names(assay_list), function(x) TRUE, simplify = FALSE),
                              reference_features_list = lapply(assay_list, rownames),
                              ncomponentsReference = 50,
                              ncomponentsSubset = 50,
                              suppressMessages = TRUE,
                              projectAll = FALSE,
                              maxFeatures = 1000,
                              plot = TRUE,
                              scale.center = FALSE,
                              scale.scale = FALSE) {

  # require packages
  require(igraph)
  require(scater)

  # defensive programming, check various things
  # the columns of each assay_list (cells) should all have different names
  # each of the assays should have rownames and colnames

  # assay_list should have names

  # if labels_list given the entries should match the ncol(assay_list)

  # remove messages from calculatePCA
  if (suppressMessages) {
    sm = function(expr) {
      suppressWarnings(suppressMessages(eval(expr)))
    }
  } else {
    sm = function(expr) {
      eval(expr)
    }
  }

  # reference_list should have names
  # if it's not a list, then should be a character vector of the names
  # of assay_list
  if (is.character(reference_list)) {
    reference_list <- sapply(names(assay_list), function(x) x %in% reference_list, simplify = FALSE)
  }

  if (plot) plotFeatureOverlaps(assay_list)

  assay_network = featureNetwork(assay_list)
  if (plot) plot(assay_network)

  # check whether the network is a connected component
  # the number of components should be 1:
  if (components(assay_network)$no != 1) {
    stop("feature network is not connected, features must overlap in some way")
  }

  # if needed, scale the data
  if (any(c(scale.center, scale.scale))) {
    assay_list <- lapply(
      assay_list,
      function(x) t(scale(t(x), center = scale.center, scale = scale.scale))
    )
  }


  all_embeddings_list = list()

  for (reference_dataset in names(assay_list)) {
    # reference_dataset = names(assay_list)[1]

    # if not a reference just go to the next one
    if (!reference_list[[reference_dataset]]) next

    message(paste0("treating \"", reference_dataset, "\" as reference"))

    # when the graph has a weight, then by default it will
    # use them
    to_nodes = names(sort(distances(assay_network, to = reference_dataset)[,reference_dataset]))
    all_paths = lapply(all_shortest_paths(assay_network,
                                          from = reference_dataset)$res, names)
    names(all_paths) <- unlist(lapply(all_paths, function(x) rev(x)[1]))
    all_paths <- all_paths[to_nodes]

    # the PC space of the reference dataset
    nPC = min(ncomponentsReference, nrow(assay_list[[reference_dataset]]))

    for (projectionType in c("PC", "LD")) {
      # projectionType = "LD"

      if (projectionType %in% "PC") {



        if (TRUE) {
          reference_scores = sm(calculatePCA(assay_list[[reference_dataset]][reference_features_list[[reference_dataset]],], ncomponents = nPC,
                                             scale = FALSE))
          attr(reference_scores, "loadings") <- list(attr(reference_scores, "loadings"),
                                                     rowMeans(assay_list[[reference_dataset]]))
        }

        # if labels are given, identify the rotation of PCs that represent LDs
        d_nPC = diag(nPC)
        colnames(d_nPC) <- paste0(reference_dataset, "_", colnames(reference_scores))

        P_0 = d_nPC

      }


      if (projectionType %in% "LD") {

        if (is.null(labels_list[[reference_dataset]])) next

        require(MASS)

        message(paste0("labels provided for \"", reference_dataset, "\", adding LD components"))

        # features = Reduce(intersect, lapply(assay_list[path_current[1:2]], rownames))
        features = Reduce(intersect, lapply(assay_list, rownames))
        # features = rownames(assay_list[[reference_dataset]])
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
        # experimental: remove features with zero variance for LDA
        vars = rowMaxs(apply(fac2sparse(labels_train), 1, function(x)
          rowWeightedVars(assay_list[[reference_dataset]][features,],x)), na.rm = TRUE)
        if (any(vars == 0)) message("removing genes with zero intra-class variance")
        features <- features[vars > 0]


        data_train = t(assay_list[[reference_dataset]][features,])

        lda.fit = sm(lda(data_train[!is.na(labels_train),],
                         grouping = labels_train[!is.na(labels_train)]))
        colnames(lda.fit$scaling) <- paste0(reference_dataset, "_", colnames(lda.fit$scaling))
        # lda.scaling = lda.fit$scaling
        # colnames(lda.scaling) <- paste0(reference_dataset, "_", colnames(lda.scaling))

        # reference_scores = t(assay_list[[reference_dataset]][features,])
        reference_scores = t(assay_list[[reference_dataset]][features,]) %projpred% lda.fit

        d_nLD = diag(ncol(reference_scores))
        colnames(d_nLD) <- colnames(reference_scores)

        P_0 = d_nLD



        # if (!is.null(labels_list[[reference_dataset]])) {
        if (FALSE) {

          require(MASS)

          message(paste0("labels provided for \"", reference_dataset, "\", adding LD components"))

          labels_train = labels_list[[reference_dataset]]
          data_train = reference_scores


          if (TRUE) {
            lda.fit = lda(data_train[!is.na(labels_train),],
                          grouping = labels_train[!is.na(labels_train)])
            # lda.scaling = lda.fit$scaling
            # colnames(lda.scaling) <- paste0(reference_dataset, "_", colnames(lda.scaling))

            # experimental:
            # colnames(lda.fit$scaling) <- paste0(reference_dataset, "_", colnames(lda.fit$scaling))
          } else {
            # experimental:
            require(e1071)
            lda.fit = svm(x = data_train[!is.na(labels_train),],
                          y = factor(labels_train[!is.na(labels_train)]),
                          kernel = "linear")
            message("fitted labels model")
          }

          if (FALSE) {
            # P_0 = cbind(d_nPC, lda.scaling)
          } # experimental:
          P_0 = list(d_nPC, lda.fit)
        } else {
          # because the reference scores has nPC components
          # P_0 = d_nPC
        }


      }


      embedding_list = list()

      for (path in all_paths) {
        # path = all_paths[[2]]

        message(paste0("generating embedding for path with reference \"",
                       reference_dataset,
                       "\": ",
                       paste0(rev(paste0("\"", path, "\"")), collapse = " -> ")))

        if (identical(path, reference_dataset)) {
          if (FALSE) {
            embedding_list[[reference_dataset]] <- reference_scores %*% P_0
          }
          # experimental:
          embedding_list[[reference_dataset]] <- reference_scores %projpred% P_0
          next
        }

        path_current = path

        P = P_0

        # experimental:
        obj = list(P)
        ops = list("%projpred%")

        while(length(path_current) > 1) {

          features_current = Reduce(intersect, lapply(assay_list[path_current[1:2]], rownames))
          if (path_current[1] == reference_dataset) {
            current_scores = as.matrix(reference_scores)
          } else {
            current_obj = obj[[length(obj)]]
            if (is.list(current_obj)) current_obj <- current_obj[[1]]
            scores_features = setdiff(rownames(current_obj), "intercept")
            # if (FALSE) {
            #   current_scores = as.matrix(t(assay_list[[path_current[1]]][rownames(P),]))
            # }
            current_scores = as.matrix(t(assay_list[[path_current[1]]][scores_features,]))
          }


          if (length(path_current) > 2) {
            nPC_sub = min(ncomponentsSubset, length(features_current))

            # if (TRUE) {
            dimred_current = sm(calculatePCA(assay_list[[path_current[1]]][features_current,],
                                             ncomponents = nPC_sub, scale = FALSE))
            attr(dimred_current, "rotation") <- list(attr(dimred_current, "rotation"),
                                                     rowMeans(assay_list[[path_current[1]]][features_current,]))
            # } else {
            #   # experimental: replace calculatePCA with irlba
            #   dimred_current = sm(calculatePCA_irlba(assay_list[[path_current[1]]][features_current,],
            #                                          nPC = nPC_sub))
            # }
            loadings_current = attr(dimred_current, "rotation")

            # if (length(path_current) > 2) {
            # if (FALSE) {
            coef = lm.fit(cbind(intercept = 1, dimred_current), current_scores)$coefficients
            coef <- na.omit(coef)

            # even more experimental:
            # obj <- c(coef, obj)
            obj[[length(obj) + 1]] <- coef

            obj[[length(obj) + 1]] <- loadings_current
            ops <- c("%*1%", ops)
            ops <- c("%**%", ops)
            # }


          } else {
            # if (FALSE) {
            # coef = lm.fit(cbind(intercept = 1, as.matrix(t(assay_list[[path_current[1]]][features_current,]))), current_scores)$coefficients
            # }


            # experimental:
            # if there are more than 5000 in features_current,
            # then restrict to HVGs
            # (set as a parameter maxFeatures ?)

            # maxFeatures = 5000
            if (length(features_current) > maxFeatures) {
              message("more input features than maxFeatures, subsetting features using variance ranking")
              require(scran)
              genevars = modelGeneVar(assay_list[[path_current[1]]][features_current,])
              genevars_sorted = genevars[order(genevars$bio, decreasing = TRUE),]
              features_current <- rownames(genevars_sorted)[seq_len(maxFeatures)]
            }

            # print("fitting linear models")

            coef = lm.fit(
              cbind(
                intercept = 1,
                as.matrix(t(assay_list[[path_current[1]]][features_current,]))
              ),
              current_scores
            )$coefficients
            coef <- na.omit(coef)

            # print("fitted linear models")

            # even more experimental:
            # obj <- c(coef, obj)
            obj[[length(obj) + 1]] <- coef

            ops <- c("%*1%", ops)


          }
          # if length(path_current) == 2 then this is the last step,
          # i.e. can use all genes in the remaining assay for coefficient estimation


          # experimental: checking for numerical overflow?
          # matList[[length(matList) + 1]] <- coef
          # matList[[length(matList) + 1]] <- cbind(1, loadings_current)

          # # projection matrix
          # if (FALSE) {
          #   P <- cbind(interecept = 1, loadings_current) %*% coef %*% P
          # }
          # # experimental:
          # if (FALSE) {
          #   P <- cbind(intercept = 1, loadings_current) %*% coef %projpred% P
          # }

          ## experimental:
          path_previous <- path_current[1]

          path_current <- path_current[-1]
        }

        # now that the path is just length 1, perform the projection
          obj[[length(obj) + 1]] <- t(assay_list[[path_current]])
          if (length(obj) - length(ops) != 1) {
            ops <- c("%**%",ops)
          }
          embedding_list[[path_current]] <- .runOps(rev(obj), ops, leftToRight = FALSE)


        ## experimental: also project the prior data too

        if (projectAll) {
          obj[[length(obj)]] <- t(assay_list[[path_previous]])

          embedding_list[[path_previous]] <- .runOps(rev(obj), ops, leftToRight = FALSE)

          for (path_previous_previous in path) {

            if (path_previous_previous == path_previous) break

            # path_previous_previous = path[which(path == path_previous) - 1]
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

      # recenter the embedding:
      embedding <- t(t(embedding) - colMeans(embedding))

      all_embeddings_list[[paste0(reference_dataset, "_", projectionType)]] <- embedding
    }

  }

  # experimental: make sure the rownames are matching when cbinding
  all_cells = rownames(all_embeddings_list[[1]])

  # all_embeddings = do.call(cbind, all_embeddings_list)
  all_embeddings = do.call(cbind, lapply(all_embeddings_list, "[", all_cells, ))

  return(all_embeddings)
}
