#' Write morpho data to file
#'
#' @description
#' Export components of a morpho object to file. Can write trees, character
#' matrices, or fossil ages in various formats.
#'
#' @param data A morpho object
#' @param file File name
#' @param type type to write: "tree", "matrix", or "ages"
#' @param reconstructed If TRUE, write the reconstructed version. Default FALSE.
#' @param uncertainty Numeric. Age uncertainty for fossil ages. Default 0.
#' @param all Logical. If TRUE, write data for all sampled specimens. Applies to
#'   `type = "matrix"` and `type = "ages"` Default FALSE.
#'
#' @return No return value, called for its side effect of writing data to a file.
#'
#' @importFrom ape write.tree write.nexus.data
#'
#' @export
#'
#' @examples
#' data(morpho_data)
#' tmp <- tempfile(fileext = ".tre")
#' write.morpho(morpho_data, file = tmp, type = "tree")
#'
write.morpho <- function(data, file, type = "tree", all = FALSE,
                         reconstructed = FALSE, uncertainty = 0) {

  if (!is.morpho(data)) stop("Error: data must be a morpho object")
  if (is.null(file)) stop("Error: No file name specified")



  if (type == "tree") {
    if (all) {
      write.recon.tree(data, file,  all, reconstructed)
    } else {
      ape::write.tree(data$trees$EvolTree, file)
    }

  } else if (type == "matrix") {
    if (reconstructed) {
      write.recon.matrix(data, file, all)
    } else if (all) {
      ape::write.nexus.data(c(data$sequences$tips, data$sequences$SA),
                            file, format = "standard")
    } else {
      ape::write.nexus.data(data$sequences$tips, file, format = "standard")
    }
  } else if (type == "ages") {
    if (reconstructed) {
      write.recon.tsv(data, file, uncertainty, all)
    } else {
      write.tsv(data, file, uncertainty, all)
    }

  } else {
    stop("Error: 'type' must be 'tree', 'matrix', or 'ages'")
  }
}



#' Write reconstructed tree to file
#'
#' @description
#' Write the reconstructed tree to Newick string
#'
#' @param data Morpho object
#' @param file File name
#'
#' @return
#' No return value, called for its side effect of writing data to a file.
#'
#' @examples
#' data(morpho_data)
#' tmp <- tempfile(fileext = ".tre")
#' write.recon.tree(data = morpho_data, file = tmp)
#'
#'
write.recon.tree <- function(data = NULL, file = NULL, all = TRUE, reconstructed = TRUE) {

  if (is.null(data) || !inherits(data, "morpho")) {
    stop("Error: `data` must be a morpho object.")
  }
  if (is.null(file)) stop("Error: No file name specified")
  if (is.null(data$fossil)) {
    stop("Error: Cannot reconstruct tree as no fossil data in morpho object")
  }

  tt  <- data$trees$TimeTree
  fos <- data$fossil

  # separate extant samples (age 0) from real fossils
  is_extant <- fos$hmax < 1e-8

  if (!reconstructed) {
    # complete time tree with all fossils attached as SAs / fossil tips
    fos_only <- if (any(is_extant)) FossilSim::as.fossils(fos[!is_extant, ]) else fos
    w_tree <- FossilSim::SAtree.from.fossils(tt, fos_only,
                                             tip_order = "youngest_first")$tree

  } else if (all) {

    if (any(is_extant)) {
      # extant sampling was simulated: keep only those extant species
      extant_sampled <- tt$tip.label[fos$sp[is_extant]]
      fos_only <- FossilSim::as.fossils(fos[!is_extant, ])

      tr <- FossilSim::SAtree.from.fossils(tt, fos_only,
                                           tip_order = "youngest_first")$tree

      d        <- ape::node.depth.edgelength(tr)
      tip_ages <- max(d) - d[1:ape::Ntip(tr)]

      extant_labels <- tr$tip.label[tip_ages < 1e-6]
      species       <- sub("_[0-9]+$", "", extant_labels)
      sampled_tips  <- extant_labels[species %in% extant_sampled]

      w_tree <- FossilSim::sampled.tree.from.combined(tr,
                                                      sampled_tips = sampled_tips)
    } else {
      # no extant sampling simulated: keep all extant tips
      tr <- FossilSim::SAtree.from.fossils(tt, fos,
                                           tip_order = "youngest_first")$tree
      w_tree <- FossilSim::sampled.tree.from.combined(tr, rho = 1)
    }

  } else {
    w_tree <- FossilSim::reconstructed.tree.fossils.objects(fossils = fos,
                                                            tree = tt,
                                                            tip_order = "youngest_first")$tree
  }

  ape::write.tree(w_tree, file = file)
}

#' Write reconstructed character matrix to file
#'
#' @description
#' Write the character matrix for the reconstructed tree to a nexus file
#'
#' @param data Morpho object
#' @param file File name
#' @param all Logical. If TRUE, the sequences of all sampled ancestors in the
#'  reconstructed tree are written, not just the first on each lineage.
#'  Default FALSE.
#'
#' @return
#' No return value, called for its side effect of writing data to a file.
#'
#' @examples
#' data(morpho_data)
#' tmp <- tempfile(fileext = ".nex")
#' write.recon.matrix(data = morpho_data, file = tmp)
#'
#' @export
#'
write.recon.matrix <- function (data, file = NULL, all) {

  if (is.null(data) || !inherits(data, "morpho")) {
    stop("Error: `data` must be a morpho object.")
  }

  if (is.null(file)) stop("Error: No file name specified")
  if(is.null(data$fossil)) stop ("Error: Cannot reconstruct tree as no fossil data in morpho object")


  mat <- reconstruct.matrix(data)

  ## all sampled ancestors not already in the reconstructed matrix
  ## remaining sampled ancestors (_3, _4, ...), renamed to match the tree
  if (all) {
    extra <- morphsim_fossilsim(data, all = TRUE)
    extra <- extra[!(extra[, "Fossilsim"] %in% names(mat)), , drop = FALSE]
    if (nrow(extra) > 0) {
      extra_seqs <- data$sequences$SA[extra[, "Morphsim"]]
      names(extra_seqs) <- extra[, "Fossilsim"]
      mat <- c(mat, extra_seqs)
    }
  }

  ape::write.nexus.data(mat, file = file, format = "standard")

}


#' Write the taxa ages
#'
#' @description
#' Writes the ages of the specimens in the true tree to a file. The tsv format used
#' here is directly compatible with RevBayes
#'
#' @param data Morpho object
#' @param file File name
#' @param uncertainty Numeric. Adds uncertainty to fossil ages in the morpho object.
#'  The ages in the object are point estimates by default; setting `uncertainty`
#'  will create an age range of ± this value (in millions of years).
#' @param all Logical. If TRUE, also write the ages of internal nodes.
#'  Internal node ages are written as point values (no uncertainty). Default FALSE.
#'
#' @return
#' No return value, called for its side effect of writing data to a file.
#'
#' @examples
#' data(morpho_data)
#' tmp <- tempfile(fileext = ".tsv")
#' write.tsv(data = morpho_data, file = tmp)
#'
#' @export
write.tsv <- function (data, file, uncertainty = 0, all) {

  if (is.null(data) || !inherits(data, "morpho")) {
    stop("Error: `data` must be a morpho object.")
  }

  if (is.null(file)) stop("Error: No file name specified")

  ## ages of full tree
  tip_depths <- ape::node.depth.edgelength(data$trees$TimeTree)[1:length(data$trees$TimeTree$tip.label)]
  tree_height <- max(ape::node.depth.edgelength(data$trees$TimeTree))
  tip_ages <-   round(abs(tree_height - tip_depths),3)
  # extant_tips <- data$trees$TimeTree$tip.label[abs(tip_depths - tree_height) < 1e-8]

  cat("taxon", "min_age", "max_age", sep = "\t", "\n", file = file)
  for ( i in 1:length(tip_ages)){
    if (tip_ages[i] == 0){
      cat(data$trees$TimeTree$tip.label[i], tip_ages[i] ,
          tip_ages[i], sep = "\t", file = file, append = T )
      cat("\n", file = file, append = TRUE)
    } else {
      cat(data$trees$TimeTree$tip.label[i], (tip_ages[i] -  uncertainty) ,
          (tip_ages[i] + uncertainty), sep = "\t", file = file, append = T )
      cat("\n", file = file, append = TRUE)
    }
  }

  ## sampled ancestors

  if (all) {
    SA_labels <- names(data$sequences$SA)

    for (i in 1:length(SA_labels)){

      parts <- as.numeric(strsplit(SA_labels[i], "_")[[1]])
      specimen_num <- parts[1]
      branch_num   <- parts[2]

      # Subset the data frame to get hmin
      hmin <- data$fossil$hmin[data$fossil$ape.branch == branch_num &
                                 data$fossil$specimen  == specimen_num]

      nm <- SA_labels[i]
      if (hmin - uncertainty < 0){
        min_age <- 0
      } else {
        min_age <- hmin - uncertainty
      }
      cat(nm, min_age, (hmin + uncertainty),
          sep = "\t", file = file, append = T )
      cat("\n", file = file, append = TRUE)
    }
  }
}

#' Write the taxa ages of reconstructed tree
#'
#' @description
#' Writes the ages of the specimen in the reconstructed tree to a file. The tsv format used
#' here is directly compatible with RevBayes
#'
#' @param data Morpho object
#' @param file File name
#' @param uncertainty Numeric. Adds uncertainty to fossil ages in the morpho object.
#'  The ages in the object are point estimates by default; setting `uncertainty`
#'  will create an age range of ± this value (in millions of years).
#' @param all Logical. If TRUE, the ages of all sampled ancestors in the
#'  reconstructed tree are written, not just the first on each lineage.
#'  Default FALSE.
#'
#' @return
#' No return value, called for its side effect of writing data to a file.
#'
#' @examples
#' data(morpho_data)
#' tmp <- tempfile(fileext = ".tsv")
#' write.recon.tsv(data = morpho_data, file = tmp)
#'
#' @export

write.recon.tsv <- function (data, file, uncertainty = 0, all){

  if (is.null(data) || !inherits(data, "morpho")) {
    stop("Error: `data` must be a morpho object.")
  }

  if (is.null(file)) stop("Error: No file name specified")

  if(is.null(data$fossil)) stop ("Error: Cannot reconstruct tree as no fossil data in morpho object")


  r_tree <- FossilSim::reconstructed.tree.fossils.objects(fossils  = data$fossil,
                                                          tree = data$trees$TimeTree,
                                                          tip_order = "youngest_first")
  transformations <- morphsim_fossilsim(data, all = all)

  cat("taxon", "min_age", "max_age", sep = "\t", "\n", file = file)


  ## true tree tips
  tps <- unname(r_tree$tree$tip.label)
  matches <- grepl("_1$", tps )
  # Extract elements that match
  reconTreeTips <- gsub("_1$", "", tps[matches])

  ## add these tip labels to the file + plus all sampled ancestor
  seq_tips <- which(names(data$sequences$tips) %in% reconTreeTips)


  for ( i in 1:length(reconTreeTips)){
    ord <-  which(data$trees$TimeTree$tip.label ==reconTreeTips[i])

    if (length(ord) == 0) {
      ed <- which(data$trees$EvolTree$edge[, 2] == as.numeric(sub("t", "", reconTreeTips[i])))
      tip_ages <- round(min(data$fossil$hmin[data$fossil$ape.branch == ed]), 3)
    } else {
      node_pos <- ape::node.depth.edgelength(data$trees$TimeTree)[ord]
      tree_height <- max(ape::node.depth.edgelength(data$trees$TimeTree))
      tip_ages <- round(abs(node_pos - tree_height), 3)
    }

    if(tip_ages == 0){
      nm <- paste0(reconTreeTips[i], "_1")
      cat(nm, tip_ages, tip_ages,
          sep = "\t", file = file, append = T )
      cat("\n", file = file, append = TRUE)
    } else {
      nm <- paste0(reconTreeTips[i], "_1")
      if (tip_ages - uncertainty < 0){
        min_age <- 0
      } else {
        min_age <- tip_ages - uncertainty
      }
      cat(nm,min_age, (tip_ages + uncertainty),
          sep = "\t", file = file, append = T )
      cat("\n", file = file, append = TRUE)

     }
    }

  ## sampled ancestors

  ## sampled ancestors (_2 only, or _2, _3, ... if all = TRUE)

  for (i in seq_len(nrow(transformations))){

    parts <- as.numeric(strsplit(transformations[i, "Morphsim"], "_")[[1]])
    specimen_num <- parts[1]
    branch_num   <- parts[2]

    # Subset the data frame to get hmin
    hmin <- data$fossil$hmin[data$fossil$ape.branch == branch_num &
                               data$fossil$specimen  == specimen_num]

    nm <- transformations[i, "Fossilsim"]
    if (hmin - uncertainty < 0){
      min_age <- 0
    } else {
      min_age <- hmin - uncertainty
    }
    cat(nm,min_age, (hmin + uncertainty),
        sep = "\t", file = file, append = T )
    cat("\n", file = file, append = TRUE)
  }
}

#' Match sampled ancestor labels
#'
#' @description
#' Match the sampled ancestor labels from \code{Morphsim} and \code{Fossilsim}
#' @param data Morpho object containing fossils
#' @param all Logical. If FALSE, match only the sampled ancestors that appear as
#'  \code{_2} tips in the reconstructed tree. If TRUE, match every fossil that
#'  is not already a tip of the reconstructed tree. Default FALSE.
#' @return
#' A character matrix mapping sampled ancestor labels between the naming
#' conventions used by \code{Morphsim} (specimen_branch) and \code{Fossilsim}
#' (lineage and sample number, e.g. "t3_2")
#'
morphsim_fossilsim <- function (data = NULL, all = FALSE){

  if(is.null(data$fossil)) stop("Error: Morpho object does not contian fossils")

  r_tree <- FossilSim::reconstructed.tree.fossils.objects(fossils  = data$fossil,
                                                          tree = data$trees$TimeTree,
                                                          tip_order =  "youngest_first")
  tps <- unname(r_tree$tree$tip.label)

  tree        <- data$trees$TimeTree
  ntips       <- length(tree$tip.label)
  depths      <- ape::node.depth.edgelength(tree)
  tree_height <- max(depths)

  ## node each fossil sits above
  fos <- data$fossil
  fos$node  <- data$trees$EvolTree$edge[fos$ape.branch, 2]
  fos$label <- NA_character_

  ## name every fossil: lineage + sample number (youngest first)
  for (nd in unique(fos$node)) {
    rows <- which(fos$node == nd)
    rows <- rows[order(fos$hmin[rows])]

    if (nd <= ntips) {
      ## lineage ending in a tip of the true tree
      lineage <- tree$tip.label[nd]
      extant  <- round(abs(tree_height - depths[nd]), 3) == 0
      ## extant tip is _1, unless it's already in the fossil table at hmin 0
      if (extant && !any(fos$hmin[rows] == 0)) {
        offset <- 1
      } else {
        offset <- 0
      }
    } else {
      ## lineage not ending in a tip: _1 is its youngest fossil
      lineage <- paste0("t", nd)
      offset  <- 0
    }
    fos$label[rows] <- paste0(lineage, "_", seq_along(rows) + offset)
  }

  if (all) {
    ## every fossil that isn't a _1 tip of the reconstructed tree
    keep <- !(fos$label %in% tps[grepl("_1$", tps)])
  } else {
    ## only _2 sampled ancestors that are in the reconstructed tree (old behaviour)
    keep <- fos$label %in% tps & grepl("_2$", fos$label)
  }

  transformation <- cbind(Morphsim  = paste0(fos$specimen, "_", fos$ape.branch)[keep],
                          Fossilsim = fos$label[keep])
  return(transformation)
}

