library(testthat)
library(MorphSim)
library(ape)

# ============================================================
# Shared setup
# ============================================================

set.seed(42)
tree <- ape::rtree(6)
time.tree <- TreeSim::sim.bd.taxa(n = 6, numbsim = 1, lambda = 0.1, mu = 0.05)[[1]]

# ============================================================
# sim.morpho
# ============================================================

test_that("sim.morpho returns a morpho object with correct structure", {
  data <- sim.morpho(tree = tree, k = 2, trait.num = 5)
  expect_true(is.morpho(data))
  expect_equal(length(data$sequences$tips[[1]]), 5)
  expect_setequal(names(data$sequences$tips), tree$tip.label)
  expect_equal(length(data$transition_history), 5)
  expect_equal(length(data$root.states), 5)
})

test_that("sim.morpho works with time tree and clock models", {
  # strict clock
  data <- sim.morpho(time.tree = time.tree, k = 2, trait.num = 5, br.rates = 0.1)
  expect_true(is.morpho(data))
  expect_s3_class(data$trees$TimeTree, "phylo")

  # relaxed clock
  br_rates <- runif(nrow(time.tree$edge), min = 0.01, max = 0.2)
  data2 <- sim.morpho(time.tree = time.tree, k = 2, trait.num = 5, br.rates = br_rates)
  expect_true(is.morpho(data2))
})

test_that("sim.morpho works with partitions", {
  data <- sim.morpho(tree = tree, k = c(2, 3), trait.num = 10, partition = c(5, 5))
  expect_true(is.morpho(data))
  expect_equal(length(data$sequences$tips[[1]]), 10)
  expect_equal(length(data$model$Specified), 2)
})

test_that("sim.morpho works with MkV and ACRV models", {
  data <- sim.morpho(tree = tree, k = 3, trait.num = 10,
                     variable = TRUE, ACRV = "gamma", alpha.gamma = 1, ACRV.ncats = 4)
  expect_true(is.morpho(data))
  expect_false(is.null(data$model$RateVar))

  # MkV: no invariant characters
  tip_matrix <- t(as.data.frame(data$sequences$tips))
  for (i in seq_len(ncol(tip_matrix))) {
    expect_gt(length(unique(tip_matrix[, i])), 1)
  }
})

test_that("sim.morpho works with custom Q matrix", {
  ord_Q <- matrix(c(-0.5, 0.5, 0.0,
                    0.3333333, -0.6666667, 0.3333333,
                    0.0, 0.5, -0.5), nrow = 3, byrow = TRUE)
  data <- sim.morpho(tree = tree, k = 3, trait.num = 5, define.Q = ord_Q)
  expect_true(is.morpho(data))
})

test_that("sim.morpho validates inputs", {
  expect_error(sim.morpho(k = 2, trait.num = 5), "Must provide a tree")
  expect_error(sim.morpho(tree = tree, k = 1, trait.num = 5), "at least 2 states")
  expect_error(sim.morpho(tree = tree, k = c(2, 3), trait.num = 10), "1 partition")
})

# ============================================================
# FossilSim integration
# ============================================================

test_that("sim.morpho generates data for tips and sampled ancestors", {
  set.seed(42)
  f <- FossilSim::sim.fossils.poisson(rate = 0.5, tree = time.tree, root.edge = FALSE)
  f2 <- FossilSim::sim.extant.samples(fossils = f, tree = time.tree, rho = 1)

  data <- sim.morpho(time.tree = time.tree, k = 2, trait.num = 5,
                     br.rates = 0.1, fossil = f2)
  expect_true(is.morpho(data))
  expect_false(is.null(data$fossil))
  if (nrow(f) > 0) {
    expect_gt(length(data$sequences$SA), 0)
    expect_equal(length(data$sequences$SA[[1]]), 5)
  }
})

# ============================================================
# sim.missing.data
# ============================================================

test_that("sim.missing.data introduces and controls missing data", {
  data <- sim.morpho(tree = tree, k = 2, trait.num = 10)

  # probability 0 = no missing
  m0 <- sim.missing.data(data = data, method = "random", seq = "tips", probability = 0)
  expect_false(any(t(as.data.frame(m0$sequences$tips)) == "?"))

  # probability 1 = all missing
  m1 <- sim.missing.data(data = data, method = "random", seq = "tips", probability = 1)
  expect_true(all(t(as.data.frame(m1$sequences$tips)) == "?"))

  # intermediate probability
  m5 <- sim.missing.data(data = data, method = "random", seq = "tips", probability = 0.5)
  expect_true(is.morpho(m5))
  expect_true(any(t(as.data.frame(m5$sequences$tips)) == "?"))
})

test_that("sim.missing.data works with structured methods", {
  # trait method
  data <- sim.morpho(tree = tree, k = 2, trait.num = 10)
  missing <- sim.missing.data(data = data, method = "trait",
                              seq = "tips", probability = 1, traits = c(1, 2))
  tip_matrix <- t(as.data.frame(missing$sequences$tips))
  expect_true(all(tip_matrix[, 1] == "?"))
  expect_true(all(tip_matrix[, 2] == "?"))
})

# ============================================================
# stats.morpho
# ============================================================

test_that("stats.morpho returns CI, RI, and convergence info", {
  set.seed(42)
  data <- sim.morpho(tree = tree, k = 2, trait.num = 20)
  stats <- stats.morpho(data)

  expect_true(is.list(stats))
  expect_true(stats$Statistics$CI >= 0 & stats$Statistics$CI <= 1)
  expect_true(stats$Statistics$RI >= 0 & stats$Statistics$RI <= 1)
  expect_true("trait" %in% colnames(stats$Convergent_Traits))
})

# ============================================================
# combine.morpho
# ============================================================

test_that("combine.morpho merges two morpho objects", {
  data1 <- sim.morpho(tree = tree, k = 2, trait.num = 5)
  data2 <- sim.morpho(tree = tree, k = 3, trait.num = 5)
  combined <- combine.morpho(data1, data2)
  expect_true(is.morpho(combined))
  expect_equal(length(combined$sequences$tips[[1]]), 10)
  expect_equal(length(combined$transition_history), 10)
})

# ============================================================
# Write / export functions
# ============================================================

test_that("export functions write valid files", {
  set.seed(42)
  f <- FossilSim::sim.fossils.poisson(rate = 0.5, tree = time.tree, root.edge = FALSE)
  f2 <- FossilSim::sim.extant.samples(fossils = f, tree = time.tree, rho = 1)
  data <- sim.morpho(time.tree = time.tree, k = 2, trait.num = 5,
                     br.rates = 0.1, fossil = f2)

  # reconstructed tree
  tmp_tree <- tempfile(fileext = ".nwk")
  write.recon.tree(data, file = tmp_tree)
  expect_s3_class(ape::read.tree(tmp_tree), "phylo")

  # reconstructed matrix
  tmp_nex <- tempfile(fileext = ".nex")
  write.recon.matrix(data, file = tmp_nex)
  expect_gt(file.info(tmp_nex)$size, 0)

  # age file
  tmp_tsv <- tempfile(fileext = ".tsv")
  write.tsv(data, file = tmp_tsv)
  ages <- read.delim(tmp_tsv, sep = "\t")
  expect_true(all(c("taxon", "min_age", "max_age") %in% names(ages)))

  unlink(c(tmp_tree, tmp_nex, tmp_tsv))
})

# ============================================================
# Reproducibility
# ============================================================

test_that("sim.morpho is deterministic with same seed", {
  set.seed(123)
  data1 <- sim.morpho(tree = tree, k = 2, trait.num = 5)
  set.seed(123)
  data2 <- sim.morpho(tree = tree, k = 2, trait.num = 5)
  expect_identical(data1$sequences$tips, data2$sequences$tips)
})
