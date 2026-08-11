library(microbenchmark)
M <- 4
true_pi <- c(0.1, 0.3, 0.6)
nbNodes <- 100
K <- length(true_pi)
true_alpha <- matrix(rev(seq(1, 10 * K^2, by = 10)) / 100, nrow = K)
sbm_sampler <- sbm::sampleSimpleSBM(nbNodes = nbNodes, blockProp = true_pi, connectParam = list(mean = true_alpha), directed = TRUE)
set.seed(1234)
sbm_realisations <- lapply(seq(M), sbm_sampler[["rNetwork"]])
A <- lapply(seq(M), function(m) sbm_realisations[[m]][["networkData"]])

# A random start of Z before optimization
fake_Z <- lapply(seq(M), function(m) {
    t(sapply(seq(nbNodes), function(idx) as.integer(sample.int(n = K, size = 1, replace = TRUE, prob = true_pi) == seq(K))))
})
bench <- microbenchmark(
    "R" = { # colSBM R logic
        R_fit <- fitSimpleSBMPop$new(A = A, Q = 3L, Z = fake_Z, directed = TRUE, free_mixture = FALSE, free_density = FALSE, init_method = "spectral", fit_opts = default_fit_opts_unipartite())
        R_fit$optimize()
    },
    "C++" = {
        # C++ logic
        ptr <- colsbm_create(A = A, Q = 3L, tau = fake_Z, directed = TRUE, distribution = "bernoulli", free_mixture = FALSE, free_density = FALSE)
        colsbm_optimize(ptr, max_step = 100L, tol = 1e-9)
    }
)
