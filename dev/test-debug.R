devtools::load_all()
set.seed(1234)
Net <- lapply(
    list(.7, .7, .7, .7),
    function(p) {
        A <- matrix(0L, 15L, 15L)
        A[lower.tri(A)][sample(15 * 14 / 2, size = round(p * 15 * 14 / 2))] <- 1L
        A <- A + t(A)
        A
    }
)
fit <- estimate_colSBM(Net,
    colsbm_model = "deltapi",
    directed = FALSE,
    distribution = "poisson",
    nb_run = 1L,
    global_opts = list(nb_cores = 1L),
    fit_opts = list(
        Q_max = 2L,
        nb_init = 1L,
        nb_models = 1L,
        depth = 1L
    )
)
fit$best_fit$Q
