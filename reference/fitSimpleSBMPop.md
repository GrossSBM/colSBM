# An R6 Class object, a fitted population of netowrks sbm once \$optimize() is done

An R6 Class object, a fitted population of netowrks sbm once
\$optimize() is done

An R6 Class object, a fitted population of netowrks sbm once
\$optimize() is done

## Public fields

- `n`:

  A list of size M with the number of nodes per network

- `M`:

  Number of networks

- `A`:

  List of incidence matrices of size `n \times n`

- `mask`:

  List of M masks, indicating NAs in the matrices. 1 for NA, 0 else

- `nb_inter`:

  A vector of length M the number of unique non NA entries

- `directed`:

  A boolean indicating if the networks are directed or not

- `Q`:

  An integer indicating the number of blocks

- `tau`:

  List of length M, variational parameters `n[m]xQ[m]` matrices

- `alpha`:

  Matrix of size QxQ, connection parameters

- `delta`:

  Vector of M, density parameters with `delta[1] = 1`

- `pi`:

  List of M vectors of size Q, the mixture parameters

- `e`:

  Vector of size M, the sum of unique entries

- `emqr`:

  List of M QxQ matrix, the sum of edges between q and r in m

- `nmqr`:

  List of M QxQ matrix, the number of entries between q and r in m

- `pim`:

  List of M vectors of size Q, the mixture parameters (pi_tilde)

- `alpham`:

  list of M QxQ matrix, the classic sbm parameters (alpha_tilde)

- `free_mixture`:

  A boolean indicating if the model is with free mixture

- `free_density`:

  A boolean indicating if the model is with free density

- `weight`:

  A vector of size M for weighted likelihood

- `distribution`:

  The emission distribution, either bernoulli or poisson

- `Cpi`:

  A list of matrices of size Q x M containing TRUE (1) or FALSE (0) if
  the cluster is represented in the network m

- `Calpha`:

  The corresponding support on the connectivity parameters computed with
  Cpi.

- `mloss`:

  Loss on the M step of the VEM

- `vloss`:

  Loss on the VE step of the VEM

- `vbound`:

  The variational bound

- `net_id`:

  A vector containing the "ids" or names of the networks (if none given,
  they are set to their number in A list)

- `df_mixture`:

  The degrees of freedom for mixture parameters pi,used to compute
  penalty

- `df_connect`:

  The degrees of freedom for connection parameters alpha,used to compute
  penalty

- `df_density`:

  The degrees of freedom for density parameters delta, used to compute
  penalty

- `logfactA`:

  A quantity used with the Poisson probability distribution

- `init_method`:

  The initialization method used for the first clustering

- `penalty`:

  The penalty computed based on the number of parameters

- `Z`:

  The clusters memberships, a list of size M of two matrices : 1 for
  rows clusters memberships and 2 for columns clusters memberships

- `map`:

  Maximum a posteriori

- `map_parameters`:

  MAP params

- `ICL`:

  Stores the ICL of the model

- `penalty_clustering`:

  Unused attribute

- `BICL`:

  Stores the BICL of the model

- `net_clustering`:

  Unused parameter

- `counter_merge`:

  A counter for the merge (backward) steps

- `counter_split`:

  A counter for the splitting (forward) steps

- `fit_opts`:

  Fit parameters, used to determine the fitting method/

## Active bindings

- `dircoef`:

  The coefficients used change if the network is directed or not

- `nb_nodes`:

  Returns n a list of the number of nodes per network

- `nb_clusters`:

  Returns Q an integer with the number of blocks

- `support`:

  Returns the Cpi, a list of M boolean matrices indicating which blocks
  are populated

- `memberships`:

  Returns the tau, the probabilities of memberships "a posteriori",
  after seeing the data

- `parameters`:

  Returns the list of parameters of the model, alpha, pi and delta

- `pred_dyads`:

  Predicted dyads from the estimated probabilities and parameters

## Methods

### Public methods

- [`fitSimpleSBMPop$new()`](#method-fitSimpleSBMPop-new)

- [`fitSimpleSBMPop$compute_map()`](#method-fitSimpleSBMPop-compute_map)

- [`fitSimpleSBMPop$objective()`](#method-fitSimpleSBMPop-objective)

- [`fitSimpleSBMPop$vb_tau_alpha()`](#method-fitSimpleSBMPop-vb_tau_alpha)

- [`fitSimpleSBMPop$vb_tau_pi()`](#method-fitSimpleSBMPop-vb_tau_pi)

- [`fitSimpleSBMPop$entropy_tau()`](#method-fitSimpleSBMPop-entropy_tau)

- [`fitSimpleSBMPop$fn_vb_alpha_delta()`](#method-fitSimpleSBMPop-fn_vb_alpha_delta)

- [`fitSimpleSBMPop$gr_vb_alpha_delta()`](#method-fitSimpleSBMPop-gr_vb_alpha_delta)

- [`fitSimpleSBMPop$eval_g0_vb_alpha_delta()`](#method-fitSimpleSBMPop-eval_g0_vb_alpha_delta)

- [`fitSimpleSBMPop$eval_jac_g0_vb_alpha_delta()`](#method-fitSimpleSBMPop-eval_jac_g0_vb_alpha_delta)

- [`fitSimpleSBMPop$update_alpha_delta()`](#method-fitSimpleSBMPop-update_alpha_delta)

- [`fitSimpleSBMPop$compute_vbound()`](#method-fitSimpleSBMPop-compute_vbound)

- [`fitSimpleSBMPop$compute_penalty()`](#method-fitSimpleSBMPop-compute_penalty)

- [`fitSimpleSBMPop$compute_icl()`](#method-fitSimpleSBMPop-compute_icl)

- [`fitSimpleSBMPop$compute_BICL()`](#method-fitSimpleSBMPop-compute_BICL)

- [`fitSimpleSBMPop$compute_exact_icl()`](#method-fitSimpleSBMPop-compute_exact_icl)

- [`fitSimpleSBMPop$compute_exact_icl_iid()`](#method-fitSimpleSBMPop-compute_exact_icl_iid)

- [`fitSimpleSBMPop$update_map_parameters()`](#method-fitSimpleSBMPop-update_map_parameters)

- [`fitSimpleSBMPop$fixed_point_tau()`](#method-fitSimpleSBMPop-fixed_point_tau)

- [`fitSimpleSBMPop$fixed_point_alpha_delta()`](#method-fitSimpleSBMPop-fixed_point_alpha_delta)

- [`fitSimpleSBMPop$update_pi()`](#method-fitSimpleSBMPop-update_pi)

- [`fitSimpleSBMPop$update_pim()`](#method-fitSimpleSBMPop-update_pim)

- [`fitSimpleSBMPop$update_alpham()`](#method-fitSimpleSBMPop-update_alpham)

- [`fitSimpleSBMPop$update_alpha()`](#method-fitSimpleSBMPop-update_alpha)

- [`fitSimpleSBMPop$init_clust()`](#method-fitSimpleSBMPop-init_clust)

- [`fitSimpleSBMPop$make_permutation()`](#method-fitSimpleSBMPop-make_permutation)

- [`fitSimpleSBMPop$m_step()`](#method-fitSimpleSBMPop-m_step)

- [`fitSimpleSBMPop$ve_step()`](#method-fitSimpleSBMPop-ve_step)

- [`fitSimpleSBMPop$update_mqr()`](#method-fitSimpleSBMPop-update_mqr)

- [`fitSimpleSBMPop$optimize()`](#method-fitSimpleSBMPop-optimize)

- [`fitSimpleSBMPop$show()`](#method-fitSimpleSBMPop-show)

- [`fitSimpleSBMPop$print()`](#method-fitSimpleSBMPop-print)

- [`fitSimpleSBMPop$clone()`](#method-fitSimpleSBMPop-clone)

------------------------------------------------------------------------

### Method `new()`

Initializes the fitBipartiteSBMPop object

#### Usage

    fitSimpleSBMPop$new(
      A = NULL,
      Q = NULL,
      Z = NULL,
      mask = NULL,
      net_id = NULL,
      distribution = "bernoulli",
      free_mixture = TRUE,
      free_density = TRUE,
      directed = NULL,
      init_method = "spectral",
      weight = NULL,
      Cpi = NULL,
      Calpha = NULL,
      logfactA = NULL,
      fit_opts = list(algo_ve = "fp", approx_pois = FALSE, minibatch = TRUE, verbosity = 1)
    )

#### Arguments

- `A`:

  List of incidence Matrix of size `n[[2]][m]xn[[2]][m]`

- `Q`:

  The number of blocks

- `Z`:

  The block memberships, a list of size M of two matrices : 1 for rows
  clusters memberships and 2 for columns clusters memberships

- `mask`:

  List of M masks, indicating NAs in the matrices. 1 for NA, 0 else

- `net_id`:

  A vector containing the "ids" or names of the networks (if none given,
  they are set to their number in A list)

- `distribution`:

  Emission distribution either : "poisson" or "bernoulli"

- `free_mixture`:

  A boolean indicating if there is a free mixture

- `free_density`:

  A boolean indicating if there is a free_density

- `directed`:

  A boolean specifying if the networks are directed or not

- `init_method`:

  The initialization method used for the first clustering

- `weight`:

  A vector of size M for weighted likelihood

- `Cpi`:

  A list of matrices of size Qd x M containing TRUE (1) or FALSE (0) if
  the d-th dimension cluster is represented in the network m

- `Calpha`:

  The corresponding support on the connectivity parameters computed with
  Cpi.

- `logfactA`:

  A quantity used with the Poisson probability distribution

- `fit_opts`:

  Fit parameters, used to determine the fitting method/ Method to
  compute the maximum a posteriori for Z clustering

- `greedy_exploration_starting_point`:

  Stores the coordinates Q1 & Q2 from the greedy exploration to keep
  track of the starting_point

------------------------------------------------------------------------

### Method `compute_map()`

#### Usage

    fitSimpleSBMPop$compute_map()

#### Returns

nothing; stores the values Objective function

------------------------------------------------------------------------

### Method `objective()`

#### Usage

    fitSimpleSBMPop$objective()

#### Returns

The evaluation of the function Computes the portion of the vbound with
tau and alpha

------------------------------------------------------------------------

### Method `vb_tau_alpha()`

#### Usage

    fitSimpleSBMPop$vb_tau_alpha(m, map = FALSE)

#### Arguments

- `m`:

  The id of the network for which to compute

- `map`:

  Wether to use the MAP parameters or not, a boolean, defaults to FALSE.

#### Returns

The computed quantity.

------------------------------------------------------------------------

### Method `vb_tau_pi()`

#### Usage

    fitSimpleSBMPop$vb_tau_pi(m, map = FALSE)

#### Arguments

- `m`:

  The id of the network for which to compute

- `map`:

  Wether to use the MAP parameters or not, a boolean, defaults to FALSE.

#### Returns

The computed quantity. Computes the entropy of the model

------------------------------------------------------------------------

### Method `entropy_tau()`

#### Usage

    fitSimpleSBMPop$entropy_tau(m)

#### Arguments

- `m`:

  The id of the network for which to compute

#### Returns

The computed quantity. Objective function for the variational bound
regarding the alpha and delta parameters.

------------------------------------------------------------------------

### Method `fn_vb_alpha_delta()`

#### Usage

    fitSimpleSBMPop$fn_vb_alpha_delta(par, emqr, nmqr)

#### Arguments

- `par`:

  The parameters, alpha and delta combined in one big vector.

- `emqr`:

  List of M QxQ matrix, the sum of edges between q and r in m

- `nmqr`:

  List of M QxQ matrix, the number of entries between q and r in m

#### Returns

The evaluation of the function Gradient of the objective function for
the variational bound regarding the alpha and delta parameters.

------------------------------------------------------------------------

### Method `gr_vb_alpha_delta()`

#### Usage

    fitSimpleSBMPop$gr_vb_alpha_delta(par, emqr, nmqr)

#### Arguments

- `par`:

  The parameters, alpha and delta combined in one big vector.

- `emqr`:

  List of M QxQ matrix, the sum of edges between q and r in m

- `nmqr`:

  List of M QxQ matrix, the number of entries between q and r in m

#### Returns

The evaluation of the function Constraint

------------------------------------------------------------------------

### Method `eval_g0_vb_alpha_delta()`

#### Usage

    fitSimpleSBMPop$eval_g0_vb_alpha_delta(par, emqr, nmqr)

#### Arguments

- `par`:

  The parameters, alpha and delta combined in one big vector.

- `emqr`:

  List of M QxQ matrix, the sum of edges between q and r in m

- `nmqr`:

  List of M QxQ matrix, the number of entries between q and r in m

#### Returns

The evaluation of the function Jacobian of the constraint

------------------------------------------------------------------------

### Method `eval_jac_g0_vb_alpha_delta()`

#### Usage

    fitSimpleSBMPop$eval_jac_g0_vb_alpha_delta(par, emqr, nmqr)

#### Arguments

- `par`:

  The parameters, alpha and delta combined in one big vector.

- `emqr`:

  List of M QxQ matrix, the sum of edges between q and r in m

- `nmqr`:

  List of M QxQ matrix, the number of entries between q and r in m

#### Returns

The evaluation of the function Updates the alpha and delta parameters

------------------------------------------------------------------------

### Method `update_alpha_delta()`

#### Usage

    fitSimpleSBMPop$update_alpha_delta(map = FALSE)

#### Arguments

- `map`:

  Wether to use the MAP parameters or not, a boolean, defaults to FALSE.

#### Returns

nothing; but stores the values Computes the variational bound (vbound)

------------------------------------------------------------------------

### Method `compute_vbound()`

#### Usage

    fitSimpleSBMPop$compute_vbound()

#### Returns

The variational bound for the model. Computes the penalty for the model

------------------------------------------------------------------------

### Method `compute_penalty()`

#### Usage

    fitSimpleSBMPop$compute_penalty()

#### Returns

the computed penalty using the formulae. Computes the ICL criterion

------------------------------------------------------------------------

### Method `compute_icl()`

#### Usage

    fitSimpleSBMPop$compute_icl(map = FALSE)

#### Arguments

- `map`:

  Wether to use the MAP parameters or not, a boolean, defaults to FALSE.

#### Returns

The ICL for the model. Computes the BICL criterion

------------------------------------------------------------------------

### Method `compute_BICL()`

#### Usage

    fitSimpleSBMPop$compute_BICL(map = TRUE)

#### Arguments

- `map`:

  Wether to use the MAP parameters or not, a boolean, defaults to FALSE.

#### Returns

The ICL for the model. Computes the exact ICL criterion

------------------------------------------------------------------------

### Method `compute_exact_icl()`

#### Usage

    fitSimpleSBMPop$compute_exact_icl()

#### Returns

The exact ICL for the model. Computes the exact ICL criterion for iid

------------------------------------------------------------------------

### Method `compute_exact_icl_iid()`

#### Usage

    fitSimpleSBMPop$compute_exact_icl_iid()

#### Returns

The exact ICL for the iid model. Updates the MAP parameters

------------------------------------------------------------------------

### Method `update_map_parameters()`

#### Usage

    fitSimpleSBMPop$update_map_parameters()

#### Returns

nothing; but stores the values Method to update tau values with a
fixed-point algorithm

------------------------------------------------------------------------

### Method `fixed_point_tau()`

#### Usage

    fitSimpleSBMPop$fixed_point_tau(m, max_iter = 1, tol = 0.01)

#### Arguments

- `m`:

  The number of the network in the netlist

- `max_iter`:

  The maximum number of iterations to perform, defaults to 1

- `tol`:

  The tolerance for which to stop iterating defaults to 1e-2

#### Returns

The new tau values Fixed point to update alpha and delta

------------------------------------------------------------------------

### Method `fixed_point_alpha_delta()`

#### Usage

    fitSimpleSBMPop$fixed_point_alpha_delta(
      map = FALSE,
      max_iter = 50,
      tol = 1e-06
    )

#### Arguments

- `map`:

  A boolean wether to use MAP parameters or not, defaults to FALSE

- `max_iter`:

  The maximum number of iterations, default to 50

- `tol`:

  The tolerance for which to stop iterating

#### Returns

nothing; stores the values Computes the pi for the whole model

------------------------------------------------------------------------

### Method `update_pi()`

#### Usage

    fitSimpleSBMPop$update_pi(m, map = FALSE)

#### Arguments

- `m`:

  The number of the network in the netlist

- `map`:

  A boolean wether to use the MAP parameters or not, defaults to FALSE

#### Returns

the pi and stores the values Computes the pi per network, known as the
pim

------------------------------------------------------------------------

### Method `update_pim()`

#### Usage

    fitSimpleSBMPop$update_pim(m, map = FALSE)

#### Arguments

- `m`:

  The number of the network in the netlist

- `map`:

  A boolean wether to use the MAP parameters or not, defaults to FALSE

#### Returns

nothing; stores the values Computes the alpha per network, known as the
alpham

------------------------------------------------------------------------

### Method `update_alpham()`

#### Usage

    fitSimpleSBMPop$update_alpham(m, map = FALSE)

#### Arguments

- `m`:

  The number of the network in the netlist

- `map`:

  A boolean wether to use the MAP parameters or not, defaults to FALSE

#### Returns

the alpham and stores the values Computes the alpha for the whole model

------------------------------------------------------------------------

### Method `update_alpha()`

#### Usage

    fitSimpleSBMPop$update_alpha(map = FALSE)

#### Arguments

- `map`:

  A boolean wether to use the MAP parameters or not, defaults to FALSE

#### Returns

the alpha and stores the values

------------------------------------------------------------------------

### Method `init_clust()`

Initialize clusters

#### Usage

    fitSimpleSBMPop$init_clust()

------------------------------------------------------------------------

### Method `make_permutation()`

TODO Remove

#### Usage

    fitSimpleSBMPop$make_permutation()

------------------------------------------------------------------------

### Method `m_step()`

The M step of the VEM

#### Usage

    fitSimpleSBMPop$m_step(map = FALSE, max_iter = 100, tol = 0.001, ...)

#### Arguments

- `map`:

  A boolean wether to use the MAP parameters or not, defaults to FALSE

- `max_iter`:

  The maximum number of iterations, default to 2

- `tol`:

  The tolerance for which to stop iterating defaults to 1e-3

- `...`:

  Other parameters

#### Returns

nothing; stores values An optimization version for the VE step of the
VEM but currently a placeholder

------------------------------------------------------------------------

### Method `ve_step()`

#### Usage

    fitSimpleSBMPop$ve_step(m, max_iter = 20, tol = 0.001, ...)

#### Arguments

- `m`:

  The number of the network in the netlist

- `max_iter`:

  The maximum number of iterations, default to 2

- `tol`:

  The tolerance for which to stop iterating defaults to 1e-3

- `...`:

  Other parameters

#### Returns

nothing; stores values Updates the mqr quantities

------------------------------------------------------------------------

### Method `update_mqr()`

Namely, it updates the emqr and nmqr.

#### Usage

    fitSimpleSBMPop$update_mqr(m)

#### Arguments

- `m`:

  The number of the network in the netlist Perform the whole
  initialization and VEM algorithm

------------------------------------------------------------------------

### Method [`optimize()`](https://rdrr.io/r/stats/optimize.html)

#### Usage

    fitSimpleSBMPop$optimize(max_step = self$fit_opts$max_step, tol = 0.001, ...)

#### Arguments

- `max_step`:

  The maximum number of steps to perform optimization

- `tol`:

  The tolerance for which to stop iterating, default to 1e-3

- `...`:

  Other parameters

#### Returns

nothing; stores values The message printed when one prints the object

------------------------------------------------------------------------

### Method `show()`

#### Usage

    fitSimpleSBMPop$show(type = "Fitted Collection of Simple SBM")

#### Arguments

- `type`:

  The title above the message. The print method

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

#### Usage

    fitSimpleSBMPop$print()

#### Returns

nothing; print to console

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    fitSimpleSBMPop$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
