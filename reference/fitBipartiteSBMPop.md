# An R6 Class object, a fitted population of netowrks sbm once \$optimize() is done

An R6 Class object, a fitted population of netowrks sbm once
\$optimize() is done

An R6 Class object, a fitted population of netowrks sbm once
\$optimize() is done

## Public fields

- `n`:

  A list with two dimensions, each of size M for the rows and cols

- `M`:

  Number of networks

- `A`:

  List of incidence Matrix of size `n[[1]][m]xn[[2]][m]`

- `mask`:

  List of M masks, indicating NAs in the matrices. 1 for NA, 0 else

- `nonNAs`:

  List of M masks, indicating non NAs in the matrices. 1 - mask, so 0
  for NA, 1 for non NA

- `nb_inter`:

  A vector of length M the number of unique non NA entries

- `Q`:

  Number of clusters, vectors of size2

- `tau`:

  List of size M of list of two variational parameters. `n[[1]][m]xQ`
  matrices and `n[[2]][m]xQ` matrices

- `alpha`:

  Matrix of size QxQ, connection parameters

- `pi`:

  List of M vectors of size Q, the mixture parameters

- `pim`:

  List of M vectors of size Q, the mixture parameters in case of
  free_mixture

- `e`:

  Vector of size M, the sum of unique entries

- `emqr`:

  List of M QxQ matrix, the sum of edges between q and r in m, ie the
  edges that are observed.

- `nmqr`:

  list of M QxQ matrix, the number of entries between q and r in m, ie
  all the possible edges.

- `alpham`:

  list of M QxQ matrix, the classic sbm parameters.

- `free_mixture_row`:

  A boolean indicating if there is a free mixture on the rows

- `free_mixture_col`:

  A boolean indicating if there is a free mixture on the columns

- `weight`:

  A vector of size M for weighted likelihood

- `distribution`:

  Emission distribution either : "poisson" or "bernoulli"

- `mloss`:

  Loss on the M step of the VEM

- `vloss`:

  Loss on the VE step of the VEM

- `vbound`:

  The variational bound

- `entropy`:

  The entropy of the variational distribution

- `net_id`:

  A vector containing the "ids" or names of the networks (if none given,
  they are set to their number in A list)

- `df_mixture`:

  The degrees of freedom for mixture parameters pi,used to compute
  penalty

- `df_connect`:

  The degrees of freedom for connection parameters alpha,used to compute
  penalty

- `Cpi`:

  A list of matrices of size Qd x M containing TRUE (1) or FALSE (0) if
  the d-th dimension cluster is represented in the network m

- `Calpha`:

  The corresponding support on the connectivity parameters computed with
  Cpi.

- `logfactA`:

  A quantity used with the Poisson probability distribution

- `init_method`:

  The initialization method used for the first clustering

- `penalty`:

  The penalty computed based on the number of parameters

- `Z`:

  The clusters memberships, a list of size M of two matrices : 1 for
  rows clusters memberships and 2 for columns clusters memberships

- `MAP`:

  Maximum a posteriori

- `MAP_parameters`:

  MAP params

- `ICL`:

  Stores the ICL of the model

- `BICL`:

  Stores the BICL of the model

- `fit_opts`:

  Fit parameters, used to determine the fitting method/

- `step_counter`:

  Counts the number of passes

- `greedy_exploration_starting_point`:

  Stores the coordinates Q1 & Q2 from the greedy exploration to keep
  track of the starting_point

- `effective_clustering_list`:

  A list of size M storing the number of the clusters that contains at
  least one point. Used for safety checks.

- `clustering_is_complete`:

  A boolean used to know if the model real blocks match the expected
  blocks.

- `tested_taus`:

  A vector of taus values for taus given by init_clust

- `tested_taus_vbound`:

  A vector of vbound values for taus given by init_clust

- `has_converged`:

  A boolean, indicating wether the current fit object VEM converged or
  not

## Active bindings

- `nb_nodes`:

  Returns n a list of the number of nodes per network

- `nb_blocks`:

  Returns Q a vector with 2 coordinates, Q1 and Q2 for the row blocks
  and the column blocks

- `support`:

  Returns the Cpi, a list of M boolean matrices indicating which blocks
  are populated

- `prob_memberships`:

  Returns the tau, the probabilities of memberships "a posteriori",
  after seeing the data

- `parameters`:

  Returns the list of parameters of the model, alpha, pi and rho

- `pred_dyads`:

  Predicted dyads from the estimated probabilities and parameters

- `memberships`:

  The block memberships

## Methods

### Public methods

- [`fitBipartiteSBMPop$new()`](#method-fitBipartiteSBMPop-new)

- [`fitBipartiteSBMPop$compute_MAP()`](#method-fitBipartiteSBMPop-compute_MAP)

- [`fitBipartiteSBMPop$vb_tau_alpha()`](#method-fitBipartiteSBMPop-vb_tau_alpha)

- [`fitBipartiteSBMPop$vb_tau_pi()`](#method-fitBipartiteSBMPop-vb_tau_pi)

- [`fitBipartiteSBMPop$entropy_tau()`](#method-fitBipartiteSBMPop-entropy_tau)

- [`fitBipartiteSBMPop$compute_vbound()`](#method-fitBipartiteSBMPop-compute_vbound)

- [`fitBipartiteSBMPop$compute_penalty()`](#method-fitBipartiteSBMPop-compute_penalty)

- [`fitBipartiteSBMPop$compute_icl()`](#method-fitBipartiteSBMPop-compute_icl)

- [`fitBipartiteSBMPop$compute_entropy()`](#method-fitBipartiteSBMPop-compute_entropy)

- [`fitBipartiteSBMPop$compute_BICL()`](#method-fitBipartiteSBMPop-compute_BICL)

- [`fitBipartiteSBMPop$update_MAP_parameters()`](#method-fitBipartiteSBMPop-update_MAP_parameters)

- [`fitBipartiteSBMPop$fixed_point_tau()`](#method-fitBipartiteSBMPop-fixed_point_tau)

- [`fitBipartiteSBMPop$update_pim()`](#method-fitBipartiteSBMPop-update_pim)

- [`fitBipartiteSBMPop$update_pi()`](#method-fitBipartiteSBMPop-update_pi)

- [`fitBipartiteSBMPop$update_alpham()`](#method-fitBipartiteSBMPop-update_alpham)

- [`fitBipartiteSBMPop$update_alpha()`](#method-fitBipartiteSBMPop-update_alpha)

- [`fitBipartiteSBMPop$taus_order()`](#method-fitBipartiteSBMPop-taus_order)

- [`fitBipartiteSBMPop$init_clust()`](#method-fitBipartiteSBMPop-init_clust)

- [`fitBipartiteSBMPop$m_step()`](#method-fitBipartiteSBMPop-m_step)

- [`fitBipartiteSBMPop$ve_step()`](#method-fitBipartiteSBMPop-ve_step)

- [`fitBipartiteSBMPop$update_mqr()`](#method-fitBipartiteSBMPop-update_mqr)

- [`fitBipartiteSBMPop$make_permutation()`](#method-fitBipartiteSBMPop-make_permutation)

- [`fitBipartiteSBMPop$compute_effective_clustering()`](#method-fitBipartiteSBMPop-compute_effective_clustering)

- [`fitBipartiteSBMPop$optimize()`](#method-fitBipartiteSBMPop-optimize)

- [`fitBipartiteSBMPop$reorder_parameters()`](#method-fitBipartiteSBMPop-reorder_parameters)

- [`fitBipartiteSBMPop$show()`](#method-fitBipartiteSBMPop-show)

- [`fitBipartiteSBMPop$print()`](#method-fitBipartiteSBMPop-print)

- [`fitBipartiteSBMPop$clone()`](#method-fitBipartiteSBMPop-clone)

------------------------------------------------------------------------

### Method `new()`

Initializes the fitBipartiteSBMPop object

#### Usage

    fitBipartiteSBMPop$new(
      A = NULL,
      Q = NULL,
      Z = NULL,
      mask = NULL,
      net_id = NULL,
      distribution = NULL,
      free_mixture_row = TRUE,
      free_mixture_col = TRUE,
      Cpi = NULL,
      Calpha = NULL,
      init_method = "spectral",
      weight = NULL,
      greedy_exploration_starting_point = NULL,
      fit_opts = list(algo_ve = "fp", minibatch = TRUE, verbosity = 1)
    )

#### Arguments

- `A`:

  List of incidence Matrix of size `n[[2]][m]xn[[2]][m]`

- `Q`:

  A vector of size 2 with the number of row blocks and column blocks

- `Z`:

  The clusters memberships, a list of size M of two matrices : 1 for
  rows clusters memberships and 2 for columns clusters memberships

- `mask`:

  List of M masks, indicating NAs in the matrices. 1 for NA, 0 else

- `net_id`:

  A vector containing the "ids" or names of the networks (if none given,
  they are set to their number in A list)

- `distribution`:

  Emission distribution either : "poisson" or "bernoulli"

- `free_mixture_row`:

  A boolean indicating if there is a free mixture on the rows

- `free_mixture_col`:

  A boolean indicating if there is a free mixture on the columns

- `Cpi`:

  A list of matrices of size Qd x M containing TRUE (1) or FALSE (0) if
  the d-th dimension cluster is represented in the network m

- `Calpha`:

  The corresponding support on the connectivity parameters computed with
  Cpi.

- `init_method`:

  The initialization method used for the first clustering

- `weight`:

  A vector of size M for weighted likelihood

- `greedy_exploration_starting_point`:

  Stores the coordinates Q1 & Q2 from the greedy exploration to keep
  track of the starting_point

- `fit_opts`:

  Fit parameters, used to determine the fitting method/ Method to
  compute the maximum a posteriori for Z clustering

------------------------------------------------------------------------

### Method `compute_MAP()`

#### Usage

    fitBipartiteSBMPop$compute_MAP()

#### Returns

nothing; stores the values Computes the portion of the vbound with tau
and alpha

------------------------------------------------------------------------

### Method `vb_tau_alpha()`

#### Usage

    fitBipartiteSBMPop$vb_tau_alpha(m, MAP = FALSE)

#### Arguments

- `m`:

  The id of the network for which to compute

- `MAP`:

  Wether to use the MAP parameters or not, a boolean, defaults to FALSE.

#### Returns

The computed quantity. Computes the portion of the vbound with tau pi
and rho

------------------------------------------------------------------------

### Method `vb_tau_pi()`

#### Usage

    fitBipartiteSBMPop$vb_tau_pi(m, MAP = FALSE)

#### Arguments

- `m`:

  The id of the network for which to compute

- `MAP`:

  Wether to use the MAP parameters or not, a boolean, defaults to FALSE.

#### Returns

The computed quantity. Computes the entropy of the model

------------------------------------------------------------------------

### Method `entropy_tau()`

#### Usage

    fitBipartiteSBMPop$entropy_tau(m)

#### Arguments

- `m`:

  The id of the network for which to compute

#### Returns

The computed quantity. Computes the variational bound (vbound)

------------------------------------------------------------------------

### Method `compute_vbound()`

#### Usage

    fitBipartiteSBMPop$compute_vbound()

#### Returns

The variational bound for the model. Computes the penalty for the model

------------------------------------------------------------------------

### Method `compute_penalty()`

#### Usage

    fitBipartiteSBMPop$compute_penalty(
      penalty_factor = self$fit_opts$penalty_factor
    )

#### Arguments

- `penalty_factor`:

  The penalty factor, a numeric, defaults to
  self\$fit_opts\$penalty_factor, which is 0.5 by default.

#### Returns

the computed penalty using the formulae. Computes the ICL criterion

------------------------------------------------------------------------

### Method `compute_icl()`

#### Usage

    fitBipartiteSBMPop$compute_icl(MAP = FALSE)

#### Arguments

- `MAP`:

  Wether to use the MAP parameters or not, a boolean, defaults to FALSE.

#### Returns

The ICL for the model. Computes the entropy of the variational
distribution

------------------------------------------------------------------------

### Method `compute_entropy()`

#### Usage

    fitBipartiteSBMPop$compute_entropy()

#### Returns

The entropy of the variational distribution computed during VEM Computes
the BICL criterion

------------------------------------------------------------------------

### Method `compute_BICL()`

#### Usage

    fitBipartiteSBMPop$compute_BICL(
      MAP = TRUE,
      penalty_factor = self$fit_opts$penalty_factor,
      store = TRUE
    )

#### Arguments

- `MAP`:

  Wether to use the MAP parameters or not, a boolean, defaults to FALSE.

- `penalty_factor`:

  The penalty factor, a numeric, defaults to
  self\$fit_opts\$penalty_factor, which is 0.5 by default.

- `store`:

  A boolean indicating if the BICL should be stored in the object or
  not, defaults to TRUE.

#### Returns

The BICL for the model. Updates the MAP parameters

------------------------------------------------------------------------

### Method `update_MAP_parameters()`

#### Usage

    fitBipartiteSBMPop$update_MAP_parameters()

#### Returns

nothing; but stores the values Method to update tau values

------------------------------------------------------------------------

### Method `fixed_point_tau()`

Not really a fixed point as tau^1 depends only tau^2.

#### Usage

    fitBipartiteSBMPop$fixed_point_tau(m, d, tol = self$fit_opts$tolerance)

#### Arguments

- `m`:

  The number of the network in the netlist

- `d`:

  The dimension to update

- `tol`:

  The tolerance for which to stop iterating. Defaults to
  self\$fit_opts\$tolerance

#### Returns

The new tau values Computes the pi per network, known as the pim

------------------------------------------------------------------------

### Method `update_pim()`

#### Usage

    fitBipartiteSBMPop$update_pim(m, MAP = FALSE)

#### Arguments

- `m`:

  The number of the network in the netlist

- `MAP`:

  A boolean wether to use the MAP parameters or not, defaults to FALSE

#### Returns

nothing; stores the values Computes the pi for the whole model

------------------------------------------------------------------------

### Method `update_pi()`

#### Usage

    fitBipartiteSBMPop$update_pi(MAP = FALSE)

#### Arguments

- `MAP`:

  A boolean wether to use the MAP parameters or not, defaults to FALSE

#### Returns

the pi and stores the values Computes the alpha per network, known as
the alpham

------------------------------------------------------------------------

### Method `update_alpham()`

#### Usage

    fitBipartiteSBMPop$update_alpham(m, MAP = FALSE)

#### Arguments

- `m`:

  The number of the network in the netlist

- `MAP`:

  A boolean wether to use the MAP parameters or not, defaults to FALSE

#### Returns

the alpham and stores the values Computes the alpha for the whole model

------------------------------------------------------------------------

### Method `update_alpha()`

#### Usage

    fitBipartiteSBMPop$update_alpha(MAP = FALSE)

#### Arguments

- `MAP`:

  A boolean wether to use the MAP parameters or not, defaults to FALSE

#### Returns

the alpha and stores the values

------------------------------------------------------------------------

### Method `taus_order()`

The goal of this function is to test different values of tau and select
the best one in the sense of the BICL (or vbound) ?

#### Usage

    fitBipartiteSBMPop$taus_order(taus_list)

#### Arguments

- `taus_list`:

  List of possible taus for which to provide a ranking

#### Returns

A vector with the order of the taus in regard of vbound

------------------------------------------------------------------------

### Method `init_clust()`

Initialize clusters

#### Usage

    fitBipartiteSBMPop$init_clust()

#### Returns

nothing; stores The M step of the VEM

------------------------------------------------------------------------

### Method `m_step()`

#### Usage

    fitBipartiteSBMPop$m_step(MAP = FALSE, tol = self$fit_opts$tolerance, ...)

#### Arguments

- `MAP`:

  A boolean wether to use the MAP parameters or not, defaults to FALSE

- `tol`:

  The tolerance for which to stop iterating defaults to
  self\$fit_opts\$tolerance

- `...`:

  Other parameters

#### Returns

nothing; stores values An optimization version for the VE step of the
VEM but currently a placeholder

------------------------------------------------------------------------

### Method `ve_step()`

#### Usage

    fitBipartiteSBMPop$ve_step(m, max_iter = 2, tol = self$fit_opts$tolerance, ...)

#### Arguments

- `m`:

  The number of the network in the netlist

- `max_iter`:

  The maximum number of iterations, default to 2

- `tol`:

  The tolerance for which to stop iterating defaults to
  self\$fit_opts\$tolerance

- `...`:

  Other parameters

#### Returns

nothing; stores values Updates the mqr quantities

------------------------------------------------------------------------

### Method `update_mqr()`

Namely, it updates the emqr and nmqr.

#### Usage

    fitBipartiteSBMPop$update_mqr(m)

#### Arguments

- `m`:

  The number of the network in the netlist TODO Investigate what its
  supposed to do

------------------------------------------------------------------------

### Method `make_permutation()`

#### Usage

    fitBipartiteSBMPop$make_permutation()

#### Returns

nothing Computes the number of blocks that are effectively populated

------------------------------------------------------------------------

### Method `compute_effective_clustering()`

#### Usage

    fitBipartiteSBMPop$compute_effective_clustering()

#### Returns

nothing; but stores the value Perform the whole initialization and VEM
algorithm

------------------------------------------------------------------------

### Method [`optimize()`](https://rdrr.io/r/stats/optimize.html)

#### Usage

    fitBipartiteSBMPop$optimize(
      max_step = self$fit_opts$max_vem_steps,
      tol = self$fit_opts$tolerance,
      ...
    )

#### Arguments

- `max_step`:

  The maximum number of steps to perform optimization

- `tol`:

  The tolerance for which to stop iterating, default to
  self\$fit_opts\$tolerance

- `...`:

  Other parameters

#### Returns

nothing Reorder the blocks putting the "strongest" ones first in order
to have a coherent ordering of blocks with SBM and LBM for
visualisation.

------------------------------------------------------------------------

### Method `reorder_parameters()`

#### Usage

    fitBipartiteSBMPop$reorder_parameters()

#### Returns

nothing; stores the new ordering The message printed when one prints the
object

------------------------------------------------------------------------

### Method `show()`

#### Usage

    fitBipartiteSBMPop$show(type = "Fitted Collection of Bipartite SBM")

#### Arguments

- `type`:

  The title above the message. The print method

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

#### Usage

    fitBipartiteSBMPop$print()

#### Returns

nothing; print to console

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    fitBipartiteSBMPop$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
