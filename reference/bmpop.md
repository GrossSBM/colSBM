# An R6 Class object, a collection of model for population of sbm netowrks

An R6 Class object, a collection of model for population of sbm netowrks

An R6 Class object, a collection of model for population of sbm netowrks

## Public fields

- `n`:

  A list of size M with the number of nodes per network

- `A`:

  List of incidence Matrix of size `n \times n`

- `M`:

  Number of networks

- `mask`:

  List of M masks, indicating NAs in the matrices. 1 for NA, 0 else

- `directed`:

  A boolean indicating if the networks are directed or not

- `distribution`:

  Emission distribution either : "poisson" or "bernoulli"

- `net_id`:

  A vector containing the "ids" or names of the networks (if none given,
  they are set to their number in A list)

- `model_list`:

  A list of size Q containing the best models

- `global_opts`:

  A list of options for the model space exploration. See details for
  more information on which options are available.

- `fit_opts`:

  A list of options specifically for fitting the models.

- `fit_sbm`:

  Pre-fitted sbm objects

- `Z_init`:

  A list of initializations for the Z memberships.

- `free_density`:

  A boolean indicating if we consider free density or if all networks
  have the same density

- `free_mixture`:

  A boolean signaling if there is free mixture

- `ICL_sbm`:

  A list storing the ICL for each sbm fitted

- `ICL`:

  A list of size Q storing the best ICL found for each Q

- `BICL`:

  A list of size Q storing the best BICL found for each Q

- `vbound`:

  A list of size Q storing the best vbound found for each Q

- `best_fit`:

  A fitSimpleSBMPop object changing regularly to store the current best
  fit.

- `logfactA`:

  A quantity used with the Poisson probability distribution

- `improved`:

  A field use at each step to check if it has improved.

## Methods

### Public methods

- [`bmpop$new()`](#method-bmpop-new)

- [`bmpop$optimize_sbm()`](#method-bmpop-optimize_sbm)

- [`bmpop$optimize_from_sbm()`](#method-bmpop-optimize_from_sbm)

- [`bmpop$optimize_spectral()`](#method-bmpop-optimize_spectral)

- [`bmpop$optimize_init()`](#method-bmpop-optimize_init)

- [`bmpop$optimize_from_zinit()`](#method-bmpop-optimize_from_zinit)

- [`bmpop$burn_in()`](#method-bmpop-burn_in)

- [`bmpop$forward_pass()`](#method-bmpop-forward_pass)

- [`bmpop$backward_pass()`](#method-bmpop-backward_pass)

- [`bmpop$optimize()`](#method-bmpop-optimize)

- [`bmpop$choose_models()`](#method-bmpop-choose_models)

- [`bmpop$show()`](#method-bmpop-show)

- [`bmpop$print()`](#method-bmpop-print)

- [`bmpop$clone()`](#method-bmpop-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new instance of the bisbmpop object

This class is generally called via the user function `estimate_colBiSBM`

#### Usage

    bmpop$new(
      netlist = NULL,
      net_id = NULL,
      directed = NULL,
      distribution = "bernoulli",
      free_density = FALSE,
      free_mixture = FALSE,
      fit_sbm = NULL,
      Z_init = NULL,
      global_opts = list(),
      fit_opts = list()
    )

#### Arguments

- `netlist`:

  The list of M networks

- `net_id`:

  A list of name for the networks, defaults to 1 to M if not provided

- `directed`:

  A boolean indicating if the networks are directed or not

- `distribution`:

  The emission distribution either "bernoulli" or "poisson"

- `free_density`:

  If we account for different density between networks

- `free_mixture`:

  A boolean indicating if there is free mixture

- `fit_sbm`:

  The pre-fitted SBM

- `Z_init`:

  A list providing a clustering of the nodes

- `global_opts`:

  A list of global options used by the algorithm. See details of the
  user function for more information.

- `fit_opts`:

  A list of fit options used by the algorithm. See details of the user
  function for more information.

#### Returns

A new 'sbmpop' object. Fit a list of SBM if fit_sbm == TRUE

------------------------------------------------------------------------

### Method `optimize_sbm()`

#### Usage

    bmpop$optimize_sbm()

#### Returns

nothing; but stores the values Fit the colSBM mode using sbm as
initializations

------------------------------------------------------------------------

### Method `optimize_from_sbm()`

#### Usage

    bmpop$optimize_from_sbm(index, Q, nb_clusters)

#### Arguments

- `index`:

  The sequence of networks number going from 1 to M

- `Q`:

  The number of clusters

- `nb_clusters`:

  A subindex for model list, in practice always 1

#### Returns

bmpop object Fit the colSBM mode using spectral decompositions as
initializations

------------------------------------------------------------------------

### Method `optimize_spectral()`

#### Usage

    bmpop$optimize_spectral(index, Q, nb_clusters)

#### Arguments

- `index`:

  The sequence of networks number going from 1 to M

- `Q`:

  The number of clusters

- `nb_clusters`:

  A subindex for model list, in practice always 1

#### Returns

bmpop object Fit a colSBM model using a given initialization provided
with Z for a specific network

------------------------------------------------------------------------

### Method `optimize_init()`

#### Usage

    bmpop$optimize_init(index, Z, Q, nb_clusters, Cpi = NULL, Calpha = NULL)

#### Arguments

- `index`:

  The sequence of networks number going from 1 to M

- `Z`:

  The provided initialization

- `Q`:

  The number of clusters

- `nb_clusters`:

  A subindex for model list, in practice always 1

- `Cpi`:

  A list of size M containing the support for the possibly absent
  blocks, defaults to NULL

- `Calpha`:

  A list of size M containing the support for the possibly absent
  interaction parameter, defaults to NULL

#### Returns

bmpop object Fit the whole colSBM mode using a given initialization
provided with a given Z_init

------------------------------------------------------------------------

### Method `optimize_from_zinit()`

#### Usage

    bmpop$optimize_from_zinit(index, Q, nb_clusters)

#### Arguments

- `index`:

  The sequence of networks number going from 1 to M

- `Q`:

  The number of clusters

- `nb_clusters`:

  A subindex for model list, in practice always 1

#### Returns

bmpop object Burn-in method that performs the initialization necessary
to next begin the search and model selection

------------------------------------------------------------------------

### Method `burn_in()`

#### Usage

    bmpop$burn_in()

#### Returns

nothing; but stores the values This forward pass split the clusters
found previously and reallocate the nodes to the new clusters

------------------------------------------------------------------------

### Method `forward_pass()`

#### Usage

    bmpop$forward_pass(
      Q_min = self$global_opts$Q_min,
      Q_max = self$global_opts$Q_max,
      index = seq(self$M),
      nb_clusters = 1L
    )

#### Arguments

- `Q_min`:

  The minimal number of clusters, defaults to the global options Q_min
  parameter

- `Q_max`:

  The max number of clusters, defaults to the global options Q_max
  parameter

- `index`:

  The sequence of networks number going from 1 to M

- `nb_clusters`:

  A subindex for model list, in practice always 1

#### Returns

Q - 1 This backward pass merges the clusters found previously

------------------------------------------------------------------------

### Method `backward_pass()`

#### Usage

    bmpop$backward_pass(
      Q_min = self$global_opts$Q_min,
      Q_max = self$global_opts$Q_max,
      index = seq(self$M),
      nb_clusters = 1L
    )

#### Arguments

- `Q_min`:

  The minimal number of clusters, defaults to the global options Q_min
  parameter

- `Q_max`:

  The max number of clusters, defaults to the global options Q_max
  parameter

- `index`:

  The sequence of networks number going from 1 to M

- `nb_clusters`:

  A subindex for model list, in practice always 1

#### Returns

Q + 1 The optimization method

------------------------------------------------------------------------

### Method [`optimize()`](https://rdrr.io/r/stats/optimize.html)

This method performs the burn in and the steps of moving window with
cluster splitting (forward pass) and merging (backward pass) around the
mode found.

#### Usage

    bmpop$optimize()

#### Returns

nothing; but stores the values Performs the model selection based on the
BICL criterion

------------------------------------------------------------------------

### Method `choose_models()`

#### Usage

    bmpop$choose_models(models, Q, index = seq(self$M), nb_clusters = 1L)

#### Arguments

- `models`:

  The list of models in which to choose the best

- `Q`:

  The value of Q for which the model selection is performed

- `index`:

  The sequence of networks number going from 1 to M

- `nb_clusters`:

  A subindex for model list, in practice always 1

#### Returns

The best models The message printed when one prints the object

------------------------------------------------------------------------

### Method `show()`

#### Usage

    bmpop$show(type = "Fitted Collection of Simple SBM")

#### Arguments

- `type`:

  The title above the message. The print method

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

#### Usage

    bmpop$print()

#### Returns

nothing; print to console

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    bmpop$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
