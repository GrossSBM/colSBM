# An R6 Class object, a collection of model for population of BiSBM networks

This object contain all the logic and methods to estimate a colBiSBM
model

## Public fields

- `n`:

  A list with two dimensions, each of size M for the rows and cols

- `A`:

  List of incidence Matrix of size `n[[1]][m]xn[[2]][m]`

- `M`:

  Number of networks

- `mask`:

  List of M masks, indicating NAs in the matrices. 1 for NA, 0 else

- `distribution`:

  Emission distribution either : "poisson" or "bernoulli"

- `net_id`:

  A vector containing the "ids" or names of the networks (if none given,
  they are set to their number in A list)

- `model_list`:

  A list of size Q1max \* Q2max containing the best models

- `discarded_model_list`:

  A list of size Q1max \* Q2max \* (nb_models - 1) containing the
  discarded models.

- `global_opts`:

  A list of options for the model space exploration. See details for
  more information on which options are available.

- `fit_opts`:

  A list of options specifically for fitting the models.

- `separated_inits`:

  A nested list : Q1 init containing Q2 init with each entry containing
  list of size M storing the separated inits for the class.

- `exploration_order_list`:

  A list used to store the path taken in the state space.

- `Z_init`:

  A list of initializations for the Z memberships.

- `free_mixture_row`:

  A boolean signaling if there is free mixture on the row blocks

- `free_mixture_col`:

  A boolean signaling if there is free mixture on the columns blocks

- `ICL`:

  A matrix of size Q1\*Q2 storing the best ICL found for each value of
  Q1, Q2.

- `sep_BiSBM`:

  A named list containing attributes for each sep BiSBM.

- `BICL`:

  A matrix of size Q1\*Q2 storing the best BICL found for each value of
  Q1, Q2.

- `vbound`:

  A matrix of size Q1\*Q2 storing the best vbound found for each value
  of Q1, Q2.

- `best_fit`:

  A fitBipartiteSBMPop object changing regularly to store the current
  best fit.

- `adjusted_fit`:

  Defining a field in case the user wants to adjust a prefitted
  bisbmpop.

- `logfactA`:

  A quantity used with the Poisson probability distribution

- `improved`:

  A field use at each step to check if it has improved.

- `moving_window_coordinates`:

  A list of size two containing the coordinates of the bottom left and
  top right points of the square

- `old_moving_window_coordinates`:

  A list containing the previous coordinates of the moving window.

- `joint_modelisation_preferred`:

  A boolean to store the preferred modelisation.

## Methods

### Public methods

- [`bisbmpop$new()`](#method-bisbmpop-initialize)

- [`bisbmpop$split_clustering()`](#method-bisbmpop-split_clustering)

- [`bisbmpop$merge_clustering()`](#method-bisbmpop-merge_clustering)

- [`bisbmpop$greedy_exploration()`](#method-bisbmpop-greedy_exploration)

- [`bisbmpop$optimize_from_zinit()`](#method-bisbmpop-optimize_from_zinit)

- [`bisbmpop$burn_in()`](#method-bisbmpop-burn_in)

- [`bisbmpop$optimize()`](#method-bisbmpop-optimize)

- [`bisbmpop$moving_window()`](#method-bisbmpop-moving_window)

- [`bisbmpop$truncate_discarded_model_list()`](#method-bisbmpop-truncate_discarded_model_list)

- [`bisbmpop$point_is_in_limits()`](#method-bisbmpop-point_is_in_limits)

- [`bisbmpop$store_criteria_and_best_fit()`](#method-bisbmpop-store_criteria_and_best_fit)

- [`bisbmpop$compute_sep_BiSBM_BICL()`](#method-bisbmpop-compute_sep_BiSBM_BICL)

- [`bisbmpop$choose_joint_or_separated()`](#method-bisbmpop-choose_joint_or_separated)

- [`bisbmpop$print_metrics()`](#method-bisbmpop-print_metrics)

- [`bisbmpop$show()`](#method-bisbmpop-show)

- [`bisbmpop$print()`](#method-bisbmpop-print)

- [`bisbmpop$clone()`](#method-bisbmpop-clone)

------------------------------------------------------------------------

### `bisbmpop$new()`

Create a new instance of the bisbmpop object

This class is generally called via the user function `estimate_colBiSBM`

#### Usage

    bisbmpop$new(
      netlist = NULL,
      net_id = NULL,
      distribution = NULL,
      free_mixture_row = FALSE,
      free_mixture_col = FALSE,
      Z_init = NULL,
      global_opts = list(),
      fit_opts = list()
    )

#### Arguments

- `netlist`:

  The list of M networks

- `net_id`:

  A list of name for the networks, defaults to 1 to M if not provided

- `distribution`:

  The emission distribution either "bernoulli" or "poisson"

- `free_mixture_row`:

  A boolean indicating if there is free mixture on the row blocks

- `free_mixture_col`:

  A boolean indicating if there is free mixture on the column blocks

- `Z_init`:

  A bidimensional list providing a clustering of the row and column
  nodes

- `global_opts`:

  A list of global options used by the algorithm. See details of the
  user function for more information.

- `fit_opts`:

  A list of fit options used by the algorithm. See details of the user
  function for more information.

#### Returns

A new 'bisbmpop' object. A method to perform the splitting of the
clusters

------------------------------------------------------------------------

### `bisbmpop$split_clustering()`

#### Usage

    bisbmpop$split_clustering(origin_model, axis = "row")

#### Arguments

- `origin_model`:

  a model (fitBipartite object) to split from.

- `axis`:

  a string to indicate if this is a column split or a row split. "row"
  or "col".

#### Returns

best of the possible models tested A method to perform the merging of
the clusters

------------------------------------------------------------------------

### `bisbmpop$merge_clustering()`

#### Usage

    bisbmpop$merge_clustering(origin_model, axis = "row")

#### Arguments

- `origin_model`:

  a model (fitBipartite object) to merge from.

- `axis`:

  a string to indicate if this is a "row", "col" or "both" merge

#### Returns

best of the possible models tested Method to greedily explore state of
space looking for the mode and storing the models discovered along the
way

------------------------------------------------------------------------

### `bisbmpop$greedy_exploration()`

#### Usage

    bisbmpop$greedy_exploration(
      starting_point,
      max_iter = self$fit_opts$greedy_exploration_max_steps,
      max_step_without_improvement =
        self$fit_opts$greedy_exploration_max_steps_without_improvement
    )

#### Arguments

- `starting_point`:

  A vector of the two coordinates c(Q1,Q2) which are the starting point.

- `max_iter`:

  The maximum number of iteration of greedy exploration to perform,
  defaults to 12.

- `max_step_without_improvement`:

  defaults to 2, the number of steps to try improving before stopping
  the search.

#### Returns

c(Q1_mode, Q2_mode) which indicates the Q1 and Q2 for which the BICL was
maximal. Optimization from a given Z_init

The method takes no parameters

------------------------------------------------------------------------

### `bisbmpop$optimize_from_zinit()`

#### Usage

    bisbmpop$optimize_from_zinit()

#### Returns

nothing but stores the models in the model list Burn-in method to start
exploring the state of space

The method takes no parameters but modify the object

------------------------------------------------------------------------

### `bisbmpop$burn_in()`

#### Usage

    bisbmpop$burn_in()

#### Returns

nothing; but stores the values The optimization method

------------------------------------------------------------------------

### `bisbmpop$optimize()`

This method performs the burn in (ie initialization + greedy
exploration) and the steps of moving window with cluster splitting and
merging around the mode found.

#### Usage

    bisbmpop$optimize()

#### Returns

nothing; but stores the values The moving window application

------------------------------------------------------------------------

### `bisbmpop$moving_window()`

This method is a moving windows over the Q1xQ2 space for the number of
clusters

#### Usage

    bisbmpop$moving_window(center, depth = 1)

#### Arguments

- `center`:

  is the coordinates (Q1, Q2) in the model list of the mode

- `depth`:

  is how far away from the center the method should be applied in a grid
  style going from center - (depth,depth) to center + (depth, depth)

#### Returns

nothing; but updates the object by adding new models To truncate the
discarder model list

------------------------------------------------------------------------

### `bisbmpop$truncate_discarded_model_list()`

This method remove the worst models regarding the BICL criterion.

#### Usage

    bisbmpop$truncate_discarded_model_list()

#### Returns

nothing Check if the point is in the limits

------------------------------------------------------------------------

### `bisbmpop$point_is_in_limits()`

#### Usage

    bisbmpop$point_is_in_limits(point)

#### Arguments

- `point`:

  A vector of size 2, containing the first coordinate Q1, the number of
  row blocks and Q2 the number of column blocks.

#### Returns

A boolean if the point if in the limits of `[0,Q1_max]x[0,Q2_max]` Store
the criteria and best fit

------------------------------------------------------------------------

### `bisbmpop$store_criteria_and_best_fit()`

This method stores the criteria (vbound, ICL, BICL) of the models in
model_list attribute and at the end choose the best model using the BICL
criterion.

#### Usage

    bisbmpop$store_criteria_and_best_fit()

#### Returns

nothing; modifies the object Computation of the separated Bi SBM

------------------------------------------------------------------------

### `bisbmpop$compute_sep_BiSBM_BICL()`

This method performs the computation of BiSBM for each of the network in
the netlist.

#### Usage

    bisbmpop$compute_sep_BiSBM_BICL()

#### Returns

nothing; stores the values Method to choose the collection modelisation
or a separated one

------------------------------------------------------------------------

### `bisbmpop$choose_joint_or_separated()`

Using the BICL criterion, the best model fitted after the procedure and
the separated Bi SBM this method choose the one that maximizes the BICL
criterion.

#### Usage

    bisbmpop$choose_joint_or_separated()

#### Returns

nothing; stores a boolean Print the vbound, the ICL and the BICL

------------------------------------------------------------------------

### `bisbmpop$print_metrics()`

#### Usage

    bisbmpop$print_metrics()

#### Returns

nothing The message printed when one prints the object

------------------------------------------------------------------------

### `bisbmpop$show()`

#### Usage

    bisbmpop$show(
      type = "Fitted Collection of Bipartite SBM\nwith all fitted models"
    )

#### Arguments

- `type`:

  The title above the message. The print method

------------------------------------------------------------------------

### `bisbmpop$print()`

#### Usage

    bisbmpop$print()

#### Returns

nothing; prints to console

------------------------------------------------------------------------

### `bisbmpop$clone()`

The objects of this class are cloneable with this method.

#### Usage

    bisbmpop$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
