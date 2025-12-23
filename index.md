# colSBM

colSBM is an R package which implements methods for clustering and
inferring the mesoscale structure for collection of networks. In
particular, it allows to:

- Find a common mesoscale structure in a collection of networks by using
  a Stochastic Block Model (SBM) extension for the joint modeling of
  collection of networks.

- Provide a clustering of the nodes of each networks. Classifies
  networks in groups of networks with similar mesoscale structures.

Mathematical details of the methods as well as simulation studies and
applications can be find in Chabert-Liddell, Barbillon, and Donnet
(2024) .

## Installation

You can install the development version of colSBM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GrossSBM/colSBM")
```

# References

Chabert-Liddell, Saint-Clair, Pierre Barbillon, and Sophie Donnet. 2024.
“Learning Common Structures in a Collection of Networks. An Application
to Food Webs.” *The Annals of Applied Statistics* 18 (2): 1213–35.
<https://doi.org/10.1214/23-AOAS1831>.
