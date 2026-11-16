# Tutorial on food webs

``` r

library(colSBM)
library(patchwork)
data("foodwebs")
```

## Estimation with colSBM

We load a list of 8 foodwebs. They are binary directed networks with
different number of species. First, we are going to model jointly the
first $`3`$ networks, using the iid-colSBM model.

``` r

# global_opts = list(nb_cores = 1L,
#                    nb_models = 5L,
#                    nb_init = 10L,
#                    depth = 2L,
#                    verbosity = 1,
#                    spectral_init = FALSE,
#                    Q_max = 8L,
#                    plot_details = 1)

set.seed(1234)
res_fw_iid <- estimate_colSBM(
  netlist = foodwebs[1:3], # A list of networks
  colsbm_model = "iid", # The name of the model
  directed = TRUE, # Foodwebs are directed networks
  net_id = names(foodwebs)[1:3], # Name of the networks
  nb_run = 1L, # Number of runs of the algorithm
  global_opts = list(
    verbosity = 0,
    plot_details = 0,
    Q_max = 8
  ) # Max number of clusters
)
```

We can look at how the variational bound and the model selection
criteria evolve with the number of clusters. Here, the BICL criterion
selects Q = 5 blocks.

``` r

plot(res_fw_iid)
```

![](tutorial_files/figure-html/unnamed-chunk-3-1.png)

``` r

best_fit <- res_fw_iid$best_fit
```

## Results and analysis

Here are some useful fields to analyze the results.

``` r

best_fit
#> Fitted Collection of Simple SBM -- bernoulli variant for 3 networks 
#> =====================================================================
#> Dimension = ( 105 58 71 ) - ( 5 ) blocks.
#> BICL =  -1970.397  -- #Empty blocks :  0  
#> =====================================================================
#> * Useful fields 
#>   $distribution, $nb_nodes, $nb_clusters, $support, $Z 
#>   $memberships, $parameters, $BICL, $vbound, $pred_dyads
```

We can get:

- the estimation of the model parameters

``` r

best_fit$parameters
#> $alpha
#>              [,1]         [,2]         [,3]        [,4]         [,5]
#> [1,] 2.008852e-06 1.427805e-09 1.414338e-02 0.008266716 1.447384e-01
#> [2,] 9.342556e-03 2.916666e-08 7.196511e-01 0.509407102 8.036370e-01
#> [3,] 7.972743e-12 8.772641e-11 5.037851e-05 0.025876863 1.183569e-03
#> [4,] 2.379962e-11 4.028603e-10 1.457957e-06 0.123033657 2.922992e-06
#> [5,] 1.294199e-04 1.073739e-09 2.807068e-02 0.223372403 1.192370e-02
#> 
#> $pi
#> $pi[[1]]
#> [1] 0.46658271 0.02564103 0.22553420 0.06645679 0.21578527
#> 
#> $pi[[2]]
#> [1] 0.46658271 0.02564103 0.22553420 0.06645679 0.21578527
#> 
#> $pi[[3]]
#> [1] 0.46658271 0.02564103 0.22553420 0.06645679 0.21578527
#> 
#> 
#> $delta
#> [1] 1 1 1
```

- The block memberships:

``` r

best_fit$Z
#> [[1]]
#>      Unidentified sp1 FW_009   Terrestrial plant material 
#>                            2                            2 
#>    Terrestrial invertebrates        Achnanthes lanceolata 
#>                            5                            1 
#>              Batrachospermum         Calothrix sp1 FW_009 
#>                            1                            1 
#>         Cocconeis placentula         Cosmarium sp1 FW_009 
#>                            1                            1 
#>        Cyclotella sp1 FW_009              Cymbella aspera 
#>                            1                            1 
#>           Cymbella cuspidata               Cymbella kappi 
#>                            1                            1 
#>              Cymbella tumida              Diatoma heimale 
#>                            1                            1 
#>              Epithemia sorex            Epithemia turgida 
#>                            1                            1 
#>           Eunotia serpentina      Unidentified sp2 FW_009 
#>                            1                            1 
#>           Euntoia pectinalis        Fragilaria sp1 FW_009 
#>                            1                            1 
#>        Fragilaria vaucheriae         Frustulia rhomboides 
#>                            1                            1 
#>        Gomphoneis herculeana       Gomphonema accuminatum 
#>                            1                            1 
#>        Gomphonema angustatum         Gomphonema truncatum 
#>                            1                            1 
#>      Unidentified sp3 FW_009             Melosira varians 
#>                            1                            1 
#>            Navicula avenacea           Navicula dicephala 
#>                            1                            1 
#>              Nitzschia dubia        Oedogonium sp1 FW_009 
#>                            1                            1 
#>        Phormidium sp1 FW_009         Pinnularia mesolepta 
#>                            1                            1 
#>           Pinnularia viridis    Pleaurotaenium sp1 FW_009 
#>                            1                            1 
#>         Rhoicospenia curvata        Rhopalodia sp1 FW_009 
#>                            1                            1 
#>       Schizothrix sp1 FW_009       Staurastrum sp1 FW_009 
#>                            1                            1 
#>            Surirella elegans             Surirella tenera 
#>                            1                            1 
#>                 Synedra ulna        Tabellaria fenestrata 
#>                            1                            1 
#>        Tabellaria flocculosa          Ulothrix sp1 FW_009 
#>                            1                            1 
#>      Unidentified sp4 FW_009      Unidentified sp5 FW_009 
#>                            1                            1 
#>        Acroneuria sp1 FW_009          Aelosoma sp1 FW_009 
#>                            4                            5 
#>         Alloperla sp1 FW_009         Anepeorus sp1 FW_009 
#>                            3                            4 
#>             Antocha saxicola                       Baetis 
#>                            5                            5 
#>               Boyeria vinosa Bryophaenocladius sp1 FW_009 
#>                            3                            5 
#>        Chauliodes sp1 FW_009            Chimarra atterima 
#>                            4                            5 
#>      Unidentified sp6 FW_009     Conchapelopia sp1 FW_009 
#>                            5                            4 
#>       Cryptolabis sp1 FW_009         Cyrnellus sp1 FW_009 
#>                            3                            3 
#>     Dicrotendipes sp1 FW_009          Diplectrona modesta 
#>                            5                            5 
#>             Ectopria nervosa   Endochironomous sp1 FW_009 
#>                            5                            5 
#>       Ephemerella sp1 FW_009  Eukieferiella pseudomontana 
#>                            3                            3 
#>   Eukiefferiella 'dark' type        Glossosoma sp1 FW_009 
#>                            5                            3 
#>          Gyraulus sp1 FW_009            Haploperla brevis 
#>                            5                            3 
#>          Hexatoma sp1 FW_009                Hydrophilidae 
#>                            4                            3 
#>  Hydropsyche sp1 FW_009arana            Larsia sp1 FW_009 
#>                            4                            3 
#>        Leucrocuta sp1 FW_009           Leuctra sp1 FW_009 
#>                            5                            5 
#>      Unidentified sp7 FW_009      Metriocnemus sp1 FW_009 
#>                            5                            5 
#>         Micrasema sp1 FW_009        Ochthebius sp1 FW_009 
#>                            3                            3 
#>      Ophiogomphus sp1 FW_009  Paraleptophlebia sp1 FW_009 
#>                            4                            3 
#>  Paranyctiophylax sp1 FW_009     Polycentropus sp1 FW_009 
#>                            3                            3 
#>         Probezzia sp1 FW_009        Promoresia sp1 FW_009 
#>                            3                            5 
#>         Psephenus sp1 FW_009 Pseudolimnolphila sp1 FW_009 
#>                            5                            3 
#>       Rhyacophila sp1 FW_009          Simulium sp1 FW_009 
#>                            3                            5 
#>        Sphaerium occidentale     Stempelinella sp1 FW_009 
#>                            3                            5 
#>            Stenelmis crenata         Stenelmis sp1 FW_009 
#>                            5                            5 
#>           Suwalia sp1 FW_009        Tallaperla sp1 FW_009 
#>                            3                            5 
#>           Tanytarsus Genus A     Tricorythodes sp1 FW_009 
#>                            3                            3 
#>        Notropis heterolepsis                  Brook trout 
#>                            4                            4 
#>           Orcnocetes virilis       Rhinichthys cataractae 
#>                            3                            4 
#>            Cambarus bartonii 
#>                            3 
#> 
#> [[2]]
#>     Unidentified sp1 FW_012_01      Terrestrial invertebrates 
#>                              2                              1 
#>                 Plant material          Achnanthes lanceolata 
#>                              2                              1 
#>         Achnanthes minutissima      Audouinella sp1 FW_012_01 
#>                              1                              1 
#>                Batrachospermum               Blue-green algae 
#>                              1                              1 
#>                      Calothrix               Cymbella cistula 
#>                              1                              1 
#>               Cymbella mulleri                Diatoma heimale 
#>                              1                              1 
#>                Epithemia sorex              Epithemia turgida 
#>                              1                              1 
#>             Eunotia pectinalis          Eunotia sp1 FW_012_01 
#>                              1                              1 
#>           Frustulia rhomboides          Gomphoneis herculeana 
#>                              1                              1 
#>          Gomphonema intricatum       Gomphonema sp1 FW_012_01 
#>                              1                              1 
#>             Meridion circulare              Navicula avenacea 
#>                              1                              1 
#>                  Pleurotaenium          Rhoicosphenia curvata 
#>                              1                              1 
#>                  Stigeoclonium                   Synedra ulna 
#>                              1                              1 
#>                       Ulothrix     Unidentified sp2 FW_012_01 
#>                              1                              1 
#>                       Aelosoma    Brachycentrus sp1 FW_012_01 
#>                              5                              5 
#>               Cambarus bartoni       Chauliodes sp1 FW_012_01 
#>                              4                              3 
#>         Cordulegaster maculata        Dicranota sp1 FW_012_01 
#>                              3                              3 
#>             Ectopria thoracica                 Epeorus dispar 
#>                              5                              5 
#>       Glossosoma sp1 FW_012_01      Homoplectra sp1 FW_012_01 
#>                              3                              5 
#>       Hudsonimya sp1 FW_012_01      Hydropsyche sp1 FW_012_01 
#>                              3                              5 
#>       Leucrocuta sp1 FW_012_01          Leuctra sp1 FW_012_01 
#>                              5                              3 
#>       Lumbriculiid oligochaete Parametriocnemus sp1 FW_012_01 
#>                              5                              5 
#>     Neureclipsis sp1 FW_012_01     Ophiogomphus sp1 FW_012_01 
#>                              3                              3 
#>        Palpomyia sp1 FW_012_01        Palpomyia sp2 FW_012_01 
#>                              3                              3 
#>       Promoresia sp1 FW_012_01        Psephenus sp1 FW_012_01 
#>                              5                              5 
#>        Soliperla sp1 FW_012_01                Stenelmis adult 
#>                              4                              3 
#>        Stenelmis sp1 FW_012_01          Suwalia sp1 FW_012_01 
#>                              3                              3 
#>               Tallaperla maria        Thaumalea sp1 FW_012_01 
#>                              5                              5 
#>             Tipula abdominalis                     Salamander 
#>                              5                              4 
#> 
#> [[3]]
#>    Unidentified sp1 FW_012_02            Terrestrial plants 
#>                             2                             2 
#>              Terrestrial bugs Achnanthes inflata var. elata 
#>                             1                             1 
#>         Achnanthes lanceolata           Achnanthes linearis 
#>                             1                             1 
#>        Achnanthes minutissima           Auodinella hermanii 
#>                             1                             1 
#>              Blue Green algae                     Calothrix 
#>                             1                             1 
#>          Cocconeis placentula                Cymbella kappi 
#>                             1                             1 
#>              Cymbella mulleri               Diatoma heimale 
#>                             1                             1 
#>             Epithemia turgida              Eunotia meisteri 
#>                             1                             1 
#>            Eunotia pectinalis         Fragilaria vaucheriae 
#>                             1                             1 
#>          Frustulia rhomboides         Gomphoneis herculeana 
#>                             1                             1 
#>        Gomphonema accuminatum         Gomphonema angustatum 
#>                             1                             1 
#>         Gomphonema intricatum         Gomphonema tennuellum 
#>                             1                             1 
#>                  Marssoniella             Navicula avenacea 
#>                             1                             1 
#>        Navicula cryptocephala               Navicula mutica 
#>                             1                             1 
#>        Navicula sp1 FW_012_02            Pinnularia viridis 
#>                             1                             1 
#>         Rhoicosphenia curvata                    Rhopalodia 
#>                             1                             1 
#>         Surirella brebbisonii             Surirella elegans 
#>                             1                             1 
#>                   Synechoccus                  Synedra ulna 
#>                             1                             1 
#>                      Ulothrix    Unidentified sp2 FW_012_02 
#>                             1                             1 
#>       Aeolosoma sp1 FW_012_02                 Ajax longipes 
#>                             5                             4 
#>               Amphinemura wui           Anchytarsus bicolor 
#>                             5                             3 
#>                        Baetis                    Hudsonimya 
#>                             5                             3 
#>                    Cricotopus                     Dicranota 
#>                             5                             3 
#>  Eukieffidrella pseudomontana           Diplectrona modesta 
#>                             3                             3 
#>                       Dixella                  Dolophilodes 
#>                             3                             4 
#>            Ectopria thoracica                Epeorus dispar 
#>                             5                             5 
#>                  Fatigia pele        Hexatoma sp1 FW_012_02 
#>                             3                             3 
#>                       Leuctra             Oligo Lumbr. Blue 
#>                             5                             5 
#>            Oligo. Lumbr. Pink                  Ophiogomphus 
#>                             5                             3 
#>             Paraleptophelebia                      Pericoma 
#>                             3                             5 
#>                       Pilaria      Pentaneuri sp1 FW_012_02 
#>                             3                             5 
#>       Polycentropus maculatus                     Stenelmis 
#>                             3                             5 
#>              Tallaperla maria                     Tanyderid 
#>                             3                             3 
#>                 Conchapelopia                        Tipula 
#>                             1                             3 
#>              Wormaldia moesta                      Crayfish 
#>                             1                             3 
#>                    Salamander 
#>                             3
```

- The prediction for each dyads in the networks, here for network
  number 3. If your goal is dyad prediction, then you should use
  `colsbm_model = "delta"`, instead of `colsbm_model = "iid"`.

``` r

best_fit$pred_dyads[[3]][1:10, 1:5]
#>                               Unidentified sp1 FW_012_02 Terrestrial plants
#> Unidentified sp1 FW_012_02                  0.000000e+00       3.120870e-08
#> Terrestrial plants                          3.120870e-08       0.000000e+00
#> Terrestrial bugs                            1.594940e-09       1.594940e-09
#> Achnanthes inflata var. elata               1.594955e-09       1.594955e-09
#> Achnanthes lanceolata                       1.594955e-09       1.594955e-09
#> Achnanthes linearis                         1.594933e-09       1.594933e-09
#> Achnanthes minutissima                      1.594956e-09       1.594956e-09
#> Auodinella hermanii                         1.594277e-09       1.594277e-09
#> Blue Green algae                            1.594955e-09       1.594955e-09
#> Calothrix                                   1.594953e-09       1.594953e-09
#>                               Terrestrial bugs Achnanthes inflata var. elata
#> Unidentified sp1 FW_012_02        9.384482e-03                  9.344968e-03
#> Terrestrial plants                9.384482e-03                  9.344968e-03
#> Terrestrial bugs                  0.000000e+00                  2.451596e-06
#> Achnanthes inflata var. elata     9.443825e-06                  0.000000e+00
#> Achnanthes lanceolata             9.443483e-06                  2.445121e-06
#> Achnanthes linearis               9.443514e-06                  2.445265e-06
#> Achnanthes minutissima            9.443463e-06                  2.445098e-06
#> Auodinella hermanii               9.452042e-06                  2.457419e-06
#> Blue Green algae                  9.443875e-06                  2.445532e-06
#> Calothrix                         9.444669e-06                  2.446368e-06
#>                               Achnanthes lanceolata
#> Unidentified sp1 FW_012_02             9.342819e-03
#> Terrestrial plants                     9.342819e-03
#> Terrestrial bugs                       2.044082e-06
#> Achnanthes inflata var. elata          2.037948e-06
#> Achnanthes lanceolata                  0.000000e+00
#> Achnanthes linearis                    2.037739e-06
#> Achnanthes minutissima                 2.037565e-06
#> Auodinella hermanii                    2.050105e-06
#> Blue Green algae                       2.038000e-06
#> Calothrix                              2.038839e-06
```

We can also plot the networks individually, with the groups reordered by
trophic levels:

``` r

p <- gtools::permutations(best_fit$Q, best_fit$Q)
ind <- which.min(
  sapply(
    seq(nrow(p)),
    function(x) {
      sum((tcrossprod(best_fit$pi[[1]]) * best_fit$alpha)[p[x, ], p[x, ]][
        upper.tri(best_fit$alpha)
      ])
    }
  )
)
ord <- p[ind, ]
plot(res_fw_iid$best_fit, type = "block", net_id = 1, ord = ord) +
  plot(res_fw_iid$best_fit, type = "block", net_id = 2, ord = ord) +
  plot(res_fw_iid$best_fit, type = "block", net_id = 3, ord = ord)
```

![](tutorial_files/figure-html/plot-block-1.png)

Or make different plots to exhibit the mesoscale structure:

``` r

plot(res_fw_iid$best_fit, type = "graphon", ord = ord)
```

![](tutorial_files/figure-html/unnamed-chunk-8-1.png)

``` r

plot(res_fw_iid$best_fit, type = "meso", mixture = TRUE, ord = ord)
```

![](tutorial_files/figure-html/unnamed-chunk-8-2.png)

## Clustering of networks

Let simulate some directed networks with a lower triangular structure
that looks alike foodwebs.

``` r

set.seed(1234)
alpha <- matrix(c(
  .05, .01, .01, .01,
  .3, .05, .01, .01,
  .5, .4, .05, .01,
  .1, .8, .1, .05
), 4, 4, byrow = TRUE)
pi <- c(.1, .2, .6, .1)
sim_net <-
  replicate(3,
    {
      X <-
        sbm::sampleSimpleSBM(100,
          blockProp = pi, connectParam = list(mean = alpha),
          directed = TRUE
        )
      X$rNetwork
      X$networkData
    },
    simplify = FALSE
  )
```

``` r

set.seed(1234)

net_clust <- clusterize_unipartite_networks(
  netlist = c(foodwebs[1:3], sim_net), # A list of networks
  colsbm_model = "iid", # The name of the model
  directed = TRUE, # Foodwebs are directed networks
  net_id = c(names(foodwebs)[1:3], "sim1", "sim2", "sim3"), # Name of the networks
  nb_run = 3L, # Nmber of runs of the algorithm
  global_opts = # List of options
    list(
      verbosity = 0, # Verbosity level of the algorithm
      plot_details = 0, # Monitoring plot of the algorithm
      Q_max = 9, # Max number of clusters
      backend = "parallel" # Backend for parallel computing
    ),
  verbose = FALSE
)
```

We can extract the best partition:

``` r

best_partition <- net_clust$partition
```

The plot of the mesoscale structure of the whole collection is the
following:

``` r

plot(best_partition[[1]])
```

![](tutorial_files/figure-html/unnamed-chunk-12-1.png)

but then we can compare the mesoscale structures of the 2 groups:

``` r

plot(best_partition[[1]],
  type = "graphon",
  ord = order(best_partition[[1]]$alpha %*% best_partition[[1]]$pi[[1]])
) +
  plot(best_partition[[2]],
    type = "graphon",
    ord = order(best_partition[[2]]$alpha %*% best_partition[[2]]$pi[[1]])
  ) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "1")
```

![](tutorial_files/figure-html/plot-part-1.png)
