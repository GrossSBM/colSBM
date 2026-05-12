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
#> BICL =  -1965.85  -- #Empty blocks :  0  
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
#>             [,1]        [,2]         [,3]         [,4]       [,5]
#> [1,] 0.001788899 0.000021912 1.043364e-10 5.504046e-12 0.02699628
#> [2,] 0.026313659 0.006637112 9.852776e-10 2.755205e-10 0.20532212
#> [3,] 0.686730891 0.798570147 2.916666e-08 8.445247e-10 0.58270643
#> [4,] 0.007863483 0.141342696 1.444459e-09 7.963034e-11 0.01108018
#> [5,] 0.004367508 0.002717288 4.651027e-10 2.410493e-11 0.12571554
#> 
#> $pi
#> $pi[[1]]
#> [1] 0.21040608 0.23172923 0.02564103 0.45604944 0.07617421
#> 
#> $pi[[2]]
#> [1] 0.21040608 0.23172923 0.02564103 0.45604944 0.07617421
#> 
#> $pi[[3]]
#> [1] 0.21040608 0.23172923 0.02564103 0.45604944 0.07617421
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
#>                            3                            3 
#>    Terrestrial invertebrates        Achnanthes lanceolata 
#>                            2                            4 
#>              Batrachospermum         Calothrix sp1 FW_009 
#>                            4                            4 
#>         Cocconeis placentula         Cosmarium sp1 FW_009 
#>                            4                            4 
#>        Cyclotella sp1 FW_009              Cymbella aspera 
#>                            4                            4 
#>           Cymbella cuspidata               Cymbella kappi 
#>                            4                            4 
#>              Cymbella tumida              Diatoma heimale 
#>                            4                            4 
#>              Epithemia sorex            Epithemia turgida 
#>                            4                            4 
#>           Eunotia serpentina      Unidentified sp2 FW_009 
#>                            4                            4 
#>           Euntoia pectinalis        Fragilaria sp1 FW_009 
#>                            4                            4 
#>        Fragilaria vaucheriae         Frustulia rhomboides 
#>                            4                            4 
#>        Gomphoneis herculeana       Gomphonema accuminatum 
#>                            4                            4 
#>        Gomphonema angustatum         Gomphonema truncatum 
#>                            4                            4 
#>      Unidentified sp3 FW_009             Melosira varians 
#>                            4                            4 
#>            Navicula avenacea           Navicula dicephala 
#>                            4                            4 
#>              Nitzschia dubia        Oedogonium sp1 FW_009 
#>                            4                            4 
#>        Phormidium sp1 FW_009         Pinnularia mesolepta 
#>                            4                            4 
#>           Pinnularia viridis    Pleaurotaenium sp1 FW_009 
#>                            4                            4 
#>         Rhoicospenia curvata        Rhopalodia sp1 FW_009 
#>                            4                            4 
#>       Schizothrix sp1 FW_009       Staurastrum sp1 FW_009 
#>                            4                            4 
#>            Surirella elegans             Surirella tenera 
#>                            4                            4 
#>                 Synedra ulna        Tabellaria fenestrata 
#>                            4                            4 
#>        Tabellaria flocculosa          Ulothrix sp1 FW_009 
#>                            4                            4 
#>      Unidentified sp4 FW_009      Unidentified sp5 FW_009 
#>                            4                            4 
#>        Acroneuria sp1 FW_009          Aelosoma sp1 FW_009 
#>                            5                            2 
#>         Alloperla sp1 FW_009         Anepeorus sp1 FW_009 
#>                            1                            5 
#>             Antocha saxicola                       Baetis 
#>                            2                            2 
#>               Boyeria vinosa Bryophaenocladius sp1 FW_009 
#>                            1                            2 
#>        Chauliodes sp1 FW_009            Chimarra atterima 
#>                            5                            2 
#>      Unidentified sp6 FW_009     Conchapelopia sp1 FW_009 
#>                            2                            5 
#>       Cryptolabis sp1 FW_009         Cyrnellus sp1 FW_009 
#>                            1                            1 
#>     Dicrotendipes sp1 FW_009          Diplectrona modesta 
#>                            2                            2 
#>             Ectopria nervosa   Endochironomous sp1 FW_009 
#>                            2                            2 
#>       Ephemerella sp1 FW_009  Eukieferiella pseudomontana 
#>                            1                            2 
#>   Eukiefferiella 'dark' type        Glossosoma sp1 FW_009 
#>                            2                            1 
#>          Gyraulus sp1 FW_009            Haploperla brevis 
#>                            2                            1 
#>          Hexatoma sp1 FW_009                Hydrophilidae 
#>                            5                            1 
#>  Hydropsyche sp1 FW_009arana            Larsia sp1 FW_009 
#>                            5                            1 
#>        Leucrocuta sp1 FW_009           Leuctra sp1 FW_009 
#>                            2                            2 
#>      Unidentified sp7 FW_009      Metriocnemus sp1 FW_009 
#>                            2                            2 
#>         Micrasema sp1 FW_009        Ochthebius sp1 FW_009 
#>                            1                            1 
#>      Ophiogomphus sp1 FW_009  Paraleptophlebia sp1 FW_009 
#>                            5                            1 
#>  Paranyctiophylax sp1 FW_009     Polycentropus sp1 FW_009 
#>                            1                            1 
#>         Probezzia sp1 FW_009        Promoresia sp1 FW_009 
#>                            1                            2 
#>         Psephenus sp1 FW_009 Pseudolimnolphila sp1 FW_009 
#>                            2                            1 
#>       Rhyacophila sp1 FW_009          Simulium sp1 FW_009 
#>                            1                            2 
#>        Sphaerium occidentale     Stempelinella sp1 FW_009 
#>                            2                            2 
#>            Stenelmis crenata         Stenelmis sp1 FW_009 
#>                            2                            2 
#>           Suwalia sp1 FW_009        Tallaperla sp1 FW_009 
#>                            1                            2 
#>           Tanytarsus Genus A     Tricorythodes sp1 FW_009 
#>                            1                            1 
#>        Notropis heterolepsis                  Brook trout 
#>                            5                            5 
#>           Orcnocetes virilis       Rhinichthys cataractae 
#>                            1                            5 
#>            Cambarus bartonii 
#>                            1 
#> 
#> [[2]]
#>     Unidentified sp1 FW_012_01      Terrestrial invertebrates 
#>                              3                              4 
#>                 Plant material          Achnanthes lanceolata 
#>                              3                              4 
#>         Achnanthes minutissima      Audouinella sp1 FW_012_01 
#>                              4                              4 
#>                Batrachospermum               Blue-green algae 
#>                              4                              4 
#>                      Calothrix               Cymbella cistula 
#>                              4                              4 
#>               Cymbella mulleri                Diatoma heimale 
#>                              4                              4 
#>                Epithemia sorex              Epithemia turgida 
#>                              4                              4 
#>             Eunotia pectinalis          Eunotia sp1 FW_012_01 
#>                              4                              4 
#>           Frustulia rhomboides          Gomphoneis herculeana 
#>                              4                              4 
#>          Gomphonema intricatum       Gomphonema sp1 FW_012_01 
#>                              4                              4 
#>             Meridion circulare              Navicula avenacea 
#>                              4                              4 
#>                  Pleurotaenium          Rhoicosphenia curvata 
#>                              4                              4 
#>                  Stigeoclonium                   Synedra ulna 
#>                              4                              4 
#>                       Ulothrix     Unidentified sp2 FW_012_01 
#>                              4                              4 
#>                       Aelosoma    Brachycentrus sp1 FW_012_01 
#>                              2                              2 
#>               Cambarus bartoni       Chauliodes sp1 FW_012_01 
#>                              1                              1 
#>         Cordulegaster maculata        Dicranota sp1 FW_012_01 
#>                              1                              1 
#>             Ectopria thoracica                 Epeorus dispar 
#>                              2                              2 
#>       Glossosoma sp1 FW_012_01      Homoplectra sp1 FW_012_01 
#>                              1                              2 
#>       Hudsonimya sp1 FW_012_01      Hydropsyche sp1 FW_012_01 
#>                              1                              5 
#>       Leucrocuta sp1 FW_012_01          Leuctra sp1 FW_012_01 
#>                              2                              2 
#>       Lumbriculiid oligochaete Parametriocnemus sp1 FW_012_01 
#>                              2                              2 
#>     Neureclipsis sp1 FW_012_01     Ophiogomphus sp1 FW_012_01 
#>                              2                              1 
#>        Palpomyia sp1 FW_012_01        Palpomyia sp2 FW_012_01 
#>                              1                              1 
#>       Promoresia sp1 FW_012_01        Psephenus sp1 FW_012_01 
#>                              2                              2 
#>        Soliperla sp1 FW_012_01                Stenelmis adult 
#>                              5                              1 
#>        Stenelmis sp1 FW_012_01          Suwalia sp1 FW_012_01 
#>                              2                              1 
#>               Tallaperla maria        Thaumalea sp1 FW_012_01 
#>                              2                              2 
#>             Tipula abdominalis                     Salamander 
#>                              2                              5 
#> 
#> [[3]]
#>    Unidentified sp1 FW_012_02            Terrestrial plants 
#>                             3                             3 
#>              Terrestrial bugs Achnanthes inflata var. elata 
#>                             4                             4 
#>         Achnanthes lanceolata           Achnanthes linearis 
#>                             4                             4 
#>        Achnanthes minutissima           Auodinella hermanii 
#>                             4                             4 
#>              Blue Green algae                     Calothrix 
#>                             4                             4 
#>          Cocconeis placentula                Cymbella kappi 
#>                             4                             4 
#>              Cymbella mulleri               Diatoma heimale 
#>                             4                             4 
#>             Epithemia turgida              Eunotia meisteri 
#>                             4                             4 
#>            Eunotia pectinalis         Fragilaria vaucheriae 
#>                             4                             4 
#>          Frustulia rhomboides         Gomphoneis herculeana 
#>                             4                             4 
#>        Gomphonema accuminatum         Gomphonema angustatum 
#>                             4                             4 
#>         Gomphonema intricatum         Gomphonema tennuellum 
#>                             4                             4 
#>                  Marssoniella             Navicula avenacea 
#>                             4                             4 
#>        Navicula cryptocephala               Navicula mutica 
#>                             4                             4 
#>        Navicula sp1 FW_012_02            Pinnularia viridis 
#>                             4                             4 
#>         Rhoicosphenia curvata                    Rhopalodia 
#>                             4                             4 
#>         Surirella brebbisonii             Surirella elegans 
#>                             4                             4 
#>                   Synechoccus                  Synedra ulna 
#>                             4                             4 
#>                      Ulothrix    Unidentified sp2 FW_012_02 
#>                             4                             4 
#>       Aeolosoma sp1 FW_012_02                 Ajax longipes 
#>                             2                             5 
#>               Amphinemura wui           Anchytarsus bicolor 
#>                             2                             2 
#>                        Baetis                    Hudsonimya 
#>                             2                             1 
#>                    Cricotopus                     Dicranota 
#>                             2                             1 
#>  Eukieffidrella pseudomontana           Diplectrona modesta 
#>                             1                             1 
#>                       Dixella                  Dolophilodes 
#>                             1                             1 
#>            Ectopria thoracica                Epeorus dispar 
#>                             2                             2 
#>                  Fatigia pele        Hexatoma sp1 FW_012_02 
#>                             1                             1 
#>                       Leuctra             Oligo Lumbr. Blue 
#>                             2                             2 
#>            Oligo. Lumbr. Pink                  Ophiogomphus 
#>                             2                             1 
#>             Paraleptophelebia                      Pericoma 
#>                             1                             2 
#>                       Pilaria      Pentaneuri sp1 FW_012_02 
#>                             1                             2 
#>       Polycentropus maculatus                     Stenelmis 
#>                             1                             2 
#>              Tallaperla maria                     Tanyderid 
#>                             1                             1 
#>                 Conchapelopia                        Tipula 
#>                             1                             1 
#>              Wormaldia moesta                      Crayfish 
#>                             1                             1 
#>                    Salamander 
#>                             5
```

- The prediction for each dyads in the networks, here for network
  number 3. If your goal is dyad prediction, then you should use
  `colsbm_model = "delta"`, instead of `colsbm_model = "iid"`.

``` r

best_fit$pred_dyads[[3]][1:10, 1:5]
#>                               Unidentified sp1 FW_012_02 Terrestrial plants
#> Unidentified sp1 FW_012_02                  0.000000e+00       3.123467e-08
#> Terrestrial plants                          3.123467e-08       0.000000e+00
#> Terrestrial bugs                            1.316430e-09       1.316430e-09
#> Achnanthes inflata var. elata               1.604745e-09       1.604745e-09
#> Achnanthes lanceolata                       1.604746e-09       1.604746e-09
#> Achnanthes linearis                         1.604745e-09       1.604745e-09
#> Achnanthes minutissima                      1.604746e-09       1.604746e-09
#> Auodinella hermanii                         1.604695e-09       1.604695e-09
#> Blue Green algae                            1.604746e-09       1.604746e-09
#> Calothrix                                   1.604746e-09       1.604746e-09
#>                               Terrestrial bugs Achnanthes inflata var. elata
#> Unidentified sp1 FW_012_02         0.142102573                  1.364604e-06
#> Terrestrial plants                 0.142102573                  1.364604e-06
#> Terrestrial bugs                   0.000000000                  1.340432e-07
#> Achnanthes inflata var. elata      0.002847148                  0.000000e+00
#> Achnanthes lanceolata              0.002847142                  1.588838e-07
#> Achnanthes linearis                0.002847143                  1.588838e-07
#> Achnanthes minutissima             0.002847141                  1.588838e-07
#> Auodinella hermanii                0.002847173                  1.588813e-07
#> Blue Green algae                   0.002847142                  1.588838e-07
#> Calothrix                          0.002847142                  1.588838e-07
#>                               Achnanthes lanceolata
#> Unidentified sp1 FW_012_02             1.668519e-08
#> Terrestrial plants                     1.668519e-08
#> Terrestrial bugs                       1.607678e-09
#> Achnanthes inflata var. elata          1.917182e-09
#> Achnanthes lanceolata                  0.000000e+00
#> Achnanthes linearis                    1.917181e-09
#> Achnanthes minutissima                 1.917182e-09
#> Auodinella hermanii                    1.917145e-09
#> Blue Green algae                       1.917182e-09
#> Calothrix                              1.917182e-09
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
