devtools::load_all()
library(here)

fit_iid <- readRDS(here("dev", "prefit_9_sim_iid.Rds"))

netlist <- fit_iid$A
names(netlist) <- fit_iid$net_id

partition_init <- readRDS(here("dev", "partition.Rds"))

# desc_res <- clusterize_bipartite_networks(
#     netlist = netlist,
#     net_id = names(netlist),
#     colsbm_model = "iid",
#     full_collection_init = fit_iid,
#     # partition_init = partition_init,
#     global_opts = list(
#         nb_cores = 1L,
#         backend = "no_mc"
#     )
# )
library(future)
options(future.globals.maxSize = 10 * 1024^3) # 100 MB
plan("multisession", workers = 2L)
res <- clusterize_bipartite_networks_d_a(
    netlist = netlist,
    net_id = names(netlist),
    colsbm_model = "iid",
    full_collection_init = fit_iid,
    # partition_init = partition_init,
    global_opts = list(
        nb_cores = 1L,
        backend = "future"
    )
)
