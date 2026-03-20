test_that("extract_nodes_groups errors on invalid input", {
  expect_error(
    extract_nodes_groups(1L),
    "must be a fitSimpleSBMPop, fitBipartiteSBMPop, bmpop or bisbmpop object"
  )
})

test_that("extract_nodes_groups works on fitSimpleSBMPop", {
  fit_simple <- structure(
    list(
      Z = list(
        c(a = 1L, b = 2L),
        c(c = 2L, d = 1L)
      ),
      net_id = c("net_1", "net_2")
    ),
    class = "fitSimpleSBMPop"
  )

  out <- extract_nodes_groups(fit_simple)

  expect_s3_class(out, "data.frame")
  expect_equal(colnames(out), c("network", "node_name", "cluster"))
  expect_equal(nrow(out), 4L)
  expect_setequal(out$network, c("net_1", "net_2"))
  expect_setequal(out$node_name, c("a", "b", "c", "d"))
  expect_setequal(as.integer(out$cluster), c(1L, 2L, 2L, 1L))
})

test_that("extract_nodes_groups works on fitBipartiteSBMPop", {
  fit_bipartite <- structure(
    list(
      memberships = list(
        list(
          row = c(r1 = 1L, r2 = 2L),
          col = c(c1 = 2L, c2 = 1L)
        )
      ),
      net_id = "net_bi"
    ),
    class = "fitBipartiteSBMPop"
  )

  out <- extract_nodes_groups(fit_bipartite)

  expect_s3_class(out, "data.frame")
  expect_equal(colnames(out), c("network", "node_name", "cluster", "node_type"))
  expect_equal(nrow(out), 4L)
  expect_true(all(out$network == "net_bi"))
  expect_setequal(out$node_type, c("row", "col"))
  expect_setequal(out$node_name[out$node_type == "row"], c("r1", "r2"))
  expect_setequal(out$node_name[out$node_type == "col"], c("c1", "c2"))
})

test_that("extract_nodes_groups extracts from best_fit for bmpop and bisbmpop", {
  fit_simple <- structure(
    list(
      Z = list(c(a = 1L, b = 2L)),
      net_id = "net_1"
    ),
    class = "fitSimpleSBMPop"
  )

  bmpop_fit <- structure(list(best_fit = fit_simple), class = "bmpop")
  bisbmpop_fit <- structure(list(best_fit = fit_simple), class = "bisbmpop")

  expect_equal(extract_nodes_groups(bmpop_fit), extract_nodes_groups(fit_simple))
  expect_equal(extract_nodes_groups(bisbmpop_fit), extract_nodes_groups(fit_simple))
})
