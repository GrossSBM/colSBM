# colSBM modeling context

This context captures the shared language for the statistical model-fitting workflow in colSBM, especially the parts that are candidates for numerical acceleration.

## Language

**Model exploration**:
The outer search over candidate block counts and model families that evaluates many fits and keeps the most promising ones.
_Avoid_: brute-force search, tuning

**Model fitting**:
The inner optimization of a single statistical model for a fixed block count and initialization.
_Avoid_: training, estimation loop

**Variational bound**:
The objective used to judge whether a model fit is improving, often summarized as a vbound or BIC-L value.
_Avoid_: loss, score

**Initialization**:
The starting assignment of nodes to blocks used to begin a model fit.
_Avoid_: seed, config

**Membership**:
A probabilistic assignment of a node or network to a block.
_Avoid_: label, cluster

**Collection**:
A set of networks treated jointly by the colSBM procedure.
_Avoid_: dataset, corpus

**Fit object**:
The structured result of a model fit, including parameters, memberships, and diagnostics.
_Avoid_: result, output
