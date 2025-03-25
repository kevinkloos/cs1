library(batchtools)

big_sim_layout <- tidyr::expand_grid(n_train = 1000,
                                     n_test = c(100, 1000),
                                     sdev_pos = c(0.5, 1.5),
                             	     sdev_neg = c(0.5, 1.5),
                              	     prev_train = 0.5,
                                     prev_test = c(0.3, 0.5, 0.9),
                                     skew = c(1, 2, 4),
                                     reps = 1000)

reg <- batchtools::makeRegistry(file.dir = "cs_paper1/sim3",
                                source = "cs_paper1/cluster-input/sim3.R",
                                packages = c("tidyverse", "rootSolve", "pracma",
                                             "sn","EnvStats", "magrittr", "ks"),
                                seed = 1)

batchtools::batchMap(one_sim_function, args = big_sim_layout)
submitJobs(resources = list(walltime = "23:59:59", memory = 1024, ncpus = 1, partition = "cpu-medium"))
batchtools::getJobStatus()
