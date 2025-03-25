library(batchtools)

big_sim_layout <- expand_grid(n_train = c(100, 1000),
			      n_test = c(100, 1000),
                              prev_train = c(0.4, 0.5, 0.6),
                              prev_test = c(0.3, 0.5, 0.9),
                              sdev_pos = c(0.5, 1, 1.5),
                              sdev_neg = c(0.5, 1, 1.5),
                              reps = 10000)

reg <- batchtools::makeRegistry(file.dir = "cs_paper1/sim2",
                                source = "cs_paper1/cluster-input/sim2.R",
                                packages = c("tidyverse", "rootSolve", "pracma",
                                             "sn","EnvStats", "magrittr", "ks"),
                                seed = 1)

batchtools::batchMap(one_sim_function, args = big_sim_layout)
submitJobs(resources = list(walltime = "23:59:59", memory = 1024, ncpus = 1, partition = "cpu-medium"))
batchtools::getJobStatus()
