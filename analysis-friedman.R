source("_R/data-load-and-clean.R")

friedman_dfr1 <- 
  expand_grid(n_obs = c(100, 1000), 
              sdev_pos = c(0.5, 1, 1.5),
              sdev_neg = c(0.5, 1, 1.5),
              prev_test = c(0.3, 0.5, 0.9)) |>
  mutate(cond = 1:54) |>
  right_join(y = sim1_clean) |>
  mutate(sqe = err^2) |>
  group_by(cond, quantifier) |>
  summarise(rmse = mean(sqrt(sqe))) |>
  bind_rows(sim1_addon |> mutate(cond = 1:54) |> select(cond, quantifier, rmse)) |>
  pivot_wider(names_from = "quantifier", values_from = "rmse")

nemenyi1 <-
  scmamp::postHocTest(data = friedman_dfr1, algorithms = c("O-CS", "T-CS", "DyS", "O-MS", "T-MS", "SLD", "AC"),
                      group.by = "cond", test = 'friedman', correct = 'holm')

scmamp::plotCD(friedman_dfr1 |> ungroup() |> select(-cond), 
               alpha = 0.05, cex = 1.2, main = "Critical Difference Plot (AIC)",  
               y.axis = 1, decreasing = FALSE)

friedman_dfr2 <- 
  sim2_clean |>
  select(iter, quantifier, err, n_train:sdev_neg) |>
  group_by(n_train, n_test, prev_train, prev_test, sdev_pos, sdev_neg, quantifier) |>
  summarise(rmse = sqrt(mean(err^2, na.rm = TRUE))) |>
  bind_rows(sim2_addon) |>
  pivot_wider(names_from = "quantifier", values_from = "rmse") |>
  ungroup() |>
  select(`O-CS` : `AC`)

nemenyi2 <-
  scmamp::postHocTest(data = friedman_dfr2, algorithms = c("O-CS", "T-CS", "DyS", "O-MS", "T-MS", "SLD", "AC"),
                      test = 'friedman', correct = 'holm')

scmamp::plotCD(friedman_dfr2, 
               alpha = 0.05, cex = 1.25, main = "Critical Difference Plot (AIC)",  
               y.axis = 1, decreasing = FALSE)

friedman_dfr3a <- 
  sim3_clean |>
  filter(skew == 1) |>
  select(iter, quantifier, err, n_train:prev_test) |>
  group_by(n_train, n_test, prev_train, prev_test, sdev_pos, sdev_neg, quantifier) |>
  summarise(rmse = sqrt(mean(err^2, na.rm = TRUE))) |>
  bind_rows(sim3_addon |> filter(skew == 1) |> select(-skew)) |>
  pivot_wider(names_from = "quantifier", values_from = "rmse") |>
  ungroup() |>
  select(`O-CS (norm)` : `AC`)

nemenyi3a <-
  scmamp::postHocTest(data = friedman_dfr3a, 
                      algorithms = c("O-CS (norm)", "O-CS (skew)", "T-CS (norm)",
                                     "T-CS (skew)", "DyS", "O-MS", "T-MS", "SLD", "AC"),
                      test = 'friedman', correct = 'holm')

scmamp::plotCD(friedman_dfr3a, 
               alpha = 0.05, cex = 1, main = "Critical Difference Plot (AIC)",  
               y.axis = 1, decreasing = FALSE)

friedman_dfr3b <- 
  sim3_clean |>
  filter(skew == 2) |>
  select(iter, quantifier, err, n_train:prev_test) |>
  group_by(n_train, n_test, prev_train, prev_test, sdev_pos, sdev_neg, quantifier) |>
  summarise(rmse = sqrt(mean(err^2, na.rm = TRUE))) |>
  bind_rows(sim3_addon |> filter(skew == 2) |> select(-skew)) |>
  pivot_wider(names_from = "quantifier", values_from = "rmse") |>
  ungroup() |>
  select(`O-CS (norm)` : `AC`)

nemenyi3b <-
  scmamp::postHocTest(data = friedman_dfr3b, 
                      algorithms = c("O-CS (norm)", "O-CS (skew)", "T-CS (norm)",
                                     "T-CS (skew)", "DyS", "O-MS", "T-MS", "SLD", "AC"),
                      test = 'friedman', correct = 'holm')

scmamp::plotCD(friedman_dfr3b, 
               alpha = 0.05, cex = 1, main = "Critical Difference Plot (AIC)",  
               y.axis = 1, decreasing = FALSE)

friedman_dfr3c <- 
  sim3_clean |>
  filter(skew == 4) |>
  select(iter, quantifier, err, n_train:prev_test) |>
  group_by(n_train, n_test, prev_train, prev_test, sdev_pos, sdev_neg, quantifier) |>
  summarise(rmse = sqrt(mean(err^2, na.rm = TRUE))) |>
  bind_rows(sim3_addon |> filter(skew == 4) |> select(-skew)) |>
  pivot_wider(names_from = "quantifier", values_from = "rmse") |>
  ungroup() |>
  select(`O-CS (norm)` : `AC`)

nemenyi3c <-
  scmamp::postHocTest(data = friedman_dfr3c, 
                      algorithms = c("O-CS (norm)", "O-CS (skew)", "T-CS (norm)",
                                     "T-CS (skew)", "DyS", "O-MS", "T-MS", "SLD", "AC"),
                      test = 'friedman', correct = 'holm')

scmamp::plotCD(friedman_dfr3c, 
               alpha = 0.05, cex = 1, main = "Critical Difference Plot (AIC)",  
               y.axis = 1, decreasing = FALSE)
