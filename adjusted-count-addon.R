# AC variance

## SIM 1: tpr and fpr are known

ac_sim1 <- function(n_obs, sdev_pos, sdev_neg, prev, reps) {
  
  n_pos <- n_obs * prev
  n_neg <- n_obs - n_pos
  
  est_pos <- map_dbl(1:reps, ~ sum(rnorm(n_pos, mean = 1, sd = sdev_pos) >= 0.5))
  est_neg <- map_dbl(1:reps, ~ sum(rnorm(n_neg, mean = 0, sd = sdev_neg) >= 0.5))
  
  tpr <- pnorm(0.5, mean = 1, sd = sdev_pos, lower.tail = FALSE)
  fpr <- pnorm(0.5, mean = 0, sd = sdev_neg, lower.tail = FALSE)
  
  cc <- (est_pos + est_neg) / n_obs
  ac <- (cc - fpr) / (tpr - fpr)
  ac <- ac |> pmax(0) |> pmin(1)
  
  out <- sqrt(mean((ac - prev)^2, na.rm = TRUE))
  return(out)
}

## SIM 2: tpr and fpr are estimated

ac_sim2 <- function(n_train, n_test, prev_train, prev_test, sdev_pos, sdev_neg, reps) {
  
  n1_train <- n_train * prev_train
  n0_train <- n_train - n1_train
  
  n1_test <- n_test * prev_test
  n0_test <- n_test - n1_test
  
  tpr_est <- map_dbl(1:reps, ~ mean(rnorm(n1_train, mean = 1, sd = sdev_pos) >= 0.5))
  fpr_est <- map_dbl(1:reps, ~ mean(rnorm(n0_train, mean = 0, sd = sdev_pos) >= 0.5))
  
  est_pos <- map_dbl(1:reps, ~ sum(rnorm(n1_test, mean = 1, sd = sdev_pos) >= 0.5))
  est_neg <- map_dbl(1:reps, ~ sum(rnorm(n0_test, mean = 0, sd = sdev_neg) >= 0.5))
  
  cc <- (est_pos + est_neg) / n_test
  ac <- (cc - fpr_est) / (tpr_est - fpr_est)
  ac <- ac |> pmax(0) |> pmin(1)
  
  out <- sqrt(mean((ac - prev_test)^2, na.rm = TRUE))
  return(out)
}

## SIM 2: tpr and fpr are estimated

rskew <- function(n_obs, mu, sdev, skew) {
  dlt <- skew / sqrt(1 + skew^2)
  b <- sqrt(pi / (pi - 2*dlt^2)) * sdev
  a <- mu - b * dlt * sqrt(2 / pi)  
  
  out <- sn::rsn(n = n_obs, xi = a, omega = b, alpha = skew) %>% as.numeric()
  return(out)
}

## SIM 3: tpr and fpr are estimated with skewness parameter
ac_sim3 <- function(n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test, skew, reps) {
  
  n1_train <- n_train * prev_train
  n0_train <- n_train - n1_train
  
  n1_test <- n_test * prev_test
  n0_test <- n_test - n1_test
  
  tpr_est <- map_dbl(1:reps, ~ mean(rskew(n1_train, mu = 1, sdev = sdev_pos, skew = skew) >= 0.5))
  fpr_est <- map_dbl(1:reps, ~ mean(rskew(n0_train, mu = 0, sdev = sdev_pos, skew = skew) >= 0.5))
  
  est_pos <- map_dbl(1:reps, ~ sum(rskew(n1_test, mu = 1, sdev = sdev_pos, skew = skew) >= 0.5))
  est_neg <- map_dbl(1:reps, ~ sum(rskew(n0_test, mu = 0, sdev = sdev_neg, skew = skew) >= 0.5))
  
  cc <- (est_pos + est_neg) / n_test
  ac <- (cc - fpr_est) / (tpr_est - fpr_est)
  ac <- ac |> pmax(0) |> pmin(1)
  
  out <- sqrt(mean((ac - prev_test)^2, na.rm = TRUE))
  return(out)
}

set.seed(123)

### SIM 1 
sim1_layout <- 
  tidyr::expand_grid(n_obs = c(100, 1000),
                     sdev_pos = c(0.5, 1, 1.5),
                     sdev_neg = c(0.5, 1, 1.5),
                     prev = c(0.3, 0.5, 0.9),
                     reps = 10000)
sim1_addon <- 
  sim1_layout |>
  select(-reps) |>
  mutate(quantifier = "AC",
         rmse = pmap_dbl(sim1_layout, ac_sim1)) |>
  rename(prev_test = prev)

sim2_layout <- 
  tidyr::expand_grid(n_train = c(100, 1000),
                     n_test = c(100, 1000),
                     prev_train = c(0.4, 0.5, 0.6),
                     prev_test = c(0.3, 0.5, 0.9),
                     sdev_pos = c(0.5, 1, 1.5),
                     sdev_neg = c(0.5, 1, 1.5),
                     reps = 10000)

sim2_addon <- 
  sim2_layout |>
  select(-reps) |>
  mutate(quantifier = "AC",
         rmse = pmap_dbl(sim2_layout, ac_sim2))

sim3_layout <- 
  tidyr::expand_grid(n_train = 1000,
                     n_test = c(100, 1000),
                     sdev_pos = c(0.5, 1.5),
                     sdev_neg = c(0.5, 1.5),
                     prev_train = 0.5,
                     prev_test = c(0.3, 0.5, 0.9),
                     skew = c(1, 2, 4),
                     reps = 10000)

sim3_addon <- 
  sim3_layout |>
  select(-reps) |>
  mutate(quantifier = "AC",
         rmse = pmap_dbl(sim3_layout, ac_sim3))


