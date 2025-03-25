######################### FUNCTIONS #########################
simulate_test_set_gauss <- function(n_obs, sdev_pos, sdev_neg, prev) {
  stopifnot("n_obs * prev should be an integer" = ((n_obs * prev) %% 1) == 0)
  
  n_pos <- n_obs * prev
  n_neg <- n_obs - n_pos
  
  pos <- rnorm(n_pos, mean = 1, sd = sdev_pos)
  neg <- rnorm(n_neg, mean = 0, sd = sdev_neg)
  
  return(sort(c(pos, neg), decreasing = TRUE))
}

tpr_gauss <- function(theta, sdev_pos) {
  out <- pnorm(theta, mean = 1, sd = sdev_pos, lower.tail = FALSE)
  return(out)
}

fpr_gauss <- function(theta, sdev_neg) {
  out <- pnorm(theta, mean = 0, sd = sdev_neg, lower.tail = FALSE)
  return(out)
}

get_threshold_barriers_gauss <- function(sdev_pos, sdev_neg, trunc_value) {
  .tpr <- function(theta) {pnorm(theta, mean = 1, sd = sdev_pos, lower.tail = FALSE)}
  .fpr <- function(theta) {pnorm(theta, mean = 0, sd = sdev_neg, lower.tail = FALSE)}
  
  p1 <- uniroot(f = function(theta) {.tpr(theta) - .fpr(theta) - trunc_value},
                lower = -50, upper = 0.5)$root
  p2 <- uniroot(f = function(theta) {.tpr(theta) - .fpr(theta) - trunc_value},
                lower = 0.5, upper = 50)$root
  
  return(c(p1, p2))
}


median_sweep_gauss <- function(test_set, sdev_pos, sdev_neg, trunc_value = 0.25) {
  tpr <- tpr_gauss(test_set, sdev_pos)
  fpr <- fpr_gauss(test_set, sdev_neg)
  test_set <- sort(test_set, decreasing = TRUE)
  
  cac_vec <- (1 : length(test_set)) / length(test_set)
  ac <- (cac_vec - fpr) / (tpr - fpr)
  if(!is.null(trunc_value)) {
    barriers <- get_threshold_barriers_gauss(sdev_pos, sdev_neg, trunc_value)
    bool_crit <- between(test_set, barriers[1], barriers[2])
    ac <- subset(ac, bool_crit)
  }
  ms <- ac %>% median()
  return(ms)
}

adjusted_count_gauss <- function(cac, theta, sdev_pos, sdev_neg) {
  .tpr <- tpr_gauss(theta, sdev_pos)
  .fpr <- fpr_gauss(theta, sdev_neg)
  return((cac - .fpr) / (.tpr - .fpr))
}

continuous_sweep_gauss <- function(test_set, sdev_pos, sdev_neg, trunc_value) {
  
  .calculate_one_area <- function(cac, theta_l, theta_r, sdev_pos, sdev_neg) {
    out <- integrate(f = function(x) {adjusted_count_gauss(cac, x, sdev_pos, sdev_neg)},
                     lower = theta_l, upper = theta_r)$value
    return(out)
  }
  
  n_obs <- length(test_set)
  # make 3 vectors containing the left bound, right bound, and cac-value
  all_theta_r <- test_set %>% head(-1)
  all_theta_l <- test_set %>% tail(-1)
  all_cac <- seq(1/n_obs, 1 - 1/n_obs, 1/n_obs) 
  
  # filter out the intervals that fall outside the "classification rate" `trunc_value`
  barriers <- get_threshold_barriers_gauss(sdev_pos, sdev_neg, trunc_value)
  diff_range <- barriers %>% diff()
  bool_crit <- all_theta_r >= barriers[1] & all_theta_l <= barriers[2]
  n_trunc <- sum(bool_crit)
  
  # replace the values slightly below and sightly about barriers with correct values.
  all_theta_r <- subset(all_theta_r, bool_crit) %>% replace(1, barriers[2])
  all_theta_l <- subset(all_theta_l, bool_crit) %>% replace(n_trunc, barriers[1])
  all_cac <- subset(all_cac, bool_crit)
  
  all_areas <- 
    pmap_dbl(.l = list(all_cac, all_theta_l, all_theta_r), 
             .f =  ~ .calculate_one_area(..1, ..2, ..3, sdev_pos, sdev_neg))
  
  area_sum <- sum(all_areas)
  return(area_sum / diff_range)
}

compare_quantifiers_gauss <- function(n_obs, prev, sdev_pos, sdev_neg, trunc_value) {
  test_set <- simulate_test_set_gauss(n_obs, sdev_pos, sdev_neg, prev)
  tbl <- tibble(cs = continuous_sweep_gauss(test_set, sdev_pos, sdev_neg, trunc_value),
                ms = median_sweep_gauss(test_set, sdev_pos, sdev_neg, trunc_value))
  return(tbl)
}

ev_adjusted_count_gauss <- function(n_obs, sdev_pos, sdev_neg, prev, theta) {
  .tpr <- tpr_gauss(theta, sdev_pos)
  .fpr <- fpr_gauss(theta, sdev_neg)
  
  dfr_pos <- tibble(n_pos = 0 : (n_obs * prev),
                    prob_pos = dbinom(n_pos, n_obs * prev, .tpr))
  dfr_neg <- tibble(n_neg = 0 : (n_obs * (1 - prev)),
                    prob_neg = dbinom(n_neg, n_obs * (1 - prev), .fpr))
  
  dfr_tot <- 
    expand_grid(dfr_pos, dfr_neg) %>%
    mutate(cac = (n_pos + n_neg) / n_obs,
           prob_tot = prob_pos * prob_neg)
  
  dfr_summ <- dfr_tot %>%
    group_by(cac) %>%
    summarise(prob = sum(prob_tot)) %>%
    mutate(ac = map_dbl(cac, ~ adjusted_count_gauss(.x, theta, sdev_pos, sdev_neg)),
           ac_ev = ac * prob) 
  
  out <- dfr_summ %>% pull(ac_ev) %>% sum()
  
  return(out)
  
}

classify_and_count_covariance_gauss <- function(theta_x, theta_y, n_obs, 
                                                sdev_pos, sdev_neg, prev) {
  
  .tpr <- function(theta) {pnorm(theta, mean = 1, sd = sdev_pos, lower.tail = FALSE)}
  .fpr <- function(theta) {pnorm(theta, mean = 0, sd = sdev_neg, lower.tail = FALSE)}
  
  n0 <- n_obs * (1 - prev)
  n1 <- n_obs * prev
  
  pt1 <- .tpr(pmax(theta_x, theta_y)) - .tpr(theta_x) * .tpr(theta_y)
  pt2 <- .fpr(pmax(theta_x, theta_y)) - .fpr(theta_x) * .fpr(theta_y)
  
  return(n1 / n_obs^2 * pt1 + n0 / n_obs^2 * pt2)
} 

adjusted_count_covariance_gauss <- function(theta_x, theta_y, n_obs, 
                                            sdev_pos, sdev_neg, prev) {
  
  .tpr <- function(theta) {pnorm(theta, mean = 1, sd = sdev_pos, lower.tail = FALSE)}
  .fpr <- function(theta) {pnorm(theta, mean = 0, sd = sdev_neg, lower.tail = FALSE)}
  .diff_rates <- function(theta) {.tpr(theta) - .fpr(theta)}
  
  cov_th <- classify_and_count_covariance_gauss(theta_x, theta_y, n_obs, 
                                                sdev_pos, sdev_neg, prev)
  diff_rates_x <- .diff_rates(theta_x)
  diff_rates_y <- .diff_rates(theta_y)
  
  out <- cov_th / (diff_rates_x * diff_rates_y) 
  return(out)
}

get_threshold_barriers_gauss <- function(sdev_pos, sdev_neg, diff_thres) {
  .tpr <- function(theta) {pnorm(theta, mean = 1, sd = sdev_pos, lower.tail = FALSE)}
  .fpr <- function(theta) {pnorm(theta, mean = 0, sd = sdev_neg, lower.tail = FALSE)}
  
  p1 <- uniroot(f = function(theta) {.tpr(theta) - .fpr(theta) - diff_thres},
                lower = -50, upper = 0.5)$root
  p2 <- uniroot(f = function(theta) {.tpr(theta) - .fpr(theta) - diff_thres},
                lower = 0.5, upper = 50)$root
  
  return(c(p1, p2))
}


continuous_sweep_integration_gauss <- function(n_obs, sdev_pos, sdev_neg, 
                                               trunc_value = NULL, prev) {
  # obtain decision boundaries
  barrier_levels <- get_threshold_barriers_gauss(sdev_pos, sdev_neg, trunc_value)
  
  # make sure that y is never larger than x in the integral
  ymax <- function(x) {x}
  # compute integral
  intgrl <- 
    pracma::integral2(fun = function(xi, yi) {
      adjusted_count_covariance_gauss(xi, yi, n_obs, sdev_pos, sdev_neg, prev)
    }, barrier_levels[1], barrier_levels[2], barrier_levels[1], ymax)$Q 
  
  scaling_factor <- diff(barrier_levels)^2
  out <- 2 * intgrl / scaling_factor
  
  return(out)
}

data_optim_thres <- function(n_obs, prev, sdev_pos, sdev_neg, step_rates) {
  
  out <- 
    tibble(diff_rates = seq(0, 1, length.out = step_rates),
           var_val = map_dbl(seq(0, 1, length.out = step_rates), function(x) {
             val_temp <- try(continuous_sweep_integration_gauss(n_obs, sdev_pos, 
                                                                sdev_neg, x, prev),
                             silent = T)
             out <- ifelse(class(val_temp) == "try-error", NaN, val_temp)
             return(out)})) %>%
    mutate(min_point = (var_val <= min(var_val, na.rm = T)),
           min_point = replace_na(min_point, FALSE),
           last_num = map2_lgl(var_val, lead(var_val), ~ !is.nan(.x) & is.nan(.y)))
  
  return(out)
}

SLD_gauss <- function(test_set, sdev_pos, sdev_neg) {
  # compute scores
  alpha_train <- c(0.5, 0.5)
  test_scores <- 
    dnorm(test_set, 1, sdev_pos) / (dnorm(test_set, 1, sdev_pos) + dnorm(test_set, 0, sdev_neg))
  
  SLD_alg <- function(prob_test, alpha_train, maxiter = 1000, crit = 0.0001) { 
    # maxiter and crit from QuaPy
    
    prob_test <- cbind(1 - prob_test, prob_test)
    alpha_train <- alpha_train %>% matrix(nrow = 1)
    alpha_new <- alpha_train
    prob_test <- as.matrix(prob_test)
    iter <- 0
    diff <- Inf
    
    .updater <- function(alpha_new, alpha_train, prob_test) {
      n <- length(alpha_train)
      # E-step
      unscaled <- 
        map(1:2, ~ alpha_new[.x] / alpha_train[.x] * prob_test[,.x]) %>%
        bind_cols(.name_repair = "unique") %>%
        suppressMessages()
      scaled <- unscaled / rowSums(unscaled)
      
      # U-step
      out <- colMeans(scaled)
      names(out) <- NULL
      return(out)
    }
    
    while(iter < maxiter && crit <= diff) {
      alpha_upd <- .updater(alpha_new, alpha_train, prob_test)
      diff <- abs(alpha_upd - alpha_new) %>% mean()
      iter <- iter + 1
      alpha_new <- alpha_upd
    }
    return(alpha_new)
  }
  
  test_predict_SLD <- SLD_alg(test_scores, alpha_train)
  return(test_predict_SLD)
}

DyS_prep <- function(sdev_pos, sdev_neg, n_bins = 8, tol = 10^(-5)) {
  .compute_histogram <- function(scores, n_bins, rng) {
    bns_selec <- cut(x = scores, breaks = seq(rng[1], rng[2], length.out = n_bins + 1))
    dns <- table(bns_selec) / length(bns_selec)
    names(dns) <- NULL
    dns <- as.vector(dns)
    return(dns)
  }
  
  rng <- c(0, 1)
  # obtain train and test scores (prob)
  rnd_p <- rnorm(1e7, 1, sdev_pos)
  rnd_n <- rnorm(1e7, 0, sdev_neg)
  
  pos_scores <- dnorm(rnd_p, 1, sdev_pos) / 
    (dnorm(rnd_p, 1, sdev_pos) + dnorm(rnd_p, 0, sdev_neg))
  neg_scores <- dnorm(rnd_n, 1, sdev_pos) / 
    (dnorm(rnd_n, 1, sdev_pos) + dnorm(rnd_n, 0, sdev_neg))
  
  # compute histogram of positive and negative scores
  h0 <- .compute_histogram(neg_scores, n_bins, rng)
  h1 <- .compute_histogram(pos_scores, n_bins, rng)
  
  return(list(h0 = h0, h1 = h1))
}

DyS_gauss <- function(test_set, sdev_pos, sdev_neg, n_bins = 8, tol = 10^(-5), prep) {
  # helper function to compute DyS
  .compute_histogram <- function(scores, n_bins, rng) {
    bns_selec <- cut(x = scores, breaks = seq(rng[1], rng[2], length.out = n_bins + 1))
    dns <- table(bns_selec) / length(bns_selec)
    names(dns) <- NULL
    dns <- as.vector(dns)
    return(dns)
  }
  
  .hellinger_distance <- function(p_vec, q_vec) {
    sqsq <- (sqrt(p_vec) - sqrt(q_vec))^2
    out <- sqrt(sum(sqsq)) / sqrt(2)
    return(out)
  }
  
  .tenary_search <- function(lft, rgt, h_dist, tol, h0, h1, h_test, mthd) {
    conv <- FALSE
    while(!conv) {
      lft3 <- lft + (rgt - lft) / 3
      rgt3 <- rgt - (rgt - lft) / 3
      
      if(h_dist(lft3, h0, h1, h_test, mthd) > h_dist(rgt3, h0, h1, h_test, mthd)) {
        lft <- lft3
      }
      else {
        rgt <- rgt3
      }
      conv <- abs(lft-rgt) < tol
    }
    out <- (lft + rgt) / 2
    return(out)
  }
  
  .histogram_distance <- function(prev_est, h0, h1, h_test, mthd) {
    h_train <- (1 - prev_est) * h0 + prev_est * h1
    dist_trte <- mthd(h_train, h_test)
    return(dist_trte)
  }
  
  .compute_dys <- function(tst_scr, h0, h1, n_bins, tol, rng, mthd) {
    h_test <- .compute_histogram(tst_scr, n_bins, rng)
    out <- .tenary_search(0, 1, .histogram_distance, tol, h0, h1, h_test, mthd)
    return(out)
  } 
  
  # obtain range and distance method, can be expanded to more
  rng <- c(0, 1)
  mthd <- get(".hellinger_distance")
  
  
  test_scores <- 
    dnorm(test_set, 1, sdev_pos) / (dnorm(test_set, 1, sdev_pos) + dnorm(test_set, 0, sdev_neg))
  
  # minimize distance mixture and test histograms with DyS
  out <- .compute_dys(test_scores, prep$h0, prep$h1, n_bins, tol, rng, mthd)
  return(out)
}

one_sim_function <- function(n_obs, sdev_pos, sdev_neg, prev, reps, row_nrs = NULL) {
  
  opt_thres <- 
    data_optim_thres(n_obs, prev, sdev_pos, sdev_neg, 10000) %>% 
    filter(min_point) %>% 
    pull(diff_rates)
  
  prep_dys <- DyS_prep(sdev_pos, sdev_neg)

  out <- 
    map(1:reps, function(x) {
      test_set <- simulate_test_set_gauss(n_obs, sdev_pos, sdev_neg, prev)
      quants <- list(cs_trad = continuous_sweep_gauss(test_set, sdev_pos, sdev_neg, 0.25),
                     ms_trad = median_sweep_gauss(test_set, sdev_pos, sdev_neg, 0.25),
                     cs_opt = continuous_sweep_gauss(test_set, sdev_pos, sdev_neg, opt_thres),
                     ms_opt = median_sweep_gauss(test_set, sdev_pos, sdev_neg, opt_thres),
                     dys = DyS_gauss(test_set, sdev_pos, sdev_neg, prep = prep_dys),
                     sld = SLD_gauss(test_set, sdev_pos, sdev_neg)[2],
                     opt_thres = opt_thres)

      return(quants)
    })

  return(out)
}
