######################### FUNCTIONS ##########################
rskew <- function(n_obs, mu, sdev, skew) {
  dlt <- skew / sqrt(1 + skew^2)
  b <- sqrt(pi / (pi - 2*dlt^2)) * sdev
  a <- mu - b * dlt * sqrt(2 / pi)  
  
  out <- sn::rsn(n = n_obs, xi = a, omega = b, alpha = skew) %>% as.numeric()
  return(out)
}

dskew <- function(theta, mu, sdev, skew) {
  dlt <- skew / sqrt(1 + skew^2)
  b <- sqrt(pi / (pi - 2*dlt^2)) * sdev
  a <- mu - b * dlt * sqrt(2 / pi)  
  
  out <- sn::dsn(x = theta, xi = a, omega = b, alpha = skew) %>% as.numeric()
  return(out)
}

pskew <- function(theta, mu, sdev, skew) {
  dlt <- skew / sqrt(1 + skew^2)
  b <- sqrt(pi / (pi - 2*dlt^2)) * sdev
  a <- mu - b * dlt * sqrt(2 / pi)  
  
  out <- 1 - sn::psn(x = theta, xi = a, omega = b, alpha = skew) %>% as.numeric()
  return(out)
}

simulate_test_skew <- function(n_obs, sdev_pos, sdev_neg, prev, skew, prev_train = NULL) {
  stopifnot("n_obs * prev should be an integer" = ((n_obs * prev) %% 1) == 0)
  n_pos <- n_obs * prev
  n_neg <- n_obs - n_pos
  
  pos <- rskew(n_pos, mu = 1, sdev = sdev_pos, skew = -skew)
  neg <- rskew(n_neg, mu = 0, sdev = sdev_neg, skew = skew)
  
  .raw_to_prob <- function(dat, prev, sdev_pos, sdev_neg, skew) {
    p1 <- prev * dskew(dat, mu = 1, sdev = sdev_pos, skew = -skew)
    p2 <- (1 - prev) * dskew(dat, mu = 0, sdev = sdev_neg, skew = skew)
    return(p1 / (p1 + p2))
  }
  
  prob_pos <- .raw_to_prob(pos, prev_train, sdev_pos, sdev_neg, skew)
  prob_neg <- .raw_to_prob(neg, prev_train, sdev_pos, sdev_neg, skew)
  
  out <- list(raw = sort(c(pos, neg), decreasing = TRUE),
              prob = sort(c(prob_pos, prob_neg), decreasing = FALSE))
  return(out)
}


simulate_train_skew <- function(n_obs, sdev_pos, sdev_neg, prev, skew) {
  stopifnot("n_obs * prev should be an integer" = ((n_obs * prev) %% 1) == 0)
  n_pos <- n_obs * prev
  n_neg <- n_obs - n_pos
  
  data_pos <- rskew(n_pos, mu = 1, sdev = sdev_pos, skew = -skew)
  data_neg <- rskew(n_neg, mu = 0, sdev = sdev_neg, skew = skew)
  
  mean_pos <- mean(data_pos)
  mean_neg <- mean(data_neg)
  
  sd_pos <- sd(data_pos)
  sd_neg <- sd(data_neg)
  
  skew_pos <- sn::msn.mle(y = data_pos)$dp$alpha
  skew_neg <- sn::msn.mle(y = data_neg)$dp$alpha
  
  cs_pos <- c("mean" = mean_pos, "sd" = sd_pos, "skew" = skew_pos)
  cs_neg <- c("mean" = mean_neg, "sd" = sd_neg, "skew" = skew_neg)
  
  ecdf_pos <- ecdf(data_pos)
  ecdf_neg <- ecdf(data_neg)
  
  ms_pos <- function(x) {1 - ecdf_pos(x)}
  ms_neg <- function(x) {1 - ecdf_neg(x)}
  
  .raw_to_prob <- function(dat, prev, sdev_pos, sdev_neg, skew) {
    p1 <- prev * dskew(dat, mu = 1, sdev = sdev_pos, skew = -skew)
    p2 <- (1 - prev) * dskew(dat, mu = 0, sdev = sdev_neg, skew = skew)
    return(p1 / (p1 + p2))
  }
  
  prob_pos <- .raw_to_prob(data_pos, prev, sdev_pos, sdev_neg, skew) |> sort()
  prob_neg <- .raw_to_prob(data_neg, prev, sdev_pos, sdev_neg, skew) |> sort()
  
  return(list(cs_pos = cs_pos, 
              cs_neg = cs_neg,
              ms_pos = ms_pos,
              ms_neg = ms_neg,
              prob_pos = prob_pos,
              prob_neg = prob_neg,
              prev_train = prev,
              rng_raw = range(c(data_pos, data_neg))))
}

get_prevalence_estimates <- function(train, test) {
  
  # helper functions for continuous sweep and median sweep
  .adjusted_count_ms <- function(cac, tpr, fpr) {
    out <- (cac - fpr) / (tpr - fpr)
    return(out)
  }
  
  .adjusted_count_cs <- function(cac, theta, mean_pos, mean_neg, sdev_pos, sdev_neg) {
    .tpr <- pnorm(theta, mean = mean_pos, sd = sdev_pos, lower.tail = FALSE)
    .fpr <- pnorm(theta, mean = mean_neg, sd = sdev_neg, lower.tail = FALSE)
    return((cac - .fpr) / (.tpr - .fpr))
  }
  
  .adjusted_count_cs_skew <- function(cac, theta, mean_pos, mean_neg, sdev_pos, sdev_neg,
                                      skew_pos, skew_neg) {
    .tpr <- pskew(theta, mu = mean_pos, sdev = sdev_pos, skew = skew_pos)
    .fpr <- pskew(theta, mu = mean_neg, sdev = sdev_neg, skew = skew_neg)
    return((cac - .fpr) / (.tpr - .fpr))
  }
  
  .calculate_one_area <- function(cac, theta_l, theta_r, mean_pos, mean_neg, 
                                  sdev_pos, sdev_neg) {
    
    out <- integrate(f = function(x) {
      .adjusted_count_cs(cac, x, mean_pos, mean_neg, sdev_pos, sdev_neg)
    }, lower = theta_l, upper = theta_r)$value
    
    return(out)
  }
  
  .calculate_one_area_skew <- function(cac, theta_l, theta_r, mean_pos, mean_neg, 
                                       sdev_pos, sdev_neg, skew_pos, skew_neg) {
    
    out <- integrate(f = function(x) {
      .adjusted_count_cs_skew(cac, x, mean_pos, mean_neg, sdev_pos, sdev_neg,
                              skew_pos, skew_neg)
    }, lower = theta_l, upper = theta_r)$value
    
    return(out)
  }
  
  .get_threshold_barriers <- function(mean_pos, mean_neg, sdev_pos, sdev_neg, pdelta, train) {
    tpr <- function(theta) {pnorm(theta, mean = mean_pos, sd = sdev_pos, lower.tail = F)}
    fpr <- function(theta) {pnorm(theta, mean = mean_neg, sd = sdev_neg, lower.tail = F)}
    diff_rates <- function(theta) {.tpr(theta) - .fpr(theta)}
    
    barriers <-
      rootSolve::uniroot.all(f = function(x) {
        out <- tpr(x) - fpr(x) - pdelta
        return(out)
      }, interval = train$rng_raw)
    
    if(length(barriers) == 0) {
      lft_range <- tpr(train$rng_raw[1]) - fpr(train$rng_raw[1])
      rgt_range <- tpr(train$rng_raw[2]) - fpr(train$rng_raw[2])
      if(lft_range >= pdelta && rgt_range >= pdelta) {
        barriers <- rng_raw
      }
      else {
        return(NA)
      }
    }
    
    if(length(barriers) %% 2 == 1) {
      lft <- train$rng_raw[1]
      rgt <- train$rng_raw[2]
      # append it with otherwise the left or right boundary 
      # (whichever is suitable)
      if(tpr(lft) - fpr(lft) > pdelta) {
        barriers <- c(lft, barriers)
      }
      else if(tpr(rgt) - fpr(rgt) > pdelta){
        barriers <- c(barriers, rgt)
      }
    }
    return(barriers)
  }
  
  .get_threshold_barriers_skew <- function(mean_pos, mean_neg, sdev_pos, sdev_neg, 
                                           skew_pos, skew_neg, pdelta, train) {
    tpr <- function(theta) {pskew(theta, mu = mean_pos, sdev = sdev_pos, skew = skew_pos)}
    fpr <- function(theta) {pskew(theta, mu = mean_neg, sdev = sdev_neg, skew = skew_neg)}
    diff_rates <- function(theta) {.tpr(theta) - .fpr(theta)}
    
    barriers <-
      rootSolve::uniroot.all(f = function(x) {
        out <- tpr(x) - fpr(x) - pdelta
        return(out)
      }, interval = train$rng_raw)
    
    if(length(barriers) == 0) {
      lft_range <- tpr(train$rng_raw[1]) - fpr(train$rng_raw[1])
      rgt_range <- tpr(train$rng_raw[2]) - fpr(train$rng_raw[2])
      if(lft_range >= pdelta && rgt_range >= pdelta) {
        barriers <- rng_raw
      }
      else {
        return(NA)
      }
    }
    
    if(length(barriers) %% 2 == 1) {
      lft <- train$rng_raw[1]
      rgt <- train$rng_raw[2]
      # append it with otherwise the left or right boundary 
      # (whichever is suitable)
      if(tpr(lft) - fpr(lft) > pdelta) {
        barriers <- c(lft, barriers)
      }
      else if(tpr(rgt) - fpr(rgt) > pdelta){
        barriers <- c(barriers, rgt)
      }
    }
    return(barriers)
  }
  
  .adjusted_count_cov <- function(theta_x, theta_y, n_obs, mean_pos, mean_neg, sdev_pos, 
                                  sdev_neg) {
    
    .tpr <- function(theta) {pnorm(theta, mean = mean_pos, sd = sdev_pos, lower.tail = F)}
    .fpr <- function(theta) {pnorm(theta, mean = mean_neg, sd = sdev_neg, lower.tail = F)}
    .diff_rates <- function(theta) {.tpr(theta) - .fpr(theta)}
    
    prev <- 0.5
    n0 <- n_obs * (1 - prev)
    n1 <- n_obs * prev
    
    pt1 <- .tpr(pmax(theta_x, theta_y)) - .tpr(theta_x) * .tpr(theta_y)
    pt2 <- .fpr(pmax(theta_x, theta_y)) - .fpr(theta_x) * .fpr(theta_y)
    
    cov_th <- n1 / n_obs^2 * pt1 + n0 / n_obs^2 * pt2
    
    diff_rates_x <- .diff_rates(theta_x)
    diff_rates_y <- .diff_rates(theta_y)
    
    out <- cov_th / (diff_rates_x * diff_rates_y) 
    return(out)
  }
  
  .adjusted_count_cov_skew <- function(theta_x, theta_y, n_obs, mean_pos, mean_neg, 
                                       sdev_pos, sdev_neg, skew_pos, skew_neg) {
    
    .tpr <- function(theta) {pskew(theta, mu = mean_pos, sdev = sdev_pos, skew = skew_pos)}
    .fpr <- function(theta) {pskew(theta, mu = mean_neg, sdev = sdev_neg, skew = skew_neg)}
    .diff_rates <- function(theta) {.tpr(theta) - .fpr(theta)}
    
    prev <- 0.5
    n0 <- n_obs * (1 - prev)
    n1 <- n_obs * prev
    
    pt1 <- .tpr(pmax(theta_x, theta_y)) - .tpr(theta_x) * .tpr(theta_y)
    pt2 <- .fpr(pmax(theta_x, theta_y)) - .fpr(theta_x) * .fpr(theta_y)
    
    cov_th <- n1 / n_obs^2 * pt1 + n0 / n_obs^2 * pt2
    
    diff_rates_x <- .diff_rates(theta_x)
    diff_rates_y <- .diff_rates(theta_y)
    
    out <- cov_th / (diff_rates_x * diff_rates_y) 
    return(out)
  }
  
  .continuous_sweep_int <- function(n_obs, mean_pos, mean_neg, sdev_pos, 
                                    sdev_neg, trunc_value = NULL) {
    # obtain decision boundaries
    barrier_levels <- .get_threshold_barriers(mean_pos, mean_neg, sdev_pos, 
                                              sdev_neg, trunc_value)
    # compute integral
    intgrl <- 
      pracma::integral2(fun = function(xi, yi) {
        .adjusted_count_cov(xi, yi, n_obs, mean_pos, mean_neg, sdev_pos, sdev_neg)
      }, barrier_levels[1], barrier_levels[2], barrier_levels[1], barrier_levels[2])$Q 
    
    scaling_factor <- diff(barrier_levels)^2
    out <- intgrl / scaling_factor
    
    return(out)
  }
  
  .continuous_sweep_int_skew <- function(n_obs, mean_pos, mean_neg, sdev_pos, sdev_neg, 
                                         skew_pos, skew_neg, trunc_value = NULL) {
    # obtain decision boundaries
    barrier_levels <- .get_threshold_barriers_skew(mean_pos, mean_neg, sdev_pos, sdev_neg, 
                                                   skew_pos, skew_neg, trunc_value)
    # compute integral
    intgrl <- 
      pracma::integral2(fun = function(xi, yi) {
        .adjusted_count_cov_skew(xi, yi, n_obs, mean_pos, mean_neg, sdev_pos, 
                                 sdev_neg, skew_pos, skew_neg)
      }, barrier_levels[1], barrier_levels[2], barrier_levels[1], barrier_levels[2])$Q 
    
    scaling_factor <- diff(barrier_levels)^2
    out <- intgrl / scaling_factor
    
    return(out)
  }
  

  .get_optim_threshold <- function(mean_pos, mean_neg, sdev_pos, sdev_neg, train) {
    tpr <- function(theta) {pnorm(theta, mean = mean_pos, sd = sdev_pos, lower.tail = F)}
    fpr <- function(theta) {pnorm(theta, mean = mean_neg, sd = sdev_neg, lower.tail = F)}
    diff_rates <- function(theta) {tpr(theta) - fpr(theta)}
    rng <- train$rng
    
    # avoid boundaries
    err <- 10^(-5)
    
    max_pdelta <- 
      optimize(f = function(x) {tpr(x) - fpr(x)}, 
               interval = rng, maximum = TRUE) |>
      purrr::pluck("objective")
    
    min_pdelta <- 
      pmax(tpr(rng[1]) - fpr(rng[1]),
           tpr(rng[2]) - fpr(rng[2]))
    
    .estimate_mse <- function(tpr, fpr, rng, pdelta) {
      # compute covariance at point x, y
      .cov_ac <- function(x_coor, y_coor, rates_cs) {
        prev <- 0.5
        diff_rates <- function(x) {tpr(x) - fpr(x)}
        
        pt1 <- purrr::map2_dbl(x_coor, y_coor, ~ tpr(max(.x, .y)) - tpr(.x) * tpr(.y))
        pt2 <- purrr::map2_dbl(x_coor, y_coor, ~ fpr(max(.x, .y)) - fpr(.x) * fpr(.y))
        
        cov_th <-  prev * (pt1 + pt2)
        
        diff_rates_x <- purrr::map_dbl(x_coor, ~ diff_rates(.x))
        diff_rates_y <- purrr::map_dbl(y_coor, ~ diff_rates(.x))
        
        out <- cov_th / (diff_rates_x * diff_rates_y) 
        return(out)
      }
      # find the roots where F^+ - F^- = pdelta
      roots_pdelta <- 
        rootSolve::uniroot.all(f = function(x) {
          out <- tpr(x) - fpr(x) - pdelta
          return(out)
        }, interval = rng)
      
      # check whether length roots pdelta is not zero, otherwise infinite variance
      if(length(roots_pdelta) == 0) {
        return(Inf)
      }
      
      # check whether roots pdelta is an odd number
      if(length(roots_pdelta) %% 2 == 1) {
        # append it with otherwise the left or right boundary (whichever is suitable)
        if(tpr(rng[1]) - fpr(rng[1]) > pdelta) {
          roots_pdelta <- c(rng[1], roots_pdelta)
        }
        else if(tpr(rng[2]) - fpr(rng[2]) > pdelta){
          roots_pdelta <- c(roots_pdelta, rng[2])
        }
      }
      # sum the integrals (also works if there is one interval)
      cov_auc <- numeric(length(roots_pdelta) / 2)
      idx_sum <- 1
      for(idx in seq(1, length(roots_pdelta), 2)) {
        theta_l <- roots_pdelta[idx]
        theta_r <- roots_pdelta[idx + 1]
        if(theta_l == theta_r) {
          cov_auc[idx_sum] <- .cov_ac(theta_l, theta_r, rates_cs)
        }
        else{
          cov_auc[idx_sum] <- pracma::integral2(f = function(x, y) {
            .cov_ac(x, y, rates_cs)
          }, theta_l, theta_r, theta_l, theta_r)$Q / (theta_r - theta_l)^2
        }
        idx_sum <- idx_sum + 1
      }
      
      return(sum(cov_auc))
    }
    
    # need estimate of mse
    optim_pdelta <- 
      optimise(f = function(x) {.estimate_mse(tpr, fpr, rng, x)},
               interval = c(min_pdelta + err, max_pdelta - err)) |>
      purrr::pluck("minimum")
    
    return(optim_pdelta)
  }
  
  .get_optim_threshold_skew <- function(mean_pos, mean_neg, sdev_pos, sdev_neg, skew_pos,
                                        skew_neg, train) {
    tpr <- function(theta) {pskew(theta, mu = mean_pos, sdev = sdev_pos, skew = skew_pos)}
    fpr <- function(theta) {pskew(theta, mu = mean_neg, sdev = sdev_neg, skew = skew_neg)}
    diff_rates <- function(theta) {tpr(theta) - fpr(theta)}
    rng <- train$rng
    
    # avoid boundaries
    err <- 10^(-5)
    
    max_pdelta <- 
      optimize(f = function(x) {tpr(x) - fpr(x)}, 
               interval = rng, maximum = TRUE) |>
      purrr::pluck("objective")
    
    min_pdelta <- 
      pmax(tpr(rng[1]) - fpr(rng[1]),
           tpr(rng[2]) - fpr(rng[2]))
    
    .estimate_mse <- function(tpr, fpr, rng, pdelta) {
      # compute covariance at point x, y
      .cov_ac <- function(x_coor, y_coor, rates_cs) {
        prev <- 0.5
        diff_rates <- function(x) {tpr(x) - fpr(x)}
        
        pt1 <- purrr::map2_dbl(x_coor, y_coor, ~ tpr(max(.x, .y)) - tpr(.x) * tpr(.y))
        pt2 <- purrr::map2_dbl(x_coor, y_coor, ~ fpr(max(.x, .y)) - fpr(.x) * fpr(.y))
        
        cov_th <-  prev * (pt1 + pt2)
        
        diff_rates_x <- purrr::map_dbl(x_coor, ~ diff_rates(.x))
        diff_rates_y <- purrr::map_dbl(y_coor, ~ diff_rates(.x))
        
        out <- cov_th / (diff_rates_x * diff_rates_y) 
        return(out)
      }
      # find the roots where F^+ - F^- = pdelta
      roots_pdelta <- 
        rootSolve::uniroot.all(f = function(x) {
          out <- tpr(x) - fpr(x) - pdelta
          return(out)
        }, interval = rng)
      
      # check whether length roots pdelta is not zero, otherwise infinite variance
      if(length(roots_pdelta) == 0) {
        return(Inf)
      }
      
      # check whether roots pdelta is an odd number
      if(length(roots_pdelta) %% 2 == 1) {
        # append it with otherwise the left or right boundary (whichever is suitable)
        if(tpr(rng[1]) - fpr(rng[1]) > pdelta) {
          roots_pdelta <- c(rng[1], roots_pdelta)
        }
        else if(tpr(rng[2]) - fpr(rng[2]) > pdelta){
          roots_pdelta <- c(roots_pdelta, rng[2])
        }
      }
      # sum the integrals (also works if there is one interval)
      cov_auc <- numeric(length(roots_pdelta) / 2)
      idx_sum <- 1
      for(idx in seq(1, length(roots_pdelta), 2)) {
        theta_l <- roots_pdelta[idx]
        theta_r <- roots_pdelta[idx + 1]
        if(theta_l == theta_r) {
          cov_auc[idx_sum] <- .cov_ac(theta_l, theta_r, rates_cs)
        }
        else{
          cov_auc[idx_sum] <- pracma::integral2(f = function(x, y) {
            .cov_ac(x, y, rates_cs)
          }, theta_l, theta_r, theta_l, theta_r)$Q / (theta_r - theta_l)^2
        }
        idx_sum <- idx_sum + 1
      }
      
      return(sum(cov_auc))
    }
    
    # need estimate of mse
    optim_pdelta <- 
      optimise(f = function(x) {.estimate_mse(tpr, fpr, rng, x)},
               interval = c(min_pdelta + err, max_pdelta - err)) |>
      purrr::pluck("minimum")
    
    return(optim_pdelta)
  }
  
  # functions to compute continuous sweep and median sweep
  .get_median_sweep <- function(train, test, trunc_value) {
    pos <- train %>% pluck("ms_pos")
    neg <- train %>% pluck("ms_neg")
    
    out <- 
      tibble(tpr = pos(test),
             fpr = neg(test),
             cac = (1 : length(test)) / length(test)) %>%
      mutate(diff = tpr - fpr) %>%
      filter(diff >= trunc_value) %>%
      mutate(ac = pmap_dbl(.l = list(cac, tpr, fpr), 
                           .f = ~ .adjusted_count_ms(..1, ..2, ..3))) %>%
      pull(ac) %>%
      median()
    return(out)
  }
  
  .get_continuous_sweep <- function(train, test, trunc_value) {
    # set variables
    n_obs <- length(test)
    mean_pos <- train$cs_pos[1]
    sdev_pos <- train$cs_pos[2]
    mean_neg <- train$cs_neg[1]
    sdev_neg <- train$cs_neg[2]
    
    # make 3 vectors containing the left bound, right bound, and cac-value
    all_theta_r <- test %>% head(-1)
    all_theta_l <- test %>% tail(-1)
    all_cac <- seq(1/n_obs, 1 - 1/n_obs, 1/n_obs) 
    
    # filter out the intervals that fall outside the "classification rate" `trunc_value`
    barriers <- .get_threshold_barriers(mean_pos, mean_neg, 
                                        sdev_pos, sdev_neg, 
                                        trunc_value, train)
    diff_range <- barriers %>% diff()
    bool_crit <- all_theta_r >= barriers[1] & all_theta_l <= barriers[2]
    n_trunc <- sum(bool_crit)
    
    # replace the values slightly below and sightly about barriers with correct values.
    all_theta_r <- subset(all_theta_r, bool_crit) %>% replace(1, barriers[2])
    all_theta_l <- subset(all_theta_l, bool_crit) %>% replace(n_trunc, barriers[1])
    all_cac <- subset(all_cac, bool_crit)
    
    all_areas <- 
      pmap_dbl(.l = list(all_cac, all_theta_l, all_theta_r), 
               .f =  ~ .calculate_one_area(..1, ..2, ..3, 
                                           mean_pos, mean_neg, sdev_pos, sdev_neg))
    
    area_sum <- sum(all_areas)
    return(area_sum / diff_range)
  }
  
  .get_continuous_sweep_skew <- function(train, test, trunc_value) {
    # set variables
    n_obs <- length(test)
    mean_pos <- train$cs_pos[1]
    sdev_pos <- train$cs_pos[2]
    skew_pos <- train$cs_pos[3]
    mean_neg <- train$cs_neg[1]
    sdev_neg <- train$cs_neg[2]
    skew_neg <- train$cs_neg[3]
    
    # make 3 vectors containing the left bound, right bound, and cac-value
    all_theta_r <- test %>% head(-1)
    all_theta_l <- test %>% tail(-1)
    all_cac <- seq(1/n_obs, 1 - 1/n_obs, 1/n_obs) 
    
    # filter out the intervals that fall outside the "classification rate" `trunc_value`
    barriers <- .get_threshold_barriers_skew(mean_pos, mean_neg, 
                                             sdev_pos, sdev_neg,
                                             skew_pos, skew_neg, 
                                             trunc_value, train)
    diff_range <- barriers %>% diff()
    bool_crit <- all_theta_r >= barriers[1] & all_theta_l <= barriers[2]
    n_trunc <- sum(bool_crit)
    
    # replace the values slightly below and sightly about barriers with correct values.
    all_theta_r <- subset(all_theta_r, bool_crit) %>% replace(1, barriers[2])
    all_theta_l <- subset(all_theta_l, bool_crit) %>% replace(n_trunc, barriers[1])
    all_cac <- subset(all_cac, bool_crit)
    
    all_areas <- 
      pmap_dbl(.l = list(all_cac, all_theta_l, all_theta_r), 
               .f =  ~ .calculate_one_area_skew(..1, ..2, ..3, 
                                                mean_pos, mean_neg, sdev_pos, sdev_neg,
                                                skew_pos, skew_neg))
    
    area_sum <- sum(all_areas)
    return(area_sum / diff_range)
  }
  
  .get_SLD <- function(train, test_scores) {
    # compute scores
    alpha_train <- c(1 - train$prev_train, train$prev_train)
    
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
    
    test_predict_SLD <- SLD_alg(test_scores, alpha_train) %>% dplyr::last()
    
    return(test_predict_SLD)
  }
  
  .get_DyS <- function(train, test_scores, n_bins = 8, tol = 10^(-5)) {
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
    
    # obtain train and test scores (prob)
    
    # compute histogram of positive and negative scores
    h0 <- .compute_histogram(train$prob_neg, n_bins, rng)
    h1 <- .compute_histogram(train$prob_pos, n_bins, rng)
    
    # minimize distance mixture and test histograms with DyS
    out <- .compute_dys(test_scores, h0, h1, n_bins, tol, rng, mthd)
    return(out)
  }
  
  
  optim_thres <- .get_optim_threshold(train$cs_pos[1], train$cs_neg[1], 
                                      train$cs_pos[2], train$cs_neg[2], 
                                      train)
  
  optim_thres_skew <- .get_optim_threshold_skew(train$cs_pos[1], train$cs_neg[1], 
                                                train$cs_pos[2], train$cs_neg[2], 
                                                train$cs_pos[3], train$cs_neg[3], 
                                                train)
  
  out <- list(cs_trad = .get_continuous_sweep(train, test$raw, 0.25),
              cs_trad_skew = .get_continuous_sweep_skew(train, test$raw, 0.25),
              ms_trad = .get_median_sweep(train, test$raw, 0.25),
              cs_opt = .get_continuous_sweep(train, test$raw, optim_thres),
              cs_opt_skew = .get_continuous_sweep_skew(train, test$raw, optim_thres_skew),
              ms_opt = .get_median_sweep(train, test$raw, optim_thres),
              dys = .get_DyS(train, test$prob),
              sld = .get_SLD(train, test$prob),
              optim_thres = optim_thres,
              optim_thres_skew = optim_thres_skew)
  return(out)
}


one_sim_function <- function(n_train, n_test, sdev_pos, sdev_neg,
                             prev_train, prev_test, skew, reps) {
  
  out <- map(1:reps, function(x) {
    train <- simulate_train_skew(n_train, sdev_pos, sdev_neg, prev_train, skew)
    test <- simulate_test_skew(n_test, sdev_pos, sdev_neg, prev_test, skew, prev_train)
    dfr_prev <- get_prevalence_estimates(train, test)
    return(dfr_prev)
  })
  
  return(out)
}
