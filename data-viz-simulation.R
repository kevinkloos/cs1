library(tidyverse)
library(gridExtra)

########################################################################################
# load the data first!!
# boxplots MSE sim1

p1 <- 
  sim1_clean |>
  group_by(quantifier, n_obs, sdev_pos, sdev_neg, prev_test) |>
  reframe(rmse = mean(sqrt(err^2), na.rm = T)) |>
  ungroup() |>
  bind_rows(sim1_addon) |>
  mutate(clr = case_when(quantifier == "O-CS" ~ 'A', 
                         quantifier == "T-MS" ~ 'C',
                         quantifier %in% c("DyS", "SLD") ~ "D",
                         .default = "B")) |>
  ggplot(aes(x = fct_relevel(quantifier, c("O-CS", "T-CS", "O-MS", "T-MS", "AC", "DyS", "SLD")),
             y = rmse,
             fill = clr)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "quantifier", y = "RMSE") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
  scale_fill_manual(values = c("A" = "#FFFFFF", "B" = "#DED4D4", 
                               "C" = "#999090", "D" = "#5D5555")) +
  theme_bw()

p2 <- 
  sim1_clean |>
  group_by(quantifier, n_obs, sdev_pos, sdev_neg, prev_test) |>
  reframe(rmse = mean(sqrt(err^2), na.rm = T)) |>
  bind_rows(sim1_addon) |>
  pivot_wider(names_from = quantifier, values_from = rmse) |>
  mutate(across("T-CS": "AC", ~ .x - `O-CS`)) |>
  select(-`O-CS`) |>
  pivot_longer(cols = "T-CS": "AC", names_to = "quantifier", values_to = "rmse_diff") |>
  mutate(clr = case_when(quantifier %in% c("DyS", "SLD") ~ "D",
                         quantifier == "T-MS" ~ 'C',
                         .default = "B")) |>
  ggplot(aes(x = rmse_diff, fill = clr)) +
  geom_histogram(boundary = 0, binwidth = 0.005, show.legend = FALSE,color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "orange", linewidth = 1) +
  facet_wrap(~ fct_relevel(quantifier, c("T-CS", "O-MS", "T-MS", "AC", "DyS", "SLD"))) +
  scale_fill_manual(values = c("B" = "#DED4D4","C" = "#999090", "D" = "#5D5555")) +
  labs(x = "difference in RMSE", y = NULL) +
  scale_x_continuous(limits = c(-0.05, 0.10)) +
  theme_bw() 


p3 <- grid.arrange(p1, p2, nrow = 1, widths = c(3, 4), heights = 3)
ggsave("_img/sim1/s1_box_hist.pdf", p3, width = 7, height = 4)

rm(p1, p2)
#########################################################################################
p1 <- 
  sim2_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(rmse = mean(sqrt(err^2), na.rm = T)) |>
  ungroup() |>
  bind_rows(sim2_addon) |>
  mutate(clr = case_when(quantifier == "O-CS" ~ 'A', 
                         quantifier == "T-MS" ~ 'C',
                         quantifier %in% c("DyS", "SLD") ~ "D",
                         .default = "B")) |>
  ggplot(aes(x = fct_relevel(quantifier, c("O-CS", "T-CS", "O-MS", "T-MS", "AC", "DyS", "SLD")),
             y = rmse,
             fill = clr)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "quantifier", y = "RMSE") +
  scale_y_continuous(limits = c(0, 0.25), breaks = seq(0, 0.25, 0.05)) +
  scale_fill_manual(values = c("A" = "#FFFFFF", "B" = "#DED4D4",
                               "C" = "#999090", "D" = "#5D5555")) +
  theme_bw()


p2 <- 
  sim2_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(rmse = mean(sqrt(err^2), na.rm = T)) |>
  bind_rows(sim2_addon) |>
  pivot_wider(names_from = quantifier, values_from = rmse) |>
  mutate(across("T-CS": "AC", ~ .x - `O-CS`)) |> 
  select(-`O-CS`) |>
  pivot_longer(cols = "T-CS": "AC", names_to = "quantifier", values_to = "rmse_diff") |>
  mutate(clr = case_when(quantifier %in% c("DyS", "SLD") ~ "D",
                         quantifier == "T-MS" ~ 'C',
                         .default = "B")) |>
  ggplot(aes(x = rmse_diff, fill = clr)) +
  geom_histogram(boundary = 0, binwidth = 0.005, show.legend = FALSE,color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "orange", linewidth = 1) +
  facet_wrap(~ fct_relevel(quantifier, c("T-CS", "O-MS", "T-MS", 'AC', "DyS", "SLD"))) +
  scale_fill_manual(values = c("B" = "#DED4D4","C" = "#999090", "D" = "#5D5555")) +
  labs(x = "difference in RMSE", y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2))
  

p3 <- gridExtra::grid.arrange(p1, p2, nrow = 1, widths = c(3, 4), heights = 4)
ggsave("_img/sim2/s2_box_hist.pdf", p3, width = 7, height = 4)

rm(p1, p2)
#######################################################################################
p1 <- 
  sim3_clean |>
  filter(skew == 1) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(rmse = mean(sqrt(err^2), na.rm = T)) |>
  bind_rows(sim3_addon |> filter(skew == 1) |> select(-skew)) |>
  ungroup() |>
  mutate(clr = case_when(quantifier == "O-CS" ~ 'A', 
                         quantifier == "T-MS" ~ 'C',
                         quantifier %in% c("DyS", "SLD") ~ "D",
                         .default = "B")) |>
  ggplot(aes(x = fct_relevel(quantifier, c("O-CS (skew)", "T-CS (skew)", 
                                           "O-CS (norm)", "T-CS (norm)",
                                           "O-MS", "T-MS", "AC", "DyS", "SLD")),
             y = rmse,
             fill = clr)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "quantifier", y = "RMSE") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
  scale_fill_manual(values = c("A" = "#FFFFFF", "B" = "#DED4D4",
                               "C" = "#999090", "D" = "#5D5555")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2)) 

p2 <- 
  sim3_clean |>
  filter(skew == 1) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(rmse = mean(sqrt(err^2), na.rm = T)) |>
  bind_rows(sim3_addon |> filter(skew == 1) |> select(-skew)) |>
  ungroup() |>
  pivot_wider(names_from = quantifier, values_from = rmse) |>
  select(n_train:prev_test, `O-CS (skew)`, `T-CS (skew)`, everything()) |>
  mutate(across("T-CS (skew)": "AC", ~ .x - `O-CS (skew)`)) |>
  select(-`O-CS (skew)`) |>
  pivot_longer(cols = "T-CS (skew)": "AC", names_to = "quantifier", values_to = "rmse_diff") |>
  mutate(clr = case_when(quantifier %in% c("DyS", "SLD") ~ "D",
                         quantifier == "T-MS" ~ 'C',
                         .default = "B")) |>
  ggplot(aes(x = rmse_diff, fill = clr)) +
  geom_histogram(boundary = 0, binwidth = 0.005, show.legend = FALSE,color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "orange", linewidth = 1) +
  facet_wrap(~ fct_relevel(quantifier, c("T-CS (skew)", "O-CS (norm)", "T-CS (norm)", "O-MS", "T-MS", 
                                         "AC", "DyS", "SLD")),
             ncol = 5) +
  scale_fill_manual(values = c("B" = "#DED4D4","C" = "#999090", "D" = "#5D5555")) +
  scale_x_continuous(limits = c(-0.06, 0.06)) +
  labs(x = "difference in RMSE", y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2)) +
  theme(axis.text = element_text(size = 7),
        strip.text = element_text(size = 7)) 


p3 <- gridExtra::grid.arrange(p1, p2, nrow = 1, widths = c(3, 4), heights = 3.5)
ggsave("_img/sim3/box_hist_skew1.pdf", p3, width = 7, height = 3)


p1 <- 
  sim3_clean |>
  filter(skew == 2) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(rmse = mean(sqrt(err^2), na.rm = T)) |>
  ungroup() |>
  bind_rows(sim3_addon |> filter(skew == 2) |> select(-skew)) |>
  mutate(clr = case_when(quantifier == "O-CS" ~ 'A', 
                         quantifier == "T-MS" ~ 'C',
                         quantifier %in% c("DyS", "SLD") ~ "D",
                         .default = "B")) |>
  ggplot(aes(x = fct_relevel(quantifier, c("O-CS (skew)", "T-CS (skew)", 
                                           "O-CS (norm)", "T-CS (norm)",
                                           "O-MS", "T-MS", "AC", "DyS", "SLD")),
             y = rmse,
             fill = clr)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "quantifier", y = "RMSE") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
  scale_fill_manual(values = c("A" = "#FFFFFF", "B" = "#DED4D4",
                               "C" = "#999090", "D" = "#5D5555")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2))

p2 <- 
  sim3_clean |>
  filter(skew == 2) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(rmse = mean(sqrt(err^2), na.rm = T)) |>
  bind_rows(sim3_addon |> filter(skew == 2) |> select(-skew)) |>
  ungroup() |>
  pivot_wider(names_from = quantifier, values_from = rmse) |>
  select(n_train:prev_test, `O-CS (skew)`, `T-CS (skew)`, everything()) |>
  mutate(across("T-CS (skew)": "AC", ~ .x - `O-CS (skew)`)) |>
  select(-`O-CS (skew)`) |>
  pivot_longer(cols = "T-CS (skew)": "AC", names_to = "quantifier", values_to = "rmse_diff") |>
  mutate(clr = case_when(quantifier %in% c("DyS", "SLD") ~ "D",
                         quantifier == "T-MS" ~ 'C',
                         .default = "B")) |>
  ggplot(aes(x = rmse_diff, fill = clr)) +
  geom_histogram(boundary = 0, binwidth = 0.005, show.legend = FALSE,color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "orange", linewidth = 1) +
  facet_wrap(~ fct_relevel(quantifier, c("T-CS (skew)", "O-CS (norm)", "T-CS (norm)", "O-MS", "T-MS", 
                                         "AC", "DyS", "SLD")),
             ncol = 5) +
  scale_fill_manual(values = c("B" = "#DED4D4", "C" = "#999090", "D" = "#5D5555")) +
  scale_x_continuous(limits = c(-0.06, 0.06)) +
  labs(x = "difference in RMSE", y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2)) +
  theme(axis.text = element_text(size = 7),
        strip.text = element_text(size = 7)) 

p3 <- gridExtra::grid.arrange(p1, p2, nrow = 1, widths = c(3, 4), heights = 3.5)
ggsave("_img/sim3/box_hist_skew2.pdf", p3, width = 7, height = 3)

p1 <- 
  sim3_clean |>
  filter(skew == 4) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(rmse = mean(sqrt(err^2), na.rm = T)) |>
  bind_rows(sim3_addon |> filter(skew == 4) |> select(-skew)) |>
  ungroup() |>
  mutate(clr = case_when(quantifier == "O-CS" ~ 'A', 
                         quantifier == "T-MS" ~ 'C',
                         quantifier %in% c("DyS", "SLD") ~ "D",
                         .default = "B")) |>
  ggplot(aes(x = fct_relevel(quantifier, c("O-CS (skew)", "T-CS (skew)", 
                                           "O-CS (norm)", "T-CS (norm)",
                                           "O-MS", "T-MS", "AC", "DyS", "SLD")),
             y = rmse,
             fill = clr)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "quantifier", y = "RMSE") +
  scale_y_continuous(limits = c(0, 0.20), breaks = seq(0, 0.20, 0.05)) +
    scale_fill_manual(values = c("A" = "#FFFFFF", "B" = "#DED4D4",
                                 "C" = "#999090", "D" = "#5D5555")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2))

p2 <- 
  sim3_clean |>
  filter(skew == 4) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(rmse = mean(sqrt(err^2), na.rm = T)) |>
  bind_rows(sim3_addon |> filter(skew == 4) |> select(-skew)) |>
  ungroup() |>
  pivot_wider(names_from = quantifier, values_from = rmse) |>
  select(n_train:prev_test, `O-CS (skew)`, `T-CS (skew)`, everything()) |>
  mutate(across("T-CS (skew)": "AC", ~ .x - `O-CS (skew)`)) |>
  select(-`O-CS (skew)`) |>
  pivot_longer(cols = "T-CS (skew)": "AC", names_to = "quantifier", values_to = "rmse_diff") |>
  mutate(clr = case_when(quantifier %in% c("DyS", "SLD") ~ "D",
                         quantifier == "T-MS" ~ 'C',
                         .default = "B")) |>
  ggplot(aes(x = rmse_diff, fill = clr)) +
  geom_histogram(boundary = 0, binwidth = 0.005, show.legend = FALSE,color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "orange", linewidth = 1) +
  facet_wrap(~ fct_relevel(quantifier, c("T-CS (skew)", "O-CS (norm)", "T-CS (norm)", 
                                         "O-MS", "T-MS", "AC", "DyS", "SLD")),
             ncol = 5) +
    scale_fill_manual(values = c("B" = "#DED4D4", "C" = "#999090", "D" = "#5D5555")) +
  scale_x_continuous(limits = c(-0.09, 0.09)) +
  labs(x = "difference in RMSE", y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2)) +
  theme(axis.text = element_text(size = 7),
        strip.text = element_text(size = 7)) 

p3 <- gridExtra::grid.arrange(p1, p2, nrow = 1, widths = c(3, 4), heights = 3)
ggsave("_img/sim3/box_hist_skew4.pdf", p3, width = 7, height = 3)


################################

sim1_aggregated <- 
  sim1_clean |>
  group_by(quantifier, n_obs, sdev_pos, sdev_neg, prev_test) |>
  reframe(rmse = mean(sqrt(err^2), na.rm = T)) |>
  bind_rows(sim1_addon)

sim2_aggregated <- 
  sim2_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(rmse = mean(sqrt(err^2), na.rm = T)) |>
  bind_rows(sim2_addon) 

sim1_clean |>
  group_by(quantifier, n_obs, sdev_pos, sdev_neg, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ungroup() |>
  arrange(ae) |>
  group_by(n_obs, sdev_pos, sdev_neg, prev_test) |>
  reframe(quantifier = quantifier,
          rank_nr = as.factor(rank(ae)),
          ae = ae) |>
  ungroup() |>
  ggplot(aes(y = fct_reorder(quantifier, as.numeric(rank_nr), mean) |> fct_rev(), 
             fill = rank_nr)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_fill_viridis_d() +
  labs(x = "number of conditions",
       y = "quantifier",
       fill = "rank")

ggsave("_img/sim1/barplot.png", height = 5, width = 7)

sim2_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ungroup() |>
  arrange(ae) |>
  group_by(n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(quantifier = quantifier,
          rank_nr = as.factor(rank(ae)),
          ae = ae) |>
  ungroup() |>
  ggplot(aes(y = fct_reorder(quantifier, as.numeric(rank_nr), mean) |> fct_rev(), 
             fill = rank_nr)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_fill_viridis_d() +
  labs(x = "number of conditions",
       y = "quantifier",
       fill = "rank")

ggsave("_img/sim2/barplot.png", height = 5, width = 7)

sim3_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test, skew) |>
  filter(skew == 1) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ungroup() |>
  arrange(ae) |>
  group_by(n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test, skew) |>
  reframe(quantifier = quantifier,
          rank_nr = as.factor(rank(ae)),
          ae = ae) |>
  ungroup() |>
  ggplot(aes(y = fct_reorder(quantifier, as.numeric(rank_nr), mean) |> fct_rev(), 
             fill = rank_nr)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_fill_viridis_d() +
  labs(x = "number of conditions",
       y = "quantifier",
       fill = "rank")

ggsave("_img/sim3/barplot_skew1.png", height = 5, width = 7)

sim3_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test, skew) |>
  filter(skew == 2) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ungroup() |>
  arrange(ae) |>
  group_by(n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test, skew) |>
  reframe(quantifier = quantifier,
          rank_nr = as.factor(rank(ae)),
          ae = ae) |>
  ungroup() |>
  ggplot(aes(y = fct_reorder(quantifier, as.numeric(rank_nr), mean) |> fct_rev(),
             fill = rank_nr)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_fill_viridis_d() +
  labs(x = "number of conditions",
       y = "quantifier",
       fill = "rank")

ggsave("_img/sim3/barplot_skew2.png", height = 5, width = 7)

sim3_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test, skew) |>
  filter(skew == 4) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ungroup() |>
  arrange(ae) |>
  group_by(n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test, skew) |>
  reframe(quantifier = quantifier,
          rank_nr = as.factor(rank(ae)),
          ae = ae) |>
  ungroup() |>
  ggplot(aes(y = fct_reorder(quantifier, as.numeric(rank_nr), mean) |> fct_rev(), 
             fill = rank_nr)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_fill_viridis_d() +
  labs(x = "number of conditions",
       y = "quantifier",
       fill = "rank")

ggsave("_img/sim3/barplot_skew4.png", height = 5, width = 7)

################################

sim1_clean |>
  group_by(quantifier, n_obs, sdev_pos, sdev_neg, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae)) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim1/boxplot_full.png", height = 5, width = 7)

sim1_clean |>
  group_by(quantifier, n_obs, sdev_pos, sdev_neg, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(n_obs))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "number of observations") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim1/boxplot_nobs.png", height = 5, width = 7)

sim1_clean |>
  group_by(quantifier, n_obs, sdev_pos, sdev_neg, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(sdev_pos))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "SD positive class") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim1/boxplot_sdevpos.png", height = 5, width = 7)

sim1_clean |>
  group_by(quantifier, n_obs, sdev_pos, sdev_neg, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(sdev_neg))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "SD negative class") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim1/boxplot_sdevneg.png", height = 5, width = 7)

sim1_clean |>
  group_by(quantifier, n_obs, sdev_pos, sdev_neg, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(prev_test))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "prevalence test data") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim1/boxplot_prevtest.png", height = 5, width = 7)

##########################
sim2_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(n_train))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "# observations train data") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim2/boxplot_ntrain.png", height = 5, width = 7)

sim2_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(n_test))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "# observations test data") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim2/boxplot_ntest.png", height = 5, width = 7)

sim2_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(sdev_pos))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "SD positive class") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim2/boxplot_sdevpos.png", height = 5, width = 7)

sim2_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(sdev_neg))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "SD negative class") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim2/boxplot_sdevneg.png", height = 5, width = 7)

sim2_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(prev_train))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "prevalence train data") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim2/boxplot_prevtrain.png", height = 5, width = 7)

sim2_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(prev_test))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "prevalence test data") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim2/boxplot_prevtest.png", height = 5, width = 7)

##########################
sim3_clean |>
  filter(skew == 1) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(n_test))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "# observations test data") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew1_ntest.png", height = 5, width = 7)

sim3_clean |>
  filter(skew == 1) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(sdev_pos))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "SD positive class") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew1_sdevpos.png", height = 5, width = 7)

sim3_clean |>
  filter(skew == 1) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(sdev_neg))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "SD negative class") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew1_sdevneg.png", height = 5, width = 7)

sim3_clean |>
  filter(skew == 1) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(prev_test))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "prevalence test data") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew1_prevtest.png", height = 5, width = 7)

##########################
sim3_clean |>
  filter(skew == 2) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(n_test))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "# observations test data") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew2_ntest.png", height = 5, width = 7)

sim3_clean |>
  filter(skew == 2) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(sdev_pos))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "SD positive class") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew2_sdevpos.png", height = 5, width = 7)

sim3_clean |>
  filter(skew == 2) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(sdev_neg))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "SD negative class") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew2_sdevneg.png", height = 5, width = 7)

sim3_clean |>
  filter(skew == 2) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(prev_test))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "prevalence test data") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew2_prevtest.png", height = 5, width = 7)
###########################

sim3_clean |>
  filter(skew == 4) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(n_test))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "# observations test data") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew4_ntest.png", height = 5, width = 7)

sim3_clean |>
  filter(skew == 4) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(sdev_pos))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "SD positive class") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew4_sdevpos.png", height = 5, width = 7)

sim3_clean |>
  filter(skew == 4) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(sdev_neg))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "SD negative class") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew4_sdevneg.png", height = 5, width = 7)

sim3_clean |>
  filter(skew == 4) |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(prev_test))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "prevalence test data") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew4_prevtest.png", height = 5, width = 7)
###########################

sim3_clean |>
  group_by(quantifier, n_train, n_test, sdev_pos, sdev_neg, prev_train, prev_test, skew) |>
  reframe(ae = mean(abs(err), na.rm = T)) |>
  ggplot(aes(x = fct_reorder(quantifier, ae, median), y = ae, fill = as.factor(skew))) +
  geom_boxplot() +
  labs(x = "quantifier", y = "MAE across conditions", fill = "skewness of the data") +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05))

ggsave("_img/sim3/boxplot_skew.png", height = 5, width = 7)

###########################
