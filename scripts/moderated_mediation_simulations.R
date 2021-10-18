# test moderated mediation model
setwd("~/Box/svenson_hf_DO/scripts")
.libPaths("~/GitHub/mediation_measurement_noise/lib/")
library(bmediatR)
library(tidyverse)
library(reshape2)
library(patchwork)
source("bmediatR_moderated_mediation.R")
rename <- dplyr::rename
select <- dplyr::select

library(doFuture)
library(doRNG)
registerDoFuture()
registerDoRNG()
plan(multisession)


theme_set(theme_classic())

colors <- c("Total Mediation" = "seagreen1", "Moderated Mediation A" = "seagreen3", "Moderated Mediation B" = "seagreen",
            "Partial Mediation" = "orchid1", "Moderated Partial A" = "orchid", "Moderated Partial B" = "orchid4", "Moderated Partial C" = "darkorchid3",
            "Co-local" = "dodgerblue", "Moderated Co-local A" = "dodgerblue3", "Moderated Co-local C" = "dodgerblue4",
            "Reactive" = "goldenrod1", "Moderated Reactive B" = "goldenrod", "Moderated Reactive C" = "goldenrod4",
            "Partial Reactive" = "chocolate4", "Other" = "gray", "Other Moderation" = "gray45"
)

transform_df <- function(df1){
  df1 <- df1 %>%
    rename(`Total Mediation` = `1,1,0`, `Reactive` = `0,*,1`, `Partial Mediation` = `1,1,1`,
           `Partial Reactive` = `1,*,1`, `Co-local` = `1,0,1`, `Moderated Mediation A` = `i,1,0`, `Moderated Mediation B` = `1,i,0`,
           `Moderated Reactive C` = `0,*,i`, `Moderated Reactive B` = `0,*i,1`, `Moderated Partial A` = `i,1,1`,
           `Moderated Partial B` = `1,i,1` , `Moderated Partial C` = `1,1,i`,
           `Moderated Co-local C` = `1,0,i`, `Moderated Co-local A` = `i,0,1`)

  df1 %>% select(sim, contains(",")) %>% melt(id = "sim") %>% mutate(moderated = grepl("i", variable)) %>%
    group_by(sim, moderated) %>% summarise(Other = sum(value[!moderated]), `Other Moderation` = sum(value[moderated])) %>%
    group_by(sim) %>% summarise(Other = sum(Other), `Other Moderation` = sum(`Other Moderation`)) %>%
    merge(df1 %>% select(-contains(",")), by = "sim") %>%
    select(-sim) %>% apply(2, mean) %>%
    melt() %>% rownames_to_column("Model") %>%
    mutate(Model = factor(Model, levels = c("Total Mediation", "Moderated Mediation A", "Moderated Mediation B",
                                            "Partial Mediation", "Moderated Partial A", "Moderated Partial B", "Moderated Partial C",
                                            "Co-local", "Moderated Co-local A", "Moderated Co-local C",
                                            "Reactive", "Moderated Reactive B", "Moderated Reactive C",
                                            "Partial Reactive", "Other", "Other Moderation")))
}

plot_bar <- function(df){
  df %>% transform_df() %>%
    ggplot(aes("M", value, fill = Model)) + geom_bar(stat= "identity", position = "fill") +
    labs(x = "", y = "Posterior Probability", fill = "Model") +
    scale_fill_manual(values = colors)
}

plot_bar2 <- function(df_list, text = T){
  if(text){
    text_function <- function() geom_text(aes(color = Model, label = round(value, 3)),
                                          fontface = 2, show.legend = F, vjust = -0.25)
  }else{
    text_function <- function() NULL
  }

  df_list %>% lapply(transform_df) %>% bind_rows(.id = "case") %>%
    ggplot(aes(Model, value, fill = Model)) +
    geom_bar(stat= "identity") +
    facet_wrap(~case, ncol = 1) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme_classic() +
    theme(axis.text.x = element_blank()) +
    ylim(0,1) +
    labs(y = "Posterior Probability", x = "") +
    text_function()
}

N = 200

temp <- function(i, n){
  # Total Mediation with moderation on X -> M
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X + t + X*t)
  Y = rnorm(n, M)

  # data.frame(Rxm = cor(X, M), Rtm = cor(t, M), Rym = cor(Y, M), Rxy = cor(X, Y), Rty = cor(t, Y))

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", 200), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

nsim = 100
df1 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")


# Total Mediation with moderation on M -> Y

temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X)
  Y = rnorm(n, M + t + M*t)
  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", 200), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df2 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

mod_med <- plot_bar2(list("Case 1" = df1, "Case 2" = df2))

ggsave("../fig/moderation_sims/moderated_mediation_prob.png",
       mod_med,
       width = 10, height = 7)


# Total Mediation with no moderation
temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X + t)
  Y = rnorm(n, M + t)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", 200), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df1 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")


temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X + t)
  Y = rnorm(n, M)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", 200), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df2 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")


temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X )
  Y = rnorm(n, M + t)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df3 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

med_modCov <- plot_bar2(list("Case 1" = df1, "Case 2" = df2, "Case 3" = df3))

ggsave("../fig/moderation_sims/mediation_mod_covariate.png",
       med_modCov,
       width = 10, height = 7)



## Partial Mediation with Moderation
temp <- function(i, n){
  # Partial Mediation with moderation on X -> M
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X + t + X*t)
  Y = rnorm(n, M + X)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df1 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")


temp <- function(i, n){
  # Partial Mediation with moderation on X -> M
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X )
  Y = rnorm(n, M + X + t + M*t)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df2 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

temp <- function(i, n){
  # Partial Mediation with moderation on X -> Y
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X )
  Y = rnorm(n, M + X + t + X*t)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df3 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

part_mod <- plot_bar2(list("Case 1" = df1, "Case 2" = df2, "Case 3" = df3))

ggsave("../fig/moderation_sims/moderated_partial_med.png",
       part_mod,
       width = 10, height = 10)


## Partial Mediation without Moderation
temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X + t)
  Y = rnorm(n, M + X + t )

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df1 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")


temp <- function(i, n){
  # Total Mediation with moderation on X -> M
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X + t)
  Y = rnorm(n, M + X )

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df2 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

temp <- function(i, n){
  # Total Mediation with moderation on X -> M
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X )
  Y = rnorm(n, M + X + t )

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df3 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

part_modCov <- plot_bar2(list("Case 1" = df1, "Case 2" = df2, "Case 3" = df3))

ggsave("../fig/moderation_sims/partial_mediation_mod_covariate.png",
       part_modCov,
       width = 10, height = 7)

## Reactive with moderation
temp <- function(i, n){
  # Reactive with moderation on X-> Y
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  Y = rnorm(n, X + t + X*t)
  M = rnorm(n, Y)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df1 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")


temp <- function(i, n){
  # Reactive with moderation on Y -> M
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  Y = rnorm(n, X )
  M = rnorm(n, Y + t + Y*t)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df2 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

reactive_mod <- plot_bar2(list("Case 1" = df1, "Case 2" = df2))

ggsave("../fig/moderation_sims/moderated_reactive_med.png",
       reactive_mod,
       width = 10, height = 7)

# Reactive no moderation
temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  Y = rnorm(n, X + t)
  M = rnorm(n, Y + t)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df1 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  Y = rnorm(n, X + t)
  M = rnorm(n, Y )

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df2 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  Y = rnorm(n, X )
  M = rnorm(n, Y + t)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df3 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

reactive_modCov <- plot_bar2(list("Case 1" = df1, "Case 2" = df2, "Case 3" = df3))

ggsave("../fig/moderation_sims/reactive_mod_covariate.png",
       reactive_modCov,
       width = 10, height = 7)


## Co-local with moderation
temp <- function(i, n){
  # Co-local with moderation on X-> Y
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  Y = rnorm(n, X + t + X*t)
  M = rnorm(n, X)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df1 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  Y = rnorm(n, X )
  M = rnorm(n, X + t + X*t)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df2 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

colocal_mod <- plot_bar2(list("Case 1" = df1, "Case 2" = df2))

ggsave("../fig/moderation_sims/moderated_colocal.png",
       colocal_mod,
       width = 10, height = 7)


# Co-local no moderation
temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  Y = rnorm(n, X + t)
  M = rnorm(n, X + t)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df1 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")


temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  Y = rnorm(n, X + t)
  M = rnorm(n, X )

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df2 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")


temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  Y = rnorm(n, X )
  M = rnorm(n, X + t)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df3 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

colocal_modCov <- plot_bar2(list("Case 1" = df1, "Case 2" = df2, "Case 3" = df3))

ggsave("../fig/moderation_sims/colocal_mod_covariate.png",
       colocal_modCov,
       width = 10, height = 7)


## Scenario for which we do not have a moderation model included

# Mediation with moderation on both edges
temp <- function(i, n){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X + t + X*t)
  Y = rnorm(n, M + t + M*t)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df1 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

medAB <- plot_bar2(list(df1)) +
  theme(strip.background = element_blank(), strip.text = element_blank())

ggsave("../fig/moderation_sims/mediation_moderated_AB.png", medAB,
       width = 10, height = 5)


## Partial Mediation with moderation on two arms
temp <- function(i, n){
  # Total Mediation with moderation on X -> M
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X + t + X*t)
  Y = rnorm(n, M + X + t + X*t)

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

df1 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, n = N, .id = "sim")

plot_bar(df1)

################################################################################
## Simulate complete mediation with error on M and Y
################################################################################

# Total Mediation with moderation on X -> M
temp <- function(i, n, error_M, error_Y){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X + t + X*t)
  Y = rnorm(n, M)

  M = sim_target_from_mediator(as.matrix(M), error_M^2)$data
  Y = sim_target_from_mediator(as.matrix(Y), error_Y^2)$data

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

# both error cors = 0.9
df1 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.9, error_M = 0.9, n = N, .id = "sim")

# both error cors = 0.8
df2 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.8, error_M = 0.8, n = N, .id = "sim")

# both error cors = 0.7
df3 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.7, error_M = 0.7, n = N, .id = "sim")

mod_medA_equal_err <- plot_bar2(
  list(df1, df2, df3) %>% setNames(paste("Case", seq(1, length(.))))
) +
  theme(text = element_text(size = 15))

ggsave("../fig/moderation_sims/moderated_mediationA_with_equal_error.png",
       mod_medA_equal_err,
       width = 10, height = 7)

# varying error cor
df4 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.9, error_M = 0.8, n = N, .id = "sim")

# varying error cor
df5 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.9, error_M = 0.7, n = N, .id = "sim")

# varying error cor
df6 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.8, error_M = 0.9, n = N, .id = "sim")

# varying error cor
df7 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.7, error_M = 0.9, n = N, .id = "sim")

mod_med_unequal_err <- plot_bar2(list(df4, df5, df6, df7) %>%
                                   setNames(paste("Case", 4:7))) +
  theme(text = element_text(size = 15))

ggsave("../fig/moderation_sims/moderated_mediationA_with_unequal_error.png",
       mod_med_unequal_err,
       width = 10, height = 10)

# Total Mediation with moderation on M -> Y
temp <- function(i, n, error_M, error_Y){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X)
  Y = rnorm(n, M + t + M*t)

  M = sim_target_from_mediator(as.matrix(M), error_M^2)$data
  Y = sim_target_from_mediator(as.matrix(Y), error_Y^2)$data

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

# both error cors = 0.9
df1 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.9, error_M = 0.9, n = N, .id = "sim")

# both error cors = 0.8
df2 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.8, error_M = 0.8, n = N, .id = "sim")

# both error cors = 0.7
df3 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.7, error_M = 0.7, n = N, .id = "sim")

mod_medB_equal_err <- plot_bar2(list(df1, df2, df3) %>%
                                  setNames(paste("Case", seq(1, length(.))))) +
  theme(text = element_text(size = 15))

ggsave("../fig/moderation_sims/moderated_mediationB_with_equal_error.png",
       mod_medB_equal_err,
       width = 10, height = 7)

# varying error cor
df4 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.9, error_M = 0.8, n = N, .id = "sim")

# varying error cor
df5 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.9, error_M = 0.7, n = N, .id = "sim")

# varying error cor
df6 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.8, error_M = 0.9, n = N, .id = "sim")

# varying error cor
df7 <- 1:100 %>% setNames(1:100) %>% plyr::ldply(temp, error_Y = 0.7, error_M = 0.9, n = N, .id = "sim")

mod_medB_unequal_err <- plot_bar2(list(df4, df5, df6, df7) %>%
                                    setNames(paste("Case", 4:7))) +
  theme(text = element_text(size = 15))

ggsave("../fig/moderation_sims/moderated_mediationB_with_unequal_error.png",
       mod_medB_unequal_err,
       width = 10, height = 10)


################################################################################
## Larger scale moderated mediation with error simulations
## We assume that X and t are measured without error
################################################################################

# Total Mediation with moderation on X -> M
temp <- function(i, n, error_M, error_Y){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X + t + X*t)
  Y = rnorm(n, M)

  M = sim_target_from_mediator(as.matrix(M), error_M^2)$data
  Y = sim_target_from_mediator(as.matrix(Y), error_Y^2)$data

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

mod_medA_error <- expand.grid(error_M = seq(0.05, 0.95, by = 0.05),
                              error_Y = seq(0.05, 0.95, by = 0.05)) %>%
  mutate(N = 200, nsim = 100) %>%
  plyr::mdply(function(error_M, error_Y, N, nsim){
    1:nsim %>% setNames(1:nsim) %>%
      plyr::ldply(temp, error_Y = error_Y, error_M = error_M, n = N, .id = "sim")
  }, .parallel = T) %>%
  select(-nsim, -N)

saveRDS(mod_medA_error, file = "../data/moderated_mediationA_with_error.rds")


# Total Mediation with moderation on M -> Y
temp <- function(i, n, error_M, error_Y){
  X = rbinom(n, 1, 0.5)
  t = rbinom(n, 1, 0.5)
  M = rnorm(n, X)
  Y = rnorm(n, M + t + M*t)

  M = sim_target_from_mediator(as.matrix(M), error_M^2)$data
  Y = sim_target_from_mediator(as.matrix(Y), error_Y^2)$data

  med_results <- bmediatR_moderation(Y, matrix(M, dimnames = list(rep("", n), "M1")),
                                     create_design_mat(X), create_design_mat(t),
                                     ln_prior_c = "reactive", verbose = F, align_data = F,
                                     options_X = list(sum_to_zero = T, center = F, scale = F),
                                     options_t = list(sum_to_zero = T, center = F, scale = F))
  as.data.frame(exp(med_results$ln_post_c))
}

mod_medB_error <- expand.grid(error_M = seq(0.05, 0.95, by = 0.05),
                              error_Y = seq(0.05, 0.95, by = 0.05)) %>%
  mutate(N = 200, nsim = 100) %>%
  plyr::mdply(function(error_M, error_Y, N, nsim){
    1:nsim %>% setNames(1:nsim) %>%
      plyr::ldply(temp, error_Y = error_Y, error_M = error_M, n = N, .id = "sim")
  }, .parallel = T) %>%
  select(-nsim, -N)

saveRDS(mod_medB_error, file = "../data/moderated_mediationB_with_error.rds")


mod_medB_error_long <- mod_medB_error %>%
  rename(`Total Mediation` = `1,1,0`, `Reactive` = `0,*,1`, `Partial Mediation` = `1,1,1`,
         `Partial Reactive` = `1,*,1`, `Co-local` = `1,0,1`, `Moderated Mediation A` = `i,1,0`, `Moderated Mediation B` = `1,i,0`,
         `Moderated Reactive C` = `0,*,i`, `Moderated Reactive B` = `0,*i,1`, `Moderated Partial A` = `i,1,1`,
         `Moderated Partial B` = `1,i,1` , `Moderated Partial C` = `1,1,i`,
         `Moderated Co-local C` = `1,0,i`, `Moderated Co-local A` = `i,0,1`) %>%
  rowwise() %>%
  mutate(`Other Moderation` = sum(c_across(c(contains("i") & contains(",")))),
         Other = sum(c_across(c(contains(",") & !contains("i"))))) %>%
  ungroup() %>%
  select(-contains(",")) %>%
  pivot_longer(!c(error_M, error_Y, sim), names_to = "Model", values_to = "post_prob") %>%
  mutate(Model = factor(Model, levels = c("Total Mediation", "Moderated Mediation A", "Moderated Mediation B",
                                          "Partial Mediation", "Moderated Partial A", "Moderated Partial B", "Moderated Partial C",
                                          "Co-local", "Moderated Co-local A", "Moderated Co-local C",
                                          "Reactive", "Moderated Reactive B", "Moderated Reactive C",
                                          "Partial Reactive", "Other", "Other Moderation")))

mod_medB_error_long %>%
  group_by(error_M, error_Y, Model) %>%
  summarise( post_prob = mean(post_prob)) %>%
  ggplot(aes(error_M, error_Y, post_prob, fill = post_prob)) +
  geom_tile() +
  facet_wrap(~Model) +
  scale_fill_viridis_c()









