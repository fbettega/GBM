##########################################################
##                    packages                          ##
##########################################################
library(tidyverse)
library(nnet)
library(twang)
library(HTLR)
library(PSweight)
library(CBPS)
# Packages GBM
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("slice", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
##########################################################
##                    function                          ##
##########################################################
source("Source_papier_prett.R",
       encoding = "UTF-8"
)

get_ps_custom <- function(twang_res, desc_str) {
  data.frame(
    PS_es.max = 1 / get.weights(twang_res,
                                stop.method = "es.max",
                                withSampW = FALSE
    ) %>%
      unlist(use.names = FALSE),
    PS_es.mean = 1 / get.weights(twang_res,
                                 stop.method = "es.mean",
                                 withSampW = FALSE
    ) %>%
      unlist(use.names = FALSE)
  ) %>% setNames(paste0(desc_str, names(.)))
}



name_expo <- paste0(quote(Treatment))
name_outcome <- paste0(quote(out_come))
Id <- paste0(quote(id))


simulation_fun <- function(n = 30000) {
  df_sim <- data.frame(
    confu_1 = rnorm(n),
    confu_2 = rnorm(n),
    confu_3 = rnorm(n),
    confu_4 = rnorm(n),
    confu_5 = rnorm(n),
    confu_6 = rnorm(n),
    confu_7 = rnorm(n),
    confu_8 = rnorm(n),
    confu_9 = rnorm(n),
    confu_10 = rnorm(n),
    confu_bool_1 = rbinom(n, 1, 0.5),
    confu_bool_2 = rbinom(n, 1, 0.5),
    confu_bool_3 = rbinom(n, 1, 0.2),
    confu_bool_4 = rbinom(n, 1, 0.2),
    confu_bool_5 = rbinom(n, 1, 0.05),
    confu_bool_6 = rbinom(n, 1, 0.2),
    confu_bool_7 = rbinom(n, 1, 0.5),
    confu_bool_8 = rbinom(n, 1, 0.05),
    confu_bool_9 = rbinom(n, 1, 0.2),
    confu_bool_10 = rbinom(n, 1, 0.5)
  )
  
  coef <- c(
    1.29, 0.932, -0.713, 0.562, # interaction 1 # confu1
    -0.838, -0.068, -0.798, 1.035, -0.998,# interaction 1 # confu2
    0.778, 1.588, 0.419, # interaction 2 # confu3
    0.641, 0.478, 0.893, # interaction 2 # confu4
    -0.164, 0.294, 0.602, 0.5, 0.417, # interaction 3 # confu7
    1.172, 0.538, -2.147, 0.702, # confu9
    -0.639, 0.462, -0.594,0.479, -0.385, # interaction 3 # confu10
    1.435, 1.295, 0.871, 
    0.767, .578, 0.266,
    0.89, -0.161, -0.242, -0.624, 0.318,
    1.432, 0.772, 0.942, 1.524,
    -0.444,-0.8, 0.33, 0.143, 0.71,
    0.514, -0.215, 0.869, 
    0.652, 0.483, 0.358,
    0.873, 0.462, 1.328, 0.106, 0.957,
    # end additive
    -0.413, -0.649, -0.212, # interaction 1 
    0.471, 0.397, 0.249,# interaction 2
    0.315, 0.654, 0.173, # interaction 3
    -1.206, 1.23, -0.3, 0.506,
    0.146, 1.543, -0.656, 0.6, 0.035, -0.78, 1.648, 0.092, 0.818, 0.92,
    0.193, 0.549, 0.967, 0.612, -0.766, -1.383,-1.67, 
    1,1.75,2.5, # outcome coef 
    0.173 , 0.476 , 1.371 ,-1.144, - 0.189 , 0.217, 0.170 , 0.602, 0.284 ,#instru coef 
    -1.344,0.299,1.187# outcome interaction coef
  ) # Treatment coef

  
  
  
  treatments_prob_add <-   gendata_MLR(n,
                                       p = ncol(df_sim) -5 ,
                                       NC = 3, 
                                       X = as.matrix(df_sim)[,-c(3,6,8,16,18)],
                                       betas = rbind(c(0,0.2,0),
                                                     matrix(coef[1:60][-c(3,6,8,16,18,
                                                                          23,26,28,36,38,
                                                                          43,46,48,56,58)],
                                                            ncol = 3,
                                                            byrow = TRUE,
                                                     dimnames = list(colnames(as.matrix(df_sim)[,-c(3,6,8,16,18)]),
                                                                     c(1:3)
                                                                     )
                                                     )
                                                     )
  )
  

  df_sim_add <- df_sim %>%
    mutate(
      id = 1:n,
      Treatment = treatments_prob_add$y,
      # Treatment = apply(
      #   treatments_prob_add, 1,
      #   function(x) {
      #     sample(1:3, 1, prob = softmax(x))
      #   }
      # ),
      out_come = rnorm(n,0,0.1) +
        coef[70] * (confu_1) +
        coef[71] * (confu_2) +
        coef[72] * (confu_3) +
        coef[73] * (confu_4) +
        coef[74] * (confu_5) +
        coef[75] * (confu_6) +
        coef[76] * (confu_7) +
        coef[77] * (confu_8) +
        coef[78] * (confu_9) +
        coef[79] * (confu_10) +
        (coef[80] * confu_bool_1) +
        (coef[81] * confu_bool_2) +
        (coef[82] * confu_bool_3) +
        (coef[83] * confu_bool_4) +
        (coef[84] * confu_bool_5) +
        (coef[85] * confu_bool_6) +
        (coef[86] * confu_bool_7) +
        (coef[87] * confu_bool_8) +
        (coef[88] * confu_bool_9) +
        (coef[89] * confu_bool_10) +
        coef[91] * as.numeric(Treatment == 1) +
        coef[92] * as.numeric(Treatment == 2) +
        coef[93] * as.numeric(Treatment == 3)
    ) %>%
    mutate_at(vars(starts_with("confu_bool_")), as.factor) %>%
    mutate(Treatment = as.factor(Treatment)
           )
  
  
  df_sim_compl_temp <- df_sim %>% mutate(confuint_1_2 = confu_1 * confu_2 ,
                                         confuint_4_5 = confu_4 * confu_5 ,
                                         confuint_7_10 = confu_7 * confu_10 )
  
  
  treatments_prob_compl <- gendata_MLR(n,p=ncol(df_sim_compl_temp) -5 ,
                                       NC = 3,
                                       X = as.matrix(df_sim_compl_temp)[,-c(3,6,8,16,18)],
                                       betas = rbind(c(0,0.2,0),
                                                     matrix(coef[1:69][-c(3 ,6 ,8 ,16,18,
                                                                          23,26,28,36,38,
                                                                          43,46,48,56,58)],
                                                            ncol = 3,
                                                            byrow = TRUE,
                                                            dimnames = list(c(
                                                              colnames(as.matrix(df_sim)[,-c(3,6,8,16,18)]),
                                                              "confuint_1_2","confuint_4_5","confuint_7_10"
                                                              ),
                                                                                        c(1:3)
                                                            )
                                                     )
                                       )
  )
  
  
  
  
  df_noise <- df_sim %>% 
    mutate(
      instru_1 = rnorm(n),
      instru_2 = rnorm(n),
      instru_bin =  rbinom(n, 1, 0.64))
  
  
  treatments_prob_noise <-   gendata_MLR(n,
                                       p=ncol(df_noise) - 5 ,
                                       NC = 3, 
                                       X = as.matrix(df_noise)[,-c(3,6,8,16,18)],
                                       betas = rbind(c(0,0.2,0),
                                                     matrix(coef[c(1:60,94:102)][-c(3,6,8,16,18,
                                                                          23,26,28,36,38,
                                                                          43,46,48,56,58)],
                                                            ncol = 3,
                                                            byrow = TRUE))
  )
  

  
  df_sim_noise <- df_noise %>%
    mutate(
      id = 1:n,
      Treatment = treatments_prob_noise$y,
      out_come = rnorm(n,0.1) +
        coef[70] * (confu_1) +
        coef[71] * (confu_2) +
        coef[72] * (confu_3) +
        coef[73] * (confu_4) +
        coef[74] * (confu_5) +
        coef[75] * (confu_6) +
        coef[76] * (confu_7) +
        coef[77] * (confu_8) +
        coef[78] * (confu_9) +
        coef[79] * (confu_10) +
        (coef[80] * confu_bool_1) +
        (coef[81] * confu_bool_2) +
        (coef[82] * confu_bool_3) +
        (coef[83] * confu_bool_4) +
        (coef[84] * confu_bool_5) +
        (coef[85] * confu_bool_6) +
        (coef[86] * confu_bool_7) +
        (coef[87] * confu_bool_8) +
        (coef[88] * confu_bool_9) +
        (coef[89] * confu_bool_10) +
        coef[91] * as.numeric(Treatment == 1) +
        coef[92] * as.numeric(Treatment == 2) +
        coef[93] * as.numeric(Treatment == 3),
        noise_1_continous = rnorm(n),
        noise_2_continous = rnorm(n),
        noise_3_bin = rbinom(n, 1, 0.42)
    ) %>%
    mutate_at(vars(starts_with("confu_bool_")), as.factor) %>%
    mutate(Treatment = as.factor(Treatment),
           noise_3_bin = as.factor(noise_3_bin))
  
  
  
  
  df_sim_compl <- df_sim %>%
    mutate(
      id = 1:n,
      Treatment = treatments_prob_compl$y,
      # Treatment = apply(
      #   treatments_prob_compl, 1,
      #   function(x) {
      #     sample(1:3, 1, prob = softmax(x))
      #   }
      # ),
      out_come = rnorm(n,0.1) +
        coef[70] * (confu_1) +
        coef[71] * (confu_2) +
        coef[72] * (confu_3) +
        coef[73] * (confu_4) +
        coef[74] * (confu_5) +
        coef[75] * (confu_6) +
        coef[76] * (confu_7) +
        coef[77] * (confu_8) +
        coef[78] * (confu_9) +
        coef[79] * (confu_10) +
        (coef[80] * confu_bool_1) +
        (coef[81] * confu_bool_2) +
        (coef[82] * confu_bool_3) +
        (coef[83] * confu_bool_4) +
        (coef[84] * confu_bool_5) +
        (coef[85] * confu_bool_6) +
        (coef[86] * confu_bool_7) +
        (coef[87] * confu_bool_8) +
        (coef[88] * confu_bool_9) +
        (coef[89] * confu_bool_10) +
        coef[103] * df_sim_compl_temp$confuint_1_2 +
        coef[104] * df_sim_compl_temp$confuint_4_5 +
        coef[105] * df_sim_compl_temp$confuint_7_10 +
        coef[91] * as.numeric(Treatment == 1) +
        coef[92] * as.numeric(Treatment == 2) +
        coef[93] * as.numeric(Treatment == 3)
    ) %>%
    mutate_at(vars(starts_with("confu_bool_")), as.factor) %>%
    mutate(Treatment = as.factor(Treatment))
  
  

  
  
  return(list(
    coefficient = coef,
    out_model_add = treatments_prob_add,
    out_model_compl = treatments_prob_compl,
    df_sim = df_sim,
    df_add = df_sim_add,
    df_compl = df_sim_compl,
    df_noise = df_sim_noise
  ))
}


PS_estimation <- function(Sim) {

  IPTW_add <- calcule_pds_stage(
    Sim$df_add,
    expo = Treatment,
    covar = colnames(Sim$df_sim),
    # var_mod_pds ,
    Id_data = id,
    weighting_function = multinomial_IPTW,
    out_come = out_come # ,
    # percentile_tronc = 0.01
  )
  
  PS_IPTW_add <-  1/IPTW_add$res_intermediaire$poids$poids_tronc$poids_trunc$`(0; 1)`
  
  #rm(IPTW_add)
  
  
  IPTW_comp <- calcule_pds_stage(
    Sim$df_compl,
    expo = Treatment,
    covar = colnames(Sim$df_sim), 
    # var_mod_pds ,
    Id_data = id,
    weighting_function = multinomial_IPTW,
    out_come = out_come # ,
    # percentile_tronc = 0.01
  )
  PS_IPTW_comp <- 1/IPTW_comp$res_intermediaire$poids$poids_tronc$poids_trunc$`(0; 1)`
  #rm(IPTW_comp)
  
  IPTW_comp_exact_formula <- calcule_pds_stage(
    Sim$df_compl,
    expo = Treatment,
    covar = colnames(Sim$df_sim),
    # var_mod_pds ,
    Id_data = id,
    weighting_function = multinomial_IPTW_modif_compl,
    out_come = out_come # ,
    # percentile_tronc = 0.01
  )
  PS_IPTW_comp_exact_formula <- 1/IPTW_comp_exact_formula$res_intermediaire$poids$poids_tronc$poids_trunc$`(0; 1)`
  #rm(IPTW_comp_exact_formula)
  
  
  IPTW_noise <- calcule_pds_stage(
    Sim$df_noise,
    expo = Treatment,
    covar = colnames(Sim$df_noise %>% 
                       select(-id,-out_come,-Treatment)),
    # var_mod_pds ,
    Id_data = id,
    weighting_function = multinomial_IPTW,
    out_come = out_come # ,
    # percentile_tronc = 0.01
  )
  
  PS_IPTW_noise <- 1/IPTW_noise$res_intermediaire$poids$poids_tronc$poids_trunc$`(0; 1)`
  #rm(IPTW_noise)
 
  IPTW_noise_exact <- calcule_pds_stage(
    Sim$df_noise[,-c(21:23,27:29)],
    expo = Treatment,
    covar = colnames(Sim$df_noise[,-c(21:23,27:29)] %>% 
                       select(-id,-out_come,-Treatment)),
    # var_mod_pds ,
    Id_data = id,
    weighting_function = multinomial_IPTW,
    out_come = out_come # ,
    # percentile_tronc = 0.01
  )
  
  PS_IPTW_noise_exact <- 1/IPTW_noise_exact$res_intermediaire$poids$poids_tronc$poids_trunc$`(0; 1)`
  #rm(IPTW_noise_exact)  
  
  
  
  IPTW_res <-
    list(
      IPTW_add = PS_IPTW_add,
      IPTW_comp = PS_IPTW_comp,
      IPTW_comp_exact_formula = PS_IPTW_comp_exact_formula,
      IPTW_noise = PS_IPTW_noise,
      IPTW_noise_exact = PS_IPTW_noise_exact
    )
  
  print("regression")
  print(Sys.time())
  formula_gbm <-
    as.formula(paste0("Treatment", " ~ ", paste0(colnames(Sim$df_sim),#[-c(3,6,8,16,18)],
                                                 collapse = " + ")))

  twang_res_ATE_add <-
    mnps(
      formula_gbm,
      data = Sim$df_add,
      estimand = "ATE",
      verbose = FALSE,
      # perm.test.iters = 500,
      stop.method = c("es.mean"),
      n.trees = 5000,
      interaction.depth = 2,
      shrinkage = 0.05,
      bag.fraction = 0.8
    )
  
  PS_GBM_add <-   1 / get.weights(
    twang_res_ATE_add,
    stop.method = "es.mean",
    withSampW = TRUE
  ) %>%
    unlist(use.names = FALSE)
  
  #rm(twang_res_ATE_add)
  print("GBM1")
  print(Sys.time())
  
  twang_res_ATE_compl <-
    mnps(
      formula_gbm,
      data = Sim$df_compl,
      estimand = "ATE",
      verbose = FALSE,
      #perm.test.iters = 500,
      stop.method = c("es.mean"),
      n.trees = 5000,
      interaction.depth = 2,
      shrinkage = 0.05,
      bag.fraction = 0.8
    )
  
  PS_GBM_comp <-  1 / get.weights(
    twang_res_ATE_compl,
    stop.method = "es.mean",
    withSampW = TRUE
  ) %>%
    unlist(use.names = FALSE)
  
  #rm(twang_res_ATE_compl)
  print("GBM2")
  print(Sys.time())
  
  
  
  formula_gbm2 <-
    as.formula(paste0("Treatment", " ~ ",
                      paste0(colnames(Sim$df_noise %>% 
                                        select(-id,-out_come,-Treatment)), 
                             collapse = " + ")
                      )
               )
  
  twang_res_ATE_noise <-
    mnps(
      formula_gbm2,
      data = Sim$df_noise,
      estimand = "ATE",
      verbose = FALSE,
      #perm.test.iters = 500,
      stop.method = c("es.mean"),
      n.trees = 5000,
      interaction.depth = 2,
      shrinkage = 0.05,
      bag.fraction = 0.8
    )
  
  
  
  PS_GBM_noise <-  1 / get.weights(
    twang_res_ATE_noise,
    stop.method = "es.mean",
    withSampW = TRUE
  ) %>%
    unlist(use.names = FALSE)
  
  #rm(twang_res_ATE_noise)
  print("GBM3")
  print(Sys.time())
  
  
  
  GBM_res <-
    list(
      twang_res_ATE_add = PS_GBM_add,
      twang_res_ATE_compl = PS_GBM_comp,
      twang_res_ATE_noise = PS_GBM_noise
      
    )
  
  
  CBPS_res_add <- CBPS(as.formula(paste0("Treatment", " ~ ",
                                         paste0(colnames(Sim$df_add %>% 
                                                           select(-id,-out_come,-Treatment)), 
                                                collapse = " + ")
  )
  ), data = Sim$df_add,
    standardize = FALSE, ATT = 0)
  
  CBPS_res_comp_misspe <- CBPS(as.formula(paste0("Treatment", " ~ ",
                                                 paste0(colnames(Sim$df_compl %>% 
                                                                   select(-id,-out_come,-Treatment)), 
                                                        collapse = " + ")
  )
  ), data = Sim$df_compl,
    standardize = FALSE, ATT = 0)
  
  CBPS_res_comp_correct <- CBPS(as.formula(paste0("Treatment", " ~ ",
                                                  paste0(c(colnames(Sim$df_compl %>% 
                                                                    select(-id,-out_come,-Treatment)),
                                                           "confu_1 * confu_2 + confu_4 * confu_5 + confu_7 * confu_10"), 
                                                         collapse = " + ")
  )
  ), data = Sim$df_compl,
    standardize = FALSE, ATT = 0)
  
  
  
  
  CBPS_res_noise_misspe <- CBPS(as.formula(paste0("Treatment", " ~ ",
                                                  paste0(colnames(Sim$df_noise %>% 
                                                                    select(-id,-out_come,-Treatment)), 
                                                         collapse = " + ")
  )
  ), data = Sim$df_noise,
  
    standardize = FALSE, ATT = 0)
  
  CBPS_res_noise_correct <- CBPS(as.formula(paste0("Treatment", " ~ ",
                                                   paste0(colnames(Sim$df_noise[,-c(21:23,27:29)] %>% 
                                                                     select(-id,-out_come,-Treatment)), 
                                                          collapse = " + ")
  )
  ), data = Sim$df_noise,
  standardize = FALSE, 
  ATT = 0)
  print("CBPS")
  print(Sys.time())

  CBPS_res <-
    list(
      CBPS_res_add = CBPS_res_add,
      CBPS_res_comp_misspe = CBPS_res_comp_misspe,
      CBPS_res_comp_correct = CBPS_res_comp_correct,
      CBPS_res_noise_misspe = CBPS_res_noise_misspe,
      CBPS_res_noise_correct = CBPS_res_noise_correct
    )
  
  
  
  return(list(
    IPTW_res = IPTW_res,
    GBM_res = GBM_res,
    CBPS_res = CBPS_res,
    model_complet = list(
      add_IPTW = IPTW_add,
      comp_IPTW = IPTW_comp,
      comp_exact_IPTW = IPTW_comp_exact_formula,
      noise_IPTW = IPTW_noise,
      noise_exact_IPTW = IPTW_noise_exact,
      add_gbm = twang_res_ATE_add,
      compl_gbm = twang_res_ATE_compl,
      noise_gbm = twang_res_ATE_noise,
      CBPS_add = CBPS_res_add,
      CBPS_comp_misspe = CBPS_res_comp_misspe,
      CBPS_comp_correct = CBPS_res_comp_correct,
      CBPS_noise_misspe = CBPS_res_noise_misspe,
      CBPS_noise_correct = CBPS_res_noise_correct
      )
  ))
}


#########################################################################################
ATE_fun <- function(exposure, outcome, weight) {
  if (!is.factor(exposure)) {
    exposure <- as.factor(exposure)
    warning("Exposure is not a factor so exposure is cast as factor with ref level : ",
            levels(expo)[1])
  }
  if (any(weight <= 0)) {
    stop("negative weights at index: ", which(weight <= 0))
  }
  
  df <- data.frame(
    expo = exposure,
    outcome = outcome,
    weights = 1/weight
  )
  
  ATE_mod <- lm(outcome ~ expo,weights = weights,data = df)

  res <- tibble(expo = factor(str_remove(names(ATE_mod$coefficients[-1]),"expo")),
                ATE = ATE_mod$coefficients[-1],
                std = coef(summary(ATE_mod))[-1, "Std. Error"] 
                )
  
  return(res)
}


#########################################################################################







# IPTW_noise_exact <- calcule_pds_stage(
#   Sim$df_noise[,-c(21:23,27:29)],
#   expo = Treatment,
#   covar = colnames(Sim$df_noise[-c(3,6,8,16,18,21:23,27:29)] %>% 
#                      select(-id,-out_come,-Treatment)),
#   # var_mod_pds ,
#   Id_data = id,
#   weighting_function = multinomial_IPTW,
#   out_come = out_come # ,
#   # percentile_tronc = 0.01
# )
# 
# PS_IPTW_noise_exact <- 1/IPTW_noise_exact$res_intermediaire$poids$poids_tronc$poids_trunc$`(0; 1)`



















#########################################################################################
# Boot strap other a function
boot_strap_fun <- function(df, fonction, nb_iter, params_fun) {
  replicate(nb_iter, 
            fonction(df,
                     sample(1:nrow(df),
                            size = nrow(df),
                            replace = TRUE), 
                     params_fun)
  )
}

ATE_boot_add <- function(df, index, weight_df) {
  df_fun <- df[index, ]
  
  
  # 
  # IPTW_add <- calcule_pds_stage(
  #   df_fun,
  #   expo = Treatment,
  #   covar = colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                      select(-id,-out_come,-Treatment)),
  #   Id_data = id,
  #   weighting_function = multinomial_IPTW,
  #   out_come = out_come 
  # )
  # 
  # weight_iptw_add <-  1/IPTW_add$res_intermediaire$poids$poids_tronc$poids_trunc$`(0; 1)`
  
  weight_iptw_add <- weight_df$IPTW_res$IPTW_add[index]
  
  weight_gbm_add <- weight_df$GBM_res$twang_res_ATE_add[index]
  weight_CBPS_add <- 1/weight_df$CBPS_res$CBPS_res_add$weights[index]
  

  ATE_add_IPTW <- ATE_fun(exposure = df_fun$Treatment,
                          outcome = df_fun$out_come,
                          weight = weight_iptw_add)
  
  ATE_add_GBM <- ATE_fun(exposure = df_fun$Treatment,
                         outcome = df_fun$out_come,
                         weight = weight_gbm_add)
  
  ATE_add_CBPS <- ATE_fun(exposure = df_fun$Treatment,
                          outcome = df_fun$out_come,
                          weight = weight_CBPS_add)
 
  
  # ATE_add_AIPW <- PSweight(
  #   as.formula(paste0("Treatment", " ~ ",
  #                     paste0(colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                                       select(-id,-out_come,-Treatment)), 
  #                            collapse = " + ")
  # )
  # ),
  # data = df_fun,
  # zname = "Treatment",
  # yname = "out_come",
  # augmentation = TRUE,
  # out.formula = 
  #   as.formula(paste0("out_come", " ~ ",
  #                     paste0(colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                                       select(-id,-out_come,-Treatment)), 
  #                            collapse = " + ")
  #                    )
  #                    )
  # )
  # 
  # 
  # ATE_add_adjustment = lm(as.formula(paste0("out_come", " ~ ","Treatment +",
  #                                           paste0(colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                                                             select(-id,-out_come,-Treatment)), 
  #                                                  collapse = " + ")
  # )
  # ), data = df_fun
  # 
  # )
  
  
  
  res <-  data.frame(
    expo = ATE_add_IPTW$expo,
    ATE_IPTW = ATE_add_IPTW$ATE,
    ATE_GBM = ATE_add_GBM$ATE,
    ATE_CBPS = ATE_add_CBPS$ATE#,
    # ATE_AIPW = t(summary(ATE_add_AIPW)$estimates[c(1,2),1]),
    # ATE_AIPW = t(summary(ATE_add_adjustment)$coefficients[c(2,3),1]) 
)
  return(res)
}


ATE_boot_compl <- function(df, index, weight_df) {
  df_fun <- df[index, ]
  
  
  # IPTW_comp <- calcule_pds_stage(
  #   df_fun,
  #   expo = Treatment,
  #   covar = colnames(df_fun)[-c(3,6,8,16,18)],
  #   # var_mod_pds ,
  #   Id_data = id,
  #   weighting_function = multinomial_IPTW,
  #   out_come = out_come # ,
  #   # percentile_tronc = 0.01
  # ) 
  # 
  # weight_iptw_msp_com <- 1/IPTW_comp$res_intermediaire$poids$poids_tronc$poids_trunc$`(0; 1)`
  # 
  # IPTW_comp_exact_formula <- calcule_pds_stage(
  #   df_fun,
  #   expo = Treatment,
  #   covar = colnames(df_fun)[-c(3,6,8,16,18)],
  #   # var_mod_pds ,
  #   Id_data = id,
  #   weighting_function = multinomial_IPTW_modif_compl,
  #   out_come = out_come # ,
  #   # percentile_tronc = 0.01
  # )
  # 
  # 
  # weight_iptw_godd_mod_com <- 1/IPTW_comp_exact_formula$res_intermediaire$poids$poids_tronc$poids_trunc$`(0; 1)`
  weight_iptw_msp_com <- weight_df$IPTW_res$IPTW_comp[index]
  weight_iptw_godd_mod_com  <- weight_df$IPTW_res$IPTW_comp_exact_formula[index]
  
  weight_gbm_comp <- weight_df$GBM_res$twang_res_ATE_compl[index]
  weight_CBPS_msp_com <- 1/weight_df$CBPS_res$CBPS_res_comp_misspe$weights[index]
  weight_CBPS_godd_mod_com <- 1/weight_df$CBPS_res$CBPS_res_comp_correct$weights[index]
  
  
  
  
  ATE_iptw_msp_com <- ATE_fun(exposure = df_fun$Treatment, 
                              outcome = df_fun$out_come, 
                              weight = weight_iptw_msp_com)
  ATE_iptw_godd_mod_com <- ATE_fun(exposure = df_fun$Treatment, 
                                   outcome = df_fun$out_come, 
                                   weight = weight_iptw_godd_mod_com)
  ATE_gbm_comp <- ATE_fun(exposure = df_fun$Treatment,
                          outcome = df_fun$out_come,
                          weight = weight_gbm_comp)
  ATE_comp_msp_CBPS <- ATE_fun(exposure = df_fun$Treatment,
                          outcome = df_fun$out_come,
                          weight = weight_CBPS_msp_com)
  ATE_comp_godd_mod_CBPS <- ATE_fun(exposure = df_fun$Treatment,
                          outcome = df_fun$out_come,
                          weight = weight_CBPS_godd_mod_com)
  

  # ATE_comp_msp_AIPW <- PSweight(
  #   as.formula(paste0("Treatment", " ~ ",
  #                     paste0(colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                                       select(-id,-out_come,-Treatment)), 
  #                            collapse = " + ")
  #   )
  #   ),
  #   data = df_fun,
  #   zname = "Treatment",
  #   yname = "out_come",
  #   augmentation = TRUE,
  #   out.formula = 
  #     as.formula(paste0("out_come", " ~ ",
  #                       paste0(colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                                         select(-id,-out_come,-Treatment)), 
  #                              collapse = " + ")
  #     )
  #     )
  # )
  # 
  # 
  # ATE_comp_msp_adjustment = lm(as.formula(paste0("out_come", " ~ ","Treatment +",
  #                                           paste0(colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                                                             select(-id,-out_come,-Treatment)), 
  #                                                  collapse = " + ")
  # )
  # ), data = df_fun
  # 
  # )
  # 
  # 
  # ATE_comp_godd_mod_AIPW <- PSweight(
  #   as.formula(paste0("Treatment", " ~ ",
  #                     paste0(c(colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                                       select(-id,-out_come,-Treatment)),
  #                              "confu_1 * confu_2 + confu_4 * confu_5 + confu_7 * confu_10"),
  #                            collapse = " + ")
  #   )
  #   ),
  #   data = df_fun,
  #   zname = "Treatment",
  #   yname = "out_come",
  #   augmentation = TRUE,
  #   out.formula = 
  #     as.formula(paste0("out_come", " ~ ",
  #                       paste0(c(colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                                         select(-id,-out_come,-Treatment)), 
  #                                "confu_1 * confu_2 + confu_4 * confu_5 + confu_7 * confu_10"),
  #                              collapse = " + ")
  #     )
  #     )
  # )
  # 
  # 
  # ATE_comp_godd_mod_adjustment = lm(
  #   as.formula(paste0("out_come", " ~ ","Treatment +",
  #                     paste0(c(colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                                         select(-id,-out_come,-Treatment)),  
  #                              "confu_1 * confu_2 + confu_4 * confu_5 + confu_7 * confu_10"),
  #                                                  collapse = " + ")
  # )
  # ), data = df_fun
  # 
  # )
  
  
  res <- data.frame(
    expo = ATE_iptw_msp_com$expo,
    ATE_IPTW_misspe = ATE_iptw_msp_com$ATE,
    ATE_IPTW_good_spe = ATE_iptw_godd_mod_com$ATE,
    ATE_GBM = ATE_gbm_comp$ATE,
    ATE_CBPS_misspe = ATE_comp_msp_CBPS$ATE,
    ATE_CBPS_good_spe = ATE_comp_godd_mod_CBPS$ATE#,
    # ATE_msp_AIPW = t(summary(ATE_comp_msp_AIPW)$estimates[c(1,2),1]),
    # ATE_msp_adjustment = t(summary(ATE_comp_msp_adjustment)$coefficients[c(2,3),1]) ,
    # ATE_godd_mod_AIPW = t(summary(ATE_comp_godd_mod_AIPW)$estimates[c(1,2),1]),
    # ATE_godd_mod_adjustment = t(summary(ATE_comp_godd_mod_adjustment)$coefficients[c(2,3),1])
  )
  
  
  return(res)
}

ATE_boot_noise <- function(df, index, weight_df) {

  df_fun <- df[index, ]
  
  
  # IPTW_noise <- calcule_pds_stage(
  #   df_fun,
  #   expo = Treatment,
  #   covar = colnames(Sim$df_fun[-c(3,6,8,16,18)] %>%
  #                      select(-id,-out_come,-Treatment)),
  #   Id_data = id,
  #   weighting_function = multinomial_IPTW,
  #   out_come = out_come
  # )
  # 
  # weight_iptw_noise <- 1/IPTW_noise$res_intermediaire$poids$poids_tronc$poids_trunc$`(0; 1)`
  # 
  # 
  # IPTW_noise_exact <- calcule_pds_stage(
  #   df_fun[,-c(21:23,27:29)],
  #   expo = Treatment,
  #   covar = colnames(df_fun[-c(3,6,8,16,18,21:23,27:29)] %>%
  #                      select(-id,-out_come,-Treatment)),
  #   # var_mod_pds ,
  #   Id_data = id,
  #   weighting_function = multinomial_IPTW,
  #   out_come = out_come # ,
  #   # percentile_tronc = 0.01
  # )
  # 
  # 
  # weight_iptw_noise_exact <-1/IPTW_noise_exact$res_intermediaire$poids$poids_tronc$poids_trunc$`(0; 1)`
  # 
  
  weight_iptw_noise <- weight_df$IPTW_res$IPTW_noise[index]
  weight_iptw_noise_exact <- weight_df$IPTW_res$IPTW_noise_exact[index]
  weight_gbm_noise <- weight_df$GBM_res$twang_res_ATE_noise[index]
  weight_CBPS_noise <- 1/weight_df$CBPS_res$CBPS_res_noise_misspe$weights[index]
  weight_CBPS_noise_exact <- 1/weight_df$CBPS_res$CBPS_res_noise_correct$weights[index]
  
  
  ATE_noise_IPTW <- ATE_fun(exposure = df_fun$Treatment,
                          outcome = df_fun$out_come,
                          weight = weight_iptw_noise)
  
  ATE_noise_IPTW_exact <- ATE_fun(exposure = df_fun$Treatment,
                            outcome = df_fun$out_come,
                            weight = weight_iptw_noise_exact)
  
  ATE_noise_GBM <- ATE_fun(exposure = df_fun$Treatment,
                         outcome = df_fun$out_come,
                         weight = weight_gbm_noise)
  
  ATE_noise_msp_CBPS <- ATE_fun(exposure = df_fun$Treatment,
                               outcome = df_fun$out_come,
                               weight = weight_CBPS_noise)
  
  ATE_noise_godd_mod_CBPS <- ATE_fun(exposure = df_fun$Treatment,
                                    outcome = df_fun$out_come,
                                    weight = weight_CBPS_noise_exact)  
  
  
  # 
  # ATE_noise_msp_AIPW <- PSweight(
  #   as.formula(paste0("Treatment", " ~ ",
  #                     paste0(colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                                       select(-id,-out_come,-Treatment)), 
  #                            collapse = " + ")
  #   )
  #   ),
  #   data = df_fun,
  #   zname = "Treatment",
  #   yname = "out_come",
  #   augmentation = TRUE,
  #   out.formula = 
  #     as.formula(paste0("out_come", " ~ ",
  #                       paste0(colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                                         select(-id,-out_come,-Treatment)), 
  #                              collapse = " + ")
  #     )
  #     )
  # )
  # 
  # 
  # ATE_noise_msp_adjustment = lm(as.formula(paste0("out_come", " ~ ","Treatment +",
  #                                                paste0(colnames(df_fun[-c(3,6,8,16,18)] %>% 
  #                                                                  select(-id,-out_come,-Treatment)), 
  #                                                       collapse = " + ")
  # )
  # ), data = df_fun
  # 
  # )
  # 
  # 
  # ATE_noise_godd_mod_AIPW <- PSweight(
  #   as.formula(paste0("Treatment", " ~ ",
  #                     paste0(colnames(df_fun[-c(3,6,8,16,18,21:23,27:29)] %>% 
  #                                         select(-id,-out_come,-Treatment)),
  #                            collapse = " + ")
  #   )
  #   ),
  #   data = df_fun,
  #   zname = "Treatment",
  #   yname = "out_come",
  #   augmentation = TRUE,
  #   out.formula = 
  #     as.formula(paste0("out_come", " ~ ",
  #                       paste0(colnames(df_fun[-c(3,6,8,16,18,21:23,27:29)] %>% 
  #                                           select(-id,-out_come,-Treatment)), 
  #                              collapse = " + ")
  #     )
  #     )
  # )
  # 
  # 
  # ATE_noise_godd_mod_adjustment = lm(
  #   as.formula(paste0("out_come", " ~ ","Treatment +",
  #                     paste0(colnames(df_fun[-c(3,6,8,16,18,21:23,27:29)] %>% 
  #                                         select(-id,-out_come,-Treatment)),  
  #                            collapse = " + ")
  #   )
  #   ), data = df_fun
  #   
  # )
  # 
  
  
  
  
  
  
  
  res <- data.frame(
    expo = ATE_noise_IPTW$expo,
    ATE_IPTW = ATE_noise_IPTW$ATE,
    ATE_IPTW_exact = ATE_noise_IPTW_exact$ATE,
    ATE_GBM = ATE_noise_GBM$ATE,
    ATE_CBPS = ATE_noise_msp_CBPS$ATE, 
    ATE_CBPS_exact = ATE_noise_godd_mod_CBPS$ATE#,
    # ATE_AIPW = t(summary(ATE_noise_msp_AIPW)$estimates[c(1,2),1]),
    # ATE_adjustment = t(summary(ATE_noise_msp_adjustment)$coefficients[c(2,3),1]) ,
    # ATE_AIPW_exact = t(summary(ATE_noise_godd_mod_AIPW)$estimates[c(1,2),1]),
    # ATE_adjustment_exact = t(summary(ATE_noise_godd_mod_adjustment)$coefficients[c(2,3),1])
    )
  return(res)
}


ATE_real <- function(df, weight_df) {
  df_fun <- df
  weight_iptw_noise <- weight_df$IPTW_res$IPTW_noise
  weight_iptw_noise_exact <- weight_df$IPTW_res$IPTW_noise_exact
  weight_gbm_noise <- weight_df$GBM_res$twang_res_ATE_noise
  weight_CBPS_noise <- 1/weight_df$CBPS_res$CBPS_res_noise_misspe$weights
  weight_CBPS_noise_exact <- 1/weight_df$CBPS_res$CBPS_res_noise_correct$weights
  
  
  
  weight_iptw_msp_com <- weight_df$IPTW_res$IPTW_comp
  weight_iptw_godd_mod_com <- weight_df$IPTW_res$IPTW_comp_exact_formula
  weight_gbm_comp <- weight_df$GBM_res$twang_res_ATE_compl
  weight_CBPS_msp_com <- 1/weight_df$CBPS_res$CBPS_res_comp_misspe$weights
  weight_CBPS_godd_mod_com <- 1/weight_df$CBPS_res$CBPS_res_comp_correct$weights
  
  
  
  
  weight_iptw_add <- weight_df$IPTW_res$IPTW_add
  weight_gbm_add <- weight_df$GBM_res$twang_res_ATE_add
  weight_CBPS_add <- 1/weight_df$CBPS_res$CBPS_res_add$weights
  

  
  
  
  add_ATE_IPTW <- ATE_fun(exposure = df_fun$df_add$Treatment,
                          outcome = df_fun$df_add$out_come,
                          weight = weight_iptw_add)
  
  add_ATE_GBM <- ATE_fun(exposure = df_fun$df_add$Treatment,
                         outcome = df_fun$df_add$out_come,
                         weight = weight_gbm_add)
  
  ATE_add_CBPS <- ATE_fun(exposure = df_fun$df_add$Treatment,
                          outcome = df_fun$df_add$out_come,
                          weight = weight_CBPS_add)
  
  
  ATE_add_AIPW <- PSweight(
    as.formula(paste0("Treatment", " ~ ",
                      paste0(colnames(df_fun$df_add[-c(3,6,8,16,18)] %>% 
                                        select(-id,-out_come,-Treatment)), 
                             collapse = " + ")
    )
    ),
    data = df_fun$df_add,
    zname = "Treatment",
    yname = "out_come",
    augmentation = TRUE,
    out.formula = 
      as.formula(paste0("out_come", " ~ ",
                        paste0(colnames(df_fun$df_add %>% 
                                          dplyr::select(-id,-out_come,-Treatment)), 
                               collapse = " + ")
      )
      )#,
    # bootstrap = TRUE,
    # R = 500
  )
  
  
  ATE_add_adjustment = lm(as.formula(paste0("out_come", " ~ ","Treatment +",
                                            paste0(colnames(df_fun$df_add %>% 
                                                              select(-id,-out_come,-Treatment)), 
                                                   collapse = " + ")
  )
  ), data = df_fun$df_add
  
  )
  
  
  
  
  noise_ATE_IPTW <- ATE_fun(exposure = df_fun$df_noise$Treatment,
                            outcome = df_fun$df_noise$out_come,
                            weight = weight_iptw_noise)

  noise_ATE_IPTW_exact <- ATE_fun(exposure = df_fun$df_noise$Treatment,
                                  outcome = df_fun$df_noise$out_come,
                                  weight = weight_iptw_noise_exact)
 
  noise_ATE_GBM <- ATE_fun(exposure = df_fun$df_noise$Treatment,
                           outcome = df_fun$df_noise$out_come,
                           weight = weight_gbm_noise)
  
  
  ATE_noise_msp_CBPS <- ATE_fun(exposure = df_fun$df_noise$Treatment,
                                outcome = df_fun$df_noise$out_come,
                                weight = weight_CBPS_noise)
  
  ATE_noise_godd_mod_CBPS <- ATE_fun(exposure = df_fun$df_noise$Treatment,
                                     outcome = df_fun$df_noise$out_come,
                                     weight = weight_CBPS_noise_exact)  
  
   
  
  
  ATE_noise_msp_AIPW <- PSweight(
    as.formula(paste0("Treatment", " ~ ",
                      paste0(colnames(df_fun$df_noise[-c(3,6,8,16,18)] %>% 
                                        select(-id,-out_come,-Treatment)), 
                             collapse = " + ")
    )
    ),
    data = df_fun$df_noise,
    zname = "Treatment",
    yname = "out_come",
    augmentation = TRUE,
    out.formula = 
      as.formula(paste0("out_come", " ~ ",
                        paste0(colnames(df_fun$df_noise %>% 
                                          select(-id,-out_come,-Treatment)), 
                               collapse = " + ")
      )
      )#,
    # bootstrap = TRUE,
    # R = 500
  )
  
  
  ATE_noise_msp_adjustment = lm(as.formula(paste0("out_come", " ~ ","Treatment +",
                                                  paste0(colnames(df_fun$df_noise %>% 
                                                                    select(-id,-out_come,-Treatment)), 
                                                         collapse = " + ")
  )
  ), data = df_fun$df_noise
  
  )
  
  
  ATE_noise_godd_mod_AIPW <- PSweight(
    as.formula(paste0("Treatment", " ~ ",
                      paste0(colnames(df_fun$df_noise[-c(3,6,8,16,18,21:23,27:29)] %>% 
                                        select(-id,-out_come,-Treatment)),
                             collapse = " + ")
    )
    ),
    data = df_fun$df_noise,
    zname = "Treatment",
    yname = "out_come",
    augmentation = TRUE,
    out.formula = 
      as.formula(paste0("out_come", " ~ ",
                        paste0(colnames(df_fun$df_noise[-c(21:23,27:29)] %>% 
                                          select(-id,-out_come,-Treatment)), 
                               collapse = " + ")
      )
      )#,
    # bootstrap = TRUE,
    # R = 500
  )
  
  
  ATE_noise_godd_mod_adjustment = lm(
    as.formula(paste0("out_come", " ~ ","Treatment +",
                      paste0(colnames(df_fun$df_noise[-c(21:23,27:29)] %>% 
                                        select(-id,-out_come,-Treatment)),  
                             collapse = " + ")
    )
    ), data = df_fun$df_noise
    
  )
  
  
  
  
  compl_ATE_IPTW_misspe <- ATE_fun(exposure = df_fun$df_compl$Treatment,
                              outcome = df_fun$df_compl$out_come,
                              weight = weight_iptw_msp_com)
  
  compl_ATE_IPTW_good_spe <- ATE_fun(exposure = df_fun$df_compl$Treatment, 
                                   outcome = df_fun$df_compl$out_come,
                                   weight = weight_iptw_godd_mod_com)
  
  compl_ATE_GBM <- ATE_fun(exposure = df_fun$df_compl$Treatment,
                          outcome = df_fun$df_compl$out_come,
                          weight = weight_gbm_comp)
  
  ATE_comp_msp_CBPS <- ATE_fun(exposure = df_fun$df_compl$Treatment,
                               outcome = df_fun$df_compl$out_come,
                               weight = weight_CBPS_msp_com)
  
  ATE_comp_godd_mod_CBPS <- ATE_fun(exposure = df_fun$df_compl$Treatment,
                                    outcome = df_fun$df_compl$out_come,
                                    weight = weight_CBPS_godd_mod_com)
  
  
  ATE_comp_msp_AIPW <- PSweight(
    as.formula(paste0("Treatment", " ~ ",
                      paste0(colnames(df_fun$df_compl[-c(3,6,8,16,18)] %>% 
                                        select(-id,-out_come,-Treatment)), 
                             collapse = " + ")
    )
    ),
    data = df_fun$df_compl,
    zname = "Treatment",
    yname = "out_come",
    augmentation = TRUE,
    out.formula = 
      as.formula(paste0("out_come", " ~ ",
                        paste0(colnames(df_fun$df_compl %>% 
                                          select(-id,-out_come,-Treatment)), 
                               collapse = " + ")
      )
      )#,
    # bootstrap = TRUE,
    # R = 500
  )
  
  
  ATE_comp_msp_adjustment = lm(as.formula(paste0("out_come", " ~ ","Treatment +",
                                                 paste0(colnames(df_fun$df_compl %>% 
                                                                   select(-id,-out_come,-Treatment)), 
                                                        collapse = " + ")
  )
  ), data = df_fun$df_compl
  
  )
  
  
  ATE_comp_godd_mod_AIPW <- PSweight(
    as.formula(paste0("Treatment", " ~ ",
                      paste0(c(colnames(df_fun$df_compl[-c(3,6,8,16,18)] %>% 
                                          select(-id,-out_come,-Treatment)),
                               "confu_1 * confu_2 + confu_4 * confu_5 + confu_7 * confu_10"),
                             collapse = " + ")
    )
    ),
    data = df_fun$df_compl,
    zname = "Treatment",
    yname = "out_come",
    augmentation = TRUE,
    out.formula = 
      as.formula(paste0("out_come", " ~ ",
                        paste0(c(colnames(df_fun$df_compl %>% 
                                            select(-id,-out_come,-Treatment)), 
                                 "confu_1 * confu_2 + confu_4 * confu_5 + confu_7 * confu_10"),
                               collapse = " + ")
      )
      )# ,
    # bootstrap = TRUE,
    # R = 500
  )
  
  
  ATE_comp_godd_mod_adjustment = lm(
    as.formula(paste0("out_come", " ~ ","Treatment +",
                      paste0(c(colnames(df_fun$df_compl %>% 
                                          select(-id,-out_come,-Treatment)),  
                               "confu_1 * confu_2 + confu_4 * confu_5 + confu_7 * confu_10"),
                             collapse = " + ")
    )
    ), data = df_fun$df_compl
    
  )
  

  res <- list(
    IPTW_add = add_ATE_IPTW,
    GBM_add = add_ATE_GBM,
    CBPS_add = ATE_add_CBPS,
    AIPW_add = ATE_add_AIPW,
    adjustment_add = ATE_add_adjustment,
    
    IPTW_msp_compl = compl_ATE_IPTW_misspe,
    IPTW_exact_compl = compl_ATE_IPTW_good_spe,
    GBM_compl = compl_ATE_GBM,
    CBPS_msp_compl = ATE_comp_msp_CBPS,
    CBPS_exact_compl = ATE_comp_godd_mod_CBPS,
    AIPW_msp_compl = ATE_comp_msp_AIPW,
    adjustment_msp_compl = ATE_comp_msp_adjustment,
    AIPW_exact_compl = ATE_comp_godd_mod_AIPW,
    adjustment_exact_compl = ATE_comp_godd_mod_adjustment,
    
    IPTW_noise = noise_ATE_IPTW,
    IPTW_exact_noise = noise_ATE_IPTW_exact,
    GBM_noise = noise_ATE_GBM,
    CBPS_msp_noise = ATE_noise_msp_CBPS,
    CBPS_exact_noise = ATE_noise_godd_mod_CBPS,
    AIPW_msp_noise = ATE_noise_msp_AIPW,
    adjustment_msp_noise = ATE_noise_msp_adjustment,
    AIPW_exact_noise = ATE_noise_godd_mod_AIPW,
    adjustment_exact_noise = ATE_noise_godd_mod_adjustment
    
  )
  return(res)
}




multiple_simu <- function(i,n,nom_sim) {
  dir.create(paste0("data/genere/simu/simu_",nom_sim,"/"), 
             showWarnings = FALSE,
             recursive = TRUE)
  
  
  
  Simulated_df <- simulation_fun(n)
  Weight_table <- PS_estimation(Simulated_df)
  
  # browser()
  # 
  # library(boot)
  # 
  # test <- function( data, indices) ATE_boot_add( data, indices,Weight_table )
  # 
  # a <- boot(data = Simulated_df$df_add,
  #            statistic = test,
  #            R =1000)
  # 
  # 
  # boot.ci( a,type ="norm" ,index =4)$norm[-1] 
  # 
  # plot(a,index=2)
  
  boot_strap_full_process_add <- boot_strap_fun(Simulated_df$df_add, 
                                                ATE_boot_add,
                                                500,
                                                Weight_table)
  print("boot_strap_add")
  print(Sys.time())
  boot_strap_full_process_comp <- boot_strap_fun(Simulated_df$df_compl, 
                                                 ATE_boot_compl,
                                                 500,
                                                 Weight_table)
  print("boot_strap_comp")
  print(Sys.time())
  boot_strap_full_process_noise <- boot_strap_fun(Simulated_df$df_noise, 
                                                 ATE_boot_noise,
                                                 500,
                                                 Weight_table)
  
  print("boot_strap_noise")
  print(Sys.time())
  

  ATE_real_res <- ATE_real(Simulated_df,Weight_table)
  
  
  
  print("boot_strap")
  print(Sys.time())
  boot_strap_full_process <- list(
    ATE_add = boot_strap_full_process_add,
    ATE_comp = boot_strap_full_process_comp,
    ATE_noise = boot_strap_full_process_noise,
    data_sim = Simulated_df,
    Weight_table = Weight_table,
    ATE_real = ATE_real_res
  )
  
  
  
  
  saveRDS(boot_strap_full_process, file = paste0("data/genere/simu/simu_",nom_sim,"/bott_strap_full_process_", i, ".rds"))
}
