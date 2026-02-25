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