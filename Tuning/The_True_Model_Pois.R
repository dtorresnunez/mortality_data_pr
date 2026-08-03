################################################
### Modelo para mortalidad en áreas pequeñas ###
################################################
# 2026_07_28

# Tareas:
# David. Repetir el análisis utilizando regiones senatoriales.
#    - Hablar con el Dr. Pericchi.
# HOLD. Optimizar el análisis de sensibilidad.

# Paquetes
library(INLA)
library(SUMMER)
library(tidyverse)
library(ggh4x)
library(dplyr)
library(tidyr)
library(tidycensus)
library(haven)
library(sf)
library(this.path)
library(DemoTools)
library(demogR)
library(patchwork) 
library(MortalityEstimate)
library(MortCast) 
library(MortalityLaws) 
library(epitools)
library(PHEindicatormethods)
library(demography)
library(forecast)
library(readxl)
library(readr)
library(parallel)

script_dir         <- this.path::this.dir()
data_dir           <- file.path(script_dir, "data")
carpeta_resultados <- file.path(script_dir, "resultados_modelos/Poisson_Model")
dir.create(carpeta_resultados, recursive = TRUE, showWarnings = FALSE)

# Cargar la base de datos de población y muerte
df         <- read_csv(file.path(data_dir, "data_frame_population_deaths.csv"),
                       col_types = cols(.default = col_guess()))

# Cargar la mtriz de adyacencia
Amat       <- as.matrix(read.csv(file.path(data_dir, "adjacency_matrix.csv"),
                                 check.names = FALSE))

# Parámetros quinquenales y grupos de edad 
Age        <- c(0, 1, seq(5, 85, by = 5))
ages       <- c(
  "0", "01-04","05-09", "10-14", "15-19", "20-24",
  "25-29", "30-34", "35-39", "40-44", "45-49",
  "50-54", "55-59", "60-64", "65-69", "70-74",
  "75-79", "80-84", "85+"
)
age_params <- tibble(
  agegroup = ages,
  n_interval = c(1, 4, rep(5, 16), NA),
  ax = c(
    0.15, 1.5, 2.5, 2.5, 2.5,
    2.5, 2.5, 2.5, 2.5, 2.5,
    2.5, 2.5, 2.5, 2.5, 2.5,
    2.5, 2.5, 2.5, NA
  )
)

# Convertir la matriz de adyacencia en una matriz para INLA
g          <- INLA::inla.read.graph(Amat)

# Definir la previa
SB2.prior <- function(p = 1, q = 1, b = 1){
  sprintf(
    "expression:
      p = %f;
      q = %f;
      b = %f;
      sigma2 = exp(-theta);
      log_dens = lgamma(p+q) - lgamma(p) - lgamma(q) - log(b);
      log_dens = log_dens + (p-1) * log(sigma2/b);
      log_dens = log_dens - (p+q) * log(1 + sigma2/b);
      log_dens = log_dens - theta;
      return(log_dens);
    ",
    p, q, b
  )
}

# Ejecutar el modelo INLA. Funciona perfecto para Mac
calcular_e0_inla     <- function(modelo_inla, df, age_params, Age, nsamples = 1000) {
  datos_directo <- df %>%
    mutate(mx = pmax(deaths / population, 1e-6),
           sex_chr = ifelse(sex == 1, "m", "f")) %>%
    left_join(age_params, by = "agegroup")
  
  e0_obs_list <- list()
  for (reg in unique(datos_directo$region)) {
    for (per in unique(datos_directo$period)) {
      for (sx in c("m", "f")) {
        sub <- datos_directo %>%
          filter(region == reg, period == per, sex_chr == sx) %>%
          arrange(age_idx)
        nMx <- sub$mx
        if (length(nMx) <= 5) next
        AgeInt <- inferAgeIntAbr(vec = nMx)
        ff <- Age[Age >= 60 & Age < max(Age) & sub$deaths > 0]
        if (length(ff) < 2) ff <- Age[Age >= 60 & Age < max(Age)]
        tb <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, Sex = sx,
                          a0rule = "ak", axmethod = "pas", mod = FALSE, extrapLaw   = "kannisto", extrapFrom  = 80, extrapFit = ff)
        e0_obs_list[[length(e0_obs_list) + 1]] <- data.frame(
          region = reg, period = per, sex = ifelse(sx == "m", 1, 2),
          e0_observado = tb$ex[1]
        )
      }
    }
  }
  e0_observado_df <- bind_rows(e0_obs_list)
  
  # Muestras posteriores del predictor
  set.seed(123)
  samples <- inla.posterior.sample(nsamples, modelo_inla, seed = 123)
  
  log_lambda_matrix_all <- inla.posterior.sample.eval(
    function(...) { Predictor }, 
    samples
  )
  log_lambda_matrix <- log_lambda_matrix_all[1:nrow(df), , drop = FALSE]
  mx_matrix <- pmax(exp(log_lambda_matrix), 1e-6)
  
  # e0 estimado por muestra 
  datos_idx <- df %>%
    mutate(sex_chr = ifelse(sex == 1, "m", "f")) %>%
    left_join(age_params, by = "agegroup")
  combinaciones <- datos_idx %>% distinct(region, period, sex_chr)
  e0_sim_list <- vector("list", nsamples * nrow(combinaciones))
  contador <- 0
  for (s in seq_len(nsamples)) {
    datos_idx$mx <- mx_matrix[, s]
    for (i in seq_len(nrow(combinaciones))) {
      reg <- combinaciones$region[i]
      per <- combinaciones$period[i]
      sx  <- combinaciones$sex_chr[i]
      sub <- datos_idx %>%
        filter(region == reg, period == per, sex_chr == sx) %>%
        arrange(age_idx)
      nMx <- sub$mx
      if (length(nMx) <= 5) next
      AgeInt <- inferAgeIntAbr(vec = nMx)
      ff <- Age[Age >= 60 & Age < max(Age) & sub$deaths > 0]
      if (length(ff) < 2) ff <- Age[Age >= 60 & Age < max(Age)]
      tb <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, Sex = sx,
                        a0rule = "ak", axmethod = "pas", mod = FALSE, extrapLaw   = "kannisto", extrapFrom  = 80, extrapFit = ff)
      contador <- contador + 1
      e0_sim_list[[contador]] <- data.frame(
        sim = s, region = reg, period = per,
        sex = ifelse(sx == "m", 1, 2), e0 = tb$ex[1]
      )
    }
    if (s %% 100 == 0) message("Muestra ", s, " de ", nsamples)
  }
  e0_sim_df <- bind_rows(e0_sim_list)
  
  # mediana + IC 95%
  e0_estimado_df <- e0_sim_df %>%
    group_by(region, period, sex) %>%
    summarise(
      e0_estimado = median(e0, na.rm = TRUE),
      e0_lower    = quantile(e0, 0.025, na.rm = TRUE),
      e0_upper    = quantile(e0, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  # e0 observado y estimado
  e0_final <- left_join(e0_observado_df, e0_estimado_df, by = c("region", "period", "sex")) %>%
    arrange(region, period, sex) %>%
    mutate(est_eval = case_when(
      between(e0_observado, e0_lower, e0_upper) ~ "Estimación adecuada",
      e0_observado < e0_lower ~ "> e0 observado",
      e0_observado > e0_upper ~ "< e0 observado"
    ))
  
  return(e0_final)
}

# Ejecutar el modelo INLA optimizado. Funciona perfecto para Windows (ajustar mc.cores)
num.cores <- detectCores(logical = T)
calcular_e0_inla_opt <- function(modelo_inla, df, age_params, Age, nsamples = 1000, mc.cores = num.cores - 1, ...){
  
  # --- preparación única ----------------------------------------------------
  datos_idx <- df %>%
    mutate(sex_chr = ifelse(sex == 1, "m", "f")) %>%
    left_join(age_params, by = "agegroup") %>%
    mutate(.fila = row_number())          # alineado con las filas de df y
  # por tanto con mx_matrix
  stopifnot(nrow(datos_idx) == nrow(df))  # el join no debe duplicar filas
  
  # esto es el "group_by en vez de los for anidados": índices por grupo,
  # ya ordenados por edad, calculados una sola vez
  grupos <- datos_idx %>%
    arrange(region, period, sex_chr, age_idx) %>%
    group_by(region, period, sex_chr) %>%
    summarise(idx = list(.fila), .groups = "drop") %>%
    filter(lengths(idx) > 5)              # mismo criterio que el original
  
  G           <- nrow(grupos)
  idx_list    <- grupos$idx
  sex_list    <- grupos$sex_chr
  mx_obs      <- pmax(df$deaths / df$population, 1e-6)
  # AgeInt solo depende del largo del vector; una vez por grupo basta
  AgeInt_list <- lapply(idx_list, function(ix) inferAgeIntAbr(vec = mx_obs[ix]))
  # edades del ajuste Kannisto: dependen solo de deaths observadas; una vez por grupo
  extrapFit_list <- lapply(idx_list, function(ix) {
    ff <- Age[Age >= 60 & Age < max(Age) & df$deaths[ix] > 0]
    if (length(ff) < 2) ff <- Age[Age >= 60 & Age < max(Age)]  # Kannisto necesita >= 2 puntos
    ff
  })
  
  # función auxiliar: e0 de un grupo dado un vector de mx (df completo)
  e0_grupo <- function(g, mx_vec) {
    lt_abridged(nMx = mx_vec[idx_list[[g]]], AgeInt = AgeInt_list[[g]],
                Age = Age, Sex = sex_list[g],
                a0rule = "ak", axmethod = "pas", mod = FALSE,
                extrapLaw = "kannisto", extrapFrom = 80,
                extrapFit = extrapFit_list[[g]])$ex[1]
  }
  
  # --- e0 observado ---------------------------------------------------------
  e0_obs <- vapply(seq_len(G), e0_grupo, numeric(1), mx_vec = mx_obs)
  
  # --- muestras posteriores del predictor -----------------------------------
  set.seed(123)
  samples <- inla.posterior.sample(nsamples, modelo_inla, seed = 123, ...)
  log_lambda_matrix <- inla.posterior.sample.eval(
    function(...) { Predictor },
    samples
  )[seq_len(nrow(df)), , drop = FALSE]
  mx_matrix <- pmax(exp(log_lambda_matrix), 1e-6)
  
  # --- e0 estimado por muestra: matriz preasignada, sin dplyr en el bucle ---
  e0_una_muestra <- function(s) {
    vapply(seq_len(G), e0_grupo, numeric(1), mx_vec = mx_matrix[, s])
  }
  
  if (mc.cores > 1 && .Platform$OS.type == "windows") {
    # Windows no tiene fork: se usa un cluster PSOCK. Cada worker recibe
    # SOLO su bloque de columnas de mx_matrix (no la matriz completa).
    bloques  <- split(seq_len(nsamples), sort(rep_len(seq_len(mc.cores), nsamples)))
    sub_mats <- lapply(bloques, function(ss) mx_matrix[, ss, drop = FALSE])
    
    trabajador <- function(subm, idx_list, AgeInt_list, sex_list, Age, extrapFit_list) {
      G <- length(idx_list)
      apply(subm, 2, function(mx_vec) {
        vapply(seq_len(G), function(g) {
          DemoTools::lt_abridged(nMx = mx_vec[idx_list[[g]]],
                                 AgeInt = AgeInt_list[[g]],
                                 Age = Age, Sex = sex_list[g],
                                 a0rule = "ak", axmethod = "pas",
                                 mod = FALSE, extrapLaw = "kannisto",
                                 extrapFrom = 80,
                                 extrapFit = extrapFit_list[[g]])$ex[1]
        }, numeric(1))
      })
    }
    # entorno limpio: evita serializar mx_matrix/samples completos a cada worker
    environment(trabajador) <- globalenv()
    
    cl <- parallel::makeCluster(mc.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    res <- parallel::parLapply(cl, sub_mats, trabajador,
                               idx_list = idx_list, AgeInt_list = AgeInt_list,
                               sex_list = sex_list, Age = Age,
                               extrapFit_list = extrapFit_list)
    e0_mat <- do.call(cbind, res)   # bloques contiguos -> orden original
    
  } else if (mc.cores > 1) {
    # Mac / Linux: fork con mclapply (sin copia de datos)
    cols <- parallel::mclapply(seq_len(nsamples), e0_una_muestra,
                               mc.cores = mc.cores)
    err <- vapply(cols, inherits, logical(1), what = "try-error")
    if (any(err)) stop("Fallaron ", sum(err), " muestras en mclapply.")
    e0_mat <- do.call(cbind, cols)
  } else {
    e0_mat <- matrix(NA_real_, nrow = G, ncol = nsamples)
    for (s in seq_len(nsamples)) {
      e0_mat[, s] <- e0_una_muestra(s)
      if (s %% 100 == 0) message("Muestra ", s, " de ", nsamples)
    }
  }
  
  # --- resumen: equivalente al group_by + summarise original ----------------
  grupos %>%
    transmute(
      region, period,
      sex          = ifelse(sex_chr == "m", 1, 2),
      e0_observado = e0_obs,
      e0_estimado  = apply(e0_mat, 1, median, na.rm = TRUE),
      e0_lower     = unname(apply(e0_mat, 1, quantile, probs = 0.025, na.rm = TRUE)),
      e0_upper     = unname(apply(e0_mat, 1, quantile, probs = 0.975, na.rm = TRUE))
    ) %>%
    arrange(region, period, sex) %>%
    mutate(est_eval = case_when(
      between(e0_observado, e0_lower, e0_upper) ~ "Estimación adecuada",
      e0_observado < e0_lower ~ "> e0 observado",
      e0_observado > e0_upper ~ "< e0 observado"
    ))
}

# Graficar el e0 observado, estimado y los IC
e0_model_plot        <- function(dat, per, col, llh) {
  d <- dat %>% filter(period == per) %>%
    mutate(region2 = fct_reorder(paste(region, sex, sep = "___"), e0_observado))   # 1
  ggplot(d,
         aes(x = region2,                                                          # 2
             y = e0_estimado, ymin = e0_lower, ymax = e0_upper)) +
    geom_pointrange(color = col, size = 0.3) +
    geom_point(aes(y = e0_observado)) +
    coord_flip() +
    scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +                  # 3
    facet_wrap(sex ~ est_eval, scales = "free",
               labeller = labeller(sex = c(`1` = "Hombres", `2` = "Mujeres"))) +
    ggh4x::facetted_pos_scales(
      y = list(
        sex == "1" ~ scale_y_continuous(
          limits = c(60, 90),
          breaks = seq(60, 90, by = 1)
        ),
        sex == "2" ~ scale_y_continuous(
          limits = c(60, 90),
          breaks = seq(60, 90, by = 1)
        )
      )
    ) + 
    theme_minimal() +
    labs(title = paste0("e0 por mun. (est. vs. obs.), ", per, ", ", llh),
         y = "e0", x = "") +
    theme(axis.title.x = element_text(size = 6))
}

#######################
# Modo de uso de INLA #
#######################

# formula_sb2 <- deaths ~
#   factor(sex) +
#   f(age_idx,    model = "rw1",  constr = TRUE,
#     hyper = list(prec = list(prior = SB2.prior(1, 1, 1)))) +
#   f(region_idx, model = "bym2", graph = g, constr = TRUE,
#     hyper = list(prec = list(prior = SB2.prior(1, 1, 0.5)),
#                  phi  = list(prior = "logitbeta", param = c(0.5, 0.5)))) +
#   f(period_idx, model = "rw2",  constr = TRUE,
#     hyper = list(prec = list(prior = SB2.prior(1, 1, 0.25)))) +
#   f(region_period_idx, model = "iid",
#     hyper = list(prec = list(prior = SB2.prior(1, 1, 0.1)))) +
#   f(celda_idx, model = "iid",
#     hyper = list(prec = list(prior = SB2.prior(1, 1, 0.1))))
# 
# fit_sb2 <- inla(formula_sb2,
#                 family = "poisson",
#                 data = df,
#                 E = population,
#                 control.compute = list(config = TRUE, dic = TRUE, waic = TRUE))

# Definición de parámetros iniciales
df_ambos   <- df
df_hombres <- df %>% filter(sex == 1)
df_mujeres <- df %>% filter(sex == 2)
familias   <- names(INLA::inla.models()$likelihood)
modelos    <- names(INLA::inla.models()$latent)

# Modelo completo. Más adelante está el ejemplo de uso
modelo_completo <- function(
    tabla_df,
    familia,
    model_age,
    par_p_age,
    par_q_age,
    par_b_age,
    model_reg,
    par_p_reg,
    par_q_reg,
    par_b_reg,
    model_per,
    par_p_per,
    par_q_per,
    par_b_per,
    model_s_t,
    par_p_s_t,
    par_q_s_t,
    par_b_s_t,
    model_cel,
    par_p_cel,
    par_q_cel,
    par_b_cel,
    nsamples,
    guardar    = TRUE)
{
  familia <- match.arg(familia, c("poisson", "nbinomial"))
  # Etiqueta usada para nombrar los archivos generados
  nombre_modelo <- paste(
    "gru", tabla_df,
    "fam", familia,
    "age", model_age, par_p_age, par_q_age, par_b_age,
    "reg", model_reg, par_p_reg, par_q_reg, par_b_reg,
    "per", model_per, par_p_per, par_q_per, par_b_per,
    "s_t", model_s_t, par_p_s_t, par_q_s_t, par_b_s_t,
    "cel", model_cel, par_p_cel, par_q_cel, par_b_cel,
    "sam", nsamples,
    sep = "_"
  )
  
  # Definir la fórmula para INLA
  
  formula_sb2 <- deaths ~
    factor(sex):period + #nuevo cambio: efecto de interacción sexo y período
    f(age_idx, model = model_age, constr = TRUE,
      hyper = list(prec = list(prior = SB2.prior(par_p_age, par_q_age, par_b_age)))) +
    f(region_idx, model = model_reg, graph = g, constr = TRUE,
      hyper = list(prec = list(prior = SB2.prior(par_p_reg , par_q_reg , par_b_reg)),
                   phi = list(prior = "logitbeta", param = c(0.5, 0.5)))) +
    f(period_idx, model = model_per, constr = TRUE,
      hyper = list(prec = list(prior = SB2.prior(par_p_per, par_q_per, par_b_per)))) +
    f(region_period_idx, model = model_s_t,
      hyper = list(prec = list(prior = SB2.prior(par_p_s_t, par_q_s_t, par_b_s_t)))) +
    f(cell_idx, model = model_cel,
      hyper = list(prec = list(prior = SB2.prior(par_p_cel, par_q_cel, par_b_cel))))

  # formula_sb2 <- deaths ~
  #   factor(sex):period + #nuevo cambio: efecto de interacción sexo y período
  #   f(age_idx, model = model_age, constr = TRUE,
  #     hyper = list(prec = list(prior = SB2.prior(par_p_age, par_q_age, par_b_age)))) +
  #   f(region_idx, model = model_reg, graph = g, constr = TRUE,
  #     hyper = list(prec = list(prior = SB2.prior(par_p_reg , par_q_reg , par_b_reg)),
  #                  phi = list(prior = "logitbeta", param = c(0.5, 0.5)))) +
  #   f(period_idx, model = model_per, constr = TRUE,
  #     hyper = list(prec = list(prior = SB2.prior(par_p_per, par_q_per, par_b_per)))) +
  #   f(region_period_idx, model = model_s_t,
  #     hyper = list(prec = list(prior = SB2.prior(par_p_s_t, par_q_s_t, par_b_s_t))))

  #No descomentar - formula anterior
  # formula_sb2 <- deaths ~
  #   factor(sex) +
  #   f(age_idx, model = model_age, constr = TRUE,
  #     hyper = list(prec = list(prior = SB2.prior(par_p_age, par_q_age, par_b_age)))) +
  #   f(region_idx, model = model_reg, graph = g, constr = TRUE,
  #     hyper = list(prec = list(prior = SB2.prior(par_p_reg , par_q_reg , par_b_reg)),
  #                  phi = list(prior = "logitbeta", param = c(0.5, 0.5)))) +
  #   f(period_idx, model = model_per, constr = TRUE,
  #     hyper = list(prec = list(prior = SB2.prior(par_p_per, par_q_per, par_b_per)))) +
  #   f(region_period_idx, model = model_s_t,
  #     hyper = list(prec = list(prior = SB2.prior(par_p_s_t, par_q_s_t, par_b_s_t)))) +
  #   f(cell_idx, model = model_cel,
  #     hyper = list(prec = list(prior = SB2.prior(par_p_cel, par_q_cel, par_b_cel))))
  
  formula_h <- deaths ~
    f(age_idx, model = model_age, constr = TRUE,
      hyper = list(prec = list(prior = SB2.prior(par_p_age, par_q_age, par_b_age)))) +
    f(region_idx, model = model_reg, graph = g, constr = TRUE,
      hyper = list(prec = list(prior = SB2.prior(par_p_reg , par_q_reg , par_b_reg)),
                   phi = list(prior = "logitbeta", param = c(0.5, 0.5)))) +
    f(period_idx, model = model_per, constr = TRUE,
      hyper = list(prec = list(prior = SB2.prior(par_p_per, par_q_per, par_b_per)))) +
    f(region_period_idx, model = model_s_t,
      hyper = list(prec = list(prior = SB2.prior(par_p_s_t, par_q_s_t, par_b_s_t)))) +
    f(cell_idx, model = model_cel,
      hyper = list(prec = list(prior = SB2.prior(par_p_cel, par_q_cel, par_b_cel))))
  
  formula_m <- formula_h
  
  nombre_formula <- c(
    ambos   = "formula_sb2",
    hombres = "formula_h",
    mujeres = "formula_m"
  )[tabla_df]
  
  datos_modelo <- get(
    paste0("df_", tabla_df)
  )
  
  formula_modelo <- get(
    unname(nombre_formula)
  )
  
  #OJO
  
  # Descomentar y Ejecutar la fórmula para INLA bajo la familia Binomial Negativa
  if(familia == "nbinomial"){
    fit_sb2_nbinom <- inla(formula_modelo,
                           family = familia,
                           control.family = list(
                             hyper = list(
                               theta = list(
                                 prior="normal",
                                 param=c(log(10),0.5)  # Genera una previa normal centrada en log(10)
                                 # c(media, precisión) = c(log(10), 2) \approx (2.30, 2) 
                               )
                             )
                           ),
                           data = datos_modelo,
                           E = datos_modelo$population,
                           control.compute = list(config = TRUE, dic = TRUE, waic = TRUE))}
  
  # Descomentar y Ejecutar la fórmula para INLA bajo la familia Poisson
  if(familia == "poisson"){
    fit_sb2_pois <- inla(formula_modelo,
                         family = familia,
                         data = datos_modelo,
                         E = datos_modelo$population,
                         control.compute = list(config = TRUE, dic = TRUE, waic = TRUE))}
  
  fit_sb2_mod <- c(
    nbinomial   = "fit_sb2_nbinom",
    poisson = "fit_sb2_pois"
  )[familia]
  
  fit_sb2 <- get(unname(fit_sb2_mod))
  
  if(familia == "nbinomial"){
    hyper_nbinomial <- fit_sb2$.args$control.family[[1]]$hyper$theta$param
    hyper_familia   <- paste(round(hyper_nbinomial, 3), collapse = "_")
    nombre_modelo <- paste(
      "gru", tabla_df,
      "fam", familia,
      "hyp", hyper_familia,
      "age", model_age, par_p_age, par_q_age, par_b_age,
      "reg", model_reg, par_p_reg, par_q_reg, par_b_reg,
      "per", model_per, par_p_per, par_q_per, par_b_per,
      "s_t", model_s_t, par_p_s_t, par_q_s_t, par_b_s_t,
      "cel", model_cel, par_p_cel, par_q_cel, par_b_cel,
      "sam", nsamples,
      sep = "_"
    )
  }
  
  # Ejecutar las muestras por cada modelo
  tiempo_ajuste <- fit_sb2$cpu.used[["Total"]]
  tiempo_e0 <- system.time(modelo_final_con <- calcular_e0_inla_opt(fit_sb2,
                                                                    datos_modelo,
                                                                    age_params,
                                                                    Age, nsamples = nsamples))[["elapsed"]]
  
  # Elegir el período
  periodo   <- unique(modelo_final_con$period)
  
  # Graficar los IC
  grafica_e0 <- periodo |>
    purrr::map(
      \(per_actual) {
        e0_model_plot(
          dat = modelo_final_con,
          per = per_actual,
          col = "purple",
          llh = nombre_modelo
        )
      }
    ) |>
    rlang::set_names(periodo)
  
  # Tabular la cobertura segmentada
  tabla_cobertura <- modelo_final_con %>%
    group_by(period, sex) %>%
    summarise(
      pct_dentro             = 100 * mean(est_eval == "Estimación adecuada", na.rm = TRUE),
      pct_dentro_sobreestima = 100 * mean(est_eval == "Estimación adecuada" &
                                            e0_estimado > e0_observado, na.rm = TRUE),
      pct_dentro_subestima   = 100 * mean(est_eval == "Estimación adecuada" &
                                            e0_estimado < e0_observado, na.rm = TRUE),
      pct_fuera_sobreestima  = 100 * mean(est_eval == "> e0 observado", na.rm = TRUE),
      pct_fuera_subestima    = 100 * mean(est_eval == "< e0 observado", na.rm = TRUE),
      .groups = "drop"
    )
  
  # Tabular la cobertura completa
  tabla_cobertura_completa <- modelo_final_con %>%
    group_by(period, sex) %>%
    summarise(cobertura = 100 * mean(est_eval == "Estimación adecuada", na.rm = TRUE),
              .groups = "drop")
  
  # Tabular la cobertura completa por sexo (totales por sexo)
  tabla_cobertura_completa_sexo <- modelo_final_con %>%
    group_by(sex) %>%
    summarise(cobertura = 100 * mean(est_eval == "Estimación adecuada", na.rm = TRUE),
              .groups = "drop")
  
  ## begin EGR
  xwalk_reg <- tibble::tibble(
    reg_code = c(72001,72003,72005,72007,72009,72011,72013,72015,72017,72019,
                 72021,72023,72025,72027,72029,72031,72033,72035,72037,72039,
                 72041,72043,72045,72047,72049,72051,72053,72054,72055,72057,
                 72059,72061,72063,72065,72067,72069,72071,72073,72075,72077,
                 72079,72081,72083,72085,72087,72089,72091,72093,72095,72097,
                 72099,72101,72103,72105,72107,72109,72111,72113,72115,72117,
                 72119,72121,72123,72125,72127,72129,72131,72133,72135,72137,
                 72139,72141,72143,72145,72147,72149,72151,72153),
    name = c("Adjuntas","Aguada","Aguadilla","Aguas Buenas","Aibonito",
             "Añasco","Arecibo","Arroyo","Barceloneta","Barranquitas",
             "Bayamón","Cabo Rojo","Caguas","Camuy","Canóvanas","Carolina",
             "Cataño","Cayey","Ceiba","Ciales","Cidra","Coamo","Comerío",
             "Corozal","Culebra","Dorado","Fajardo","Florida","Guánica",
             "Guayama","Guayanilla","Guaynabo","Gurabo","Hatillo",
             "Hormigueros","Humacao","Isabela","Jayuya","Juana Díaz","Juncos",
             "Lajas","Lares","Las Marías","Las Piedras","Loíza","Luquillo",
             "Manatí","Maricao","Maunabo","Mayagüez","Moca","Morovis",
             "Naguabo","Naranjito","Orocovis","Patillas","Peñuelas","Ponce",
             "Quebradillas","Rincón","Río Grande","Sabana Grande","Salinas",
             "San Germán","San Juan","San Lorenzo","San Sebastián",
             "Santa Isabel","Toa Alta","Toa Baja","Trujillo Alto","Utuado",
             "Vega Alta","Vega Baja","Vieques","Villalba","Yabucoa","Yauco")
  ) %>%
    mutate(region = chartr("áéíóúüñ", "aeiouun", name),
           country_code = 630, include_code = 2, orden = 0)
  
  stopifnot(all(unique(datos_modelo$region) %in% xwalk_reg$region))
  
  plantilla_reg <- bind_rows(
    xwalk_reg,
    tibble::tibble(reg_code = 630, name = "Puerto Rico", region = "Puerto Rico",
                   country_code = 630, include_code = 0, orden = 1)
  ) %>% arrange(orden, reg_code)
  
  pred_sb2 <- datos_modelo %>%
    mutate(
      mx = pmax(fit_sb2$summary.fitted.values$mean, 1e-6)
    ) %>%
    left_join(age_params, by = "agegroup")
  pred_sb2$sex <- ifelse(pred_sb2$sex == 1, "m", "f")
  
  pred_sb2 <- bind_rows(
    pred_sb2,
    pred_sb2 %>%
      group_by(period, sex, agegroup, age_idx) %>%
      summarise(mx = sum(mx * population) / sum(population),
                population = sum(population), .groups = "drop") %>%
      mutate(region = "Puerto Rico")
  ) %>%
    arrange(region, period, sex, age_idx)
  
  municipios <- sort(unique(pred_sb2$region))
  periodos   <- sort(unique(pred_sb2$period))
  sexos      <- intersect(c("m", "f"), unique(pred_sb2$sex))
  
  tablas <- list()
  for (muni in municipios) {
    for (per in periodos) {
      for (sx in sexos) {
        pred_sub_sb2 <- pred_sb2 %>%
          filter(region == muni, period == per, sex == sx)
        nMx    <- pred_sub_sb2$mx
        AgeInt <- inferAgeIntAbr(vec = nMx)
        tablas[[muni]][[per]][[sx]] <- lt_abridged(nMx = nMx, AgeInt = AgeInt,
                                                   Age = Age, a0rule = "ak",
                                                   axmethod = "pas",
                                                   Sex = sx, mod = FALSE,
                                                   extrapLaw = "kannisto", extrapFrom = 80,
                                                   extrapFit = Age[Age >= 60 & Age < max(Age)])
      }
    }
  }
  
  e0_resumen_sb2 <- data.frame()
  tablas_vida    <- data.frame()
  for (m in names(tablas)) {
    for (p in names(tablas[[m]])) {
      for (s in sexos) {
        tb     <- tablas[[m]][[p]][[s]]
        sexnum <- if (s == "m") 1 else 2
        e0_resumen_sb2 <- rbind(
          e0_resumen_sb2,
          data.frame(period = p, region = m, sex = sexnum, e0 = tb$ex[1])
        )
        tablas_vida <- rbind(
          tablas_vida,
          data.frame(region = m, period = p, sex = sexnum, tb)
        )
      }
    }
  }
  e0_resumen <- e0_resumen_sb2 %>% arrange(region, period, sex)
  
  ages18    <- c(paste(seq(0, 80, 5), seq(4, 84, 5), sep = "-"), "85+")
  map_age18 <- setNames(c("0-4", "0-4", ages18[-1]), ages)
  anios     <- as.character(seq(1980, 2020, by = 5))
  
  mx18 <- pred_sb2 %>%
    mutate(age = unname(map_age18[agegroup])) %>%
    group_by(region, period, sex, age) %>%
    summarise(mx = sum(mx * population) / sum(population), .groups = "drop") %>%
    mutate(sex = ifelse(sex == "m", 1L, 2L))
  
  armar_e0 <- function(sx) {
    w <- e0_resumen %>% filter(sex == sx) %>%
      mutate(year = substr(period, 1, 4)) %>%
      select(region, year, e0) %>%
      tidyr::pivot_wider(names_from = year, values_from = e0)
    for (a in setdiff(anios, names(w))) w[[a]] <- NA_real_
    plantilla_reg %>% left_join(w, by = "region") %>%
      select(country_code, reg_code, name, include_code, all_of(anios))
  }
  
  armar_mx <- function(sx) {
    w <- mx18 %>% filter(sex == sx) %>%
      mutate(year = substr(period, 1, 4)) %>%
      select(region, age, year, mx) %>%
      tidyr::pivot_wider(names_from = year, values_from = mx)
    for (a in setdiff(anios, names(w))) w[[a]] <- NA_real_
    tidyr::expand_grid(plantilla_reg %>% select(region, reg_code, include_code),
                       age = ages18) %>%
      left_join(w, by = c("region", "age")) %>%
      select(reg_code, age, include_code, all_of(anios))
  }
  ## end EG
  
  # begin EGR
  archivo_summary <- archivo_pdf <- archivo_cobertura <- archivo_cobertura_completa <- archivo_cobertura_completa_sexo <- metadatos_modelo <- archivo_tablas_vida <- archivo_e0_resumen <- archivo_e0F <- archivo_e0M <- archivo_mxF <- archivo_mxM <- archivo_e0_IC <- NULL
  # end EGR
  
  if (guardar) {
    fecha_hora      <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
    carpeta_corrida <- file.path(carpeta_resultados, paste(fecha_hora, tabla_df, sep = "_"))
    dir.create(carpeta_corrida, recursive = TRUE, showWarnings = FALSE)
    
    archivo_summary <- file.path(carpeta_corrida, paste0(fecha_hora, "_02_summary.txt"))
    writeLines(capture.output(summary(fit_sb2)), archivo_summary)
    
    archivo_pdf <- file.path(carpeta_corrida, paste0(fecha_hora, "_22_e0_graficas.pdf"))
    grDevices::pdf(file = archivo_pdf, width = 8 * 2, height = 12 * 2)
    purrr::walk(grafica_e0, print)
    grDevices::dev.off()
    
    archivo_cobertura <- file.path(carpeta_corrida, paste0(fecha_hora, "_11_cobertura.csv"))
    readr::write_csv(as.data.frame(tabla_cobertura), archivo_cobertura)
    
    archivo_cobertura_completa <- file.path(carpeta_corrida, paste0(fecha_hora, "_12_cobertura_completa.csv"))
    readr::write_csv(as.data.frame(tabla_cobertura_completa), archivo_cobertura_completa)
    
    archivo_cobertura_completa_sexo <- file.path(carpeta_corrida, paste0(fecha_hora, "_13_cobertura_completa_sexo.csv"))
    readr::write_csv(as.data.frame(tabla_cobertura_completa_sexo), archivo_cobertura_completa_sexo)
    
    metadatos_modelo <- file.path(carpeta_corrida, paste0(fecha_hora, "_01_metadatos.txt"))
    writeLines(nombre_modelo, metadatos_modelo)
    
    ## begin EGR
    archivo_tablas_vida <- file.path(carpeta_corrida, paste0(fecha_hora, "_24_tablas_vida.csv"))
    readr::write_csv(as.data.frame(tablas_vida), archivo_tablas_vida)
    
    archivo_e0_IC <- file.path(carpeta_corrida, paste0(fecha_hora, "_21_e0_IC.csv"))
    readr::write_csv(as.data.frame(modelo_final_con), archivo_e0_IC)
    
    archivo_e0_resumen <- file.path(carpeta_corrida, paste0(fecha_hora, "_23_e0_resumen.csv"))
    readr::write_csv(as.data.frame(e0_resumen), archivo_e0_resumen)
    
    esc <- function(d, f) write.table(d, f, sep = "\t", row.names = FALSE,
                                      quote = FALSE, na = "")
    if (2 %in% e0_resumen$sex) {
      archivo_e0F <- file.path(carpeta_corrida, paste0(fecha_hora, "_31_e0F.txt"))
      archivo_mxF <- file.path(carpeta_corrida, paste0(fecha_hora, "_33_mxF.txt"))
      esc(armar_e0(2), archivo_e0F); esc(armar_mx(2), archivo_mxF)
    }
    if (1 %in% e0_resumen$sex) {
      archivo_e0M <- file.path(carpeta_corrida, paste0(fecha_hora, "_32_e0M.txt"))
      archivo_mxM <- file.path(carpeta_corrida, paste0(fecha_hora, "_34_mxM.txt"))
      esc(armar_e0(1), archivo_e0M); esc(armar_mx(1), archivo_mxM)
    }
    ## end EGR
  }
  
  invisible(list(
    nombre_modelo    = nombre_modelo,
    datos            = datos_modelo,
    formula          = formula_modelo,
    fit              = fit_sb2,
    modelo_final_con = modelo_final_con,
    periodos         = periodo,
    graficas         = grafica_e0,
    cobertura        = tabla_cobertura,
    # begin EGR
    tablas_vida      = tablas_vida,
    e0_resumen       = e0_resumen,
    cobertura_completa      = tabla_cobertura_completa,
    cobertura_completa_sexo = tabla_cobertura_completa_sexo,
    # end EGR
    tiempo_ajuste    = tiempo_ajuste,
    tiempo_e0        = tiempo_e0,
    archivos         = list(summary   = archivo_summary,
                            graficas  = archivo_pdf,
                            cobertura = archivo_cobertura,
                            cobertura_completa = archivo_cobertura_completa,
                            cobertura_completa_sexo = archivo_cobertura_completa_sexo,
                            metadatos = metadatos_modelo,
                            # begin EGR
                            tablas_vida = archivo_tablas_vida,
                            e0_resumen  = archivo_e0_resumen,
                            e0_IC = archivo_e0_IC,
                            # end EGR
                            e0F = archivo_e0F, e0M = archivo_e0M,
                            mxF = archivo_mxF, mxM = archivo_mxM)
  ))
  
}

#Resultado de muestras para la familia Poisson sb2(1,1,10)
resultado_ambos_poisson <- modelo_completo(
  tabla_df   = "ambos",  
  familia    = "poisson", 
  model_age  = "rw2",     #Mejora de RW1 a RW2
  par_p_age  = 1,
  par_q_age  = 0.5,
  par_b_age  = 25,
  model_per  = "rw2",
  par_p_per  = 0.5,
  par_q_per  = 0.5,
  par_b_per  = 25,
  model_reg  = "bym2",
  par_p_reg  = 0.5,
  par_q_reg  = 0.5,
  par_b_reg  = 25,
  model_s_t  = "iid",
  par_p_s_t  = 0.5,
  par_q_s_t  = 0.5,
  par_b_s_t  = 25,
  model_cel  = "iid",
  par_p_cel  = 1,
  par_q_cel  = 1,
  par_b_cel  = 25,
  nsamples   = 100
)

# Todos los resultados de "resultados_ambos_poisson"
summary(resultado_ambos_poisson$fit)                  # summary del fit  
View(resultado_ambos_poisson$cobertura)               # tabla de cobertura
View(resultado_ambos_poisson$cobertura_completa)      # tabla de cobertura completa
View(resultado_ambos_poisson$cobertura_completa_sexo) # tabla de cobertura completa por sexo
View(resultado_ambos_poisson$modelo_final_con)        # tabla de vida e0_observado, e0_estimado e IC
View(resultado_ambos_poisson$e0_resumen)              # tabla de e0 por region-periodo-sexo OJO: (incluye PR)
View(resultado_ambos_poisson$tablas_vida)             # tablas de vida completas (Age, nMx, nqx, lx, ndx, nLx, Tx, ex) OJO: (incluye PR)
resultado_ambos_poisson$archivos$metadatos            # parámetros del modelo
resultado_ambos_poisson$graficas[["2020-2024"]]       # gráfica de IC para un período

# Visualizando coeficientes del fit de "resultados_ambos_poisson"
resultado_ambos_poisson$fit$summary.fixed
resultado_ambos_poisson$fit$summary.hyperpar
resultado_ambos_poisson$fit$summary.random
resultado_ambos_poisson$fit$summary.fitted.values
resultado_ambos_poisson$fit$dic$dic
resultado_ambos_poisson$fit$waic$waic
resultado_ambos_poisson$fit$mlik
resultado_ambos_poisson$fit$cpu.used
resultado_ambos_poisson$fit$.args$data


