################################################
### Modelo para mortalidad en áreas pequeñas ###
################################################
# 2026_07_24

# Tareas:
# 1. Poner en la misma escala los valores de e0 observado,
#    e0 estimado y sus intervalos de credibilidad.
# 2. Segmentar la cobertura en:
#    - Estimación adecuada.
#    - Por debajo del intervalo: e0 observado < e0_lower.
#    - Por encima del intervalo: e0 observado > e0_upper.
# 3. Realizar el tuning de las previas para TODAS las etiquetas:
#    age, region, period, region_period y cell.
# 4. Incorporar previas diferentes por sexo:
#    - Hombres: p = q = 0.5.
#    - Mujeres: p = q = 1.
# 5. Escalar el número de muestras posteriores:
#    100, 500, 1000 o más muestras.
# 6. Extender el modelo de Poisson a Binomial Negativa
#    y comparar DIC, WAIC y cobertura.
# 7. Repetir el análisis utilizando regiones senatoriales.
# 8. Optimizar el análisis de sensibilidad.

# Paquetes
library(INLA)
library(SUMMER)
library(tidyverse)
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

script_dir         <- this.path::this.dir()
data_dir           <- file.path(script_dir, "data")
carpeta_resultados <- "resultados_modelos"

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
        tb <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, Sex = sx,
                          a0rule = "ak", axmethod = "pas", mod = FALSE)
        e0_obs_list[[length(e0_obs_list) + 1]] <- data.frame(
          region = reg, period = per, sex = ifelse(sx == "m", 1, 2),
          e0_observado = tb$ex[1]
        )
      }
    }
  }
  e0_observado_df <- bind_rows(e0_obs_list)
  
  # Muestras posteriores del predictor
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
      tb <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, Sex = sx,
                        a0rule = "ak", axmethod = "pas", mod = FALSE)
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
calcular_e0_inla_opt <- function(modelo_inla, df, age_params, Age, nsamples = 1000, mc.cores = 19, ...){
  
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
  
  # función auxiliar: e0 de un grupo dado un vector de mx (df completo)
  e0_grupo <- function(g, mx_vec) {
    lt_abridged(nMx = mx_vec[idx_list[[g]]], AgeInt = AgeInt_list[[g]],
                Age = Age, Sex = sex_list[g],
                a0rule = "ak", axmethod = "pas", mod = FALSE)$ex[1]
  }
  
  # --- e0 observado ---------------------------------------------------------
  e0_obs <- vapply(seq_len(G), e0_grupo, numeric(1), mx_vec = mx_obs)
  
  # --- muestras posteriores del predictor -----------------------------------
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
    
    trabajador <- function(subm, idx_list, AgeInt_list, sex_list, Age) {
      G <- length(idx_list)
      apply(subm, 2, function(mx_vec) {
        vapply(seq_len(G), function(g) {
          DemoTools::lt_abridged(nMx = mx_vec[idx_list[[g]]],
                                 AgeInt = AgeInt_list[[g]],
                                 Age = Age, Sex = sex_list[g],
                                 a0rule = "ak", axmethod = "pas",
                                 mod = FALSE)$ex[1]
        }, numeric(1))
      })
    }
    # entorno limpio: evita serializar mx_matrix/samples completos a cada worker
    environment(trabajador) <- globalenv()
    
    cl <- parallel::makeCluster(mc.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    res <- parallel::parLapply(cl, sub_mats, trabajador,
                               idx_list = idx_list, AgeInt_list = AgeInt_list,
                               sex_list = sex_list, Age = Age)
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
    mutate(region2 = fct_reorder(paste(region, sex, sep = "___"), e0_estimado))   # 1
  ggplot(d,
         aes(x = region2,                                                          # 2
             y = e0_estimado, ymin = e0_lower, ymax = e0_upper)) +
    geom_pointrange(color = col, size = 0.3) +
    geom_point(aes(y = e0_observado)) +
    coord_flip() +
    scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +                  # 3
    facet_wrap(sex ~ est_eval, scales = "free",
               labeller = labeller(sex = c(`1` = "Hombres", `2` = "Mujeres"))) +
    theme_minimal() +
    labs(title = paste0("e0 por municipio (estimada vs. observada), ", per, ", ", llh),
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

# Definición de parámetros para INLA
familias  <- names(INLA::inla.models()$likelihood)
familia   <- "poisson"
modelos   <- names(INLA::inla.models()$latent)
model_age <- "rw1"
par_p_age <- 1
par_q_age <- 1
par_b_age <- 1
model_reg <- "bym2"
par_p_reg <- 1
par_q_reg <- 1
par_b_reg <- 1
model_per <- "rw2"
par_p_per <- 1
par_q_per <- 1
par_b_per <- 1
model_s_t <- "iid"
par_p_s_t <- 1
par_q_s_t <- 1
par_b_s_t <- 1
model_cel <- "iid"
par_p_cel <- 1
par_q_cel <- 1
par_b_cel <- 1
nsamples  <- 100

# Etiqueta usada para nombrar los archivos generados
nombre_modelo <- paste(
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
  factor(sex) +
  f(age_idx, model = model_age, constr = TRUE,
    hyper = SB2.prior(par_p_age, par_q_age, par_b_age)) +
  f(region_idx, model = model_reg, graph = g, constr = TRUE,
    hyper = SB2.prior(par_p_reg , par_q_reg , par_b_reg )) +
  f(period_idx, model = model_per,  constr = TRUE,
    hyper = SB2.prior(par_p_per, par_q_per, par_b_per)) +
  f(region_period_idx, model = model_s_t,
    hyper = SB2.prior(par_p_s_t, par_q_s_t, par_b_s_t)) +
  f(cell_idx, model = model_cel,
    hyper = SB2.prior(par_p_cel, par_b_cel, par_b_cel))

# Ejecutar la fórmula para INLA
fit_sb2 <- inla(formula_sb2,
                family = familia,
                data = df,
                E = population,
                control.compute = list(config = TRUE, dic = TRUE, waic = TRUE))

# Revisar el resumen
summary(fit_sb2)

# Ejecutar las muestras por cada modelo
system.time(modelo_final_con <- calcular_e0_inla_opt(fit_sb2, df, age_params,
                                                     Age, nsamples = nsamples))

# Elegir el período
periodo   <- "2020-2024"

# Graficar los IC
grafica_e0 <- e0_model_plot(modelo_final_con, periodo, "purple", nombre_modelo)

# Tabular la cobertura
tabla <- modelo_final_con %>%
  group_by(period, sex) %>%
  summarise(cobertura = 100 * mean(est_eval == "Estimación adecuada"),
            .groups = "drop")

# Salvar el gráfico usando en el nombre una marca temporal inicial con el formato
# AAAA-MM-DD-HH-MM-SS (la hora está en formato 24 horas para un orden automático)
# junto con la etiqueta del "nombre_modelo"
ggplot2::ggsave(
  filename  = file.path(carpeta_resultados,
                        paste0(format(Sys.time(), "%Y-%m-%d-%I-%M-%S"),"-",
                               nombre_modelo, "_", periodo,"_e0.pdf")),
  plot      = grafica_e0,
  device    = grDevices::cairo_pdf,
  width     = 8,
  height    = 12,
  units     = "in",
  scale     = 2,
  limitsize = FALSE
)

# Salvar la tabla de cobertura
readr::write_csv(as.data.frame(tabla),
                 file.path(carpeta_resultados, 
                           paste0(format(Sys.time(),"%Y-%m-%d-%I-%M-%S"), "-",
                                  nombre_modelo, "_cobertura.csv")))
