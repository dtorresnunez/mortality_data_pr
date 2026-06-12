library(SUMMER)
library(INLA)
library(tidyverse)
library(tidycensus)
library(here)
library(rio)
library(sf)
library(DemoTools)
# -----------------------------------------------------------------

# Obtener nombres de municipios limpios, tidycensus
municipio <- fips_codes %>%
  filter(state == "PR") %>%
  select(county) %>%
  pull() %>%
  str_remove(" Municipio")

# ==============================================================================
# 1. PARÁMETROS BASE Y CREACIÓN DEL GRID QUINQUENAL
# ==============================================================================
regions <- municipio
periods_list <- as.character(1980:2024)
ages <- c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")

df <- expand.grid(region = regions,
                  period = periods_list, 
                  agegroup = ages)

# ==============================================================================
# 2. SIMULACIÓN DE POBLACIÓN (Escala de ~1 millón de mujeres en el país por año)
# ==============================================================================
# 1,000,000 / (78 municipios * 7 grupos de edad) ≈ 1,831 mujeres por celda.
df$women <- round(runif(nrow(df), min = 1000, max = 2660))

# ==============================================================================
# 3. GENERACIÓN DE CURVAS DINÁMICAS POR AÑO (Transición Histórica)
# ==============================================================================
# Asignamos directamente las tasas a los 7 grupos quinquenales (sin spline)

# Curva A: 1980 (Pico temprano en 20-24, escalado para que 5 * sum(tasas) = 2.70)
efectos_1980 <- c(0.06, 0.155, 0.11, 0.06, 0.03, 0.01, 0.001)
curve_1980   <- (efectos_1980 / (5 * sum(efectos_1980))) * 2.70

# Curva B: 2024 (Baja fecundidad, pico temprano, escalado para un TFR de ~1.20)
efectos_2024 <- c(0.02, 0.070, 0.045, 0.025, 0.01, 0.002, 0.0002)
curve_2024   <- (efectos_2024 / (5 * sum(efectos_2024))) * 0.9

# Crear la matriz de transición temporal (ahora de 7 filas correspondiente a los grupos)
years_numeric <- 1980:2024
num_years <- length(years_numeric)
matrix_rates <- matrix(NA, nrow = 7, ncol = num_years)
rownames(matrix_rates) <- ages
colnames(matrix_rates) <- periods_list

for(i in 1:num_years){
  alpha <- (2024 - years_numeric[i]) / (2024 - 1980)
  matrix_rates[, i] <- alpha * curve_1980 + (1 - alpha) * curve_2024
}

# ==============================================================================
# 4. MAPEADO DE PROBABILIDADES Y SIMULACIÓN DE NACIMIENTOS
# ==============================================================================
df$prob_fecundidad <- matrix_rates[cbind(df$agegroup, df$period)]

# Simulación binomial
df$births <- rbinom(nrow(df), size = df$women, prob = df$prob_fecundidad)

df$prob_fecundidad <- NULL

# ==============================================================================
# 5. DIAGNÓSTICO Y VERIFICACIÓN DE LA TRANSICIÓN (Factor multiplicador de 5)
# ==============================================================================
verificacion <- df %>%
  group_by(period) %>%
  summarise(
    Poblacion_Mujeres_Pais = sum(women),
    Nacimientos_Pais = sum(births),
    # Multiplicamos por 5 porque cada tasa representa un bloque de 5 años de edad
    TFR_Simulado = 5 * (sum(births / women) / n_distinct(region))
  )

print("--- COMPORTAMIENTO GLOBAL NACIONAL ---")
print(head(verificacion, 3))
print(tail(verificacion, 3))

# TFR desagregado por municipio y año
tfr_por_municipio <- df %>%
  group_by(period, region) %>%
  summarise(
    Poblacion_Mujeres = sum(women),
    Nacimientos = sum(births),
    # Multiplicamos por 5 para obtener el TFR municipal real
    TFR_Municipal = 5 * sum(births / women),
    .groups = "drop"
  )

print("--- MUESTRA DEL TFR POR MUNICIPIO ---")
print(head(tfr_por_municipio))
tail(tfr_por_municipio)


shapefile_sf <- st_read(here("municipios_shp", "g03_legales_municipios_edicion_octubre2015.shp"))

# Convertir a formato  SUMMER
#shapefile_sp <- as_Spatial(shapefile_sf)

Amat <- getAmat(geo = shapefile_sf$geometry, names = municipio)


df$region_idx <- as.numeric(factor(df$region, levels = regions))
df$period_idx <- as.numeric(factor(df$period, levels = periods_list))

# NUEVO: Índice numérico para el grupo de edad
df$age_idx    <- as.numeric(factor(df$agegroup, levels = ages))

# -----------------------------------------------------------------
# 3. Ajustar modelo Bayesiano Global (Un solo FIT)
# -----------------------------------------------------------------

# Añadimos f(age_idx, model = "rw2") 

formula_global <- births ~ 1 + 
  # Penaliza cambios bruscos entre edades consecutivas
  f(age_idx, model = "rw2", 
    hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))) + 
  f(region_idx, model = "bym2", graph = Amat,
    hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))) +
  f(period_idx, model = "rw2", 
    hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01))))

system.time(
  fit_global <- inla(formula_global, 
                     family = "poisson", 
                     data = df, 
                     E = df$women, 
                     control.predictor = list(compute = TRUE))
  
)

# Extraer las ASFR estimadas directamente del modelo
df$asfr_suavizada       <- fit_global$summary.fitted.values$mean
df$asfr_suavizada_lower <- fit_global$summary.fitted.values$`0.025quant`
df$asfr_suavizada_upper <- fit_global$summary.fitted.values$`0.975quant`

# -----------------------------------------------------------------
# 4. Consolidar y Calcular TFR = 5 * sum(ASFR)
# -----------------------------------------------------------------
# NOTA SOBRE LOS INTERVALOS DE LA TFR: 
# Sumar los cuantiles directos (lower y upper) asume correlación perfecta.
# Para el ejemplo mantenemos tu lógica, pero la TFR puntual ahora es mucho más robusta.

TFR_global <- df %>%
  group_by(region, period) %>%
  summarise(
    TFR_lower = 5 * sum(asfr_suavizada_lower),
    TFR       = 5 * sum(asfr_suavizada),
    TFR_upper = 5 * sum(asfr_suavizada_upper), 
    .groups = "drop"
  )

# -----------------------------------------------------------------
# 5. Imprimir Resultados Finales
# -----------------------------------------------------------------
summary(fit_global)
View(TFR_global)

TFR_global %>%
  filter(region == "San Juan")
