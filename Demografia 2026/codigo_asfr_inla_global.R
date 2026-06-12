# Instalar INLA si no está presente en tu entorno:
# install.packages("INLA", repos=c(getOption("repos"), INLA="https://r-inla-download.org"), dep=TRUE)

library(SUMMER)
library(INLA)
library(dplyr)
library(tidyr)

set.seed(123)

# -----------------------------------------------------------------
# 1. Crear datos inventados con Intervalos Quinquenales (3 periodos)
# -----------------------------------------------------------------
regions <- c("Norte", "Centro", "Sur")
periods_list <- c("2005-2009", "2010-2014", "2015-2019") 
ages <- c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")

df <- expand.grid(region = regions,
                  period = periods_list, 
                  agegroup = ages)

df$women <- round(runif(nrow(df), 1500, 4000))
age_effect <- c(0.08, 0.15, 0.12, 0.08, 0.04, 0.01, 0.002)

df$births <- rbinom(nrow(df), size = df$women,
                    prob = rep(age_effect, each = length(regions) * length(periods_list)))

# -----------------------------------------------------------------
# 2. Matriz de adyacencia espacial e Índices
# -----------------------------------------------------------------
adj <- matrix(0, 3, 3)
rownames(adj) <- colnames(adj) <- regions
adj["Norte", "Centro"] <- 1
adj["Centro", "Norte"] <- 1
adj["Centro", "Sur"] <- 1
adj["Sur", "Centro"] <- 1
Amat <- adj

df$region_idx <- as.numeric(factor(df$region, levels = regions))
df$period_idx <- as.numeric(factor(df$period, levels = periods_list))

# NUEVO: Índice numérico para el grupo de edad
df$age_idx    <- as.numeric(factor(df$agegroup, levels = ages))

# -----------------------------------------------------------------
# 3. Ajustar modelo Bayesiano Global (Un solo FIT)
# -----------------------------------------------------------------

# Añadimos f(age_idx, model = "rw2") 
formula_global <- births ~ 1 + 
  f(age_idx, model = "rw2") + 
  f(region_idx, model = "bym2", graph = Amat) + 
  f(period_idx, model = "rw2")

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
print(TFR_global)
