library(INLA)
library(SUMMER)
library(tidyverse)
library(tidycensus)
library(sf)
library(here)
set.seed(123)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. Datos con grupos qu in qu en al es del censo
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Obtener nombres de municipios limpios, tidycensus
municipio <- fips_codes %>%
  filter(state == "PR") %>%
  select(county) %>%
  pull() %>%
  str_remove(" Municipio")

regions <- municipio
periods <- 1980:2024
  
ages <- c("0-4", "5-9", "10-14", "15-19" , "20-24" ,
          "25-29" , "30-34" , "35-39" , "40-44" , "45-49" ,
          "50-54", "55-59" , "60-64" , "65-69" , "70-74" ,
          "75-79" , "80+" )
age_effect <- c(0.0150, 0.0005, 0.0003, 0.0008, 0.0012,
                0.0015, 0.0018, 0.0022, 0.0030, 0.0045,
                0.0070, 0.0110, 0.0170, 0.0260, 0.0400,
                0.0600, 0.1000)

df <- expand.grid(region = regions,
                  period = periods,
                  agegroup = ages)
df$population <- round(runif(nrow(df), 800, 2000))
df$deaths <- rbinom(nrow(df),
                    size = df$population,
                    prob = rep(age_effect, each = length(regions) * length(periods))
)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. Parametros de tabla de vida
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
age_params <- tibble(agegroup = ages,
                     n_interval = c(rep(5, length(ages) - 1), NA),
                     ax = c(2.0, 2.5, 2.5, 2.5, 2.5,
                            2.5, 2.5, 2.5, 2.5, 2.5,
                            2.5, 2.5, 2.5, 2.5, 2.5,
                            2.5, NA)
)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. Matriz de adyacencia
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
shapefile_sf <- st_read(here("municipios_shp", "g03_legales_municipios_edicion_octubre2015.shp"))
# shapefile_sf <- st_read("C:/Users/Estudiante/Desktop/Demografia 2026/municipios_shp/g03_legales_municipios_edicion_octubre2015.shp")

Amat <- getAmat(geo = shapefile_sf$geometry, names = municipio)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 3. Indices para INLA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
df <- df %>%
  mutate(
    region_idx = as.integer(factor(region, levels = regions)),
    period_idx = as.integer(factor(period, levels = periods)),
    age_idx = as.integer(factor(agegroup, levels = ages)),
    region_period_idx = as.integer(factor(paste(region, period)))
  )

g <- INLA::inla.read.graph(Amat)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 4. Modelo INLA Poisson
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
formula_inla <- deaths ~
  f(age_idx, model = "rw1", constr = TRUE) +
  f(region_idx, model = "bym2", graph = g, constr = TRUE) +
  f(period_idx, model = "rw2", constr = TRUE) +
  f(region_period_idx, model = "iid")

system.time(
  fit <- inla(
    formula = formula_inla,
    family = "poisson",
    data = df,
    E = population,
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE)
  )
)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 5. Extraer tasas suavizadas
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pred <- df %>%
  mutate(mx = pmax(fit$summary.fitted.values$mean, 1e-6)) %>%
  left_join(age_params, by = "agegroup")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 6. Calcular qx
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pred <- pred %>%
  mutate(
    qx = case_when(
      agegroup == "80+" ~ 1,
      TRUE ~ (n_interval * mx) / (1 + (n_interval - ax) * mx)
    ),
    qx = pmin(pmax(qx, 0), 1)
  )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 7. Tabla de vida
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
life_tables <- pred %>%
  arrange(period, region, age_idx) %>%
  group_by(period, region) %>%
  group_modify(~ {
    df_lt <- .x
    k <- nrow(df_lt)
    lx <- numeric(k + 1)
    dx <- numeric(k)
    Lx <- numeric(k)
    Tx <- numeric(k)
    lx[1] <- 100000
    
    for (i in seq_len(k)) {
      dx[i] <- lx[i] * df_lt$qx[i]
      lx[i + 1] <- lx[i] - dx[i]
      
      if (df_lt$agegroup[i] == "80+") {
        Lx[i] <- ifelse(df_lt$mx[i] > 0, lx[i] / df_lt$mx[i], 0)
      } else {
        Lx[i] <- df_lt$n_interval[i] * lx[i] - (df_lt$n_interval[i] - df_lt$ax[i]) * dx[i]
      }
    }
    
    Tx[k] <- Lx[k]
    for (i in (k - 1):1) Tx[i] <- Tx[i + 1] + Lx[i]
    
    tibble(
      agegroup = df_lt$agegroup,
      mx = df_lt$mx,
      qx = df_lt$qx,
      lx = lx[1:k],
      dx = dx,
      Lx = Lx,
      Tx = Tx,
      ex = Tx / lx[1:k],
      e0 = Tx[1] / lx[1]
    )
  }) %>%
  ungroup()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 8. Esperanza de vida
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
e0_resumen <- life_tables %>%
  distinct(period, region, e0) %>%
  arrange(region, period)

print(e0_resumen)
