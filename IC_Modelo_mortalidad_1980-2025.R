################################################
### Modelo para mortalidad en áreas pequeñas ###
################################################
#2026_07_11

#Paquetes
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
library(readxl)# EGR

set.seed(123)

script_dir <- this.path::this.dir()
data_dir <- file.path(script_dir, "data")
shp_dir <- file.path(data_dir, "municipios_shp")

set.seed(123)

# -----------------------------
# 1. Datos reales
# -----------------------------
muni_xwalk <- fips_codes %>%
  filter(state == "PR") %>%
  transmute(
    fips3 = county_code,
    region = str_remove(county, " Municipio")
  )

municipio <- muni_xwalk$region
regions <- municipio

#Periodo para años sencillos
periods <- 1980:2025 #2024

#"00-04"
ages <- c(
  "0", "01-04","05-09", "10-14", "15-19", "20-24",
  "25-29", "30-34", "35-39", "40-44", "45-49",
  "50-54", "55-59", "60-64", "65-69", "70-74",
  "75-79", "80-84", "85+"
)

#Periodos por años quinquenales 2000-2024
#period_breaks <- c(seq(2000,2020, by = 5), 2024)
#period_labels <- paste0(seq(2000, 2020, by = 5), "-", seq(2004, 2024, by = 5))

#Periodos por años quinquenales 1980-2024
period_breaks <- c(seq(1980, 2020, by = 5), 2025)
period_labels <- paste0(seq(1980, 2020, by = 5), "-", c(seq(1984, 2020, by = 5), 2025)) #paste0(seq(1980, 2020, by = 5), "-", c(seq(1985, 2020, by = 5), 2025))
print(period_labels) #confirmar que si está por quinquenio

#Estricto cinco años, el último quinquenio de 2000-2024
period_labels <- paste0(seq(1980, 2020, by = 5), "-", seq(1984, 2024, by = 5))


period_quinquenal <- function(year) {
  as.character(
    cut(
      year,
      breaks = period_breaks,
      labels = period_labels,
      right = FALSE,          
      include.lowest = TRUE
    )
  )
}

poblacion <- read_csv(
  file.path(data_dir, "municipios_population_1980_2025_redondeo.csv"),
  col_types = cols(
    fips3 = col_character(),
    .default = col_guess()
  )
) %>%
  filter(
    year >= 1980, #2000
    year <= 2025, #2024
    agegrp != 0,
    #sex !=0 #Esto es si solo queremos trabajar en grupo solo con dos sexos, no el total (0).
  ) %>%
  mutate(
    fips3 = str_pad(fips3, width = 3, side = "left", pad = "0"),
    period = period_quinquenal(year),
    agegroup = case_when(
      agegrp == 1 ~ "0",
      agegrp == 2  ~ "01-04",
      agegrp == 3  ~ "05-09",
      agegrp == 4  ~ "10-14",
      agegrp == 5  ~ "15-19",
      agegrp == 6  ~ "20-24",
      agegrp == 7  ~ "25-29",
      agegrp == 8  ~ "30-34",
      agegrp == 9  ~ "35-39",
      agegrp == 10  ~ "40-44",
      agegrp == 11 ~ "45-49",
      agegrp == 12 ~ "50-54",
      agegrp == 13 ~ "55-59",
      agegrp == 14 ~ "60-64",
      agegrp == 15 ~ "65-69",
      agegrp == 16 ~ "70-74",
      agegrp == 17 ~ "75-79",
      agegrp == 18 ~ "80-84",
      agegrp %in% c(19, 20) ~ "85+" #OJO: 5-9 vs 05-09  #HOLD >= vs %in%
    ) 
  ) %>%
  group_by(fips3, agegroup, sex, period) %>% #year por period, para años sencillos
  summarise(
    #population = sum(population, na.rm = TRUE), #años sencillos
    population_mean = mean(population, na.rm = TRUE),
    n_years = n(), 
    .groups = "drop"
  ) %>%
  mutate(
    population = population_mean * n_years          
  ) %>%
  transmute(
    fips3,
    #period = as.integer(year), #años sencillos
    period,
    agegroup,
    sex, 
    population = as.numeric(population)
  ) %>%
  left_join(muni_xwalk, by = "fips3")

#Nota: solo aparece la poblacion de hombre y mujeres, no la población (0).
# begin EGR
ANIOS_CORREGIDOS <- 2015:2020

defunciones_orig <- read_dta(
  file.path(data_dir, "defunciones_municipios_long_1979_2023.dta")
) %>% rename(sex = sexo) %>%
  filter(
    year >= min(periods),
    year <= max(periods),
    !(year %in% ANIOS_CORREGIDOS), # EGR: excluimos los anios a corregir
    !is.na(fips3),
    fips3 != "",
    !is.na(edad),
  ) %>%
  mutate(
    period = period_quinquenal(year),
    #period = as.integer(year), #años sencillos
    agegroup = as.character(
      cut(
        edad,
        breaks = c(0, 1, seq(5, 85, by = 5), Inf), #HOLD c(seq(0, 85, by = 5), Inf),
        labels = ages,
        right = FALSE
      )
    )
  ) %>%
  mutate(sex = as.integer(sex)) %>% # EGR: por si sexo con la nueva base no se lee como numerico
  count(fips3, period, agegroup, sex, name = "deaths")

# EGR: ahora hay que poner las defunciones de la base nueva en el formato que ya se tiene
# EGR: se deben corregir FIPS3
pedazo <- read_excel(
  file.path(data_dir, "2026-07-16_corregidas_defunciones_wide_2015-2023.xlsx"),
  sheet = "Sheet1"
)
pedazo <- unique(data.frame(
  muni       = pedazo$muni,
  fips3_xlsx = str_pad(as.character(pedazo$fips3), 3, "left", "0")
))
pedazo <- pedazo[!is.na(pedazo$muni), ]
pedazo$fips3_ok <- muni_xwalk$fips3[
  match(tolower(chartr("áéíóúüñ", "aeiouun", pedazo$muni)),
        tolower(chartr("áéíóúüñ", "aeiouun", muni_xwalk$region)))
]
sum(is.na(pedazo$fips3_ok))
pedazo[pedazo$fips3_xlsx != pedazo$fips3_ok & !is.na(pedazo$fips3_ok), ]
# EGR: correcion hecha sobre FIPS3
defunciones_corr <- read_excel(
  file.path(data_dir, "2026-07-16_corregidas_defunciones_wide_2015-2023.xlsx"),
  sheet = "Sheet1"
) %>%
  filter(year %in% ANIOS_CORREGIDOS) %>%
  select(-edad_NA) %>%                        # muertes sin edad: se descartan
  rename(sex = sexo) %>%
  pivot_longer(
    cols = starts_with("edad_"),
    names_to = "edad",
    names_prefix = "edad_",
    values_to = "deaths"
  ) %>%
  mutate(
    edad  = as.numeric(edad),
    sex   = as.integer(sex),
    fips3 = str_pad(as.character(fips3), 3, "left", "0"),   # el xlsx trae 1, 3, 5...
    fips3 = recode(fips3,
                   "011" = "013", "013" = "015", "015" = "011",
                   "055" = "057", "057" = "059", "059" = "061",
                   "061" = "063", "063" = "055") # EGR: se re-mapea FIPS3 para los pedazos indentificados en pedazo[pedazo$fips3_xlsx != pedazo$fips3_ok & !is.na(pedazo$fips3_ok), ] 
      ) %>%
  filter(
    !is.na(fips3),
    fips3 != "",
    !is.na(sex),
    deaths > 0
  ) %>%
  mutate(
    period = period_quinquenal(year),
    agegroup = as.character(
      cut(
        edad,
        breaks = c(0, 1, seq(5, 85, by = 5), Inf),
        labels = ages,
        right = FALSE
      )
    )
  ) %>%
  group_by(fips3, period, agegroup, sex) %>%
  summarise(deaths = sum(deaths), .groups = "drop")

# EGR: juntamos las defunciones de la base vieja con la nueva solo
# EGR: solo en los anios 2015 - 2020. De 2021 a 2023 nos quedamos con la
# EGR: base inicial de defunciones

defunciones <- bind_rows(defunciones_orig, defunciones_corr) %>%
  group_by(fips3, period, agegroup, sex) %>%
  summarise(deaths = sum(deaths), .groups = "drop")

# EGR: comprobemos las muertes de edad 0 por periodo
crudo_orig <- read_dta(
  file.path(data_dir, "defunciones_municipios_long_1979_2023.dta")
) %>%
  filter(
    year >= 2010, year <= 2023,
    !(year %in% ANIOS_CORREGIDOS),
    !is.na(edad)
  ) %>%
  mutate(year = as.integer(year), edad = as.numeric(edad)) %>%
  count(year, edad, name = "deaths")

crudo_corr <- read_excel(
  file.path(data_dir, "2026-07-16_corregidas_defunciones_wide_2015-2023.xlsx"),
  sheet = "Sheet1"
) %>%
  filter(year %in% ANIOS_CORREGIDOS) %>%
  select(-edad_NA) %>%
  pivot_longer(
    cols = starts_with("edad_"),
    names_to = "edad", names_prefix = "edad_", values_to = "deaths"
  ) %>%
  mutate(edad = as.numeric(edad)) %>%
  filter(deaths > 0) %>%
  count(year, edad, wt = deaths, name = "deaths")

tabla_cruda <- bind_rows(crudo_orig, crudo_corr) %>%
  filter(edad <= 12) %>%
  group_by(year, edad) %>%
  summarise(deaths = sum(deaths), .groups = "drop") %>%
  pivot_wider(names_from = edad, values_from = deaths, values_fill = 0) %>%
  arrange(year)

print(tabla_cruda, n = Inf)
# defunciones <- read_dta(
#   file.path(data_dir, "defunciones_municipios_long_1979_2023.dta")
# ) %>% rename(sex = sexo) %>%
#   filter(
#     year >= min(periods),
#     year <= max(periods),
#     !is.na(fips3),
#     fips3 != "",
#     !is.na(edad),
#   ) %>%
#   mutate(
#     period = period_quinquenal(year),
#     #period = as.integer(year), #años sencillos
#     agegroup = as.character(
#       cut(
#         edad,
#         breaks = c(0, 1, seq(5, 85, by = 5), Inf), #HOLD c(seq(0, 85, by = 5), Inf),
#         labels = ages,
#         right = FALSE
#       )
#     )
#   ) %>%
#   count(fips3, period, agegroup, sex, name = "deaths")
# end EGR

df <- poblacion %>%
  left_join(defunciones, by = c("fips3", "period", "agegroup","sex")) %>%
  mutate(
    deaths = replace_na(deaths, 0L)
  ) %>%
  select(region, period, agegroup, sex, population, deaths)

# ------------------------------
# 2. Parámetros de tabla de vida
# ------------------------------
age_params <- tibble(
  agegroup = ages,
  n_interval = c(1, 4, rep(5, 16), NA), #HOLD n_interval = c(rep(5, length(ages) - 1), NA),
  ax = c(
    0.15, 1.5, 2.5, 2.5, 2.5, #00-04 para 2.0, pero DemoTools 0: 0.15 y 1-4: 1.5
    2.5, 2.5, 2.5, 2.5, 2.5,
    2.5, 2.5, 2.5, 2.5, 2.5,
    2.5, 2.5, 2.5, NA
  )
)

#Para confirmar que si hay 19 grupos
age_params$n_interval
# -----------------------------
# 3. Matriz de adyacencia
# -----------------------------
shapefile_sf <- st_read(
  file.path(shp_dir, "g03_legales_municipios_edicion_octubre2015.shp")
)

Amat <- getAmat(
  geo = shapefile_sf$geometry,
  names = municipio
)

# -------------------------------
# 3.1. Mapa del grafo de vecindad
# -------------------------------
shapefile_sf <- shapefile_sf %>%
  mutate(.id_muni = row_number())

shapefile_largest <- shapefile_sf %>%
  st_cast("POLYGON", warn = FALSE) %>%
  mutate(.area = as.numeric(st_area(geometry))) %>%
  group_by(.id_muni) %>%
  slice_max(.area, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(.id_muni)

cent <- suppressWarnings(
  st_point_on_surface(st_geometry(shapefile_largest))
)

cXY <- st_coordinates(cent)

ij <- which(
  Amat != 0 & upper.tri(Amat),
  arr.ind = TRUE
)

edges_sf <- st_sf(
  geometry = st_sfc(
    lapply(
      seq_len(nrow(ij)),
      function(k) {
        st_linestring(cXY[c(ij[k, 1], ij[k, 2]), ])
      }
    ),
    crs = st_crs(shapefile_sf)
  )
)

nodes_sf <- st_sf(
  region = municipio,
  geometry = cent
)

map_vecindad <- ggplot() +
  geom_sf(
    data = shapefile_sf,
    fill = "grey95",
    color = "grey75",
    linewidth = 0.2
  ) +
  geom_sf(
    data = edges_sf,
    color = "#2c7fb8",
    linewidth = 0.4,
    alpha = 0.85
  ) +
  geom_sf(
    data = nodes_sf,
    color = "#d7301f",
    size = 1.3
  ) +
  labs(
    title = "Grafo de vecindad de los municipios de Puerto Rico",
    subtitle = paste0(
      nrow(Amat), " municipios y ",
      nrow(ij), " conexiones en la matriz de adyacencia"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

print(map_vecindad)

# -----------------------------
# 4. Índices para INLA
# -----------------------------
df <- df %>% filter(sex %in% c(1, 2)) %>%
  mutate(
    region_idx = as.integer(factor(region, levels = regions)),
    period_idx = as.integer(factor(period, levels = period_labels)),
    #period_idx = as.integer(factor(period, levels = periods)), #años sencillos
    age_idx = as.integer(factor(agegroup, levels = ages)),
    region_period_idx = as.integer(factor(paste(region, period)))
  )

g <- INLA::inla.read.graph(Amat)

# ---------------------------------
# 5. Modelo INLA - Familia Poisson
# ---------------------------------
formula_inla <- deaths ~
  factor(sex) +
  f(age_idx, model = "rw1", constr = TRUE) +
  f(region_idx, model = "bym2", graph = g, constr = TRUE) +
  f(period_idx, model = "rw2", constr = TRUE) +
  f(region_period_idx, model = "iid")

fit <- inla(
  formula = formula_inla,
  family = "poisson",
  data = df,
  E = population,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

# -----------------------------
# 6. Extraer tasas suavizadas
# -----------------------------
pred <- df %>%
  mutate(
    mx = pmax(fit$summary.fitted.values$mean, 1e-6)
  ) %>%
  left_join(age_params, by = "agegroup")

#Extraer de otra forma las tasas de mortalidad #HOLD
fit$summary.fitted.values$mean 
fit$summary.fitted.values$`0.025quant`
fit$summary.fitted.values$`0.975quant` 

# -------------------------------------------------
# 7. Epitools - Estandarización directa e indirecta
# -------------------------------------------------
pred_epi <- df %>%
  mutate(
    mx = pmax(fit$summary.fitted.values$mean, 1e-6)
  ) %>%
  left_join(age_params, by = "agegroup")

stdpop_pr <- pred_epi %>%
  group_by(agegroup, sex) %>%
  summarise(stdpop = sum(population))

stdcount_pr <- pred_epi %>%
  group_by(agegroup, sex) %>%
  summarise(stdcount = sum(deaths))

ageadjust.direct(
  count  = pred_epi$deaths,
  pop    = pred_epi$population,
  stdpop = stdpop_pr$stdpop
)

ageadjust.indirect(
  count    = pred_epi$deaths,
  pop      = pred_epi$population,
  stdcount = stdcount_pr$stdcount,
  stdpop   = stdpop_pr$stdpop
)

# ------------------------------
# 8. DemoTools - Tablas de vida
# ------------------------------
pred_demo <- df %>%
  mutate(
    mx = pmax(fit$summary.fitted.values$mean, 1e-6) #Extrae el mx
  ) %>%
  left_join(age_params, by = "agegroup")

# Tabla de vida para un solo pueblo (Adjuntas)
pred_demo <- pred # para no perder pred
pred_demo$sex <- ifelse(pred$sex == 1, "m", "f")
pred_demo <- pred_demo %>%
  filter(region == "Adjuntas", period == "1985-1989", sex == "f")     #ANTES 1980-1985
nMx <- pred_demo$mx
Age <- c(0, 1, (seq(5, 85, by=5)))                #HOLD Age <- c(0, 1, (seq(5, 80, by=5)))
AgeInt <- inferAgeIntAbr(vec = nMx)
PR.lifetable <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, Sex = "f", a0rule = "ak", axmethod = "pas", mod = FALSE)
PR.lifetable

# Tabla de vida para todos los municipios
pred_demo <- pred
pred_demo$sex <- ifelse(pred_demo$sex == 1, "m", "f")
municipios_demo <- sort(unique(pred_demo$region))
periodos_demo   <- sort(unique(pred_demo$period))
sexos      <- c("m", "f")
Age <- c(0, 1, seq(5, 85, by = 5)) 
#Age <- c(0, 1, seq(5, 80, by = 5)) 
#Age <- c(0, seq(5, 85, by = 5))
tablas <- list()
for (muni in municipios_demo) {
  for (per in periodos_demo) {
    for (sx in sexos) {
      pred_sub <- pred_demo %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      tablas[[muni]][[per]][[sx]] <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, Sex = sx, a0rule = "ak", axmethod = "pas", mod = FALSE)
    }
  }
}

e0_resumen_demotools <- data.frame()

for (m in names(tablas)) {
  for (p in names(tablas[[m]])) {
    for (s in c("m", "f")) {
      
      tb <- tablas[[m]][[p]][[s]]
      
      sexnum <- if (s == "m") 1 else 2
      
      e0_resumen_demotools <- rbind(
        e0_resumen_demotools,
        data.frame(period = p, region = m, sex = sexnum, e0 = tb$ex[1])
      )
    }
  }
}


muj <- e0_resumen_demotools[e0_resumen_demotools$sex == 2, ] 
fila_muj <- muj[which.max(muj$e0), ]

hom <- e0_resumen_demotools[e0_resumen_demotools$sex == 1, ] 
fila_hom <- hom[which.max(hom$e0), ]

fila_muj
fila_hom

#Periodo 2015-2019, para comparar con el e0 nacional 
muj <- e0_resumen_demotools[e0_resumen_demotools$sex == 2, ] %>% filter(period == "2015-2019") # EGR: se corrige la etiqueta
fila_muj <- muj[which.max(muj$e0), ]

hom <- e0_resumen_demotools[e0_resumen_demotools$sex == 1, ] %>% filter(period == "2015-2019") # EGR: se corrige la etiqueta
fila_hom <- hom[which.max(hom$e0), ]

fila_muj
fila_hom

head(tablas)
tb

# ---------------------------------
# 9. Tabla de vida con Toy Example 
# --------------------------------

# --------------------------------
# 9.1. Calcular qx
# --------------------------------
pred <- pred %>%
  mutate(
    qx = case_when(
      agegroup == "85+" ~ 1,
      TRUE ~ (n_interval * mx) /
        (1 + (n_interval - ax) * mx)
    ),
    qx = pmin(pmax(qx, 0), 1)
  )

# ---------------------------------
# 9.2 Construcción de Tabla de vida
# ---------------------------------
life_tables <- pred %>%
  arrange(period, region, age_idx, sex) %>%
  group_by(period, region, sex) %>%
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
      
      if (df_lt$agegroup[i] == "85+") { #80
        Lx[i] <- ifelse(
          df_lt$mx[i] > 0,
          lx[i] / df_lt$mx[i],
          0
        )
      } else {
        Lx[i] <- df_lt$n_interval[i] * lx[i] -
          (df_lt$n_interval[i] - df_lt$ax[i]) * dx[i]
      }
    }
    
    Tx[k] <- Lx[k]
    
    for (i in (k - 1):1) {
      Tx[i] <- Tx[i + 1] + Lx[i]
    }
    
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

life_tables

# -----------------------------
# 9.3. Esperanza de vida
# -----------------------------
e0_resumen <- life_tables %>%
  distinct(period, region, sex, e0) %>%
  arrange(region, period, sex)

print(e0_resumen)

# -----------------------------
# 10. Previas in INLA
# -----------------------------

#Las distribuciones a priori se establecen en la representación interna del parámetro, que puede ser diferente de la escala del parámetro en el modelo. 
#Por ejemplo, la precisión se representa en la escala interna en la escala logarítmica.
#Resulta conveniente desde el punto de vista computacional, ya que en la escala interna el parámetro no está acotado.

#Previas a considerar

#Betaprime escalada = Scale_beta2
#Half Cauchy
#PC priors 
#Por defecto = loggamma como la previa para log-precision

#Nota: Scale_beta2 y Half Cauchy se implementan en INLA mediante prior = "expression", ya que no existen como opciones por defecto.

# --------------------------------
# 10.1. Funciones para las previas
# --------------------------------

# Convertir precisión a desviación estándar
sigma <- 1 / sqrt(fit$summary.hyperpar$mean)

# Ver todas las previas disponibles en paquete INLA
inla.list.models("prior")

#Abre la documentación de una prior específica
#inla.doc("loggamma") #log-Gamma/loggamma/parameters = shape and rate
#inla.doc("gaussian")
#inla.doc("pc")

#Lista todos los nombres de priors disponibles en INLA
names(inla.models()$prior)

#Explora hiperparámetros de modelos latentes y lista los nombres de los hiperparámetros del modelo IID
names(inla.models()$latent$iid$hyper)
#El IID solo tiene un hiperparametro 
#theta = representa la precision t

#Nombre completo del hiperparámetro theta en el modelo IID
inla.models()$latent$iid$hyper$theta$name

#Nombre corto del mismo hiperparámetro
inla.models()$latent$iid$hyper$theta$short.name

#Toda la especificación completa del hiperparámetro theta del modelo IID
inla.models()$latent$iid$hyper$theta

# Reporte de las previas para cada uno de los hiperparametros
fit$summary.hyperpar

#Muestra todos los hiperparámetros del modelo ajustado con sus priors
fit$all.hyper

#Nombre de la prior asignada al predictor lineal
fit$predictor$hyper$theta$prior

#Los parámetros de esa prior del predictor
fit$predictor$hyper$theta$param

#Función de transformación que INLA aplica internamente al primer efecto aleatorio
fit$random[[1]]$group.hyper$theta$to.theta

#Resumen posterior de todos los hiperparámetros
fit$summary.hyperpar

#Realizar una prior desde cero, programarla con los hiperparametros y todo
prec.prior <- list(prec = list(prior = "loggamma", param = c(0.01, 0.01)),
                   initial = 4, fixed = FALSE)

# --------------------------------
# 10.1.1. Previas en el toy model
# --------------------------------
#Toy_example
# formula_inla <- deaths ~
#   f(age_idx, model = "rw1", constr = TRUE,
#     hyper = prec.prior) +          
#   f(region_idx, model = "bym2", graph = g, constr = TRUE) +
#   f(period_idx, model = "rw2", constr = TRUE,
#     hyper = prec.prior) +          
#   f(region_period_idx, model = "iid",
#     hyper = prec.prior)  

# -----------------------------
# 10.2. Previas a utilizar
# -----------------------------

# --------------------------------------
# 10.2.1. Modelo con previa Half-Cauchy 
# --------------------------------------
#Half_Cauchy en INLA bajo la función de expression (Sección 5.3) #gamma = 25;
HC.prior  = "expression:
  sigma = exp(-theta/2);
  gamma = 1;
  log_dens = log(2) - log(pi) - log(gamma);
  log_dens = log_dens - log(1 + (sigma / gamma)^2);
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);
"

cat(HC.prior)

#Implementada
formula_hc <- deaths ~
  factor(sex) +
  f(age_idx,    model = "rw1",  constr = TRUE,
    hyper = list(prec = list(prior = HC.prior))) +
  f(region_idx, model = "bym2", graph = g, constr = TRUE,
    hyper = list(
      prec = list(prior = HC.prior),
      phi  = list(prior = "logitbeta", param = c(0.5, 0.5))))+  # Beta(0.5, 0.5) en escala logarítmica
  f(period_idx, model = "rw2",  constr = TRUE,
    hyper = list(prec = list(prior = HC.prior))) +
  f(region_period_idx, model = "iid",
    hyper = list(prec = list(prior = HC.prior)))

fit_hc <- inla(formula_hc,
               family  = "poisson",
               data    = df,
               E       = population,
               control.compute = list(config = TRUE, dic = TRUE, waic = TRUE))

fit_hc$all.hyper
fit_hc$predictor$hyper$theta$prior

# e0 para la previa Half Cauchy
pred_hc <- df %>%
  mutate(
    mx = pmax(fit_hc$summary.fitted.values$mean, 1e-6)
  ) %>%
  left_join(age_params, by = "agegroup")

pred_hc$sex <- ifelse(pred_hc$sex == 1, "m", "f")
municipios <- sort(unique(pred_hc$region))
periodos   <- sort(unique(pred_hc$period))
sexos      <- c("m", "f")
Age <- c(0, 1, seq(5, 85, by = 5)) #Cambio 
tablas <- list()
for (muni in municipios) {
  for (per in periodos) {
    for (sx in sexos) {
      pred_sub_hc <- pred_hc %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub_hc$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      tablas[[muni]][[per]][[sx]] <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, a0rule = "ak", axmethod = "pas", Sex = sx, mod = FALSE)
    }
  }
}

e0_resumen_hc <- data.frame()

for (m in names(tablas)) {
  for (p in names(tablas[[m]])) {
    for (s in c("m", "f")) {
      
      tb <- tablas[[m]][[p]][[s]]
      
      sexnum <- if (s == "m") 1 else 2
      
      e0_resumen_hc <- rbind(
        e0_resumen_hc,
        data.frame(period = p, region = m, sex = sexnum, e0 = tb$ex[1])
      )
    }
  }
}

# -------------------------------------------------------------
# 10.2.2. Modelo con previa Scale Beta 2 = Beta prime escalada
# -------------------------------------------------------------
SB2.prior = "expression:
  p = 1;
  q = 1;
  b = 1;
  sigma2 = exp(-theta);
 sigma2 = exp(-theta);
  log_dens = lgamma(p+q) - lgamma(p) - lgamma(q) - log(b);
  log_dens = log_dens + (p-1) * log(sigma2/b);
  log_dens = log_dens - (p+q) * log(1 + sigma2/b);
  log_dens = log_dens - theta;
  return(log_dens);
"
cat(SB2.prior)

formula_sb2 <- deaths ~
  factor(sex) +
  f(age_idx,    model = "rw1",  constr = TRUE,
    hyper = SB2.prior) +
  f(region_idx, model = "bym2", graph = g, constr = TRUE,
    hyper = SB2.prior) +
  f(period_idx, model = "rw2",  constr = TRUE,
    hyper = SB2.prior) +
  f(region_period_idx, model = "iid",
    hyper = SB2.prior)

#Explícita 
#Implementada
# formula_sb2 <- deaths ~
#   factor(sex) +
#   f(age_idx,    model = "rw1",  constr = TRUE,
#     hyper = list(prec = list(prior = SB2.prior))) +
#   f(region_idx, model = "bym2", graph = g, constr = TRUE,
#     hyper = list(
#       prec = list(prior = SB2.prior),
#       phi  = list(prior = "logitbeta", param = c(0.5, 0.5))))+  # Beta(0.5, 0.5) en escala logarítmica
#   f(period_idx, model = "rw2",  constr = TRUE,
#     hyper = list(prec = list(prior = SB2.prior))) +
#   f(region_period_idx, model = "iid",
#     hyper = list(prec = list(prior = SB2.prior)))
# 

fit_sb2 <- inla(formula_sb2,
                family  = "poisson",
                data    = df,
                E       = population,
                control.compute = list(config = TRUE, dic = TRUE, waic = TRUE))

fit_sb2$all.hyper
fit_sb2$predictor$hyper$theta$prior #NULL

# e0 para la previa Scale beta 2
pred_sb2 <- df %>%
  mutate(
    mx = pmax(fit_sb2$summary.fitted.values$mean, 1e-6)
  ) %>%
  left_join(age_params, by = "agegroup")

pred_sb2$sex <- ifelse(pred_sb2$sex == 1, "m", "f")
municipios <- sort(unique(pred_sb2$region))
periodos   <- sort(unique(pred_sb2$period))
sexos      <- c("m", "f")
Age <- c(0, 1, seq(5, 85, by = 5)) #80 A 85
tablas <- list()
for (muni in municipios) {
  for (per in periodos) {
    for (sx in sexos) {
      pred_sub_sb2 <- pred_sb2 %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub_sb2$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      tablas[[muni]][[per]][[sx]] <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, a0rule = "ak", axmethod = "pas", Sex = sx, mod = FALSE)
    }
  }
}

e0_resumen_sb2 <- data.frame()

for (m in names(tablas)) {
  for (p in names(tablas[[m]])) {
    for (s in c("m", "f")) {
      
      tb <- tablas[[m]][[p]][[s]]
      
      sexnum <- if (s == "m") 1 else 2
      
      e0_resumen_sb2 <- rbind(
        e0_resumen_sb2,
        data.frame(period = p, region = m, sex = sexnum, e0 = tb$ex[1])
      )
    }
  }
}

# ------------------------------------------------
# 10.2.3. Modelo con previa previa Half-t Student
# ------------------------------------------------
HT.prior = "expression:
  sigma = exp(-theta/2);
  nu = 3;
  log_dens = 0 - 0.5 * log(nu * pi) - (-0.1207822);
  log_dens = log_dens - 0.5 * (nu + 1) * log(1 + sigma * sigma);
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);
"
cat(HT.prior)

formula_ht <- deaths ~
  factor(sex) +
  f(age_idx,    model = "rw1",  constr = TRUE,
    hyper = HT.prior) +
  f(region_idx, model = "bym2", graph = g, constr = TRUE, hyper = HT.prior) +  
  f(period_idx, model = "rw2",  constr = TRUE,
    hyper = HT.prior) +
  f(region_period_idx, model = "iid",
    hyper = HT.prior)

#Explícita 
#Implementada 
formula_ht <- deaths ~
  factor(sex) +
  f(age_idx,    model = "rw1",  constr = TRUE,
    hyper = list(prec = list(prior = HT.prior))) +
  f(region_idx, model = "bym2", graph = g, constr = TRUE,
    hyper = list(
      prec = list(prior = HT.prior)))+  
  f(period_idx, model = "rw2",  constr = TRUE,
    hyper = list(prec = list(prior = HT.prior))) +
  f(region_period_idx, model = "iid",
    hyper = list(prec = list(prior = HT.prior)))


fit_ht <- inla(formula_ht,
               family  = "poisson",
               data    = df,
               E       = population,
               control.compute = list(config = TRUE, return.marginals.predictor=TRUE, cpo = TRUE, dic = TRUE, waic = TRUE)) 

# e0 para la previa Half T
pred_ht <- df %>%
  mutate(
    mx = pmax(fit_ht$summary.fitted.values$mean, 1e-6)
  ) %>%
  left_join(age_params, by = "agegroup")

pred_ht$sex <- ifelse(pred_ht$sex == 1, "m", "f")
municipios <- sort(unique(pred_ht$region))
periodos   <- sort(unique(pred_ht$period))
sexos      <- c("m", "f")
Age <- c(0, 1, seq(5, 85, by = 5)) #80 a 85 
tablas <- list()
for (muni in municipios) {
  for (per in periodos) {
    for (sx in sexos) {
      pred_sub_ht <- pred_ht %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub_ht$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      tablas[[muni]][[per]][[sx]] <- lt_abridged(nMx = nMx, AgeInt = AgeInt, a0rule = "ak", axmethod = "pas", Age = Age, Sex = sx, mod = FALSE)
    }
  }
}

e0_resumen_ht <- data.frame()

for (m in names(tablas)) {
  for (p in names(tablas[[m]])) {
    for (s in c("m", "f")) {
      
      tb <- tablas[[m]][[p]][[s]]
      
      sexnum <- if (s == "m") 1 else 2
      
      e0_resumen_ht <- rbind(
        e0_resumen_ht,
        data.frame(period = p, region = m, sex = sexnum, e0 = tb$ex[1])
      )
    }
  }
}

# ----------------------------------------
# 10.2.4. Modelo con previa Inverse Gamma
# ----------------------------------------
IG.prior = "expression:
  a = 1;
  b = 0.00005;
  log_dens = a * log(b) - lgamma(a) + a * theta - b * exp(theta);
  return(log_dens);
"

#Implementada
formula_ig <- deaths ~
  factor(sex) +
  f(age_idx,    model = "rw1",  constr = TRUE,
    hyper = list(prec = list(prior = IG.prior))) +
  f(region_idx, model = "bym2", graph = g, constr = TRUE,
    hyper = list(prec = list(prior = IG.prior))) +
  f(period_idx, model = "rw2",  constr = TRUE,
    hyper = list(prec = list(prior = IG.prior))) +
  f(region_period_idx, model = "iid",
    hyper = list(prec = list(prior = IG.prior)))

fit_ig <- inla(formula_ig,
               family  = "poisson",
               data    = df,
               E       = population,
               control.compute = list(config = TRUE, cpo = TRUE, dic = TRUE, waic = TRUE))


# e0 para la previa Inverse Gamma
pred_ig <- df %>%
  mutate(
    mx = pmax(fit_ig$summary.fitted.values$mean, 1e-6)
  ) %>%
  left_join(age_params, by = "agegroup")

pred_ig$sex <- ifelse(pred_ig$sex == 1, "m", "f")
municipios <- sort(unique(pred_ig$region))
periodos   <- sort(unique(pred_ig$period))
sexos      <- c("m", "f")
Age <- c(0, 1, seq(5, 85, by = 5)) #80 a 85
tablas <- list()
for (muni in municipios) {
  for (per in periodos) {
    for (sx in sexos) {
      pred_sub_ig <- pred_ig %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub_ig$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      tablas[[muni]][[per]][[sx]] <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, a0rule = "ak", axmethod = "pas", Sex = sx, mod = FALSE)
    }
  }
}

e0_resumen_ig <- data.frame()

for (m in names(tablas)) {
  for (p in names(tablas[[m]])) {
    for (s in c("m", "f")) {
      
      tb <- tablas[[m]][[p]][[s]]
      
      sexnum <- if (s == "m") 1 else 2
      
      e0_resumen_ig <- rbind(
        e0_resumen_ig,
        data.frame(period = p, region = m, sex = sexnum, e0 = tb$ex[1])
      )
    }
  }
}

# ----------------------------------------
# 10.2.5. Modelo con previa PC
# ----------------------------------------
#Previas explícitas por cada parámetro e hiperparámetro de PC priors 
formula_pc <- deaths ~
  factor(sex) +
  f(age_idx, model = "rw1", constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))) +
  f(region_idx, model = "bym2", graph = g, constr = TRUE, 
    hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))) +
  f(period_idx, model = "rw2", constr = TRUE,
    hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01)))) +
  f(region_period_idx, model = "iid",
    hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01))))

fit_pc <- inla(
  formula = formula_pc,
  family = "poisson",
  data = df,
  E = population,
  control.predictor = list(compute = TRUE),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE)
)

pred_pc <- df %>%
  mutate(
    mx = pmax(fit_pc$summary.fitted.values$mean, 1e-6)
  ) %>%
  left_join(age_params, by = "agegroup")

pred_pc$sex <- ifelse(pred_pc$sex == 1, "m", "f")
municipios <- sort(unique(pred_pc$region))
periodos   <- sort(unique(pred_pc$period))
sexos      <- c("m", "f")
Age <- c(0, 1, seq(5, 85, by = 5)) #80 a 85
tablas <- list()
for (muni in municipios) {
  for (per in periodos) {
    for (sx in sexos) {
      pred_sub_pc <- pred_pc %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub_pc$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      tablas[[muni]][[per]][[sx]] <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, a0rule = "ak", axmethod = "pas", Sex = sx, mod = FALSE)
    }
  }
}

e0_resumen_pc <- data.frame()

for (m in names(tablas)) {
  for (p in names(tablas[[m]])) {
    for (s in c("m", "f")) {
      
      tb <- tablas[[m]][[p]][[s]]
      
      sexnum <- if (s == "m") 1 else 2
      
      e0_resumen_pc <- rbind(
        e0_resumen_pc,
        data.frame(period = p, region = m, sex = sexnum, e0 = tb$ex[1])
      )
    }
  }
}

##########################################################################
##########################################################################
                #Tuning 11.3 (HOLD: Trabajo para futuro)
##########################################################################
##########################################################################

###
# Parametros de prueba. No tienen sustento en la literatura
param_gamma_age           <- list(2)
param_gamma_region        <- list(1)
param_gamma_period        <- list(0.5)
param_gamma_region_period <- list(0.25)

param_pqb_age             <- list(c(1, 1, 1))
param_pqb_region          <- list(c(1, 1, 0.5))
param_pqb_period          <- list(c(1, 1, 0.25))
param_pqb_region_period   <- list(c(1, 1, 0.1))

param_ua_age              <- list(c(2, 0.01))
param_ua_region           <- list(c(1, 0.01))
param_ua_period           <- list(c(0.5, 0.01))
param_ua_region_period    <- list(c(0.25, 0.01))

make_hc_prior <- function(gamma) {
  sprintf(
    "expression:
      sigma = exp(-theta/2);
      gamma = %f;
      log_dens = log(2) - log(pi) - log(gamma);
      log_dens = log_dens - log(1 + (sigma/gamma)^2);
      log_dens = log_dens - log(2) - theta/2;
      return(log_dens);
    ",
    gamma
  )
}

make_sb2_prior <- function(p, q, b) {
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

construct_prior <- function(familia, param, k) {
  switch(familia,
         "half.cauchy" = list(prior = make_hc_prior(gamma = param[[k]][1])),
         "scale.beta2" = list(prior = make_sb2_prior(p = param[[k]][1],
                                                     q = param[[k]][2],
                                                     b = param[[k]][3])),
         "pc.prior"    = list(prior = "pc.prec", param = c(param[[k]][1], param[[k]][2]))
  )
}

# Ejemplos de uso
construct_prior("scale.beta2", param_pqb_age, 1)
construct_prior("half.cauchy", param_gamma_region, 1)
construct_prior("pc.prior", param_ua_period, 1)
construct_prior("pc.prior", list(c(1, 0.1)), 1)

ajuste_modelo_inla <- function(defun, 
                               adya,
                               k, l, m, n,
                               previa_age,
                               previa_region,
                               previa_period,
                               previa_region_period,
                               param_age,
                               param_region,
                               param_period,
                               param_region_period){
  formula_dinamic <- deaths ~
    factor(sex) +
    f(age_idx, model = "rw1", constr = TRUE,
      hyper = list(prec = construct_prior(previa_age, param_age, k))) +
    f(region_idx, model = "bym2", graph = adya, constr = TRUE,
      hyper = list(prec = construct_prior(previa_region, param_region, l),
                   phi  = list(prior = "logitbeta", param = c(0.5, 0.5)))) + # EGR: ¿Puedo agregar esa previa para phi?
    f(period_idx, model = "rw2", constr = TRUE,
      hyper = list(prec = construct_prior(previa_period, param_period, m))) +
    f(region_period_idx, model = "iid",
      hyper = list(prec = construct_prior(previa_region_period, param_region_period, n)))
  
  fit_dynamic <- inla(
    formula           = formula_dinamic, 
    family            = "poisson", 
    data              = defun, 
    E                 = population, 
    control.predictor = list(compute = TRUE),
    control.compute   = list(config = TRUE, dic = TRUE, waic = TRUE)
  )
}

# Ejemplo de uso
fit_dynamic <- ajuste_modelo_inla(defun = df, 
                                  adya = g,
                                  k = 1, l = 1, m = 1, n = 1,
                                  previa_age = "half.cauchy",
                                  previa_region = "scale.beta2",                   
                                  previa_period = "pc.prior", 
                                  previa_region_period = "scale.beta2",
                                  param_age             = param_gamma_age,
                                  param_region          = param_pqb_region,
                                  param_period          = param_ua_period,                   
                                  param_region_period   = param_pqb_region_period
)
# devtools::install_github("julianfaraway/brinla")
# library(brinla)
# bri.hyperpar.plot(fit_dynamic, together = T)

# -----------------------------
# 11. Análisis de sensitividad 
# -----------------------------

# Corre muy lento. Para 81 modelos tardo cerca de 15 minutos. Esto se puede optimizar usando 
# AI. Terminando este chunk hay dos optimizaciones que tardan mucho menos
# EGR. El comentario de arriba fue en la noche. Ahora en la manania tengo esta noticia.
# Al borrar el environment y ejecutar de nuevo el csize.models(), tardo 4.6 minutos.
# En conlusion, el codigo funciona bien y su tiempo de ejecucion es comparable con el dado por AI.
# Se tendria que probar con mas juegos de parametros para ver si el do.call() lo sigue 
# haciendo tan bien como los lapply(). Aun asi, para juegos de previas
# individuales y parametros individuales (es decir, si solo se quiere ver el DIC o WAIC de una 
# previas y su parametros para cada efecto), se puede correr solo la funcion de ajuste_modelo_inla().
# Espero que esto les sea de utilidad.
fam <- c(half.cauchy="param_gamma_", scale.beta2="param_pqb_", pc.prior="param_ua_")
efectos <- c(age = "age", region = "region", period = "period", region_period = "region_period")

combinaciones <- expand.grid(
  age           = names(fam),
  region        = names(fam),
  period        = names(fam),
  region_period = names(fam)
)

nombres_previas <- apply(combinaciones, 1, function(x) {
  paste(sapply(names(efectos), function(efecto) {
    previa <- x[efecto]
    param  <- get(paste0(fam[previa], efecto))
    paste0(previa, " ", efectos[efecto], " (", paste(unlist(param), collapse = ", "), ")")
  }),collapse = " | ")
})

csize.models <- do.call(c, lapply(seq_len(nrow(combinaciones)), function(i) {
  a <- as.character(combinaciones$age[i])
  r <- as.character(combinaciones$region[i])
  p <- as.character(combinaciones$period[i])
  s <- as.character(combinaciones$region_period[i])
  param_age <- get(paste0(fam[a], "age"))
  param_region <- get(paste0(fam[r], "region"))
  param_period <- get(paste0(fam[p], "period"))
  param_region_period <- get(paste0(fam[s], "region_period"))
  indices <- expand.grid(k = seq_along(param_age),
                         l = seq_along(param_region),
                         m = seq_along(param_period),
                         n = seq_along(param_region_period))
  modelos <- lapply(seq_len(nrow(indices)), function(j) {
    ajuste_modelo_inla(defun = df,
                       adya = g,
                       k = indices$k[j], l = indices$l[j], m = indices$m[j], n = indices$n[j],
                       previa_age = a,
                       previa_region = r,
                       previa_period = p,
                       previa_region_period = s,
                       param_age = param_age, 
                       param_region = param_region,
                       param_period = param_period,
                       param_region_period = param_region_period)
  })
  names(modelos) <- apply(indices, 1, function(x) {
    paste(paste0(a, " age (", paste(param_age[[x["k"]]], collapse = ", "), ")"),
          paste0(r, " region (", paste(param_region[[x["l"]]], collapse = ", "), ")"),
          paste0(p, " period (", paste(param_period[[x["m"]]], collapse = ", "), ")"),
          paste0(s, " region_period (", paste(param_region_period[[x["n"]]], collapse = ", "), ")"),
          sep = " | ")
  })
  modelos
}))

sensitivity_analysis_summary<- tibble(
  prior = names(csize.models),
  DIC   = sapply(csize.models, function(m) m$dic$dic),
  WAIC  = sapply(csize.models, function(m) m$waic$waic),
  sigma_age_idx = sapply(csize.models, function(m) {
    1 / sqrt(m$summary.hyperpar["Precision for age_idx", "mean"])
  }),
  sigma_region_idx = sapply(csize.models, function(m) {
    1 / sqrt(m$summary.hyperpar["Precision for region_idx", "mean"])
  }),
  sigma_period_idx = sapply(csize.models, function(m) {
    1 / sqrt(m$summary.hyperpar["Precision for period_idx", "mean"])
  }),
  sigma_region_period_idx = sapply(csize.models, function(m) {
    1 / sqrt(m$summary.hyperpar["Precision for region_period_idx", "mean"])
  })
)

print(sensitivity_analysis_summary)

###########################################
# Optimizacion 1 con AI. Tarda 5.5 minutos
###########################################
fam <- c(half.cauchy="param_gamma_", scale.beta2="param_pqb_", pc.prior="param_ua_")
op <- function(ef) unlist(lapply(names(fam), \(f) paste(f, seq_along(get(paste0(fam[f], ef))))))
G <- expand.grid(age=op("age"), region=op("region"), period=op("period"),
                 region_period=op("region_period"), stringsAsFactors=FALSE)

FA <- function(t) sub(" .*","",t)
KK <- function(t) as.integer(sub(".* ","",t))
PA <- function(ef,t) get(paste0(fam[FA(t)], ef))
tx <- function(ef,t) paste0(FA(t)," ",ef," (",paste(PA(ef,t)[[KK(t)]],collapse=", "),")")

csize.models <- lapply(seq_len(nrow(G)), \(i) {
  a<-G$age[i]; r<-G$region[i]; p<-G$period[i]; s<-G$region_period[i]
  ajuste_modelo_inla(defun=df, adya=g, k=KK(a), l=KK(r), m=KK(p), n=KK(s),
                     previa_age=FA(a), previa_region=FA(r), previa_period=FA(p), previa_region_period=FA(s),
                     param_age=PA("age",a), param_region=PA("region",r),
                     param_period=PA("period",p), param_region_period=PA("region_period",s))
})

names(csize.models) <- sapply(seq_len(nrow(G)), \(i)
                              paste(tx("age",G$age[i]), tx("region",G$region[i]),
                                    tx("period",G$period[i]), tx("region_period",G$region_period[i]), sep=" | "))

sig <- function(m,p) 1/sqrt(m$summary.hyperpar[p,"mean"])
sensitivity_analysis_summary <- tibble(
  prior=names(csize.models),
  DIC =sapply(csize.models,\(m) m$dic$dic),
  WAIC=sapply(csize.models,\(m) m$waic$waic),
  s_age=sapply(csize.models,sig,"Precision for age_idx"),
  s_reg=sapply(csize.models,sig,"Precision for region_idx"),
  s_per=sapply(csize.models,sig,"Precision for period_idx"),
  s_int=sapply(csize.models,sig,"Precision for region_period_idx"))
print(sensitivity_analysis_summary)

###########################################
# Optimizacion 2 con AI. Tarda 4.5 minutos
###########################################
fam <- c(half.cauchy="param_gamma_", scale.beta2="param_pqb_", pc.prior="param_ua_")
op <- function(ef) unlist(lapply(names(fam), \(f) paste(f, seq_along(get(paste0(fam[f], ef))))))
G <- expand.grid(age=op("age"), region=op("region"), period=op("period"),
                 region_period=op("region_period"), stringsAsFactors=FALSE)
pr <- function(ef, t) {
  f <- sub(" .*","",t); p <- get(paste0(fam[f], ef))[[as.integer(sub(".* ","",t))]]
  if (f=="half.cauchy") list(prior=make_hc_prior(p[1]))
  else if (f=="scale.beta2") list(prior=make_sb2_prior(p[1],p[2],p[3]))
  else list(prior="pc.prec", param=p)
}
tx <- function(ef, t) { f<-sub(" .*","",t)
paste0(f," ",ef," (",paste(get(paste0(fam[f],ef))[[as.integer(sub(".* ","",t))]],collapse=", "),")") }

csize.models <- lapply(seq_len(nrow(G)), \(i) inla(
  deaths ~ factor(sex) +
    f(age_idx, model="rw1", constr=TRUE, hyper=list(prec=pr("age", G$age[i]))) +
    f(region_idx, model="bym2", graph=g, constr=TRUE,
      hyper=list(prec=pr("region", G$region[i]), phi=list(prior="logitbeta", param=c(.5,.5)))) +
    f(period_idx, model="rw2", constr=TRUE, hyper=list(prec=pr("period", G$period[i]))) +
    f(region_period_idx, model="iid", hyper=list(prec=pr("region_period", G$region_period[i]))),
  family="poisson", data=df, E=population,
  control.compute=list(dic=TRUE, waic=TRUE)))

names(csize.models) <- sapply(seq_len(nrow(G)), \(i)
                              paste(tx("age",G$age[i]), tx("region",G$region[i]),
                                    tx("period",G$period[i]), tx("region_period",G$region_period[i]), sep=" | "))

sig <- function(m,p) 1/sqrt(m$summary.hyperpar[p,"mean"])
sensitivity_analysis_summary <- tibble(
  prior=names(csize.models),
  DIC =sapply(csize.models,\(m) m$dic$dic),
  WAIC=sapply(csize.models,\(m) m$waic$waic),
  s_age=sapply(csize.models,sig,"Precision for age_idx"),
  s_reg=sapply(csize.models,sig,"Precision for region_idx"),
  s_per=sapply(csize.models,sig,"Precision for period_idx"),
  s_int=sapply(csize.models,sig,"Precision for region_period_idx"))
print(sensitivity_analysis_summary)

###











































#OTRA FORMA: 
combinaciones_sb2 <- list(
  c(1, 1, 1),
  c(0.5, 0.5, 1),
  c(2, 2, 1),
  c(1, 1, 5)
)

modelos_sb2 <- lapply(combinaciones_sb2, function(param) {
  ajuste_modelo_inla(df, g,
                     familia_age = "scale.beta2", param_age = param)
})
names(modelos_sb2) <- sapply(combinaciones_sb2, paste, collapse = "_")

####

































































# -------------
# 11.3 Tuning
# -------------
#Half Cauchy
make_hc_prior <- function(gamma = 1) {
  sprintf(
    "expression:
      sigma = exp(-theta/2);
      gamma = %f;
      log_dens = log(2) - log(pi) - log(gamma);
      log_dens = log_dens - log(1 + (sigma/gamma)^2);
      log_dens = log_dens - log(2) - theta/2;
      return(log_dens);
    ",
    gamma
  )
}

#Cambio de parámetros HC
HC.prior_gamma1   <- make_hc_prior(gamma = 1)
HC.prior_gamma2.5 <- make_hc_prior(gamma = 2.5)
HC.prior_gamma0.5 <- make_hc_prior(gamma = 0.5)


#Scale Beta2 
#Nota: con p=q=1, la SBeta2(1,1,b) es prácticamente equivalente a la Half-Cauchy (ver Sección 4.3 del artículo de Pérez, Pericchi y Ramírez 2017)
make_sb2_prior <- function(p = 1, q = 1, b = 1) {
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

#Scale Beta 2 equivalente aprox. a Half-Cauchy
SB2.prior_p1_q1_b1     <- make_sb2_prior(p = 1, q = 1, b = 1) #p=q=1

#Variantes con colas más ligeras o más pesadas (p, q controlan la forma)
SB2.prior_p0.5_q0.5_b1 <- make_sb2_prior(p = 0.5, q = 0.5, b = 1) 
SB2.prior_p2_q2_b1     <- make_sb2_prior(p = 2, q = 2, b = 1)

#Variantes con distinta escala (b controla dónde se concentra la masa)
SB2.prior_p1_q1_b0.5   <- make_sb2_prior(p = 1, q = 1, b = 0.5)
SB2.prior_p1_q1_b5     <- make_sb2_prior(p = 1, q = 1, b = 5)


#Half-t 
make_halft_prior <- function(nu = 3) {
  const <- lgamma((nu+1)/2) - lgamma(nu/2)
  sprintf(
    "expression:
      sigma = exp(-theta/2);
      nu = %f;
      log_dens = %f - 0.5 * log(nu * pi);
      log_dens = log_dens - 0.5 * (nu + 1) * log(1 + sigma * sigma);
      log_dens = log_dens - log(2) - theta / 2;
      return(log_dens);
    ",
    nu, const
  )
}

# Cambio de parámetros HT
HT.prior_nu3 <- make_halft_prior(nu = 3)
HT.prior_nu5 <- make_halft_prior(nu = 5)
HT.prior_nu7 <- make_halft_prior(nu = 7)
HT.prior_nu9 <- make_halft_prior(nu = 9)

#Inverse Gamma
make_ig_prior <- function(a = 1, b = 0.00005) {
  sprintf(
    "expression:
      a = %f;
      b = %.10f;
      log_dens = a * log(b) - lgamma(a);
      log_dens = log_dens + a * theta;
      log_dens = log_dens - b * exp(theta);
      return(log_dens);
    ",
    a, b
  )
}

#Cambio de parámetros IG
IG.prior_a1_b0.00005 <- make_ig_prior(a = 1, b = 0.00005)
IG.prior_a1_b0.0001  <- make_ig_prior(a = 1, b = 0.0001)
IG.prior_a0.5_b0.0001 <- make_ig_prior(a = 0.5, b = 0.0001)
IG.prior_a0.5_b0.00001 <- make_ig_prior(a = 0.5, b = 0.00001)

#PC prior
#Cambio de parámetros PC 
#No necesita crear la expresion ni make_prior() porque INLA lo acepta nativo con param=c(u, alpha), ver inla.doc("pc.prec")
#hyper = list(<theta> = list(prior="pc.prec", param=c(<u>,<alpha>)))
PC.prior_u0.1_a0.01 <- list(prior = "pc.prec", param = c(0.1, 0.01))
PC.prior_u5_a0.01   <- list(prior = "pc.prec", param = c(5, 0.01))
PC.prior_u0.5_a0.01 <- list(prior = "pc.prec", param = c(0.5, 0.01))
PC.prior_u0.5_a0.1   <- list(prior = "pc.prec", param = c(0.5, 0.1))
PC.prior_u5_a0.1   <- list(prior = "pc.prec", param = c(5, 0.01))
PC.prior_u10_a0.1   <- list(prior = "pc.prec", param = c(10, 0.1))
PC.prior_u10_a0.01   <- list(prior = "pc.prec", param = c(10, 0.01))

class(PC.prior_u0.1_a0.01) #Solo cambia los argumentos de la expresión, no crea un objeto inla (modelo)

pc_age    <- c(0.1, 0.01)
pc_region <- c(0.1, 0.01)
pc_period <- c(0.1, 0.01)

formula <- deaths ~
  f(age_idx, model = "rw1",
    hyper = list(prec = list(prior = "pc.prec", param = pc_age))) +
  f(region_idx, model = "bym2", graph = matriz_adj,
    hyper = list(prec = list(prior = "pc.prec", param = pc_region))) +
  f(period_idx, model = "rw2",
    hyper = list(prec = list(prior = "pc.prec", param = pc_period)))


#Nota: Se deberá crear una formula_inla y un modelo (fit_inla) para cada una.

###########################################
#OTRA FORMA: Cambiar los parámetros dentro 
###########################################
# Half-Cauchy
make_hc_prior <- function(gamma = 1) {
  sprintf(
    "expression:
      sigma = exp(-theta/2);
      gamma = %f;
      log_dens = log(2) - log(pi) - log(gamma);
      log_dens = log_dens - log(1 + (sigma/gamma)^2);
      log_dens = log_dens - log(2) - theta/2;
      return(log_dens);
    ",
    gamma
  )
}

#Scale Beta2
make_sb2_prior <- function(p = 1, q = 1, b = 1) {
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

#Half-t
make_halft_prior <- function(nu = 3) {
  const <- lgamma((nu+1)/2) - lgamma(nu/2)
  sprintf(
    "expression:
      sigma = exp(-theta/2);
      nu = %f;
      log_dens = %f - 0.5 * log(nu * pi);
      log_dens = log_dens - 0.5 * (nu + 1) * log(1 + sigma * sigma);
      log_dens = log_dens - log(2) - theta / 2;
      return(log_dens);
    ",
    nu, const
  )
}

#Inverse Gamma
make_ig_prior <- function(a = 1, b = 0.00005) {
  sprintf(
    "expression:
      a = %f;
      b = %.10f;
      log_dens = a * log(b) - lgamma(a);
      log_dens = log_dens + a * theta;
      log_dens = log_dens - b * exp(theta);
      return(log_dens);
    ",
    a, b
  )
}

#Cambiar previas para el modelo ajustado en INLA

#Argumento del string familia implica el switch de la previa que se quiera emplear
construct_prior <- function(familia, param) {
  switch(familia,
         "pc.prec"       = list(prior = "pc.prec", param = param),
         "half.cauchy"   = list(prior = make_hc_prior(gamma = param[1])),
         "scale.beta2"   = list(prior = make_sb2_prior(p = param[1], q = param[2], b = param[3])),
         "half.t"        = list(prior = make_halft_prior(nu = param[1])),
         "inverse.gamma" = list(prior = make_ig_prior(a = param[1], b = param[2])),
         stop("Familia de previa no reconocida: ", familia)
  )
}
ajuste_modelo_inla <- function(df, 
                               g,
                               familia_age           = "pc.prec",
                               familia_region         = "pc.prec",
                               familia_period         = "pc.prec",
                               familia_region_period  = "pc.prec",
                               param_age              = c(0.1, 0.01),
                               param_region           = c(0.1, 0.01),
                               param_period           = c(0.1, 0.01),
                               param_region_period    = c(0.1, 0.01)) {
  
  formula_dinamic <- deaths ~
    factor(sex) +
    f(age_idx, model = "rw1", constr = TRUE,
      hyper = list(prec = construct_prior(familia_age, param_age))) +
    f(region_idx, model = "bym2", graph = g, constr = TRUE,
      hyper = list(prec = construct_prior(familia_region, param_region))) +
    f(period_idx, model = "rw2", constr = TRUE,
      hyper = list(prec = construct_prior(familia_period, param_period))) +
    f(region_period_idx, model = "iid",
      hyper = list(prec = construct_prior(familia_region_period, param_region_period)))
  
  inla(
    formula           = formula_dinamic, 
    family            = "poisson", 
    data              = df, 
    E                 = population, 
    control.predictor = list(compute = TRUE),
    control.compute   = list(config = TRUE, dic = TRUE, waic = TRUE)
  )
}

##########################################################################                       
##########################################################################
#HOLD: Para considerar (cambiar previas para cada efecto)
#Usar diferentes previas en un solo modelo
formula_mixed_priors <- deaths ~
  factor(sex) +
  f(age_idx, model = "rw1", constr = TRUE,
    hyper = list(prec = list(prior = HC.prior))) +
  f(region_idx, model = "bym2", graph = g, constr = TRUE,
    hyper = list(prec = list(prior = HT.prior))) +
  f(period_idx, model = "rw2", constr = TRUE,
    hyper = list(prec = list(prior = SB2.prior))) +
  f(region_period_idx, model = "iid",
    hyper = list(prec = list(prior = IG.prior)))

fit_mixed_priors <- inla(
  formula = formula_mixed_priors,
  family = "poisson",
  data = df,
  E = population,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

# #Otra forma de ver un modelo con lo empleado arriba
# modelo <- ajuste_modelo_inla(
#   df, g,
#   familia_age = "half.cauchy", param_age = c(1),
#   familia_region = "pc.prec", param_region = c(0.1, 0.01),
#   familia_period = "half.t", param_period = c(3),
#   familia_region_period = "scale.beta2", param_region_period = c(1, 1, 1)
# )

##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
# ------------------------------------------------------------------
# 12. Análisis de sensitividad con previas y parámetros por defecto
# ------------------------------------------------------------------
#Ejemplo: Libro de INLA - Sección 5.5

prior.list = list(
  default = list(prec = list(prior = "loggamma", param = c(0.1, 0.1))), #a = 0.1, b = 0.1
  pc.prec = list(prec = list(prior = "pc.prec", param = c(0.1, 0.01))), #u = 0.1, alpha = 0.01
  half.cauchy = list(prec = list(prior = HC.prior)),                    #gamma = 1;  a = 0.5, b = 0.5
  half.t = list(prec = list(prior = HT.prior)),                         #nu = 3
  scale.beta2 = list(prec = list(prior = SB2.prior)),                   #a = 0.5, b = 0.5
  inverse.gamma = list(prec = list(prior = IG.prior))                   #a = 1, b = 0.00005
)

csize.models <- lapply(prior.list, function(tau.prior) {
  inla(deaths ~
         factor(sex) +
         f(age_idx, model = "rw1", constr = TRUE,
           hyper = tau.prior) +
         f(region_idx, model = "bym2", graph = g, constr = TRUE) +
         f(period_idx, model = "rw2", constr = TRUE,
           hyper = tau.prior) +
         f(region_period_idx, model = "iid",
           hyper = tau.prior),
       family = "poisson",
       data = df,
       E = population,
       control.predictor = list(compute = TRUE),
       control.compute = list(dic = TRUE, waic = TRUE),
      factor(sex) +
      f(age_idx, model = "rw1", constr = TRUE,
        hyper = tau.prior) +
      f(region_idx, model = "bym2", graph = g, constr = TRUE) +
      f(period_idx, model = "rw2", constr = TRUE,
        hyper = tau.prior) +
      f(region_period_idx, model = "iid",
        hyper = tau.prior),
    family = "poisson",
    data = df,
    E = population,
    control.predictor = list(compute = TRUE),
    control.compute = list(config = TRUE, dic = TRUE, waic = TRUE)
  )
})

sensitivity_analysis_summary<- tibble(
  prior = names(csize.models),
  DIC   = sapply(csize.models, function(m) m$dic$dic),
  WAIC  = sapply(csize.models, function(m) m$waic$waic),
  sigma_age_idx = sapply(csize.models, function(m) {
    1 / sqrt(m$summary.hyperpar["Precision for age_idx", "mean"])
  }),
  sigma_period_idx = sapply(csize.models, function(m) {
    1 / sqrt(m$summary.hyperpar["Precision for period_idx", "mean"])
  }),
  sigma_region_period_idx = sapply(csize.models, function(m) {
    1 / sqrt(m$summary.hyperpar["Precision for region_period_idx", "mean"])
  })
)

print(sensitivity_analysis_summary)




#################################################################################
#################################################################################
#e0 observado (directo,crudo, no se suavizan las tasas)
pred_directo <- df %>%
  mutate(
    mx = pmax(deaths / population, 1e-6)   # <- mx crudo, no suavizado
  ) %>%
  left_join(age_params, by = "agegroup")

pred_directo$sex <- ifelse(pred_directo$sex == 1, "m", "f")
municipios <- sort(unique(pred_directo$region))
periodos   <- sort(unique(pred_directo$period))
sexos      <- c("m", "f")
Age <- c(0, 1, seq(5, 85, by = 5))

tablas_directo <- list()
for (muni in municipios) {
  for (per in periodos) {
    for (sx in sexos) {
      pred_sub <- pred_directo %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      tablas_directo[[muni]][[per]][[sx]] <- lt_abridged(
        nMx = nMx, AgeInt = AgeInt, Age = Age, 
        a0rule = "ak", axmethod = "pas", Sex = sx, mod = FALSE
      )
    }
  }
}

e0_resumen_directo <- data.frame()
for (m in names(tablas_directo)) {
  for (p in names(tablas_directo[[m]])) {
    for (s in c("m", "f")) {
      tb <- tablas_directo[[m]][[p]][[s]]
      sexnum <- if (s == "m") 1 else 2
      e0_resumen_directo <- rbind(
        e0_resumen_directo,
        data.frame(period = p, region = m, sex = sexnum, e0_observada = tb$ex[1])
      )
    }
  }
}
e0_resumen_directo
#################################################################################
#################################################################################

# --------------------------------
# 13. Sampling from the posterior
# --------------------------------
#Ejemplo: Libro de INLA - Sección 2,7; Figura 2.8

#Histrogramas de hiperparámetros desde 100, 1,000 hasta 10,000 muestras.

#PC prior
#Half Cauchy
#Half T-Student
#Scale beta 2
#Inverse gamma

class(fit_pc)
class(fit_hc)
class(fit_ht)
class(fit_sb2)
class(fit_ig)

head(inla.hyperpar.sample(100,fit_pc))
head(inla.hyperpar.sample(100,fit_hc))
head(inla.hyperpar.sample(100,fit_ht))
head(inla.hyperpar.sample(100,fit_sb2))
head(inla.hyperpar.sample(100,fit_ig))

hist(inla.hyperpar.sample(100,fit_pc))
hist(inla.hyperpar.sample(100,fit_hc))
hist(inla.hyperpar.sample(100,fit_ht))
hist(inla.hyperpar.sample(100,fit_sb2))
hist(inla.hyperpar.sample(100,fit_ig))

hist(inla.hyperpar.sample(1000,fit_pc))
hist(inla.hyperpar.sample(1000,fit_hc))
hist(inla.hyperpar.sample(1000,fit_ht))
hist(inla.hyperpar.sample(1000,fit_sb2))
hist(inla.hyperpar.sample(1000,fit_ig))

hist(inla.hyperpar.sample(10000,fit_pc))
hist(inla.hyperpar.sample(10000,fit_hc))
hist(inla.hyperpar.sample(10000,fit_ht))
hist(inla.hyperpar.sample(10000,fit_sb2))
hist(inla.hyperpar.sample(10000,fit_ig))


#Calcula el intervalo de mayor densidad posterior (HDP)
inla.hpdmarginal()

#Ver qué marginales están disponibles 
names(fit_pc$marginals.hyperpar)

#Ejemplo: HPD al 95% para la precisión del efecto de edad
hpd_age <- inla.hpdmarginal(0.95, fit_pc$marginals.hyperpar$`Precision for age_idx`)
print(hpd_age)

#Transformar el marginal de precisión a sigma usando inla.tmarginal()
marginal_sigma_age <- inla.tmarginal(function(x) 1/sqrt(x), 
                                     fit_pc$marginals.hyperpar$`Precision for age_idx`)
hpd_sigma_age <- inla.hpdmarginal(0.95, marginal_sigma_age)
print(hpd_sigma_age)

#Comparar las HPD de cada previa
modelos_previa <- list(
  "PC prior"      = fit_pc,
  "Half-Cauchy"   = fit_hc,
  "Half-t"        = fit_ht,
  "Scale Beta2"   = fit_sb2,
  "Inverse Gamma" = fit_ig
)

hpd_comparacion <- lapply(names(modelos_previa), function(nombre) {
  m <- modelos_previa[[nombre]]
  marginal_sigma <- inla.tmarginal(function(x) 1/sqrt(x), 
                                   m$marginals.hyperpar$`Precision for age_idx`)
  hpd <- inla.hpdmarginal(0.95, marginal_sigma)
  data.frame(prior = nombre, hpd_low = hpd[1], hpd_high = hpd[2])
}) %>% bind_rows()
print(hpd_comparacion)

# -----------------------
# 14. Gráficos
# -----------------------

###Población completa
population_plot <- function(df, per) {
  plot <- ggplot(df %>%
                   filter(period == per), aes(x = fct_reorder(region, population), y = population)) +
    geom_col(fill = "firebrick") + 
    scale_y_continuous() +
    coord_flip() +
    labs(title = "Población expuesta al riesgo por municipio",
         subtitle = paste0("Puerto Rico, ", per),
         x = "",
         y = "Población") + theme_minimal() +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 8))
  return(plot)
}

population_plot(df, "2020-2024")



###Tasa de mortalidad agrupada por períodos quinquenales
#OJO:Cambiar la escala de mx

mx_period_plot <- function(dat, muni) {
  
  plot <- ggplot(dat %>% filter(region == muni), aes(x = factor(period),
                                                     y = mx * 1000,
                                                     group = factor(agegroup),
                                                     color = factor(agegroup))) +
    geom_line() + geom_point() + theme_bw() +
    scale_x_discrete() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_color_discrete("Grupo de edad") +
    labs(x = "Periodo", y = "mx",
         title = paste0("Age Specific Mortality Rate (per 1,000)", ", ", muni))
  return(plot)
  
}

mx_period_plot(pred_pc, "San Juan")

###OJO: Intenté cambiarla, corroborar si está correcto!
mx_period_plot <- function(dat, muni) {
  
  plot <- ggplot(dat %>% filter(region == muni), aes(x = factor(period),
                                                     y = mx,
                                                     group = factor(agegroup),
                                                     color = factor(agegroup))) +
    geom_line() + geom_point() + theme_bw() +
    scale_x_discrete() +
    scale_y_log10(labels = scales::comma) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_color_discrete("Grupo de edad") +
    labs(x = "Periodo", y = "mx (escala log10)",
         title = paste0("Age Specific Mortality Rate (escala log), ", muni))
  return(plot)
  
}
mx_period_plot(pred_pc, "San Juan")


####Tasa de mortalidad por grupos quinquenales
mx_municipio <- function(dat, muni) {
  plot <- ggplot(dat %>% filter(region == muni), aes(x = factor(agegroup), 
                                                     y = mx * 1000, 
                                                     group = factor(period), 
                                                     color = factor(period))) +
    geom_line() + geom_point() + theme_bw() + 
    scale_x_discrete() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
      legend.text = element_text(size = 12),
      axis.text.y  = element_text(size = 11),
      legend.title = element_text(size = 13)  ) +
    scale_color_discrete("Periodo") +
    guides(color = guide_legend(ncol = 1)) + 
    labs(x = "Grupo de edad", y = "mx",
         title = paste0("Age Specific Mortality Rate (per 1,000)", ", ", muni))
  return(plot)
}

mx_municipio(pred_pc, "San Juan")



###Previas para observar los efectos de cada modelo
plot(fit_hc)

plot(fit_sb2)

plot(fit_ht)

plot(fit_ig)

plot(fit_pc)

###Comparaciones con e0 observado vs e0 estimado
#Combina e0 observado con e0 estimado de cualquiera de tus modelos
#e0_resumen_pc , e0_resumen_hc, e0_resumen_sb2, e0_resumen_ht, e0_resumen_ig
construir_comparacion_e0 <- function(e0_modelado) {
  e0_resumen_directo %>%
    rename(e0_observada = e0_observada) %>%
    left_join(
      e0_modelado %>% rename(e0_estimada = e0),
      by = c("period", "region", "sex")
    )
}

comparacion_pc <- construir_comparacion_e0(e0_resumen_pc)
comparacion_hc  <- construir_comparacion_e0(e0_resumen_hc)
comparacion_sb2 <- construir_comparacion_e0(e0_resumen_sb2)
comparacion_ht  <- construir_comparacion_e0(e0_resumen_ht)
comparacion_ig  <- construir_comparacion_e0(e0_resumen_ig)


###Comparación e0 observado vs e0 estimado por municipio 
e0_model_plot <- function(dat, per, col, llh) {
  plot <- ggplot(dat %>%
                   filter(period == per), aes(x = fct_reorder(region, e0_observada), 
                                              y = e0_estimada)) +
    geom_point(color = col, size = 1.8) +
    geom_point(aes(y = e0_observada), shape = 1) +
    coord_flip() +
    facet_wrap(~ sex, labeller = as_labeller(c(`1` = "Hombres", `2` = "Mujeres"))) +
    theme_minimal() +
    labs(title = paste0("e0 por municipio (estimada vs. observada), ", per, ", ", llh), 
         y = "e0", x = "") +
    theme(axis.title.x = element_text(size = 6))
  return(plot)
}

e0_model_plot(comparacion_pc, "2020-2024", "firebrick", "PC prior")
e0_model_plot(comparacion_hc, "2020-2024", "firebrick", "Half-Cauchy")
e0_model_plot(comparacion_sb2, "2020-2024", "firebrick", "Scale Beta2")
e0_model_plot(comparacion_ht, "2020-2024", "firebrick", "Half-t")
e0_model_plot(comparacion_ig, "2020-2024", "firebrick", "Inverse Gamma")


#Comentario: Todas se comportan de la misma forma, exceptuando la PC prior,
#quien presenta diferencia pero no tan significativa ante la comparación de las previas.
#Mirar Culebra de PC y Half Cauchy (mínimo cambio)

###Comparación e0 observado vs e0 estimado por periodo
e0_municipio_plot <- function(dat, muni, llh) {
  plot <- ggplot(dat %>%
                   filter(region == muni), aes(x = period, y = e0_estimada)) +
    geom_point(color = "red", size = 1.8) +
    geom_point(aes(y = e0_observada), shape = 1) +
    facet_wrap(~ sex, labeller = as_labeller(c(`1` = "Hombres", `2` = "Mujeres"))) +
    theme_minimal() +
    labs(title = paste0("e0 estimada vs. observada para ", muni, ", 1980-2025, ", llh), 
         y = "e0", x = "") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  return(plot)
}

e0_municipio_plot(comparacion_pc, "San Juan", "PC prior")


###HOLD
e0_dir_vs_sae_plot <- function(dat, per, llh) {
  plot <- ggplot(dat %>% filter(period == per), 
                 aes(x = e0_observada, y = e0_estimada, color = factor(sex))) +
    geom_point(shape = 1) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
    scale_color_manual(values = c("1" = "steelblue", "2" = "orange"),
                       labels = c("Hombres", "Mujeres"), name = "Sexo") +
    theme_minimal() +
    labs(title = paste0("e0: directo vs. modelado, ", per, ", ", llh),
         x = "e0 directo", y = "e0 modelado")
  return(plot)
}

e0_dir_vs_sae_plot(comparacion_pc, "2020-2024", "PC prior")


###HOLD
###e0 estimada por periodo
e0_model_plot <- function(dat, per, col, llh) {
  plot <- ggplot(dat %>%
                   filter(period == per), 
                 aes(x = fct_reorder(region, e0_estimada), 
                     y = e0_estimada,
                     ymin = e0_lower,
                     ymax = e0_upper)) +
    geom_pointrange(color = col, size = 0.4, linewidth = 0.6) +
    coord_flip() +
    facet_wrap(~ sex, labeller = as_labeller(c(`1` = "Hombres", `2` = "Mujeres"))) +
    theme_minimal() +
    labs(title = paste0("e0 por municipio (mediana e IC 95%), ", per, ", ", llh), 
         y = "e0", x = "") +
    theme(axis.title.x = element_text(size = 6),
          axis.text.y = element_text(size = 7))
  
  return(plot)
}

e0_model_plot(comparacion_pc, "2020-2024", "purple", "PC prior")

###Comparando todas las previas para sus correspondientes e0 estimados
comparacion_todas <- bind_rows(
  comparacion_pc  %>% mutate(previa = "PC prior"),
  comparacion_hc  %>% mutate(previa = "Half-Cauchy"),
  comparacion_sb2 %>% mutate(previa = "Scale Beta2"),
  comparacion_ht  %>% mutate(previa = "Half-t"),
  comparacion_ig  %>% mutate(previa = "Inverse Gamma")
)

#Mujeres, 2020-2024
ggplot(comparacion_todas %>% filter(period == "2020-2024", sex == 2),
       aes(x = fct_reorder(region, e0_estimada), y = e0_estimada, 
           ymin = e0_lower, ymax = e0_upper, color = previa)) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title = "e0 municipal (mujeres, 2020-2024) por previa",
       y = "e0", x = "")

#Hombres, 2020-2024
ggplot(comparacion_todas %>% filter(period == "2020-2024", sex == 1),
       aes(x = fct_reorder(region, e0_estimada), y = e0_estimada, 
           ymin = e0_lower, ymax = e0_upper, color = previa)) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title = "e0 municipal (hombres, 2020-2024) por previa",
       y = "e0", x = "")



####Cálculo de e0 con intervalos de credibilidad

calcular_e0_inla <- function(modelo_inla, df, age_params, Age, nsamples = 1000) {
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
          e0_observada = tb$ex[1]
        )
      }
    }
  }
  e0_observada_df <- bind_rows(e0_obs_list)
  
  # Muestras posteriores del predictor
  samples <- inla.posterior.sample(nsamples, modelo_inla)
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
  e0_estimada_df <- e0_sim_df %>%
    group_by(region, period, sex) %>%
    summarise(
      e0_estimada = median(e0, na.rm = TRUE),
      e0_lower    = quantile(e0, 0.025, na.rm = TRUE),
      e0_upper    = quantile(e0, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  # e0 observada y estimada
  e0_final <- left_join(e0_observada_df, e0_estimada_df, by = c("region", "period", "sex")) %>%
    arrange(region, period, sex) %>%
    mutate(est_eval = case_when(
      between(e0_observada, e0_lower, e0_upper) ~ "Estimación adecuada",
      e0_observada < e0_lower ~ "> e0 observada",
      e0_observada > e0_upper ~ "< e0 observada"
    ))
  
  return(e0_final)
}
# 10,100 y 1000
e0_pc_IC  <- calcular_e0_inla(fit_pc,  df, age_params, Age, nsamples = 10)
e0_hc_IC  <- calcular_e0_inla(fit_hc,  df, age_params, Age, nsamples = 10)
e0_sb2_IC <- calcular_e0_inla(fit_sb2, df, age_params, Age, nsamples = 10)
e0_ht_IC  <- calcular_e0_inla(fit_ht,  df, age_params, Age, nsamples = 10)
e0_ig_IC  <- calcular_e0_inla(fit_ig,  df, age_params, Age, nsamples = 10)

#Hombres, 2020-2024
e0_forest_plot(e0_pc_IC, "2020-2024", 1, "purple", "PC prior")
e0_forest_plot(e0_hc_IC, "2020-2024", 1, "purple", "Half-Cauchy")
e0_forest_plot(e0_sb2_IC, "2020-2024", 1, "purple", "Scale-Beta2")
e0_forest_plot(e0_ht_IC, "2020-2024", 1, "purple", "Half-t")
e0_forest_plot(e0_ig_IC, "2020-2024", 1, "purple", "Inverse-Gamma")

#Mujeres, 2020-2024
e0_forest_plot(e0_pc_IC, "2020-2024", 2, "purple", "PC prior")
e0_forest_plot(e0_hc_IC, "2020-2024", 2, "purple", "Half-Cauchy")
e0_forest_plot(e0_sb2_IC, "2020-2024", 2, "purple", "Scale-Beta2")
e0_forest_plot(e0_ht_IC, "2020-2024", 2, "purple", "Half-t")
e0_forest_plot(e0_ig_IC, "2020-2024", 2, "purple", "Inverse-Gamma")


e0_model_plot(e0_pc_IC, "2020-2024", "purple", "PC prior")
e0_model_plot(e0_hc_IC, "2020-2024", "purple", "Half-Cauchy")
e0_model_plot(e0_sb2_IC, "2020-2024", "purple", "Scale-Beta2")
e0_model_plot(e0_ht_IC, "2020-2024", "purple", "Half-t")
e0_model_plot(e0_ig_IC, "2020-2024", "purple", "Inverse-Gamma")


#Comentario:El patrón que vemos es exactamente el comportamiento esperado del shrinkage bayesiano
#que hemos venido discutiendo. En municipios pequeños con pocos eventos, el modelo se aleja del observado
#y encoge hacia el promedio. Eso no es un error del modelo, sino que es la corrección que queremos para eliminar el ruido. 
#Donde hay más datos, menos shrinkage.

########################################################################################
########################################################################################
########################################################################################
#Para el futuro - efectos de cada modelo 
n_samp <- 10000 #Comenzando en 10,000 muestras

samples <- bind_rows(
  as.data.frame(inla.hyperpar.sample(n_samp, fit_pc))  %>% mutate(prior = "PC prior"),
  as.data.frame(inla.hyperpar.sample(n_samp, fit_hc))  %>% mutate(prior = "Half-Cauchy"),
  as.data.frame(inla.hyperpar.sample(n_samp, fit_ht))  %>% mutate(prior = "Half-t"),
  as.data.frame(inla.hyperpar.sample(n_samp, fit_sb2)) %>% mutate(prior = "Scale Beta2"),
  as.data.frame(inla.hyperpar.sample(n_samp, fit_ig))  %>% mutate(prior = "Inverse Gamma")
)

#Leyenda de la gráfica
samples$prior <- factor(samples$prior,
                        levels = c("PC prior", "Half-Cauchy",
                                   "Half-t", "Scale Beta2",
                                   "Inverse Gamma"))

samples_long <- samples %>%
  pivot_longer(-prior, names_to = "hyperpar", values_to = "value")

#Para ver cuáles son los hiperparámetros empleados
unique(samples_long$hyperpar)

p1 <- ggplot(samples_long, aes(x = value, fill = prior)) +
  geom_histogram(bins = 80, alpha = 0.6, position = "identity",
                 color = "white", linewidth = 0.1) +
  facet_wrap(~ hyperpar, scales = "free", ncol = 2) +
  scale_fill_manual(values = c(
    "PC prior"      = "#2C7BB6",
    "Half-Cauchy"   = "#D7191C",
    "Half-t"        = "#1A9641",
    "Scale Beta2"   = "#FF7F00",
    "Inverse Gamma" = "#984EA3"
  )) +
  labs(
    title = "Distribución posterior de hiperparámetros por prior",
    x     = "Valor muestral",
    y     = "Frecuencia",
    fill  = "Prior"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(size = 9),
    legend.position  = "bottom"
  )

p1

p2 <- ggplot(samples_long, aes(x = value, color = prior)) +
  geom_density(linewidth = 0.8) +
  facet_wrap(~ hyperpar, scales = "free", ncol = 2) +
  scale_color_manual(values = c(
    "PC prior"      = "#2C7BB6",
    "Half-Cauchy"   = "#D7191C",
    "Half-t"        = "#1A9641",
    "Scale Beta2"   = "#FF7F00",
    "Inverse Gamma" = "#984EA3"
  )) +
  labs(
    title = "Densidad posterior de hiperparámetros por prior",
    x     = "Valor muestral",
    y     = "Densidad",
    color = "Prior"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(size = 9),
    legend.position  = "bottom"
  )

p2

hiperpar_names <- unique(samples_long$hyperpar)

plots_list <- lapply(hiperpar_names, function(hp) {
  samples_long %>%
    filter(hyperpar == hp) %>%
    ggplot(aes(x = value, fill = prior)) +
    geom_histogram(bins = 80, alpha = 0.6, position = "identity",
                   color = "white", linewidth = 0.1) +
    scale_fill_manual(values = c(
      "PC prior"      = "#2C7BB6",
      "Half-Cauchy"   = "#D7191C",
      "Half-t"        = "#1A9641",
      "Scale Beta2"   = "#FF7F00",
      "Inverse Gamma" = "#984EA3"
    )) +
    labs(title = hp, x = "", y = "Frecuencia", fill = "Prior") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 9))
})

plots_list

