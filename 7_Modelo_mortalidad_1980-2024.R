################################################
### Modelo para mortalidad en áreas pequeñas ###
################################################
#7 - 2026_30_06

#Trabajo de Prueba 7-3-2026 

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
library(epitools)
library(PHEindicatormethods)

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
periods <- 1980:2024

# ages <- c(
#   "00-04", "05-09", "10-14", "15-19", "20-24",
#   "25-29", "30-34", "35-39", "40-44", "45-49",
#   "50-54", "55-59", "60-64", "65-69", "70-74",
#   "75-79", "80-84", "85+"
# )
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
period_breaks <- c(seq(1980, 2020, by = 5), 2024)
period_labels <- paste0(seq(1980, 2020, by = 5), "-", c(seq(1985, 2020, by = 5), 2024))

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
  file.path(data_dir, "municipios_population_1980_2024.csv"),
  col_types = cols(
    fips3 = col_character(),
    .default = col_guess()
  )
) %>%
  filter(
    year >= 1980, #2000
    year <= 2023,
    agegrp != 0,
    #sex !=0 #Esto es si solo queremos trabajar en grupo solo con dos sexos, no el total (0).
  ) %>%
  mutate(
    fips3 = str_pad(fips3, width = 3, side = "left", pad = "0"),
    period = period_quinquenal(year),
    agegroup = case_when(
      agegrp == 1  ~ "00-04",
      agegrp == 2  ~ "05-09",
      agegrp == 3  ~ "10-14",
      agegrp == 4  ~ "15-19",
      agegrp == 5  ~ "20-24",
      agegrp == 6  ~ "25-29",
      agegrp == 7  ~ "30-34",
      agegrp == 8  ~ "35-39",
      agegrp == 9  ~ "40-44",
      agegrp == 10 ~ "45-49",
      agegrp == 11 ~ "50-54",
      agegrp == 12 ~ "55-59",
      agegrp == 13 ~ "60-64",
      agegrp == 14 ~ "65-69",
      agegrp == 15 ~ "70-74",
      agegrp == 16 ~ "75-79",
      agegrp == 17 ~ "80-84",
      agegrp %in% c(18, 19) ~ "85+"
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
defunciones <- read_dta(
  file.path(data_dir, "defunciones_municipios_long_1979_2023.dta")
) %>% rename(sex = sexo) %>%
  filter(
    year >= min(periods),
    year <= max(periods),
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
        breaks = c(seq(0, 85, by = 5), Inf), #85?
        labels = ages,
        right = FALSE
      )
    )
  ) %>%
  count(fips3, period, agegroup, sex, name = "deaths")

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
  n_interval = c(rep(5, length(ages) - 1), NA),
  ax = c(
    2.0, 2.5, 2.5, 2.5, 2.5,
    2.5, 2.5, 2.5, 2.5, 2.5,
    2.5, 2.5, 2.5, 2.5, 2.5,
    2.5, 2.5, NA
  )
)

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
fit_pc$summary.fitted.values$mean 
fit_pc$summary.fitted.values$`0.025quant`
fit_pc$summary.fitted.values$`0.975quant` 

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
  filter(region == "Adjuntas", period == "1985-1990", sex == "f")
nMx <- pred_demo$mx
Age <- c(0, 1, (seq(5,80,by=5)))
AgeInt <- inferAgeIntAbr(vec = nMx)
PR.lifetable <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, Sex = "f",mod = FALSE)
PR.lifetable

# Tabla de vida para todos los municipios
pred_demo <- pred
pred_demo$sex <- ifelse(pred_demo$sex == 1, "m", "f")
municipios_demo <- sort(unique(pred_demo$region))
periodos_demo   <- sort(unique(pred_demo$period))
sexos      <- c("m", "f")
Age <- c(0, 1, seq(5, 80, by = 5)) 
#Age <- c(0, seq(5, 85, by = 5))
Age_nuevo <- c(0, seq(5, 85, by = 5)) #85
tablas <- list()
for (muni in municipios_demo) {
  for (per in periodos_demo) {
    for (sx in sexos) {
      pred_sub <- pred_demo %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      tablas[[muni]][[per]][[sx]] <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, Sex = sx, mod = FALSE)
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

head(tablas)
tb

# ------------------------------------------
# 8.1. DemoTools - Función graduate_pclm()
# ------------------------------------------

# HOLD
pred_demo <- pred
pred_demo$sex <- ifelse(pred_demo$sex == 1, "m", "f")
municipios_demo <- sort(unique(pred_demo$region))
periodos_demo   <- sort(unique(pred_demo$period))
sexos      <- c("m", "f")
Age <- c(0, 1, seq(5, 80, by = 5)) 
#Age <- c(0, seq(5, 85, by = 5))
Age_nuevo <- c(0, seq(5, 85, by = 5)) #85
tablas <- list()
for (muni in municipios_demo) {
  for (per in periodos_demo) {
    for (sx in sexos) {
      pred_sub <- pred_demo %>%
        filter(region == muni, period == per, sex == sx)
      mx_pclm  <- graduate_pclm(Value = nMx, Age = Age, AgeInt = AgeInt, OAG = TRUE) #Age = Age_nuevo
      Age_pclm <- as.integer(names(mx_pclm))
      tablas[[muni]][[per]][[sx]] <- lt_single_mx(nMx = mx_pclm, AgeInt = AgeInt, Age = Age_pclm, Sex = sx, mod = FALSE)
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

head(tablas)
tb






nMx    <- pred$mx
Age    <- c(0, 1, seq(5, 80, by = 5))
AgeInt <- inferAgeIntAbr(vec = nMx)

mx_pclm <- graduate_pclm(Value = nMx, Age = Age, AgeInt = AgeInt, OAnew = 115, OAG = TRUE)
mx_pclm

#Help
#OAG = TRUE	::logical, default = TRUE is the final age group open?
#OAnew = 115 ::integer, optional new open age, higher than max(Age). See details.

pred_pclm <- pred
pred_pclm$sex <- ifelse(pred_pclm$sex == 1, "m", "f")
municipios_pclm <- sort(unique(pred_pclm$region))
periodos_pclm   <- sort(unique(pred_pclm$period))
sexos      <- c("m", "f")
Age <- c(0, 1, seq(5, 80, by = 5))
tablas_pclm <- list()
for (muni in municipios_pclm) {
  for (per in periodos_pclm) {
    for (sx in sexos) {
      pred_sub_pclm <- pred_pclm %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub_pclm$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      mx_pclm  <- graduate_pclm(Value = nMx, Age = Age, AgeInt = AgeInt, OAG = TRUE) #Age = Age_nuevo
      Age_pclm <- as.integer(names(mx_pclm))
      tablas_pclm[[muni]][[per]][[sx]] <- lt_single_mx(nMx = as.numeric(mx_pclm), Age = Age_pclm, Sex = sx, mod = FALSE)
    }
  }
}
e0_resumen_pclm <- data.frame()
for (m in names(tablas_pclm)) {
  for (p in names(tablas_pclm[[m]])) {
    for (s in c("m", "f")) {
      
      tb_pclm <- tablas_pclm[[m]][[p]][[s]]
      
      sexnum <- if (s == "m") 1 else 2
      
      e0_resumen_pclm <- rbind(
        e0_resumen_pclm,
        data.frame(period = p, region = m, sex = sexnum, e0 = tb$ex[1])
      )
    }
  }
}
muj_pclm <- e0_resumen_pclm[e0_resumen_pclm$sex == 2, ] 
fila_muj_pclm <- muj_pclm[which.max(muj_pclm$e0), ]
hom_pclm <- e0_resumen_pclm[e0_resumen_pclm$sex == 1, ]
fila_hom_pclm <- hom_pclm[which.max(hom_pclm$e0), ]
fila_muj_pclm
fila_hom_pclm
head(tablas_pclm)

# --------------------------------------------
# 8.2. DemoTools - Función graduate_sprague()
# --------------------------------------------
# nMx    <- pred_sub$mx
# Age    <- c(0, 1, seq(5, 80, by = 5))
# AgeInt <- inferAgeIntAbr(vec = nMx)
# 
# mx_sprague <- graduate_sprague(Value = nMx, Age = Age, AgeInt = AgeInt, OAG = TRUE)
# mx_sprague

pred_sprague <- pred
pred_sprague$sex <- ifelse(pred_sprague$sex == 1, "m", "f")
municipios_sprague <- sort(unique(pred_sprague$region))
periodos_sprague   <- sort(unique(pred_sprague$period))
sexos      <- c("m", "f")
Age <- c(0, 1, seq(5, 80, by = 5))
tablas_sprague <- list()
for (muni in municipios_sprague) {
  for (per in periodos_sprague) {
    for (sx in sexos) {
      pred_sub_sprague <- pred_sprague %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub_sprague$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      mx_sprague  <- graduate_sprague(Value = nMx, Age = Age, AgeInt = AgeInt, OAG = TRUE)
      Age_sprague <- as.integer(names(mx_sprague))
      tablas_sprague[[muni]][[per]][[sx]] <- lt_single_mx(nMx = as.numeric(mx_sprague), Age = Age_sprague, Sex = sx, mod = FALSE)
    }
  }
}
e0_resumen_sprague <- data.frame()
for (m in names(tablas_sprague)) {
  for (p in names(tablas_sprague[[m]])) {
    for (s in c("m", "f")) {
      
      tb <- tablas_sprague[[m]][[p]][[s]]
      
      sexnum <- if (s == "m") 1 else 2
      
      e0_resumen_sprague <- rbind(
        e0_resumen_sprague,
        data.frame(period = p, region = m, sex = sexnum, e0 = tb$ex[1])
      )
    }
  }
}
muj_sprague <- e0_resumen_sprague[e0_resumen_sprague$sex == 2, ]
fila_muj_sprague<- muj_sprague[which.max(muj_sprague$e0), ]
hom_sprague <- e0_resumen_sprague[e0_resumen_sprague$sex == 1, ]
fila_hom_sprague <- hom_sprague[which.max(hom_sprague$e0), ]
fila_muj_sprague
fila_hom_sprague
head(tablas_sprague)

# ------------------------------------------
# 8.3. DemoTools - Función graduate_beers()
# ------------------------------------------
# nMx    <- pred_sub$mx
# Age    <- c(0, 1, seq(5, 80, by = 5))
# AgeInt <- inferAgeIntAbr(vec = nMx)
# 
# mx_beers <- graduate_beers(Value = nMx, Age = Age, AgeInt = AgeInt, OAG = TRUE, method = "ord")
# mx_beers

#Help
#method = "ord" :: character. Valid values are "ord" or "mod". Default "ord".
#HOLD - johnson = FALSE :: logical. Whether or not to adjust young ages according to the Demographic Analysis & Population Projection System Software method. 
#Default FALSE.

pred_beers <- pred
pred_beers$sex <- ifelse(pred_beers$sex == 1, "m", "f")
municipios_beers <- sort(unique(pred_beers$region))
periodos_beers   <- sort(unique(pred_beers$period))
sexos      <- c("m", "f")
Age <- c(0, 1, seq(5, 80, by = 5))
tablas_beers <- list()
for (muni in municipios_beers) {
  for (per in periodos_beers) {
    for (sx in sexos) {
      pred_sub_beers <- pred_beers %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub_beers$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      mx_beers  <- graduate_beers(Value = nMx, Age = Age, AgeInt = AgeInt, OAG = TRUE, method = "ord")
      Age_beers <- as.integer(names(mx_beers))
      tablas_beers[[muni]][[per]][[sx]] <- lt_single_mx(nMx = as.numeric(mx_beers), Age = Age_beers, Sex = sx, mod = FALSE)
    }
  }
}
e0_resumen_beers <- data.frame()
for (m in names(tablas_beers)) {
  for (p in names(tablas_beers[[m]])) {
    for (s in c("m", "f")) {
      
      tb <- tablas_beers[[m]][[p]][[s]]
      
      sexnum <- if (s == "m") 1 else 2
      
      e0_resumen_beers <- rbind(
        e0_resumen_beers,
        data.frame(period = p, region = m, sex = sexnum, e0 = tb$ex[1])
      )
    }
  }
}
muj_beers <- e0_resumen_beers[e0_resumen_beers$sex == 2, ]
fila_muj_beers <- muj_beers[which.max(muj_beers$e0), ]
hom_beers <- e0_resumen_beers[e0_resumen_beers$sex == 1, ]
fila_hom_beers<- hom_beers[which.max(hom_beers$e0), ]
fila_muj_beers
fila_hom_beers
head(tablas_beers)

# --------------------------------
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
# inla.doc("loggamma") #log-Gamma/loggamma/parameters = shape and rate
# inla.doc("gaussian")
# inla.doc("pc")

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
# formula_hc <- deaths ~
#   factor(sex) +
#   f(age_idx,    model = "rw1",  constr = TRUE,
#     hyper = list(prec = list(prior = HC.prior))) +
#   f(region_idx, model = "bym2", graph = g, constr = TRUE,
#     hyper = list(
#       prec = list(prior = HC.prior),
#       phi  = list(prior = "logitbeta", param = c(0.5, 0.5))))+  # Beta(0.5, 0.5) en escala logarítmica
#   f(period_idx, model = "rw2",  constr = TRUE,
#     hyper = list(prec = list(prior = HC.prior))) +
#   f(region_period_idx, model = "iid",
#     hyper = list(prec = list(prior = HC.prior)))

#Cambio de parámetros
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
               control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE)) #HOLD

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
Age <- c(0, 1, seq(5, 80, by = 5))
tablas <- list()
for (muni in municipios) {
  for (per in periodos) {
    for (sx in sexos) {
      pred_sub_hc <- pred_hc %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub_hc$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      tablas[[muni]][[per]][[sx]] <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, Sex = sx, mod = FALSE)
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
  a = 0.5;
  b = 0.5;
  sigma2 = exp(-theta);
  log_dens = lgamma(a+b) - lgamma(a) - lgamma(b);
  log_dens = log_dens + (a-1) * log(sigma2);
  log_dens = log_dens - (a+b) * log(1 + sigma2);
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

#Cambio de parámetros
formula_sb2 <- deaths ~
  factor(sex) +
  f(age_idx,    model = "rw1",  constr = TRUE,
    hyper = list(prec = list(prior = SB2.prior))) +
  f(region_idx, model = "bym2", graph = g, constr = TRUE,
    hyper = list(
      prec = list(prior = SB2.prior),
      phi  = list(prior = "logitbeta", param = c(0.5, 0.5))))+  # Beta(0.5, 0.5) en escala logarítmica
  f(period_idx, model = "rw2",  constr = TRUE,
    hyper = list(prec = list(prior = SB2.prior))) +
  f(region_period_idx, model = "iid",
    hyper = list(prec = list(prior = SB2.prior)))

fit_sb2 <- inla(formula_sb2,
                family  = "poisson",
                data    = df,
                E       = population,
                control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))


fit_sb2$all.hyper
fit_sb2$predictor$hyper$theta$prior

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
Age <- c(0, 1, seq(5, 80, by = 5))
tablas <- list()
for (muni in municipios) {
  for (per in periodos) {
    for (sx in sexos) {
      pred_sub_sb2 <- pred_sb2 %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub_sb2$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      tablas[[muni]][[per]][[sx]] <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, Sex = sx, mod = FALSE)
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

#Cambio de parámetros
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
               control.compute = list(config = TRUE,return.marginals.predictor=TRUE, cpo = TRUE, dic = TRUE, waic = TRUE)) 

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
Age <- c(0, 1, seq(5, 80, by = 5))
tablas <- list()
for (muni in municipios) {
  for (per in periodos) {
    for (sx in sexos) {
      pred_sub_ht <- pred_ht %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub_ht$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      tablas[[muni]][[per]][[sx]] <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, Sex = sx, mod = FALSE)
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

#Cambio de parámetros
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
               control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE))


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
Age <- c(0, 1, seq(5, 80, by = 5))
tablas <- list()
for (muni in municipios) {
  for (per in periodos) {
    for (sx in sexos) {
      pred_sub_ig <- pred_ig %>%
        filter(region == muni, period == per, sex == sx)
      nMx    <- pred_sub_ig$mx
      AgeInt <- inferAgeIntAbr(vec = nMx)
      tablas[[muni]][[per]][[sx]] <- lt_abridged(nMx = nMx, AgeInt = AgeInt, Age = Age, Sex = sx, mod = FALSE)
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

# -----------------------------
# 11. Análisis de sensitividad 
# -----------------------------
#Ejemplo: Libro de INLA - Sección 5.5

prior.list = list(
  default = list(prec = list(prior = "loggamma", param = c(1, 0.00005))),
  pc.prec = list(prec = list(prior = "pc.prec", param = c(5, 0.01))),
  half.cauchy = list(prec = list(prior = HC.prior)),
  half.t = list(prec = list(prior = HT.prior)),
  scale.beta2 = list(prec = list(prior = SB2.prior)),
  inverse.gamma = list(prec = list(prior = IG.prior))
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
    control.compute = list(dic = TRUE, waic = TRUE)
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

# --------------------------------
# 12. Sampling from the posterior
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

# -----------------------
# 12. Gráficos - Previas
# -----------------------
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




