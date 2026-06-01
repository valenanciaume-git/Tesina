# ==============================================================================
# 07_graficos_verificacion.R
# Generación de todos los gráficos y cuadros de la tesina
#
# Gráficos PNG y cuadros CSV guardados en output/cuadros_y_graficos/
#   _GD  → Serie GD (Ley Nueva York)
#   _AL  → Serie AL (Ley Argentina)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(vars)
library(forecast)

select     <- dplyr::select
filter     <- dplyr::filter
mutate     <- dplyr::mutate
summarise  <- dplyr::summarise
arrange    <- dplyr::arrange
inner_join <- dplyr::inner_join

source("scripts/00_funciones.R")

CARPETA <- "output/cuadros_y_graficos"
if (!dir.exists(CARPETA)) dir.create(CARPETA, recursive = TRUE)

guardar_grafico <- function(p, nombre, ancho = 12, alto = 7) {
  ruta <- file.path(CARPETA, paste0(nombre, ".png"))
  ggsave(ruta, p, width = ancho, height = alto, dpi = 150, bg = "white")
  cat(sprintf("✓ %s.png\n", nombre))
}

guardar_cuadro <- function(df, nombre) {
  ruta <- file.path(CARPETA, paste0(nombre, ".csv"))
  write.csv(df, ruta, row.names = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("✓ %s.csv\n", nombre))
}

# ── Constantes de estilo ──────────────────────────────────────────────────────
# Notación matemática con expression() para subíndices correctos en gráficos

LABELS_FACTORES <- c(
  beta0 = expression(beta[0]~"— Nivel"),
  beta1 = expression(beta[1]~"— Pendiente"),
  beta2 = expression(beta[2]~"— Curvatura"))

COLORES_FACTORES <- c(
  beta0 = "#2166ac",
  beta1 = "#d6604d",
  beta2 = "#4dac26")

LABELS_SPREAD <- c(
  spread_nivel     = expression(beta[0]~"(nivel)"),
  spread_pendiente = expression(beta[1]~"(pendiente)"),
  spread_curvatura = expression(beta[2]~"(curvatura)"))

COLORES_SPREAD <- c(
  spread_nivel     = "#2166ac",
  spread_pendiente = "#d6604d",
  spread_curvatura = "#4dac26")

# Labeller para facet_wrap con notación matemática
labeller_factores <- as_labeller(c(
  beta0 = "beta[0]~'— Nivel'",
  beta1 = "beta[1]~'— Pendiente'",
  beta2 = "beta[2]~'— Curvatura'"),
  default = label_parsed)

# ── Cargar datos ──────────────────────────────────────────────────────────────

panel_diario <- readRDS("data/panel_diario.rds")
ancla_df     <- readRDS("data/ancla_usd.rds")
cf_list      <- readRDS("data/cashflows.rds")
factores_GD  <- readRDS("output/factores_GD.rds")
factores_AL  <- readRDS("output/factores_AL.rds")
cal_GD_sp    <- readRDS("output/calibracion_GD_sp.rds")
cal_GD_A     <- readRDS("output/calibracion_GD_A.rds")
cal_AL_sp    <- readRDS("output/calibracion_AL_sp.rds")
cal_AL_A     <- readRDS("output/calibracion_AL_A.rds")
spread_GD_AL <- readRDS("output/spread_GD_AL.rds")
pred_res     <- readRDS("output/prediccion_resultados.rds")
comp_pond    <- readRDS("output/comparacion_ponderadores.rds")

fB_GD <- factores_GD %>%
  filter(fecha >= SUBPERIODOS$B$inicio,
         fecha <= SUBPERIODOS$B$fin,
         !is.na(beta0), beta0 > 0.01)

fB_AL <- factores_AL %>%
  filter(fecha >= SUBPERIODOS$B$inicio,
         fecha <= SUBPERIODOS$B$fin,
         !is.na(beta0), beta0 > 0.01)

cat("=== Generando gráficos y cuadros ===\n\n")

# ── Función: serie precio obs vs teórico ─────────────────────────────────────

serie_precios <- function(resultados, bonos) {
  dplyr::bind_rows(Filter(Negate(is.null), lapply(resultados, function(r) {
    if (is.null(r) || is.na(r$beta0)) return(NULL)
    dplyr::bind_rows(lapply(seq_along(bonos), function(i) {
      if (is.na(r$p_obs[i]) || is.na(r$p_teo[i])) return(NULL)
      data.frame(fecha     = r$fecha,
                 bono      = bonos[i],
                 observado = r$p_obs[i],
                 teorico   = r$p_teo[i],
                 error_pct = r$error_pct[i],
                 excupon   = isTRUE(r$excupon))
    }))
  })))
}

# ──  Loadings NS ────────────────────────────────────────────────────────

plazos <- seq(0.08, 12, by = 0.1)
ldg <- data.frame(
  plazo = plazos,
  beta0 = 1,
  beta1 = (1 - exp(-plazos/TAU_FIJO)) / (plazos/TAU_FIJO),
  beta2 = (1 - exp(-plazos/TAU_FIJO)) / (plazos/TAU_FIJO) -
    exp(-plazos/TAU_FIJO)) %>%
  pivot_longer(c(beta0, beta1, beta2),
               names_to  = "factor",
               values_to = "loading")

p_ldg <- ggplot(ldg, aes(x = plazo, y = loading, color = factor)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = TAU_FIJO, linetype = "dashed",
             color = "gray50", linewidth = 0.7) +
  annotate("text", x = TAU_FIJO + 0.1, y = 0.35,
           label = sprintf("tau == %.2f~años", TAU_FIJO),
           parse = TRUE, hjust = 0, size = 3.5, color = "gray40") +
  scale_color_manual(values = COLORES_FACTORES,
                     labels  = LABELS_FACTORES) +
  scale_x_continuous(breaks = 0:12) +
  labs(title   = "Funciones de carga del modelo de Nelson-Siegel",
       subtitle = sprintf("τ = %.2f años (λ = %.4f)", TAU_FIJO, LAMBDA_FIJO),
       x = "Plazo (años)", y = "Carga del factor", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

guardar_grafico(p_ldg, "G4.1_loadings_NS")

# ──  Precios GD y AL ─────────────────────────────────────────────

for (cfg in list(
  list(bonos = BONOS_GD, titulo = "GD (Ley Nueva York)",
       pal = "Set1", archivo = "G6.1_precios_GD"),
  list(bonos = BONOS_AL, titulo = "AL (Ley Argentina)",
       pal = "Set2", archivo = "G6.2_precios_AL")
)) {
  bonos_disp <- cfg$bonos[cfg$bonos %in% names(panel_diario)]
  df_long <- panel_diario[!panel_diario$outlier,
                          c("fecha", bonos_disp)] %>%
    pivot_longer(all_of(bonos_disp),
                 names_to = "bono", values_to = "precio") %>%
    filter(!is.na(precio))
  
  p <- ggplot(df_long, aes(x = fecha, y = precio, color = bono)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    geom_vline(xintercept = as.Date("2024-04-01"),
               linetype = "dashed", color = "gray40", linewidth = 0.8) +
    annotate("rect", xmin = as.Date("2021-01-01"),
             xmax = as.Date("2024-03-31"), ymin = -Inf, ymax = Inf,
             alpha = 0.04, fill = "#d73027") +
    annotate("rect", xmin = as.Date("2024-04-01"),
             xmax = max(panel_diario$fecha), ymin = -Inf, ymax = Inf,
             alpha = 0.04, fill = "#2166ac") +
    scale_color_brewer(palette = cfg$pal) +
    labs(title   = sprintf("Evolución de precios — Bonos %s", cfg$titulo),
         subtitle = "Enero 2021 — abril 2026",
         x = NULL, y = "Precio (USD por 100 VN)", color = "Bono",
         caption = "Fuente: Precios 1816") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  guardar_grafico(p, cfg$archivo)
}

# ──  Ancla ───────────────────────────────────────────────────────────────

p_ancla <- ggplot(ancla_df, aes(x = fecha, y = tasa * 100)) +
  geom_line(color = "#2166ac", linewidth = 0.9) +
  geom_point(size = 1.5, color = "#2166ac") +
  geom_vline(xintercept = as.Date("2024-04-01"),
             linetype = "dashed", color = "gray40") +
  labs(title   = "Tasa de plazos fijos en dólares a 30 días — Ancla nominal",
       subtitle = "Bancos privados, tasa anual nominal (%)",
       x = NULL, y = "Tasa anual (%)",
       caption = "Fuente: Bloomberg (BADLARU). Frecuencia mensual.") +
  theme_minimal(base_size = 12)

guardar_grafico(p_ancla, "G6.3_ancla_pf_usd", ancho = 10, alto = 5)

# ──  Paneles precio obs vs teórico ───────────────────────────────

for (cfg in list(
  list(res = cal_GD_A,  bonos = BONOS_GD, sub = "A",
       titulo = "Sub.A: distressed (ene 2021 — mar 2024) — Serie GD",
       archivo = "G6.4_panel_subA_GD"),
  list(res = cal_GD_sp, bonos = BONOS_GD, sub = "B",
       titulo = "Sub.B: normalización (abr 2024 — dic 2025) — Serie GD",
       archivo = "G6.7_panel_subB_GD"),
  list(res = cal_AL_A,  bonos = BONOS_AL, sub = "A",
       titulo = "Sub.A: distressed — Serie AL",
       archivo = "G6.4b_panel_subA_AL"),
  list(res = cal_AL_sp, bonos = BONOS_AL, sub = "B",
       titulo = "Sub.B: normalización — Serie AL",
       archivo = "G6.7b_panel_subB_AL")
)) {
  s <- serie_precios(cfg$res, cfg$bonos)
  if (is.null(s) || nrow(s) == 0) {
    cat(sprintf("— %s: sin datos\n", cfg$archivo)); next
  }
  
  p <- ggplot(s, aes(x = fecha)) +
    geom_line(aes(y = observado, color = "Observado"), linewidth = 0.7) +
    geom_line(aes(y = teorico,   color = "Teórico NS"),
              linewidth = 0.65, linetype = "dashed") +
    geom_point(data = s[s$excupon, ], aes(y = observado),
               color = "orange", size = 1.5, shape = 17, alpha = 0.7) +
    scale_color_manual(values = c(
      "Observado"  = "#2166ac",
      "Teórico NS" = "#d73027")) +
    facet_wrap(~bono, scales = "free_y", ncol = 2) +
    labs(title   = cfg$titulo, x = NULL,
         y = "Precio (USD por 100 VN)", color = NULL,
         caption = "Triángulos naranjas: días ex-cupón | τ = 1.13 años") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  
  guardar_grafico(p, cfg$archivo, ancho = 12, alto = 10)
}

# ──  Precios por duration oct-2022 ──────────────────────────────────────

fecha_oct22 <- as.Date("2022-10-31")
obs_oct22   <- as.numeric(panel_diario[
  panel_diario$fecha == fecha_oct22, BONOS_GD])
dur_oct22   <- c(2.52, 2.74, 5.67, 5.26, 6.52, 5.95)

df_oct22 <- data.frame(
  bono = BONOS_GD, precio = obs_oct22, duration = dur_oct22) %>%
  filter(!is.na(precio))

p_oct22 <- ggplot(df_oct22, aes(x = duration, y = precio, label = bono)) +
  geom_point(size = 4, color = "#d73027") +
  geom_text(vjust = -0.8, size = 3.5, color = "#d73027") +
  geom_smooth(method = "lm", se = FALSE,
              color = "gray50", linetype = "dashed", linewidth = 0.8) +
  labs(title   = "Precios por duration modificada — Octubre 2022",
       subtitle = "Estructura no convencional: bonos cortos más baratos que los largos",
       x = "Duration modificada (años)",
       y = "Precio (USD por 100 VN)",
       caption = "Fuente: Precios 1816") +
  theme_minimal(base_size = 12)

guardar_grafico(p_oct22, "G6.5_precios_duration_oct2022",
                ancho = 10, alto = 6)

# ──  RMSE temporal ──────────────────────────────────────────────────────

factores_todos_GD <- factores_GD %>% filter(!is.na(beta0), beta0 > 0.01)

p_rmse <- ggplot(factores_todos_GD, aes(x = fecha, y = rmse)) +
  geom_line(color = "gray60", linewidth = 0.5) +
  geom_smooth(method = "loess", span = 0.15,
              color = "#2166ac", se = FALSE, linewidth = 1) +
  geom_hline(yintercept = 1.5, linetype = "dashed",
             color = "#d73027", linewidth = 0.8) +
  geom_vline(data = data.frame(
    fecha = FECHAS_EXCUPON[FECHAS_EXCUPON >= min(factores_GD$fecha)]),
    aes(xintercept = fecha),
    color = "orange", linewidth = 0.5, alpha = 0.6, linetype = "dotted") +
  geom_vline(xintercept = as.Date("2024-04-01"),
             linetype = "dotted", color = "gray40") +
  annotate("rect", xmin = min(factores_todos_GD$fecha),
           xmax = as.Date("2024-03-31"), ymin = -Inf, ymax = Inf,
           alpha = 0.04, fill = "#d73027") +
  annotate("rect", xmin = as.Date("2024-04-01"),
           xmax = max(factores_todos_GD$fecha), ymin = -Inf, ymax = Inf,
           alpha = 0.04, fill = "#2166ac") +
  labs(title   = "Evolución del RMSE — Serie GD",
       subtitle = "Línea azul: tendencia LOESS | Líneas naranjas: días ex-cupón",
       x = NULL, y = "RMSE (USD por 100 VN)",
       caption = "τ = 1.13 años | Calibración diaria por precio") +
  theme_minimal(base_size = 12)

guardar_grafico(p_rmse, "G6.6_rmse_temporal_GD")

# ── Gráficos individuales por bono ─────────────────────────────────────

sB_GD <- serie_precios(cal_GD_sp, BONOS_GD)

if (!is.null(sB_GD) && nrow(sB_GD) > 0) {
  for (b in BONOS_GD) {
    df_b   <- sB_GD %>% filter(bono == b)
    mape_b <- mean(abs(df_b$error_pct), na.rm = TRUE)
    rmse_b <- sqrt(mean((df_b$observado - df_b$teorico)^2, na.rm = TRUE))
    
    p_b <- ggplot(df_b, aes(x = fecha)) +
      geom_line(aes(y = observado, color = "Observado"), linewidth = 0.9) +
      geom_line(aes(y = teorico,   color = "Teórico NS"),
                linewidth = 0.8, linetype = "dashed") +
      geom_point(data = df_b[df_b$excupon, ], aes(y = observado),
                 color = "orange", size = 2.5, shape = 17) +
      scale_color_manual(values = c(
        "Observado"  = "#2166ac",
        "Teórico NS" = "#d73027")) +
      annotate("text",
               x = min(df_b$fecha) + (max(df_b$fecha)-min(df_b$fecha))*0.02,
               y = max(df_b$observado, na.rm = TRUE) * 0.97,
               label = sprintf("MAPE=%.3f%%\nRMSE=%.4f", mape_b, rmse_b),
               hjust = 0, vjust = 1, size = 3.5, color = "gray30") +
      labs(title   = sprintf("Precio obs vs teórico NS — %s (Serie GD)", b),
           subtitle = "Subperíodo B: normalización (abr 2024 — dic 2025)",
           x = NULL, y = "Precio (USD por 100 VN)", color = NULL,
           caption = "Triángulos: días ex-cupón | τ = 1.13 años") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    
    guardar_grafico(p_b, sprintf("G6.8_%s_obs_vs_teo", b),
                    ancho = 10, alto = 5)
  }
}

# ──  Boxplot errores GD y AL ────────────────────────────────────────────

sB_AL <- serie_precios(cal_AL_sp, BONOS_AL)

for (cfg in list(
  list(s = sB_GD, titulo = "Serie GD",  archivo = "G6.9_boxplot_GD"),
  list(s = sB_AL, titulo = "Serie AL",  archivo = "G6.9b_boxplot_AL")
)) {
  if (is.null(cfg$s) || nrow(cfg$s) == 0) next
  
  p_box <- ggplot(cfg$s, aes(x = bono, y = error_pct, fill = bono)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1.5, outlier.alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    geom_hline(yintercept = c(-1, 1), linetype = "dotted",
               color = "gray50", linewidth = 0.5) +
    scale_fill_brewer(palette = "Set2") +
    labs(title   = sprintf("Distribución del error de valuación — %s",
                           cfg$titulo),
         subtitle = "Subperíodo B: normalización",
         x = NULL, y = "Error porcentual (%)",
         caption = "Cajas: IQR | Bigotes: Q1−1.5×IQR / Q3+1.5×IQR") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  
  guardar_grafico(p_box, cfg$archivo, ancho = 10, alto = 6)
}

# ── Evolución de parámetros ───────────────────────────────────────────

for (cfg in list(
  list(df = fB_GD, serie = "GD", archivo = "G6.10_params_GD"),
  list(df = fB_AL, serie = "AL", archivo = "G6.10b_params_AL")
)) {
  df_long <- cfg$df %>%
    select(fecha, beta0, beta1, beta2) %>%
    pivot_longer(c(beta0, beta1, beta2),
                 names_to  = "factor",
                 values_to = "valor")
  
  p <- ggplot(df_long, aes(x = fecha, y = valor * 100, color = factor)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    facet_wrap(~factor, scales = "free_y", ncol = 1,
               labeller = labeller_factores) +
    scale_color_manual(values  = COLORES_FACTORES,
                       labels  = LABELS_FACTORES,
                       guide   = "none") +
    labs(title   = sprintf("Evolución de parámetros NS — Serie %s",
                           cfg$serie),
         subtitle = "Subperíodo B | τ = 1.13 años fijo",
         x = NULL, y = "Valor (%)") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(size = 10, face = "bold"))
  
  guardar_grafico(p, cfg$archivo, ancho = 12, alto = 9)
}

# ── Curvas NS en fechas clave ─────────────────────────────────────────

plazos_curva <- seq(0.08, 12, by = 0.1)
fechas_clave <- list(
  list(ym = "2022-10", label = "Oct-2022 (crisis)"),
  list(ym = "2024-04", label = "Abr-2024 (inicio norm.)"),
  list(ym = "2024-09", label = "Sep-2024"),
  list(ym = "2025-01", label = "Ene-2025"),
  list(ym = "2025-12", label = "Dic-2025"))

colores_curvas <- c(
  "Oct-2022 (crisis)"      = "#d73027",
  "Abr-2024 (inicio norm.)"= "#abd9e9",
  "Sep-2024"               = "#74add1",
  "Ene-2025"               = "#4575b4",
  "Dic-2025"               = "#08306b")

construir_curvas <- function(resultados) {
  dplyr::bind_rows(lapply(fechas_clave, function(fc) {
    res <- Filter(function(r)
      !is.null(r) && !is.na(r$beta0) &&
        format(r$fecha, "%Y-%m") == fc$ym,
      resultados)
    if (length(res) == 0) return(NULL)
    r <- res[[1]]
    data.frame(
      plazo  = plazos_curva,
      spot   = ns_spot(r$beta0, r$beta1, r$beta2,
                       TAU_FIJO, plazos_curva) * 100,
      label  = fc$label,
      valida = r$beta0 > 0.01 & r$rmse < 5)
  }))
}

for (cfg in list(
  list(res  = c(cal_GD_A, cal_GD_sp),
       titulo = "Curvas spot NS — Serie GD",
       archivo = "G6.11_curvas_GD"),
  list(res  = c(cal_AL_A, cal_AL_sp),
       titulo = "Curvas spot NS — Serie AL",
       archivo = "G6.11b_curvas_AL")
)) {
  curvas <- construir_curvas(cfg$res)
  if (is.null(curvas) || nrow(curvas) == 0) next
  
  p <- ggplot(curvas,
              aes(x = plazo, y = spot, color = label,
                  linetype = ifelse(valida, "Válida", "Inválida"))) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = colores_curvas) +
    scale_linetype_manual(
      values = c("Válida" = "solid", "Inválida" = "dashed")) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    scale_x_continuous(breaks = 0:12) +
    labs(title    = cfg$titulo,
         subtitle = "Línea punteada: calibración inválida (sub.A)",
         x = "Plazo (años)", y = "Tasa spot (%)",
         color = "Fecha", linetype = "Calibración",
         caption = "τ = 1.13 años | Ancla: PF USD 30d") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
  
  guardar_grafico(p, cfg$archivo)
}

# ── Spread GD-AL ──────────────────────────────────────────────────────

if (!is.null(spread_GD_AL) && nrow(spread_GD_AL) > 0) {
  sp_long <- spread_GD_AL %>%
    select(fecha, spread_nivel, spread_pendiente, spread_curvatura) %>%
    pivot_longer(-fecha, names_to = "factor", values_to = "spread")
  
  p_sp <- ggplot(sp_long, aes(x = fecha, y = spread * 100,
                              color = factor)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~factor, scales = "free_y", ncol = 1,
               labeller = as_labeller(c(
                 spread_nivel     = "beta[0]~'(nivel)'",
                 spread_pendiente = "beta[1]~'(pendiente)'",
                 spread_curvatura = "beta[2]~'(curvatura)'"),
                 default = label_parsed)) +
    scale_color_manual(values = COLORES_SPREAD, guide = "none") +
    labs(title   = "Spread GD-AL por jurisdicción legal",
         subtitle = "Diferencia en parámetros NS: curva GD menos curva AL",
         x = NULL, y = "Spread (puntos porcentuales)",
         caption = "Subperíodo B | τ = 1.13 años") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
  
  guardar_grafico(p_sp, "G6.12_spread_GD_AL", ancho = 12, alto = 9)
  
  # Spread en tasas spot por plazos seleccionados
  plazos_spread <- c(1, 2, 5, 10)
  
  sp_spot <- dplyr::bind_rows(lapply(seq_len(nrow(spread_GD_AL)), function(i) {
    r <- spread_GD_AL[i, ]
    vals <- sapply(plazos_spread, function(p) {
      s_gd <- ns_spot(r$beta0_GD, r$beta1_GD, r$beta2_GD, TAU_FIJO, p)
      s_al <- ns_spot(r$beta0_AL, r$beta1_AL, r$beta2_AL, TAU_FIJO, p)
      (s_gd - s_al) * 100
    })
    as.data.frame(t(vals)) %>%
      setNames(paste0("p", plazos_spread, "y")) %>%
      mutate(fecha = r$fecha)
  })) %>%
    pivot_longer(-fecha, names_to = "plazo", values_to = "spread") %>%
    mutate(plazo = recode(plazo,
                          p1y = "1 año", p2y = "2 años",
                          p5y = "5 años", p10y = "10 años"))
  
  p_sp_spot <- ggplot(sp_spot,
                      aes(x = fecha, y = spread, color = plazo)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_brewer(palette = "Set1") +
    labs(title   = "Spread de tasas spot GD-AL por plazo",
         subtitle = "Prima por jurisdicción legal (pp)",
         x = NULL, y = "Spread (pp)", color = "Plazo",
         caption = "Positivo: GD rinde más que AL al mismo plazo") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  guardar_grafico(p_sp_spot, "G6.13_spread_spot_GD_AL")
}

# ──  ACF y PACF ────────────────────────────────────────────────────────

for (cfg in list(
  list(df = fB_GD, serie = "GD", archivo = "G6.14_ACF_PACF_GD"),
  list(df = fB_AL, serie = "AL", archivo = "G6.14b_ACF_PACF_AL")
)) {
  if (nrow(cfg$df) == 0) next
  ruta <- file.path(CARPETA, paste0(cfg$archivo, ".png"))
  while (!is.null(grDevices::dev.list())) grDevices::dev.off()
  grDevices::png(ruta, width = 12, height = 8,
                 units = "in", res = 150, bg = "white")
  graphics::par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
  for (nm in c("beta0","beta1","beta2")) {
    lbl <- switch(nm,
                  beta0 = expression(beta[0]),
                  beta1 = expression(beta[1]),
                  beta2 = expression(beta[2]))
    serie <- cfg$df[[nm]]
    stats::acf(serie,
               main = bquote("ACF —" ~ .(lbl) ~ .(paste("(Serie", cfg$serie, ")"))),
               lag.max = 30)
    stats::pacf(serie,
                main = bquote("PACF —" ~ .(lbl) ~ .(paste("(Serie", cfg$serie, ")"))),
                lag.max = 30)
  }
  grDevices::dev.off()
  cat(sprintf("✓ %s.png\n", cfg$archivo))
}

# ──  IRF ───────────────────────────────────────────────────────────────

if (file.exists("output/diagnostico_dinamica.rds")) {
  diag_din <- readRDS("output/diagnostico_dinamica.rds")
  if (diag_din$n_granger_signif > 0 &&
      file.exists("output/var_modelo.rds")) {
    var_mod <- readRDS("output/var_modelo.rds")
    irf_res <- vars::irf(var_mod, n.ahead = 21,
                         boot = TRUE, runs = 200, ci = 0.90)
    grDevices::png(file.path(CARPETA, "G6.15_IRF.png"),
                   width = 12, height = 8,
                   units = "in", res = 150, bg = "white")
    plot(irf_res)
    grDevices::dev.off()
    cat("✓ G6.15_IRF.png\n")
  } else {
    cat("— G6.15 omitido: Granger no significativo\n")
  }
}

# ── RMSE predicción ───────────────────────────────────────────────────

rmse_df <- dplyr::bind_rows(lapply(pred_res$etapa1, function(r) {
  data.frame(h = r$h, Modelo = c("RW","AR","VAR"),
             RMSE = c(r$RW, r$AR, r$VAR))
}))

p_rmse_pred <- ggplot(rmse_df,
                      aes(x = h, y = RMSE, color = Modelo,
                          shape = Modelo)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    RW = "gray50", AR = "#2166ac", VAR = "#d73027")) +
  scale_x_continuous(breaks = c(1, 5, 10, 21)) +
  labs(title   = "RMSE de predicción de factores por horizonte — Serie GD",
       subtitle = "Validación interna — Subperíodo B (últimos 20%)",
       x = "Horizonte h (días hábiles)", y = "RMSE de factores",
       caption = "RW: caminata aleatoria | AR: autorregresivo | VAR: vectorial") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

guardar_grafico(p_rmse_pred, "G6.16_RMSE_prediccion", ancho = 8, alto = 5)

# ── Fan chart ─────────────────────────────────────────────────────────

if (file.exists("output/ar_modelos.rds")) {
  ar_modelos    <- readRDS("output/ar_modelos.rds")
  ultimo_factor <- fB_GD %>% tail(1)
  plazos_fan    <- seq(0.08, 12, by = 0.25)
  
  fc_b0 <- forecast::forecast(ar_modelos$beta0$modelo, h = 21)
  fc_b1 <- forecast::forecast(ar_modelos$beta1$modelo, h = 21)
  fc_b2 <- forecast::forecast(ar_modelos$beta2$modelo, h = 21)
  
  curvas_fan <- dplyr::bind_rows(
    data.frame(
      plazo = plazos_fan,
      spot  = ns_spot(ultimo_factor$beta0, ultimo_factor$beta1,
                      ultimo_factor$beta2, TAU_FIJO, plazos_fan) * 100,
      tipo  = sprintf("Estimada (%s)", ultimo_factor$fecha)),
    data.frame(
      plazo = plazos_fan,
      spot  = ns_spot(as.numeric(fc_b0$mean)[21],
                      as.numeric(fc_b1$mean)[21],
                      as.numeric(fc_b2$mean)[21],
                      TAU_FIJO, plazos_fan) * 100,
      tipo  = "Predicción AR (h=21 días)"))
  
  etiqueta_actual <- sprintf("Estimada (%s)", ultimo_factor$fecha)
  
  colores_fan   <- c("#d73027", "#2166ac")
  names(colores_fan)   <- c("Predicción AR (h=21 días)", etiqueta_actual)
  linetypes_fan <- c("dashed", "solid")
  names(linetypes_fan) <- c("Predicción AR (h=21 días)", etiqueta_actual)
  
  p_fan <- ggplot(curvas_fan,
                  aes(x = plazo, y = spot, color = tipo, linetype = tipo)) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = colores_fan) +
    scale_linetype_manual(values = linetypes_fan) +
    scale_x_continuous(breaks = 0:12) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(title   = "Curva NS estimada y predicción a 21 días — Serie GD",
         subtitle = "Predicción mediante modelos AR univariados",
         x = "Plazo (años)", y = "Tasa spot (%)",
         color = NULL, linetype = NULL,
         caption = "tau = 1.13 años | AR por AIC") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  guardar_grafico(p_fan, "G6.17_fan_chart_prediccion")
}

# ══════════════════════════════════════════════════════════════════════════════
# CUADROS CSV
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Generando cuadros CSV ---\n\n")

#  Estadísticas descriptivas de precios
c61 <- dplyr::bind_rows(lapply(c(BONOS_GD, BONOS_AL), function(b) {
  if (!b %in% names(panel_diario)) return(NULL)
  p_A <- panel_diario[panel_diario$fecha < as.Date("2024-04-01") &
                        !panel_diario$outlier, ]
  p_B <- panel_diario[panel_diario$fecha >= as.Date("2024-04-01") &
                        !panel_diario$outlier, ]
  x_A <- na.omit(as.numeric(p_A[[b]]))
  x_B <- na.omit(as.numeric(p_B[[b]]))
  if (length(x_A) == 0 || length(x_B) == 0) return(NULL)
  data.frame(
    Bono      = b,
    Serie     = ifelse(b %in% BONOS_GD, "GD", "AL"),
    Sub       = c("A (2021-mar2024)", "B (abr2024-dic2025)"),
    Media     = round(c(mean(x_A),   mean(x_B)),   2),
    Mediana   = round(c(median(x_A), median(x_B)), 2),
    Minimo    = round(c(min(x_A),    min(x_B)),    2),
    Maximo    = round(c(max(x_A),    max(x_B)),    2),
    SD        = round(c(sd(x_A),     sd(x_B)),     2))
}))
guardar_cuadro(c61, "C6.1_estadisticas_precios")

# C6.2: RMSE y MAPE por año subperíodo A — GD
factores_A_GD <- dplyr::bind_rows(Filter(Negate(is.null),
                                         lapply(cal_GD_A, function(r) {
                                           if (is.null(r)) return(NULL)
                                           data.frame(fecha=r$fecha, rmse=r$rmse, mape=r$mape,
                                                      beta2=r$beta2)
                                         })))

if (nrow(factores_A_GD) > 0) {
  c62 <- factores_A_GD %>%
    mutate(anio = format(fecha, "%Y")) %>%
    group_by(anio) %>%
    summarise(N_obs        = n(),
              RMSE_promedio = round(mean(rmse, na.rm=TRUE), 4),
              MAPE_promedio = round(mean(mape, na.rm=TRUE), 3),
              RMSE_max      = round(max(rmse,  na.rm=TRUE), 4),
              Beta2_bound   = round(mean(beta2 >= 0.999,
                                         na.rm=TRUE)*100, 1),
              .groups = "drop")
  guardar_cuadro(c62, "C6.2_RMSE_MAPE_subA_GD")
}

#  Precios oct-2022
c63 <- data.frame(
  Bono           = BONOS_GD,
  Precio_USD     = round(obs_oct22, 4),
  Duration_aprox = c(2.52, 2.74, 5.67, 5.26, 6.52, 5.95))
guardar_cuadro(c63, "C6.3_precios_oct2022_GD")

#  Comparación ponderadores
c64 <- data.frame(
  Serie  = c("GD","GD","AL","AL"),
  Modelo = rep(c("Sin ponderadores","Pond. 1/D*"), 2),
  MAPE   = round(c(mean(comp_pond$GD_sin_pond$mape),
                   mean(comp_pond$GD_pond_dur$mape),
                   mean(comp_pond$AL_sin_pond$mape),
                   mean(comp_pond$AL_pond_dur$mape)), 4),
  RMSE   = round(c(mean(comp_pond$GD_sin_pond$rmse),
                   mean(comp_pond$GD_pond_dur$rmse),
                   mean(comp_pond$AL_sin_pond$rmse),
                   mean(comp_pond$AL_pond_dur$rmse)), 4),
  Mejor  = c(
    ifelse(comp_pond$mejor_GD == "Sin ponderadores","✓",""),
    ifelse(comp_pond$mejor_GD == "Ponderado 1/D*","✓",""),
    ifelse(comp_pond$mejor_AL == "Sin ponderadores","✓",""),
    ifelse(comp_pond$mejor_AL == "Ponderado 1/D*","✓","")))
guardar_cuadro(c64, "C6.4_comparacion_ponderadores")

# MAPE y RMSE por bono subperíodo B
for (cfg in list(
  list(s = sB_GD, archivo = "C6.5_MAPE_RMSE_GD_subB"),
  list(s = sB_AL, archivo = "C6.5b_MAPE_RMSE_AL_subB")
)) {
  if (is.null(cfg$s) || nrow(cfg$s) == 0) next
  c65 <- cfg$s %>%
    group_by(bono) %>%
    summarise(
      MAPE_pct  = round(mean(abs(error_pct), na.rm=TRUE), 3),
      RMSE_USD  = round(sqrt(mean((observado-teorico)^2, na.rm=TRUE)), 4),
      MAE_USD   = round(mean(abs(observado-teorico), na.rm=TRUE), 4),
      Sesgo_pct = round(mean(error_pct, na.rm=TRUE), 3),
      .groups   = "drop")
  guardar_cuadro(c65, cfg$archivo)
}

# Estadísticas de factores β — GD y AL
c66 <- dplyr::bind_rows(lapply(list(
  list(df = fB_GD, serie = "GD"),
  list(df = fB_AL, serie = "AL")
), function(cfg) {
  cfg$df %>%
    summarise(across(c(beta0, beta1, beta2), list(
      Media  = ~round(mean(., na.rm=TRUE), 4),
      SD     = ~round(sd(.,  na.rm=TRUE), 4),
      Minimo = ~round(min(., na.rm=TRUE), 4),
      Maximo = ~round(max(., na.rm=TRUE), 4)))) %>%
    mutate(Serie = cfg$serie)
}))
guardar_cuadro(c66, "C6.6_estadisticas_factores")

#  Spread GD-AL promedio mensual
if (!is.null(spread_GD_AL) && nrow(spread_GD_AL) > 0) {
  c67 <- spread_GD_AL %>%
    mutate(mes = floor_date(fecha, "month")) %>%
    group_by(mes) %>%
    summarise(
      Spread_nivel_pp     = round(mean(spread_nivel,     na.rm=TRUE)*100, 3),
      Spread_pendiente_pp = round(mean(spread_pendiente, na.rm=TRUE)*100, 3),
      Spread_curvatura_pp = round(mean(spread_curvatura, na.rm=TRUE)*100, 3),
      .groups = "drop")
  guardar_cuadro(c67, "C6.7_spread_GD_AL_mensual")
}

# Tests estadísticos
if (file.exists("output/diagnostico_dinamica.rds")) {
  diag_din <- readRDS("output/diagnostico_dinamica.rds")
  
  c68 <- dplyr::bind_rows(lapply(names(diag_din$adf), function(nm) {
    r <- diag_din$adf[[nm]]
    data.frame(Factor       = nm,
               Estadistico  = round(r$estadistico, 4),
               p_valor      = round(r$p_valor_niveles, 4),
               Orden        = r$orden,
               Accion       = ifelse(r$necesita_diff,
                                     "Diferenciar","Niveles"))
  }))
  guardar_cuadro(c68, "C6.8_tests_ADF")
  
  c69 <- dplyr::bind_rows(lapply(names(diag_din$granger), function(nm) {
    r <- diag_din$granger[[nm]]
    data.frame(Par           = nm,
               F_stat        = round(r$f_stat, 4),
               p_valor       = round(r$p_valor, 4),
               Significativo = ifelse(r$significativo,
                                      "Sí (p<0.05)","No"))
  }))
  guardar_cuadro(c69, "C6.9_test_Granger")
  
  if (file.exists("output/ar_modelos.rds")) {
    ar_modelos <- readRDS("output/ar_modelos.rds")
    c610 <- dplyr::bind_rows(lapply(names(ar_modelos), function(nm) {
      mod <- ar_modelos[[nm]]$modelo
      data.frame(
        Factor        = nm,
        Orden_AR      = mod$arma[1],
        AIC           = round(mod$aic, 2),
        Coeficientes  = paste(round(coef(mod), 4), collapse=" | "),
        Necesita_diff = ar_modelos[[nm]]$necesita_diff)
    }))
    guardar_cuadro(c610, "C6.10_coeficientes_AR")
  }
}

# RMSE predicción
c611 <- dplyr::bind_rows(lapply(pred_res$etapa1, function(r) {
  data.frame(Horizonte_dias = r$h,
             RW  = round(r$RW,  4),
             AR  = round(r$AR,  4),
             VAR = round(r$VAR, 4),
             AR_supera_RW = ifelse(r$AR < r$RW, "Sí", "No"))
}))
guardar_cuadro(c611, "C6.11_RMSE_prediccion")

# ══════════════════════════════════════════════════════════════════════════════
# GRÁFICOS COMPARATIVOS GD vs AL
# ══════════════════════════════════════════════════════════════════════════════

cat("\n--- Generando gráficos comparativos GD vs AL ---\n\n")

factores_AL  <- readRDS("output/factores_AL.rds")

fB_AL <- factores_AL %>%
  filter(fecha >= SUBPERIODOS$B$inicio,
         fecha <= SUBPERIODOS$B$fin,
         !is.na(beta0), beta0 > 0.01)
# ── Curvas spot GD vs AL en fechas representativas ────────────────────

plazos_comp <- seq(0.08, 12, by = 0.1)
fechas_comp <- list(
  list(ym = "2024-06", label = "Jun-2024"),
  list(ym = "2024-09", label = "Sep-2024"),
  list(ym = "2025-01", label = "Ene-2025"),
  list(ym = "2025-06", label = "Jun-2025"),
  list(ym = "2025-12", label = "Dic-2025"))

curvas_comp <- dplyr::bind_rows(lapply(fechas_comp, function(fc) {
  r_gd <- Filter(function(r)
    !is.null(r) && !is.na(r$beta0) &&
      format(r$fecha, "%Y-%m") == fc$ym,
    cal_GD_sp)
  r_al <- Filter(function(r)
    !is.null(r) && !is.na(r$beta0) &&
      format(r$fecha, "%Y-%m") == fc$ym,
    cal_AL_sp)
  if (length(r_gd)==0 || length(r_al)==0) return(NULL)
  r_gd <- r_gd[[1]]; r_al <- r_al[[1]]
  
  dplyr::bind_rows(
    data.frame(
      plazo  = plazos_comp,
      spot   = ns_spot(r_gd$beta0, r_gd$beta1, r_gd$beta2,
                       TAU_FIJO, plazos_comp)*100,
      serie  = "GD (Ley NY)",
      fecha_label = fc$label),
    data.frame(
      plazo  = plazos_comp,
      spot   = ns_spot(r_al$beta0, r_al$beta1, r_al$beta2,
                       TAU_FIJO, plazos_comp)*100,
      serie  = "AL (Ley Arg.)",
      fecha_label = fc$label))
}))

p_curvas_comp <- ggplot(curvas_comp,
                        aes(x=plazo, y=spot, color=serie, linetype=serie)) +
  geom_line(linewidth=0.9) +
  facet_wrap(~fecha_label, ncol=3, scales="free_y") +
  scale_color_manual(values=c(
    "GD (Ley NY)"  = "#2166ac",
    "AL (Ley Arg.)"= "#d73027")) +
  scale_linetype_manual(values=c(
    "GD (Ley NY)"  = "solid",
    "AL (Ley Arg.)"= "dashed")) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) +
  scale_y_continuous(labels=function(x) paste0(x,"%")) +
  labs(title   = "Curvas spot NS — Comparación GD vs AL",
       subtitle = "Fechas representativas del subperíodo B",
       x="Plazo (años)", y="Tasa spot (%)",
       color=NULL, linetype=NULL,
       caption="τ = 1.13 años | AL incluye AE38") +
  theme_minimal(base_size=11) +
  theme(legend.position="bottom",
        strip.text=element_text(face="bold"))

guardar_grafico(p_curvas_comp, "G6.20_curvas_GD_vs_AL", ancho=14, alto=9)

# ── Curva spot promedio GD vs AL — subperíodo B completo ──────────────

spots_prom <- dplyr::bind_rows(lapply(plazos_comp, function(p) {
  sp_gd <- mean(sapply(seq_len(nrow(spread_df)), function(i) {
    r <- spread_df[i,]
    ns_spot(r$beta0_GD, r$beta1_GD, r$beta2_GD, TAU_FIJO, p)*100
  }))
  sp_al <- mean(sapply(seq_len(nrow(spread_df)), function(i) {
    r <- spread_df[i,]
    ns_spot(r$beta0_AL, r$beta1_AL, r$beta2_AL, TAU_FIJO, p)*100
  }))
  data.frame(
    plazo = p,
    spot  = c(sp_gd, sp_al),
    serie = c("GD (Ley NY)","AL (Ley Arg.)"))
}))

# Agregar línea de cruce
cruce_plazo <- 9.6

p_prom <- ggplot(spots_prom,
                 aes(x=plazo, y=spot, color=serie, linetype=serie)) +
  geom_line(linewidth=1.1) +
  geom_vline(xintercept=cruce_plazo, linetype="dotted",
             color="gray40", linewidth=0.8) +
  annotate("text", x=cruce_plazo+0.15, y=min(spots_prom$spot)*1.02,
           label=sprintf("Cruce\n%.1f años", cruce_plazo),
           hjust=0, size=3.2, color="gray40") +
  annotate("rect", xmin=0, xmax=cruce_plazo,
           ymin=-Inf, ymax=Inf, alpha=0.03, fill="#d73027") +
  annotate("rect", xmin=cruce_plazo, xmax=12,
           ymin=-Inf, ymax=Inf, alpha=0.03, fill="#2166ac") +
  annotate("text", x=2, y=max(spots_prom$spot)*0.97,
           label="AL rinde más", size=3.2, color="#d73027") +
  annotate("text", x=11, y=max(spots_prom$spot)*0.97,
           label="GD rinde más", size=3.2, color="#2166ac") +
  scale_color_manual(values=c(
    "GD (Ley NY)"  = "#2166ac",
    "AL (Ley Arg.)"= "#d73027")) +
  scale_linetype_manual(values=c(
    "GD (Ley NY)"  = "solid",
    "AL (Ley Arg.)"= "dashed")) +
  scale_x_continuous(breaks=0:12) +
  scale_y_continuous(labels=function(x) paste0(x,"%")) +
  labs(title   = "Curva spot promedio GD vs AL — Subperíodo B",
       subtitle = "Zona roja: AL rinde más | Zona azul: GD rinde más",
       x="Plazo (años)", y="Tasa spot promedio (%)",
       color=NULL, linetype=NULL,
       caption="τ = 1.13 años | Promedio sobre 421 días hábiles") +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom")

guardar_grafico(p_prom, "G6.21_curva_promedio_GD_vs_AL", ancho=10, alto=6)

# ── Evolución comparada de cada parámetro ─────────────────────────────

params_comp <- dplyr::bind_rows(
  fB_GD %>% select(fecha, beta0, beta1, beta2) %>% mutate(serie="GD"),
  fB_AL %>% select(fecha, beta0, beta1, beta2) %>% mutate(serie="AL"))

for (nm in c("beta0","beta1","beta2")) {
  lbl <- switch(nm,
                beta0 = expression(beta[0]~"— Nivel"),
                beta1 = expression(beta[1]~"— Pendiente"),
                beta2 = expression(beta[2]~"— Curvatura"))
  
  archivo <- switch(nm,
                    beta0 = "G6.22a_beta0_GD_vs_AL",
                    beta1 = "G6.22b_beta1_GD_vs_AL",
                    beta2 = "G6.22c_beta2_GD_vs_AL")
  
  df_nm <- params_comp %>%
    select(fecha, serie, valor = !!sym(nm))
  
  # Medias por serie
  med_GD <- mean(df_nm$valor[df_nm$serie=="GD"], na.rm=TRUE)
  med_AL <- mean(df_nm$valor[df_nm$serie=="AL"], na.rm=TRUE)
  
  p_nm <- ggplot(df_nm, aes(x=fecha, y=valor*100, color=serie)) +
    geom_line(linewidth=0.8, alpha=0.85) +
    geom_hline(yintercept=med_GD*100, linetype="dashed",
               color="#2166ac", linewidth=0.6, alpha=0.7) +
    geom_hline(yintercept=med_AL*100, linetype="dashed",
               color="#d73027", linewidth=0.6, alpha=0.7) +
    annotate("text",
             x=max(df_nm$fecha),
             y=med_GD*100,
             label=sprintf("GD: %.2f%%", med_GD*100),
             hjust=1, vjust=-0.5, size=3.2, color="#2166ac") +
    annotate("text",
             x=max(df_nm$fecha),
             y=med_AL*100,
             label=sprintf("AL: %.2f%%", med_AL*100),
             hjust=1, vjust=1.2, size=3.2, color="#d73027") +
    scale_color_manual(values=c(
      "GD"="#2166ac", "AL"="#d73027")) +
    labs(title   = bquote(.(lbl)~"— GD vs AL"),
         subtitle = "Subperíodo B: normalización (abr 2024 — dic 2025)",
         x=NULL, y="Valor (%)",
         color=NULL,
         caption="Líneas punteadas: medias del subperíodo") +
    theme_minimal(base_size=12) +
    theme(legend.position="bottom")
  
  guardar_grafico(p_nm, archivo, ancho=12, alto=5)
}

# ──  Panel triple — los tres parámetros juntos ────────────────────────

params_long <- params_comp %>%
  pivot_longer(c(beta0, beta1, beta2),
               names_to="factor", values_to="valor")

p_panel_comp <- ggplot(params_long,
                       aes(x=fecha, y=valor*100, color=serie)) +
  geom_line(linewidth=0.75, alpha=0.85) +
  facet_wrap(~factor, scales="free_y", ncol=1,
             labeller=labeller_factores) +
  scale_color_manual(values=c(
    "GD"="#2166ac", "AL"="#d73027"),
    labels=c("GD"="GD (Ley NY)", "AL"="AL (Ley Arg.)")) +
  labs(title   = "Evolución comparada de parámetros NS — GD vs AL",
       subtitle = "Subperíodo B | τ = 1.13 años",
       x=NULL, y="Valor (%)",
       color=NULL,
       caption="AL incluye AE38") +
  theme_minimal(base_size=11) +
  theme(legend.position="bottom",
        strip.text=element_text(size=10, face="bold"))

guardar_grafico(p_panel_comp, "G6.22d_params_GD_vs_AL_panel",
                ancho=12, alto=10)

# ── Resumen final ─────────────────────────────────────────────────────────────

cat("\n=== Archivos generados en", CARPETA, "===\n")
archivos <- list.files(CARPETA)
cat(sprintf("Gráficos PNG: %d\n", sum(grepl("\\.png$", archivos))))
cat(sprintf("Cuadros  CSV: %d\n", sum(grepl("\\.csv$", archivos))))
cat("\nLista completa:\n")
for (f in sort(archivos)) cat(sprintf("  %s\n", f))


cat("✓ G6.20_curvas_GD_vs_AL.png\n")
cat("✓ G6.21_curva_promedio_GD_vs_AL.png\n")
cat("✓ G6.22a_beta0_GD_vs_AL.png\n")
cat("✓ G6.22b_beta1_GD_vs_AL.png\n")
cat("✓ G6.22c_beta2_GD_vs_AL.png\n")
cat("✓ G6.22d_params_GD_vs_AL_panel.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# GRÁFICOS ANÁLISIS DINÁMICO — ADF, AR, GRANGER, FEVD
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(dplyr)
library(tidyr)
library(forecast)
library(vars)
library(grid)
library(gridExtra)

ar_modelos  <- readRDS("output/ar_modelos.rds")
var_mod     <- readRDS("output/var_modelo.rds")
diag_din    <- readRDS("output/diagnostico_dinamica.rds")
factores_GD <- readRDS("output/factores_GD.rds")

fB <- factores_GD %>%
  filter(fecha >= SUBPERIODOS$B$inicio,
         fecha <= SUBPERIODOS$B$fin,
         !is.na(beta0), beta0 > 0.01) %>%
  arrange(fecha)

# ── ADF — niveles vs diferencias  ───────────

cat("Generando G6.23 — ADF niveles vs diferencias...\n")

# Para que label_parsed funcione, los nombres deben ser expresiones plotmath válidas de R
adf_data <- dplyr::bind_rows(
  # Niveles (Columna Izquierda)
  data.frame(fecha  = fB$fecha,
             valor  = fB$beta0 * 100,
             factor = "beta[0]~~'— Nivel'",
             tipo   = "Serie en niveles",
             orden  = 1),
  data.frame(fecha  = fB$fecha,
             valor  = fB$beta1 * 100,
             factor = "beta[1]~~'— Pendiente'",
             tipo   = "Serie en niveles",
             orden  = 2),
  data.frame(fecha  = fB$fecha,
             valor  = fB$beta2 * 100,
             factor = "beta[2]~~'— Curvatura'",
             tipo   = "Serie en niveles",
             orden  = 3),
  # Diferencias (Columna Derecha)
  data.frame(fecha  = fB$fecha[-1],
             valor  = diff(fB$beta0) * 100,
             factor = "Delta*beta[0]~~'— Diferencia'",
             tipo   = "Primera diferencia",
             orden  = 4),
  data.frame(fecha  = fB$fecha,
             valor  = fB$beta1 * 100,   # β₁ ya es I(0)
             factor = "beta[1]~~'— Ya estacionaria'",
             tipo   = "Primera diferencia",
             orden  = 5),
  data.frame(fecha  = fB$fecha[-1],
             valor  = diff(fB$beta2) * 100,
             factor = "Delta*beta[2]~~'— Diferencia'",
             tipo   = "Primera diferencia",
             orden  = 6))

# Forzamos los niveles en el orden exacto para mantener la grilla de 3x2
adf_data$factor <- factor(adf_data$factor, levels = c(
  "beta[0]~~'— Nivel'", "beta[1]~~'— Pendiente'", "beta[2]~~'— Curvatura'",
  "Delta*beta[0]~~'— Diferencia'", "beta[1]~~'— Ya estacionaria'", "Delta*beta[2]~~'— Diferencia'"
))

adf_data$columna <- ifelse(adf_data$orden <= 3,
                           "Serie original", "Transformada para el modelo")

p_adf <- ggplot(adf_data, aes(x=fecha, y=valor)) +
  geom_line(aes(color=columna), linewidth=0.7) +
  geom_hline(yintercept=0, linetype="dashed",
             color="gray50", linewidth=0.4) +
  # CLAVE: agregamos labeller = label_parsed para procesar los subíndices y letras griegas
  facet_wrap(~factor, ncol=2, scales="free_y", labeller = label_parsed) +
  scale_color_manual(values=c(
    "Serie original"              = "#d73027",
    "Transformada para el modelo" = "#2166ac")) +
  scale_x_date(date_breaks="6 months", date_labels="%b\n%Y") +
  labs(title    = "Tests ADF — Series originales vs transformadas",
       subtitle = paste0(
         "β₀: ADF=−3.24, p=0.082 → I(1) | ",
         "β₁: ADF=−4.20, p=0.010 → I(0) | ",
         "β₂: ADF=−2.22, p=0.483 → I(1)"),
       x=NULL, y="Valor (%)",
       color=NULL,
       caption=paste0(
         "Columna izquierda: series sin transformar. ",
         "Columna derecha: series usadas en el modelo (en diferencias si I(1)).\n",
         "Las series diferenciadas oscilan en torno a cero, ",
         "confirmando estacionariedad.")) +
  theme_minimal(base_size=11) +
  theme(legend.position  = "bottom",
        strip.text        = element_text(face="bold", size=9),
        plot.caption      = element_text(size=8, color="gray40"),
        panel.grid.minor  = element_blank())

guardar_grafico(p_adf, "G6.23_ADF_niveles_vs_diferencias",
                ancho=14, alto=10)

# ──  AR — diagnóstico de residuos 
cat("Generando G6.24 — Diagnóstico de residuos AR...\n")

for (nm in c("beta0","beta1","beta2")) {
  mod  <- ar_modelos[[nm]]$modelo
  nd   <- ar_modelos[[nm]]$necesita_diff
  res  <- as.numeric(residuals(mod))
  n_obs <- length(res)
  p_orden <- mod$arma[1]
  
  archivo <- sprintf("G6.24_%s_residuos_AR", nm)
  
  # Abrir dispositivo gráfico
  png(file.path("output/cuadros_y_graficos", paste0(archivo,".png")),
      width=1400, height=700, res=110)
  
  # Configurar layout: Arriba gráfico temporal, abajo ACF e Histograma
  layout(matrix(c(1,1,2,3), nrow=2, byrow=TRUE))
  par(mar=c(4, 4.5, 3, 2), oma=c(0, 0, 3, 0))
  
  # 1. Gráfico de Residuos en el Tiempo
  plot(res, type="l", col="black", lwd=0.8,
       xlab="Índice de observaciones", ylab="Residuos",
       main="Residuos del modelo en el tiempo")
  abline(h=0, col="gray50", lty="dashed")
  
  # 2. Gráfico de Autocorrelación (ACF) — CORRECCIÓN DE ESCALA Y CENTRADO
  acf_res <- acf(res, lag.max=40, plot=FALSE)
  
  # Extraemos los valores omitiendo el lag 0
  lags_efectivos <- acf_res$lag[-1]
  acfs_efectivos <- acf_res$acf[-1]
  
  # Calculamos el límite dinámico para centrar el eje Y perfectamente
  max_val <- max(abs(acfs_efectivos))
  # Nos aseguramos un margen mínimo por si la serie es extremadamente limpia
  limite_y <- max(max_val * 1.3, 0.15) 
  
  plot(lags_efectivos, acfs_efectivos, type="h", lwd=2,
       xlab="Lag (Rezagos)", ylab="ACF", main="Función de Autocorrelación (ACF)",
       ylim=c(-limite_y, limite_y)) # Eje Y estirado de forma simétrica
  abline(h=0, col="black")
  
  # Líneas de significancia estadística al 95%
  limite_acf <- 1.96 / sqrt(n_obs)
  abline(h=c(-limite_acf, limite_acf), col="blue", lty="dashed")
  
  # 3. Histograma y Curva Normal
  h_res <- hist(res, breaks="FD", plot=FALSE)
  plot(h_res, col="gray70", border="white", 
       xlab="Valor del residuo", ylab="Densidad", main="Distribución de los residuos",
       freq=FALSE)
  x_axis <- seq(min(res), max(res), length=100)
  y_axis <- dnorm(x_axis, mean=mean(res), sd=sd(res))
  lines(x_axis, y_axis, col="#d73027", lwd=2)
  
  # 4. TÍTULO GENERAL — CORRECCIÓN DE SUBÍNDICES MATEMÁTICOS (SIN TEXTO PLANO)
  # Construimos la expresión matemática dinámicamente según el factor
  titulo_parsed <- if(nm == "beta0") {
    bquote(bold("Diagnóstico de Residuos — Proceso AR"*(.(p_orden))*" en diferencias ("*Delta*beta[0]*")"))
  } else if(nm == "beta1") {
    bquote(bold("Diagnóstico de Residuos — Proceso AR"*(.(p_orden))*" en niveles ("*beta[1]*")"))
  } else {
    bquote(bold("Diagnóstico de Residuos — Proceso AR"*(.(p_orden))*" en diferencias ("*Delta*beta[2]*")"))
  }
  
  # Renderizar el título exterior
  title(main=titulo_parsed, outer=TRUE, cex.main=1.3)
  
  dev.off()
  
  cat(sprintf("   ✓ %s.png\n", archivo))
}
# ──  VAR — Valores Observados vs. Predichos (Ajuste R²) ───────────────

cat("Generando G6.24b — Ajuste dentro de muestra del VAR...\n")

# 1. Extraer valores observados (ajustados a la muestra del VAR)
# var_mod$y contiene las series utilizadas en la estimación
obs_mat <- var_mod$y
n_var <- nrow(obs_mat)

# 2. Calcular los valores predichos (Fitted = Observado - Residuo)
res_mat <- residuals(var_mod)
fit_mat <- obs_mat - res_mat

# 3. Armar un dataframe largo para ggplot
var_ajuste_df <- data.frame()

# Mapeo de nombres a expresiones matemáticas para los paneles
nombres_vars <- c("beta0" = "Delta*beta[0]~~'(R'^2*' = 7.0%)'", 
                  "beta1" = "beta[1]~~'(R'^2*' = 81.6%)'", 
                  "beta2" = "Delta*beta[2]~~'(R'^2*' = 5.1%)'")

for(nm in c("beta0", "beta1", "beta2")) {
  df_temp <- data.frame(
    fecha = tail(fB$fecha, n_var), # Alinear las fechas al tamaño del VAR
    factor = nombres_vars[nm],
    observado = obs_mat[, nm] * 100, # Convertir a porcentaje si aplica
    predicho  = fit_mat[, nm] * 100
  )
  var_ajuste_df <- dplyr::bind_rows(var_ajuste_df, df_temp)
}

# Forzar el orden de los paneles
var_ajuste_df$factor <- factor(var_ajuste_df$factor, levels = nombres_vars)

# Pivotar a lo largo para mapear el color (Observado vs Predicho)
var_ajuste_long <- var_ajuste_df %>%
  tidyr::pivot_longer(c(observado, predicho), names_to = "metrica", values_to = "valor")

# 4. Graficar
p_var_fit <- ggplot(var_ajuste_long, aes(x = fecha, y = valor, color = metrica)) +
  geom_line(aes(linewidth = metrica, alpha = metrica)) +
  facet_wrap(~factor, ncol = 1, scales = "free_y", labeller = label_parsed) +
  scale_color_manual(values = c("observado" = "gray60", "predicho" = "#2166ac"),
                     labels = c("Observado (Real)", "Predicho por VAR(1)")) +
  scale_linewidth_manual(values = c("observado" = 0.6, "predicho" = 0.9), guide = "none") +
  scale_alpha_manual(values = c("observado" = 0.7, "predicho" = 1.0), guide = "none") +
  scale_x_date(date_breaks = "4 months", date_labels = "%b-%y") +
  labs(title    = "Ajuste dentro de muestra del Modelo VAR(1) por Ecuación",
       subtitle = "Comparación entre valores observados y predicciones contemporáneas del modelo estructural",
       x = NULL, y = "Valor del Factor / Cambio (%)",
       color = NULL,
       caption = "Muestra del modelo ajustada por rezagos (n = 421 observaciones).\nLos paneles en diferencias evidencian la naturaleza estocástica e impredecible de los cambios diarios.") +
  theme_minimal(base_size = 11) +
  theme(legend.position   = "bottom",
        strip.text        = element_text(face = "bold", size = 11),
        panel.grid.minor  = element_blank(),
        plot.caption      = element_text(size = 8, color = "gray40"))

guardar_grafico(p_var_fit, "G6.24b_VAR_ajuste_muestra", ancho = 12, alto = 10)
# ──  Granger — grafico de causalidad ─────────────────────────────────────

cat("Generando G6.25 — Grafo de causalidad de Granger...\n")

# Datos de Granger
granger_df <- data.frame(
  desde = c("β₂","β₀","β₁",  "β₀","β₂","β₁"),
  hasta = c("β₀","β₂","β₂",  "β₁","β₁","β₀"),
  F     = c(17.25, 9.99, 9.43, 3.20, 3.71, 0.72),
  pval  = c(0.000, 0.002, 0.002, 0.074, 0.055, 0.397),
  sig   = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))

# Posiciones de los nodos en triángulo
nodos <- data.frame(
  factor = c("β₀","β₁","β₂"),
  x      = c(0.5, 0.0, 1.0),
  y      = c(1.0, 0.0, 0.0),
  label  = c("β₀\nNivel\n(8.76%)", "β₁\nPendiente\n(−8.34%)",
             "β₂\nCurvatura\n(37.04%)"))

# Función para calcular punto en segmento con offset
offset_punto <- function(x1, y1, x2, y2, frac=0.18) {
  dx <- x2-x1; dy <- y2-y1
  dist <- sqrt(dx^2+dy^2)
  c(x1 + frac*dx/dist, y1 + frac*dy/dist)
}

# Solo flechas significativas
sig_df <- granger_df[granger_df$sig==TRUE,]

p_granger <- ggplot() +
  # Flechas NO significativas (grises, punteadas)
  geom_segment(
    data = granger_df[!granger_df$sig,],
    aes(x     = nodos$x[match(desde, nodos$factor)],
        y     = nodos$y[match(desde, nodos$factor)],
        xend  = nodos$x[match(hasta, nodos$factor)],
        yend  = nodos$y[match(hasta, nodos$factor)]),
    color    = "gray80", linewidth = 0.5,
    linetype = "dashed",
    arrow    = grid::arrow(length=unit(0.2,"cm"), type="open")) +
  # Flechas significativas — grosor proporcional a F
  geom_segment(
    data = sig_df,
    aes(x        = nodos$x[match(desde, nodos$factor)],
        y        = nodos$y[match(desde, nodos$factor)],
        xend     = nodos$x[match(hasta, nodos$factor)],
        yend     = nodos$y[match(hasta, nodos$factor)],
        linewidth= F/8),
    color = "#2166ac",
    arrow = grid::arrow(length=unit(0.3,"cm"), type="closed")) +
  # F-estadísticos sobre flechas significativas
  geom_label(
    data = sig_df,
    aes(x = (nodos$x[match(desde, nodos$factor)] +
               nodos$x[match(hasta, nodos$factor)]) / 2,
        y = (nodos$y[match(desde, nodos$factor)] +
               nodos$y[match(hasta, nodos$factor)]) / 2,
        label = sprintf("F=%.2f\np=%.3f", F, pval)),
    size=3, fill="white", color="#2166ac",
    label.size=0.3, label.padding=unit(0.15,"cm")) +
  # Nodos
  geom_point(data=nodos, aes(x=x, y=y),
             size=28, color="#2166ac", alpha=0.15) +
  geom_point(data=nodos, aes(x=x, y=y),
             size=28, color="#2166ac", shape=1, stroke=1.5) +
  geom_text(data=nodos, aes(x=x, y=y, label=label),
            size=3.5, fontface="bold", lineheight=1.1) +
  # Leyenda manual
  annotate("segment", x=0.05, xend=0.20, y=1.12, yend=1.12,
           color="#2166ac", linewidth=2,
           arrow=grid::arrow(length=unit(0.2,"cm"), type="closed")) +
  annotate("text", x=0.22, y=1.12,
           label="Causa Granger (p < 0.05)", hjust=0, size=3.2) +
  annotate("segment", x=0.55, xend=0.70, y=1.12, yend=1.12,
           color="gray70", linewidth=0.5, linetype="dashed",
           arrow=grid::arrow(length=unit(0.2,"cm"), type="open")) +
  annotate("text", x=0.72, y=1.12,
           label="No causa (p ≥ 0.05)", hjust=0, size=3.2) +
  scale_linewidth_continuous(range=c(0.5, 3), guide="none") +
  coord_cartesian(xlim=c(-0.25,1.25), ylim=c(-0.20,1.20)) +
  labs(title   = "Estructura de causalidad de Granger — Subperíodo B",
       subtitle = paste0(
         "Flechas azules: causalidad significativa al 5% | ",
         "Grosor proporcional al F-estadístico\n",
         "Flechas grises: no significativas"),
       caption = paste0(
         "Test: lmtest::grangertest(), orden VAR(1), α=5%.\n",
         "Series transformadas según ADF: Δβ₀, β₁ en niveles, Δβ₂."),
       x=NULL, y=NULL) +
  theme_void(base_size=12) +
  theme(plot.title    = element_text(face="bold", hjust=0.5),
        plot.subtitle = element_text(hjust=0.5, color="gray40"),
        plot.caption  = element_text(hjust=0.5, size=8, color="gray50"),
        plot.margin   = margin(10,10,10,10))

guardar_grafico(p_granger, "G6.25_Granger_grafo_causalidad",
                ancho=10, alto=9)

# ── FEVD — áreas apiladas ─────────────────────────────────────────────

cat("Generando G6.26 — FEVD áreas apiladas...\n")

fevd_res <- diag_din$fevd

# Convertir a data frame largo
fevd_df <- dplyr::bind_rows(lapply(c("beta0","beta1","beta2"), function(nm) {
  mat <- as.data.frame(fevd_res[[nm]])
  mat$h      <- seq_len(nrow(mat))
  mat$factor <- nm
  tidyr::pivot_longer(mat, -c(h, factor),
                      names_to="fuente", values_to="prop")
}))

# Etiquetas
fevd_df$factor_label <- factor(fevd_df$factor,
                               levels = c("beta0","beta1","beta2"),
                               labels = c("Varianza de Δβ₀ (Nivel)",
                                          "Varianza de β₁ (Pendiente)",
                                          "Varianza de Δβ₂ (Curvatura)"))

fevd_df$fuente_label <- factor(fevd_df$fuente,
                               levels = c("beta0","beta1","beta2"),
                               labels = c("Shocks de β₀ (Nivel)",
                                          "Shocks de β₁ (Pendiente)",
                                          "Shocks de β₂ (Curvatura)"))

colores_fevd <- c(
  "Shocks de β₀ (Nivel)"      = "#2166ac",
  "Shocks de β₁ (Pendiente)"  = "#fdae61",
  "Shocks de β₂ (Curvatura)"  = "#d73027")

p_fevd <- ggplot(fevd_df,
                 aes(x=h, y=prop*100, fill=fuente_label)) +
  geom_area(alpha=0.85, color="white", linewidth=0.3) +
  geom_vline(xintercept=21, linetype="dashed",
             color="gray30", linewidth=0.6) +
  annotate("text", x=21.5, y=50, label="h=21",
           hjust=0, size=3, color="gray30") +
  facet_wrap(~factor_label, ncol=1) +
  scale_fill_manual(values=colores_fevd) +
  scale_x_continuous(breaks=c(1,5,10,15,21)) +
  scale_y_continuous(labels=function(x) paste0(x,"%"),
                     breaks=seq(0,100,25)) +
  labs(title   = "Descomposición de varianza del error de predicción (FEVD)",
       subtitle = paste0(
         "Subperíodo B | VAR(1) con ordering Cholesky: β₀ → β₁ → β₂\n",
         "A h=21 días: β₀ propia=95.9% | β₁ propia=86.1% | ",
         "β₂ de β₀=62.5%"),
       x="Horizonte de predicción (días hábiles)",
       y="Proporción de la varianza (%)",
       fill=NULL,
       caption=paste0(
         "Cada panel muestra cómo se descompone la varianza del error ",
         "de predicción de cada factor entre sus propios shocks\n",
         "y los shocks de los otros factores del sistema, ",
         "a medida que aumenta el horizonte de predicción.")) +
  theme_minimal(base_size=11) +
  theme(legend.position  = "bottom",
        strip.text        = element_text(face="bold"),
        panel.grid.minor  = element_blank(),
        plot.caption      = element_text(size=8, color="gray40"))

guardar_grafico(p_fevd, "G6.26_FEVD_areas_apiladas",
                ancho=12, alto=12)

cat("\n✓ G6.23_ADF_niveles_vs_diferencias.png\n")
cat("✓ G6.24_beta0_residuos_AR.png\n")
cat("✓ G6.24_beta1_residuos_AR.png\n")
cat("✓ G6.24_beta2_residuos_AR.png\n")
cat("✓ G6.25_Granger_grafo_causalidad.png\n")
cat("✓ G6.26_FEVD_areas_apiladas.png\n")

# ──  VAR — valores reales vs ajustados ─────────────────────────────────

cat("Generando G6.27 — VAR fitted vs observed...\n")

var_mod     <- readRDS("output/var_modelo.rds")
diag_din    <- readRDS("output/diagnostico_dinamica.rds")
factores_GD <- readRDS("output/factores_GD.rds")

fB <- factores_GD %>%
  filter(fecha >= SUBPERIODOS$B$inicio,
         fecha <= SUBPERIODOS$B$fin,
         !is.na(beta0), beta0 > 0.01) %>%
  arrange(fecha)

b0 <- fB$beta0; b1 <- fB$beta1; b2 <- fB$beta2

n_min <- min(
  length(if(diag_din$adf$beta0$necesita_diff) diff(b0) else b0),
  length(if(diag_din$adf$beta1$necesita_diff) diff(b1) else b1),
  length(if(diag_din$adf$beta2$necesita_diff) diff(b2) else b2))

fitted_var <- as.data.frame(fitted(var_mod))
n_fit      <- nrow(fitted_var)
fechas_var <- tail(fB$fecha, n_min)
fechas_fit <- tail(fechas_var, n_fit)

real_b0 <- tail(tail(diff(b0), n_min), n_fit)
real_b1 <- tail(tail(b1,       n_min), n_fit)
real_b2 <- tail(tail(diff(b2), n_min), n_fit)

r2 <- sapply(c("beta0","beta1","beta2"), function(nm)
  summary(var_mod$varresult[[nm]])$r.squared)


configs <- list(
  list(nm   = "beta0",
       real = real_b0 * 100,
       fit  = fitted_var$beta0 * 100,
       # expresión para Delta β₀ con R²
       lab  = sprintf("Delta*beta[0]~~(R^2==%.1f*'%%')",
                      r2["beta0"]*100)),
  list(nm   = "beta1",
       real = real_b1 * 100,
       fit  = fitted_var$beta1 * 100,
       lab  = sprintf("beta[1]~~(R^2==%.1f*'%%')",
                      r2["beta1"]*100)),
  list(nm   = "beta2",
       real = real_b2 * 100,
       fit  = fitted_var$beta2 * 100,
       lab  = sprintf("Delta*beta[2]~~(R^2==%.1f*'%%')",
                      r2["beta2"]*100)))

df_var <- dplyr::bind_rows(lapply(configs, function(cfg) {
  data.frame(
    fecha = rep(fechas_fit, 2),
    valor = c(cfg$real, cfg$fit),
    tipo  = rep(c("Valor real","Ajustado VAR(1)"), each=n_fit),
    panel = cfg$lab)
}))

df_var$panel <- factor(df_var$panel,
                       levels = sapply(configs, `[[`, "lab"))
df_var$tipo  <- factor(df_var$tipo,
                       levels = c("Valor real","Ajustado VAR(1)"))

resid_df <- dplyr::bind_rows(lapply(configs, function(cfg) {
  res <- cfg$real - cfg$fit
  data.frame(
    fecha = fechas_fit,
    ymin  = cfg$fit - sd(res),
    ymax  = cfg$fit + sd(res),
    panel = cfg$lab)
}))
resid_df$panel <- factor(resid_df$panel,
                         levels = sapply(configs, `[[`, "lab"))

p_var_fit <- ggplot() +
  geom_ribbon(data = resid_df,
              aes(x=fecha, ymin=ymin, ymax=ymax),
              fill="#2166ac", alpha=0.08) +
  geom_line(data = df_var[df_var$tipo=="Valor real",],
            aes(x=fecha, y=valor),
            color="gray30", linewidth=0.55, alpha=0.9) +
  geom_line(data = df_var[df_var$tipo=="Ajustado VAR(1)",],
            aes(x=fecha, y=valor),
            color="#d73027", linewidth=0.9, alpha=0.85) +
  # label_parsed convierte las expresiones en notación matemática
  facet_wrap(~panel, ncol=1, scales="free_y",
             labeller = label_parsed) +
  scale_x_date(date_breaks="3 months", date_labels="%b\n%Y") +
  labs(title   = "VAR(1) \u2014 Valores reales vs ajustados por ecuaci\u00f3n",
       subtitle = paste0(
         "Subper\u00edodo B | ",
         "Gris: valor real | Rojo: ajustado VAR(1) | ",
         "Banda azul: \u00b11 DE del error"),
       x = NULL, y = NULL,
       caption = paste0(
         "\u03b2\u2081 en niveles: el VAR captura casi toda la variaci\u00f3n",
         " (R\u00b2=81.6%) \u2014 las dos l\u00edneas van pr\u00e1cticamente de la mano.\n",
         "\u0394\u03b2\u2080 y \u0394\u03b2\u2082 en diferencias: el VAR predice la tendencia",
         " suavizada mientras la serie real muestra saltos diarios",
         " impredecibles (R\u00b2\u22485-7%).\n",
         "Los R\u00b2 bajos no indican mal ajuste \u2014 las variaciones",
         " de las series I(1) son inherentemente impredecibles.")) +
  theme_minimal(base_size=11) +
  theme(strip.text       = element_text(face="bold", size=11),
        panel.grid.minor = element_blank(),
        plot.caption     = element_text(size=8, color="gray40",
                                        lineheight=1.3))

guardar_grafico(p_var_fit, "G6.27_VAR_fitted_vs_observed",
                ancho=13, alto=11)

# ── Scatter ───────────────────────────────────────────────────────────

scatter_df <- dplyr::bind_rows(lapply(configs, function(cfg) {
  data.frame(real=cfg$real, fit=cfg$fit, panel=cfg$lab)
}))
scatter_df$panel <- factor(scatter_df$panel,
                           levels = sapply(configs, `[[`, "lab"))

p_scatter <- ggplot(scatter_df, aes(x=fit, y=real)) +
  geom_point(color="#2166ac", alpha=0.25, size=0.9) +
  geom_abline(slope=1, intercept=0,
              color="#d73027", linewidth=0.8) +
  geom_smooth(method="lm", se=FALSE,
              color="gray40", linewidth=0.6, linetype="dashed") +
  facet_wrap(~panel, ncol=3, scales="free",
             labeller = label_parsed) +
  labs(title    = "VAR(1) \u2014 Dispersi\u00f3n real vs ajustado",
       subtitle = paste0(
         "L\u00ednea roja: ajuste perfecto (pendiente=1) | ",
         "L\u00ednea gris: regresi\u00f3n observada"),
       x       = "Valor ajustado (%)",
       y       = "Valor real (%)",
       caption = paste0(
         "Para \u03b2\u2081 los puntos se agrupan cerca de la diagonal.\n",
         "Para \u0394\u03b2\u2080 y \u0394\u03b2\u2082 la nube es dispersa \u2014 ",
         "las diferencias diarias son esencialmente ruido.")) +
  theme_minimal(base_size=11) +
  theme(strip.text   = element_text(face="bold", size=11),
        plot.caption = element_text(size=8, color="gray40"))

guardar_grafico(p_scatter, "G6.27b_VAR_scatter_real_vs_ajustado",
                ancho=13, alto=5)

cat("\u2713 G6.27_VAR_fitted_vs_observed.png\n")
cat("\u2713 G6.27b_VAR_scatter_real_vs_ajustado.png\n")