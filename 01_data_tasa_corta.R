# ==============================================================================
# 01_data_tasa_corta.R
# Importación y construcción del panel mensual de la tasa ancla
#
# Fuente: Bloomberg — BADLARU Index
# Descripción: tasa promedio de depósitos a plazo fijo en dólares
#              a 30-35 días, bancos privados (BCRA)
# Uso: ancla nominal del tramo corto en la función de pérdida NS
# 
#
# Output: data/ancla_usd.rds
# ==============================================================================

library(readxl)
library(dplyr)
library(lubridate)

RUTA_ANCLA   <- "data/PF_USD_30-35d.xlsx"   #Excel con los valores del ancla
FECHA_INICIO <- as.Date("2021-01-01")
FECHA_FIN    <- as.Date("2026-04-30")

# ── Carga ─────────────────────────────────────────────────────────────────────
# El archivo tiene 6 filas de encabezado 
# Las fechas vienen en orden descendente
# Los valores son tasas anuales en porcentaje (ej: 0.6875 = 0.6875%)

ancla_raw <- read_excel(RUTA_ANCLA,
                        skip      = 6,
                        col_names = c("fecha", "tasa_pct"),
                        col_types = c("date", "numeric"))

# ── Limpieza ──────────────────────────────────────────────────────────────────

ancla_raw <- ancla_raw %>%
  filter(!is.na(fecha), !is.na(tasa_pct)) %>%
  mutate(
    fecha = as.Date(fecha),
    tasa  = tasa_pct / 100    # convertir a decimal: 0.6875% -> 0.006875
  ) %>%
  filter(fecha >= FECHA_INICIO, fecha <= FECHA_FIN) %>%
  arrange(fecha)

# ── Panel mensual: último día hábil de cada mes ───────────────────────────────


ancla_mensual <- ancla_raw %>%
  mutate(anio_mes = floor_date(fecha, "month")) %>%
  group_by(anio_mes) %>%
  filter(fecha == max(fecha)) %>%
  ungroup() %>%
  select(fecha, tasa)

# ── Reporte ───────────────────────────────────────────────────────────────────

cat("=== Ancla nominal — PF USD 30 días ===\n")
cat(sprintf("Observaciones diarias:   %d\n", nrow(ancla_raw)))
cat(sprintf("Observaciones mensuales: %d\n", nrow(ancla_mensual)))
cat(sprintf("Rango: %s a %s\n",
            min(ancla_mensual$fecha), max(ancla_mensual$fecha)))
cat(sprintf("Tasa mínima: %.4f%%\n", min(ancla_mensual$tasa) * 100))
cat(sprintf("Tasa máxima: %.4f%%\n", max(ancla_mensual$tasa) * 100))
cat(sprintf("Tasa promedio: %.4f%%\n", mean(ancla_mensual$tasa) * 100))

cat("\nPanel mensual completo:\n")
print(ancla_mensual, n = Inf)

# ── Guardar ───────────────────────────────────────────────────────────────────

if (!dir.exists("data")) dir.create("data")
saveRDS(ancla_mensual, "data/ancla_usd.rds")
cat("\n✓ data/ancla_usd.rds\n")
