# ==============================================================================
# 01_data.R
# Carga y limpieza del panel de precios diarios
#
# Output: data/panel_diario.rds, data/panel_mensual.rds
# ==============================================================================

library(readxl)
library(dplyr)
library(lubridate)


# ── Parámetros ────────────────────────────────────────────────────────────────

RUTA_PRECIOS <- "data/precios_y_tir_bonos.xlsx"    #Excel con precios 1816
FECHA_INICIO <- as.Date("2021-01-01")
FECHA_FIN    <- as.Date("2026-04-30")   



# Columnas de precio en el Excel (1-based)
# Estructura: fecha | AL29 | AL30 | AE38 | AL41 | AL35 |
#                     GD29 | GD30 | GD38 | GD46 | GD41 | GD35
# Las columnas 12-22 son TIR  — no se usan
COL_PRECIOS <- c(
  AL29 = 2, AL30 = 3, AE38 = 4, AL41 = 5, AL35 = 6,
  GD29 = 7, GD30 = 8, GD38 = 9, GD46 = 10, GD41 = 11, GD35 = 12
)
# Columnas de TIR (en % anual — hay que dividir por 100)
COL_TIR <- c(
  AL29 = 13, AL30 = 14, AE38 = 15, AL41 = 16, AL35 = 17,
  GD29 = 18, GD30 = 19, GD38 = 20, GD46 = 21, GD41 = 22, GD35 = 23
)



# ── Carga ─────────────────────────────────────────────────────────────────────

raw <- read_excel(RUTA_PRECIOS, col_names = FALSE)

# Fila 1 = nombres de bonos, fila 2 en adelante = datos
# Columna 1 = fecha (número serial Excel)
fechas <- as.Date(as.POSIXct(unlist(raw[-1, 1]),
                             origin = "1970-01-01", tz = "UTC"))

# Construir data frame de precios
precios_df <- as.data.frame(
  lapply(names(COL_PRECIOS), function(bono) {
    suppressWarnings(as.numeric(unlist(raw[-1, COL_PRECIOS[[bono]]])))
  })
)
names(precios_df) <- names(COL_PRECIOS)


tir_df <- as.data.frame(
  lapply(names(COL_TIR), function(bono) {
    suppressWarnings(as.numeric(unlist(raw[-1, COL_TIR[[bono]]])))
  })
)
names(tir_df) <- paste0("tir_", names(COL_TIR))

# Unir precios y TIR en el panel
panel_diario <- cbind(fecha = fechas, precios_df, tir_df)

# ── Limpieza ──────────────────────────────────────────────────────────────────

panel_diario <- panel_diario %>%
  filter(!is.na(fecha)) %>%
  filter(fecha >= FECHA_INICIO, fecha <= FECHA_FIN) %>%
  mutate(outlier = fecha %in% OUTLIERS) %>%
  # Precios de 0 o menores a 1 son datos faltantes en el archivo
  mutate(across(all_of(names(COL_PRECIOS)),
                ~ ifelse(is.na(.), NA_real_, .))) %>%
  arrange(fecha)

# ── Panel mensual: último día hábil de cada mes ───────────────────────────────

panel_mensual <- panel_diario %>%
  mutate(anio_mes = floor_date(fecha, "month")) %>%
  group_by(anio_mes) %>%
  filter(fecha == max(fecha)) %>%
  ungroup()

# Quitar columnas que no necesitamos en el mensual
panel_mensual <- panel_mensual[, !grepl("^tir_|^anio_mes$",
                                        names(panel_mensual))]

# ── Reporte ───────────────────────────────────────────────────────────────────

cat("=== Panel diario ===\n")
cat(sprintf("Rango:          %s -> %s\n",
            min(panel_diario$fecha), max(panel_diario$fecha)))
cat(sprintf("Observaciones:  %d\n", nrow(panel_diario)))
cat(sprintf("Outliers marcados: %d\n", sum(panel_diario$outlier)))

cat("\nCompletitud por bono (% días con precio):\n")
for (b in names(COL_PRECIOS)) {
  cat(sprintf("  %-6s: %5.1f%%\n", b,
              mean(!is.na(panel_diario[[b]])) * 100))
}

cat("\n=== Panel mensual ===\n")
cat(sprintf("Observaciones:  %d meses\n", nrow(panel_mensual)))
cat(sprintf("Rango:          %s -> %s\n",
            min(panel_mensual$fecha), max(panel_mensual$fecha)))

cat("\nEstadística descriptiva — precios mensuales:\n")
cat(sprintf("  %-6s %7s %7s %7s %4s\n", "Bono", "Min", "Med", "Max", "NAs"))
for (b in names(COL_PRECIOS)) {
  x <- panel_mensual[[b]]
  x <- as.numeric(x)   # forzar numérico por las dudas
  cat(sprintf("  %-6s %7.2f %7.2f %7.2f %4d\n",
              b,
              min(x, na.rm = TRUE),
              median(x, na.rm = TRUE),
              max(x, na.rm = TRUE),
              sum(is.na(x))))
}

# ── Guardar ───────────────────────────────────────────────────────────────────

if (!dir.exists("data")) dir.create("data")
saveRDS(panel_diario,  "data/panel_diario.rds")
saveRDS(panel_mensual, "data/panel_mensual.rds")

cat("\n✓ data/panel_diario.rds\n")
cat("✓ data/panel_mensual.rds\n")
