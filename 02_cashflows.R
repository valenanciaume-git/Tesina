# ==============================================================================
# 02_cashflows.R
# Parseo de cash flows desde Excel
#
# Estructura del Excel por bono (6 columnas):
#   col + 0: fecha del flujo
#   col + 1: amortización parcial de capital
#   col + 2: saldo del período anterior (VR antes del pago)
#   col + 3: FF total = interés + amortización  ← lo que se usa para descontar
#   col + 4: saldo actual (VR después del pago)
#   col + 5: interés sobre el saldo
#
# Nota terminológica:
#   - "Amortización": devolución parcial de capital en cada fecha de pago
#   - "Interés": pago de renta sobre el saldo pendiente
#   - "FF": flujo total que recibe el tenedor = interés + amortización
#   - El período de gracia es el tramo donde la amortización = 0
#     (el tenedor solo cobra intereses, no recupera capital)
#
# Output: data/cashflows.rds
# ==============================================================================

library(readxl)
source("scripts/00_funciones.R")

RUTA_CF <- "data/cash_flows_tesis.xlsx"

BONOS <- c("AL29","AL30","AL35","AE38","AL41",
           "GD29","GD30","GD35","GD38","GD41","GD46")


raw <- read_excel(RUTA_CF, sheet = "Cash flows",
                  col_names = FALSE, col_types = "text")

parsear_bono <- function(bono_idx) {
  col_fecha <- bono_idx * 6 + 1
  col_amort <- bono_idx * 6 + 2   # amortización parcial de capital
  col_saldo_ant <- bono_idx * 6 + 3   # saldo período anterior
  col_ff    <- bono_idx * 6 + 4   # FF total = interés + amortización
  col_saldo_act <- bono_idx * 6 + 5   # saldo actual
  col_int   <- bono_idx * 6 + 6   # interés
  
  fechas_num <- suppressWarnings(
    as.numeric(unlist(raw[5:nrow(raw), col_fecha])))
  amorts     <- suppressWarnings(
    as.numeric(unlist(raw[5:nrow(raw), col_amort])))
  ffs        <- suppressWarnings(
    as.numeric(unlist(raw[5:nrow(raw), col_ff])))
  ints       <- suppressWarnings(
    as.numeric(unlist(raw[5:nrow(raw), col_int])))
  saldos_ant <- suppressWarnings(
    as.numeric(unlist(raw[5:nrow(raw), col_saldo_ant])))
  
  validos <- !is.na(fechas_num) & !is.na(ffs)
  
  data.frame(
    fecha     = as.Date(fechas_num[validos], origin = "1899-12-30"),
    amort     = amorts[validos],
    saldo_ant = saldos_ant[validos],
    interes   = ints[validos],
    ff        = ffs[validos]
  )
}

cf_list <- lapply(seq_along(BONOS) - 1, parsear_bono)
names(cf_list) <- BONOS

# ── Verificación ──────────────────────────────────────────────────────────────

cat("=== Cash flows cargados ===\n\n")

cat(sprintf("  %-6s %6s %14s %16s %10s %14s\n",
            "Bono","Flujos","Primer flujo",
            "Primera amort","Suma FF","Último flujo"))
cat(strrep("-", 75), "\n")

for (b in BONOS) {
  cf <- cf_list[[b]]
  
  # Fin del período de gracia = último flujo donde amortización = 0
  # = fecha anterior al primer flujo con amortización > 0
  idx_amort <- which(!is.na(cf$amort) & cf$amort > 0)
  primera_amort <- if (length(idx_amort) > 0)
    as.character(cf$fecha[min(idx_amort)])
  else "sin amortización"
  
  cat(sprintf("  %-6s %6d %14s %16s %10.2f %14s\n",
              b, nrow(cf),
              as.character(min(cf$fecha)),
              primera_amort,
              sum(cf$ff),
              as.character(max(cf$fecha))))
}

# Detalle por bono: primeros 3 y últimos 2 flujos
cat("\n=== Detalle de flujos por bono ===\n")
for (b in BONOS) {
  cf <- cf_list[[b]]
  
  idx_amort <- which(!is.na(cf$amort) & cf$amort > 0)
  primera_amort <- if (length(idx_amort) > 0)
    as.character(cf$fecha[min(idx_amort)])
  else "sin amortización"
  
  cat(sprintf("\n%s (%d flujos | período de gracia hasta: %s)\n",
              b, nrow(cf), primera_amort))
  cat(sprintf("  %-12s %10s %12s %10s %10s\n",
              "Fecha","Interés","Amortización","FF total","Saldo ant."))
  
  mostrar <- unique(c(1:3, (nrow(cf)-1):nrow(cf)))
  mostrar  <- mostrar[mostrar >= 1 & mostrar <= nrow(cf)]
  
  prev <- 0
  for (i in mostrar) {
    if (i > prev + 1 && prev > 0) cat("  ...\n")
    cat(sprintf("  %-12s %10.4f %12.4f %10.4f %10.4f\n",
                as.character(cf$fecha[i]),
                ifelse(is.na(cf$interes[i]), 0, cf$interes[i]),
                ifelse(is.na(cf$amort[i]),   0, cf$amort[i]),
                cf$ff[i],
                ifelse(is.na(cf$saldo_ant[i]), 0, cf$saldo_ant[i])))
    prev <- i
  }
}

# ── Guardar ───────────────────────────────────────────────────────────────────

if (!dir.exists("data")) dir.create("data")
saveRDS(cf_list, "data/cashflows.rds")
cat("\n✓ data/cashflows.rds\n")
