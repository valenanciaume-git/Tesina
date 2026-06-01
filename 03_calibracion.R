# ==============================================================================
# 03_calibracion.R
# Calibración NS por precio — Series GD (Ley NY) y AL (Ley Argentina)
#
#
# Outputs:
#   output/calibracion_GD_sp.rds   subperíodo B, sin ponderadores
#   output/calibracion_GD_pd.rds   subperíodo B, pond. 1/D*
#   output/calibracion_GD_A.rds    subperíodo A
#   output/calibracion_AL_sp.rds   subperíodo B, sin ponderadores
#   output/calibracion_AL_pd.rds   subperíodo B, pond. 1/D*
#   output/calibracion_AL_A.rds    subperíodo A
#   output/factores_GD.rds
#   output/factores_AL.rds
#   output/spread_GD_AL.rds
#   output/comparacion_ponderadores.rds
# ==============================================================================

library(dplyr)
source("scripts/00_funciones.R")

MIN_BONOS_GD <- 4
MIN_BONOS_AL <- 3

panel_diario <- readRDS("data/panel_diario.rds")
cf_list      <- readRDS("data/cashflows.rds")
ancla_df     <- readRDS("data/ancla_usd.rds")

if (!dir.exists("output")) dir.create("output")

# ── Función de calibración por día ────────────────────────────────────────────

calibrar_dia <- function(val_date, bonos, precios_obs, cf_list,
                         min_bonos = 4, usar_ponderadores = FALSE) {
  idx_disp <- !is.na(precios_obs)
  if (sum(idx_disp) < min_bonos) return(NULL)
  
  ancla        <- get_ancla(val_date, ancla_df)
  ponderadores <- NULL
  if (usar_ponderadores) {
    row_tirs     <- panel_diario[panel_diario$fecha == val_date, ]
    tirs_dia     <- as.list(row_tirs)
    ponderadores <- calcular_ponderadores(bonos, val_date, cf_list,
                                        tirs_obs = tirs_dia)}
  
  best_sse <- Inf; best_par <- NULL
  
  for (i in seq_len(nrow(GRILLA_INICIO))) {
    p0  <- as.numeric(GRILLA_INICIO[i, ])
    opt <- tryCatch(
      optim(p0, loss_precio_ancla,
            method   = "L-BFGS-B",
            lower    = BOUNDS_LOWER,
            upper    = BOUNDS_UPPER,
            val_date     = val_date,
            bonos        = bonos,
            precios_obs  = precios_obs,
            cf_list      = cf_list,
            ancla        = ancla,
            ponderadores = ponderadores,
            control  = list(maxit = 1000, factr = 1e8)),
      error = function(e) NULL)
    if (!is.null(opt) && opt$value < best_sse) {
      best_sse <- opt$value; best_par <- opt$par
    }
  }
  
  if (is.null(best_par)) return(NULL)
  
  p_teo     <- sapply(bonos, function(b)
    precio_teo(best_par, val_date, cf_list[[b]]))
  # Errores — solo sobre bonos disponibles
  obs_disp <- precios_obs[idx_disp]
  teo_disp <- p_teo[idx_disp]
  err      <- teo_disp - obs_disp          # error absoluto en USD
  err_pct  <- err / obs_disp * 100         # error porcentual (%)
  
  # ── Métricas calculadas directamente desde precios ──────────────────────────
  # Todas las métricas se calculan sobre los bonos disponibles (idx_disp)
  # RMSE  [USD]:   raíz del promedio de errores cuadráticos
  # RMSPE [%]:     raíz del promedio de errores porcentuales cuadráticos
  # MAPE  [%]:     promedio de errores porcentuales absolutos
  # MAE   [USD]:   promedio de errores absolutos
  # ME    [USD]:   sesgo medio (positivo = sobreestima)
  # MPE   [%]:     sesgo porcentual medio
  rmse_val  <- sqrt(mean(err^2))
  rmspe_val <- sqrt(mean(err_pct^2))
  mape_val  <- mean(abs(err_pct))
  mae_val   <- mean(abs(err))
  me_val    <- mean(err)
  mpe_val   <- mean(err_pct)
  
  # Vector completo de error_pct 
  error_pct_full <- rep(NA_real_, length(bonos))
  error_pct_full[idx_disp] <- err_pct
  
  list(
    fecha     = as.Date(val_date),
    beta0     = best_par[1], beta1 = best_par[2],
    beta2     = best_par[3], tau   = best_par[4],
    lambda    = 1 / best_par[4],
    ancla     = ancla,
    rmse      = rmse_val,
    rmspe     = rmspe_val,
    mape      = mape_val,
    mae       = mae_val,
    me        = me_val,
    mpe       = mpe_val,
    n_bonos   = sum(idx_disp),
    ponderada = usar_ponderadores,
    p_obs     = precios_obs,
    p_teo     = p_teo,
    error_pct = error_pct_full,
    excupon   = as.Date(val_date) %in% FECHAS_EXCUPON
  )
}
# ── Loop de calibración por subperíodo ────────────────────────────────────────

calibrar_subperiodo <- function(sub_key, bonos, min_bonos,
                                usar_ponderadores = FALSE,
                                etiqueta = "") {
  sub   <- SUBPERIODOS[[sub_key]]
  panel <- panel_diario[
    panel_diario$fecha >= sub$inicio &
      panel_diario$fecha <= sub$fin &
      !panel_diario$outlier, ]
  
  tipo <- ifelse(usar_ponderadores, "Pond. 1/D*", "Sin pond.")
  cat(sprintf("\n== Sub.%s — %s — %s (%d días) ==\n",
              sub_key, etiqueta, tipo, nrow(panel)))
  cat(sprintf("%-12s %8s %8s %8s %8s %7s\n",
              "Fecha","b0","b1","b2","RMSE","MAPE%"))
  cat(strrep("-", 57), "\n")
  
  resultados <- vector("list", nrow(panel))
  
  for (i in seq_len(nrow(panel))) {
    fecha <- panel$fecha[i]
    obs   <- as.numeric(panel[i, bonos])
    res   <- calibrar_dia(fecha, bonos, obs, cf_list,
                          min_bonos        = min_bonos,
                          usar_ponderadores = usar_ponderadores)
    resultados[[i]] <- res
    
    if (!is.null(res))
      cat(sprintf("%-12s %8.4f %8.4f %8.4f %8.4f %7.3f%s\n",
                  as.character(fecha),
                  res$beta0, res$beta1, res$beta2,
                  res$rmse,  res$mape,
                  ifelse(res$excupon, " *", "")))
  }
  resultados
}

extraer_factores <- function(resultados, serie) {
  dplyr::bind_rows(Filter(Negate(is.null), lapply(resultados, function(r) {
    if (is.null(r)) return(NULL)
    data.frame(
      fecha   = r$fecha,   serie   = serie,
      beta0   = r$beta0,   beta1   = r$beta1,
      beta2   = r$beta2,   tau     = r$tau,
      lambda  = r$lambda,  rmse    = r$rmse,
      mape    = r$mape,    n_bonos = r$n_bonos,
      excupon = r$excupon
    )
  })))
}

# ══════════════════════════════════════════════════════════════════════════════
cat("=== Calibración NS — τ =", TAU_FIJO, "años ===\n")
cat(sprintf("Grilla: %d puntos | Ancla: PF USD 30d BCRA\n\n",
            nrow(GRILLA_INICIO)))
cat("BONOS_GD: ", paste(BONOS_GD, collapse=", "), "\n")
cat("BONOS_AL: ", paste(BONOS_AL, collapse=", "), "\n\n")

# ── 1. Serie GD ───────────────────────────────────────────────────────────────

cat("══════════════════════════════════════════\n")
cat("SERIE GD — Ley Nueva York\n")
cat("══════════════════════════════════════════\n")

res_GD_B_sp <- calibrar_subperiodo("B", BONOS_GD, MIN_BONOS_GD, FALSE, "GD")
res_GD_B_pd <- calibrar_subperiodo("B", BONOS_GD, MIN_BONOS_GD, TRUE,  "GD")
res_GD_A_sp <- calibrar_subperiodo("A", BONOS_GD, MIN_BONOS_GD, FALSE, "GD")

df_GD_sp <- extraer_factores(res_GD_B_sp, "GD")
df_GD_pd <- extraer_factores(res_GD_B_pd, "GD")

mejor_GD <- ifelse(mean(df_GD_sp$mape) <= mean(df_GD_pd$mape),
                   "Sin ponderadores", "Ponderado 1/D*")
cat(sprintf("\nMejor modelo GD: %s\n", mejor_GD))
cat(sprintf("  MAPE sin pond: %.4f%%  |  MAPE pond: %.4f%%\n",
            mean(df_GD_sp$mape), mean(df_GD_pd$mape)))

factores_GD <- if (mejor_GD == "Sin ponderadores") df_GD_sp else df_GD_pd

# ── 2. Serie AL ───────────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════\n")
cat("SERIE AL — Ley Argentina\n")
cat("══════════════════════════════════════════\n")

res_AL_B_sp <- calibrar_subperiodo("B", BONOS_AL, MIN_BONOS_AL, FALSE, "AL")
res_AL_B_pd <- calibrar_subperiodo("B", BONOS_AL, MIN_BONOS_AL, TRUE,  "AL")
res_AL_A_sp <- calibrar_subperiodo("A", BONOS_AL, MIN_BONOS_AL, FALSE, "AL")

df_AL_sp <- extraer_factores(res_AL_B_sp, "AL")
df_AL_pd <- extraer_factores(res_AL_B_pd, "AL")

mejor_AL <- ifelse(mean(df_AL_sp$mape) <= mean(df_AL_pd$mape),
                   "Sin ponderadores", "Ponderado 1/D*")
cat(sprintf("\nMejor modelo AL: %s\n", mejor_AL))
cat(sprintf("  MAPE sin pond: %.4f%%  |  MAPE pond: %.4f%%\n",
            mean(df_AL_sp$mape), mean(df_AL_pd$mape)))

factores_AL  <- if (mejor_AL == "Sin ponderadores") df_AL_sp else df_AL_pd

# ── 3. Spread GD-AL ───────────────────────────────────────────────────────────

spread_GD_AL <- dplyr::inner_join(
  factores_GD[, c("fecha","beta0","beta1","beta2","rmse","mape")],
  factores_AL[, c("fecha","beta0","beta1","beta2","rmse","mape")],
  by = "fecha", suffix = c("_GD","_AL")
) %>%
  dplyr::mutate(
    spread_nivel     = beta0_GD - beta0_AL,
    spread_pendiente = beta1_GD - beta1_AL,
    spread_curvatura = beta2_GD - beta2_AL
  )

cat(sprintf("\n=== Spread GD-AL (subperíodo B) ===\n"))
cat(sprintf("  Spread nivel promedio:    %+.4f (%+.2f pp)\n",
            mean(spread_GD_AL$spread_nivel),
            mean(spread_GD_AL$spread_nivel)*100))

# ── 4. Comparación ponderadores ───────────────────────────────────────────────

cat("\n=== Comparación ponderadores — Subperíodo B ===\n")
cat(sprintf("%-6s %-16s %8s %8s\n", "Serie","Modelo","MAPE%","RMSE"))
cat(strrep("-", 42), "\n")
for (nm in list(
  list("GD", "Sin pond.",  df_GD_sp),
  list("GD", "Pond. 1/D*", df_GD_pd),
  list("AL", "Sin pond.",  df_AL_sp),
  list("AL", "Pond. 1/D*", df_AL_pd)
)) {
  cat(sprintf("%-6s %-16s %8.4f %8.4f\n",
              nm[[1]], nm[[2]],
              mean(nm[[3]]$mape), mean(nm[[3]]$rmse)))
}

# ── 5. Guardar ────────────────────────────────────────────────────────────────

saveRDS(res_GD_B_sp,  "output/calibracion_GD_sp.rds")
saveRDS(res_GD_B_pd,  "output/calibracion_GD_pd.rds")
saveRDS(res_GD_A_sp,  "output/calibracion_GD_A.rds")
saveRDS(res_AL_B_sp,  "output/calibracion_AL_sp.rds")
saveRDS(res_AL_B_pd,  "output/calibracion_AL_pd.rds")
saveRDS(res_AL_A_sp,  "output/calibracion_AL_A.rds")
saveRDS(factores_GD,  "output/factores_GD.rds")
saveRDS(factores_AL,  "output/factores_AL.rds")
saveRDS(spread_GD_AL, "output/spread_GD_AL.rds")
saveRDS(list(
  GD_sin_pond = df_GD_sp,
  GD_pond_dur = df_GD_pd,
  AL_sin_pond = df_AL_sp,
  AL_pond_dur = df_AL_pd,
  mejor_GD    = mejor_GD,
  mejor_AL    = mejor_AL
), "output/comparacion_ponderadores.rds")

cat("\n✓ output/calibracion_GD_sp.rds\n")
cat("✓ output/calibracion_GD_pd.rds\n")
cat("✓ output/calibracion_GD_A.rds\n")
cat("✓ output/calibracion_AL_sp.rds\n")
cat("✓ output/calibracion_AL_pd.rds\n")
cat("✓ output/calibracion_AL_A.rds\n")
cat("✓ output/factores_GD.rds\n")
cat("✓ output/factores_AL.rds\n")
cat("✓ output/spread_GD_AL.rds\n")
cat("✓ output/comparacion_ponderadores.rds\n")
cat("\nNota: fechas ex-cupón marcadas con (*)\n")