# ==============================================================================
# ==============================================================================
# 06_prediccion.R
# Evaluación del poder predictivo fuera de muestra
#
# ETAPA 1: Validación interna (subperíodo B) — AR vs VAR vs RW
#   - Ventana expandida: entrena en primeras 80% observaciones,
#     predice h pasos adelante con predicción dinámica iterativa
#   - AR: predicción dinámica correcta — propaga h pasos usando
#     valores predichos anteriores como inputs
#   - VAR: especificación mixta (Δβ₀, β₁, Δβ₂) consistente con 05
#   - RW: benchmark — el mejor predictor para series I(1)
#
# ETAPA 2: Validación externa (datos reales 2026)
#
# ETAPA 2b: Valuación de cartera hipotética SSN
#   V_NS (curva NS calibrada) vs V_Obs (precio de mercado)
#
# Requiere: output/factores_GD.rds, output/ar_modelos.rds,
#           output/var_modelo.rds, output/diagnostico_dinamica.rds,
#           data/panel_diario.rds, data/cashflows.rds,
#           output/calibracion_GD_sp.rds
# Output:   output/prediccion_resultados.rds
# ==============================================================================

set.seed(42)

library(dplyr)
library(vars)
library(forecast)

source("scripts/00_funciones.R")

factores_full <- readRDS("output/factores_GD.rds")
ar_modelos    <- readRDS("output/ar_modelos.rds")
var_mod       <- readRDS("output/var_modelo.rds")
diag_din      <- readRDS("output/diagnostico_dinamica.rds")
panel_diario  <- readRDS("data/panel_diario.rds")
cf_list       <- readRDS("data/cashflows.rds")
cal_GD        <- readRDS("output/calibracion_GD_sp.rds")

# Solo subperíodo B
factores_B <- factores_full %>%
  filter(fecha >= SUBPERIODOS$B$inicio,
         fecha <= SUBPERIODOS$B$fin,
         !is.na(beta0), beta0 > 0.01) %>%
  arrange(fecha)

n_total <- nrow(factores_B)
n_train <- floor(n_total * 0.80)
n_test  <- n_total - n_train

cat("=== Evaluación predictiva fuera de muestra ===\n\n")
cat(sprintf("Subperíodo B: %d observaciones\n", n_total))
cat(sprintf("Entrenamiento (80%%): %d obs hasta %s\n",
            n_train, factores_B$fecha[n_train]))
cat(sprintf("Validación (20%%):   %d obs desde %s\n\n",
            n_test, factores_B$fecha[n_train+1]))

# ── Función: predicción AR dinámica h pasos adelante ─────────────────────────
# Propaga iterativamente: usa predicciones anteriores como inputs
# para los pasos siguientes — predicción genuinamente dinámica

predecir_AR_h <- function(serie, orden_ar, h, necesita_diff, ultimo_nivel) {
  # Estimar modelo AR sobre la serie
  mod <- tryCatch(
    Arima(serie, order=c(orden_ar, 0, 0), include.mean=TRUE),
    error=function(e) NULL)
  if (is.null(mod)) return(NA)
  
  # forecast() de {forecast} ya hace predicción dinámica correctamente
  fc  <- forecast(mod, h=h)
  val <- as.numeric(fc$mean)[h]
  
  # Reconstruir nivel si la serie estaba diferenciada
  if (necesita_diff) val <- ultimo_nivel + val
  
  val
}

# ── ETAPA 1: Validación interna ───────────────────────────────────────────────

cat("=== ETAPA 1: Validación interna (subperíodo B) ===\n")
cat("Método: ventana expandida, predicción dinámica h pasos\n\n")

rmse_resultados <- list()
horizontes <- c(1, 5, 10, 21)

for (h in horizontes) {
  pred_AR  <- list()
  pred_VAR <- list()
  pred_RW  <- list()
  
  n_valid <- n_test - h + 1
  if (n_valid <= 0) next
  
  for (i in seq_len(n_valid)) {
    idx_fin <- n_train + i - 1
    datos_i <- factores_B[seq_len(idx_fin), ]
    
    # ── Caminata aleatoria ──
    pred_RW[[i]] <- c(
      beta0 = tail(datos_i$beta0, 1),
      beta1 = tail(datos_i$beta1, 1),
      beta2 = tail(datos_i$beta2, 1))
    
    # ── AR univariado — predicción dinámica correcta ──
    pred_AR_i <- setNames(numeric(3), c("beta0","beta1","beta2"))
    for (nm in c("beta0","beta1","beta2")) {
      nd  <- diag_din$adf[[nm]]$necesita_diff
      ul  <- tail(datos_i[[nm]], 1)   # último nivel de ESTE subconjunto
      serie_i <- if(nd) diff(datos_i[[nm]]) else datos_i[[nm]]
      ord <- ar_modelos[[nm]]$modelo$arma[1]
      pred_AR_i[nm] <- predecir_AR_h(serie_i, ord, h, nd, ul)
    }
    pred_AR[[i]] <- pred_AR_i
    
    # ── VAR — especificación mixta consistente con script 05 ──
    # Δβ₀, β₁ en niveles, Δβ₂
    mat_i <- cbind(
      beta0 = if(diag_din$adf$beta0$necesita_diff)
        diff(datos_i$beta0) else datos_i$beta0,
      beta1 = if(diag_din$adf$beta1$necesita_diff)
        diff(datos_i$beta1) else datos_i$beta1,
      beta2 = if(diag_din$adf$beta2$necesita_diff)
        diff(datos_i$beta2) else datos_i$beta2)
    
    var_i <- tryCatch(
      VAR(mat_i, p=diag_din$orden_var, type="const"),
      error=function(e) NULL)
    
    if (!is.null(var_i)) {
      fc_var <- predict(var_i, n.ahead=h)
      # Reconstruir niveles para Δβ₀ y Δβ₂
      ul_b0 <- tail(datos_i$beta0, 1)
      ul_b2 <- tail(datos_i$beta2, 1)
      pred_VAR[[i]] <- c(
        beta0 = if(diag_din$adf$beta0$necesita_diff)
          ul_b0 + fc_var$fcst$beta0[h,"fcst"]
        else fc_var$fcst$beta0[h,"fcst"],
        beta1 = fc_var$fcst$beta1[h,"fcst"],
        beta2 = if(diag_din$adf$beta2$necesita_diff)
          ul_b2 + fc_var$fcst$beta2[h,"fcst"]
        else fc_var$fcst$beta2[h,"fcst"])
    } else {
      pred_VAR[[i]] <- pred_RW[[i]]
    }
  }
  
  # Observaciones realizadas h pasos adelante
  obs_test <- factores_B[(n_train+h):n_total,
                         c("beta0","beta1","beta2")]
  
  rmse_factor <- function(preds, obs) {
    sqrt(mean(sapply(seq_along(preds), function(i) {
      p <- preds[[i]]
      o <- as.numeric(obs[i, ])
      if (any(is.na(p)) || any(is.na(o))) return(NA)
      sum((p - o)^2)
    }), na.rm=TRUE))
  }
  
  rmse_RW  <- rmse_factor(pred_RW,  obs_test[seq_len(n_valid), ])
  rmse_AR  <- rmse_factor(pred_AR,  obs_test[seq_len(n_valid), ])
  rmse_VAR <- rmse_factor(pred_VAR, obs_test[seq_len(n_valid), ])
  
  rmse_resultados[[as.character(h)]] <- list(
    h=h, RW=rmse_RW, AR=rmse_AR, VAR=rmse_VAR)
  
  cat(sprintf(
    "h=%2d días: RW=%.4f  AR=%.4f  VAR=%.4f  |  AR %s RW  |  VAR %s RW\n",
    h, rmse_RW, rmse_AR, rmse_VAR,
    ifelse(rmse_AR  < rmse_RW, "mejora", "no mejora"),
    ifelse(rmse_VAR < rmse_RW, "mejora", "no mejora")))
}

# ── ETAPA 2: Validación externa 2026 ─────────────────────────────────────────

cat("\n=== ETAPA 2: Validación externa — datos reales 2026 ===\n\n")

ultimo_factor <- factores_B %>%
  filter(fecha <= as.Date("2025-12-31")) %>%
  tail(1)

cat(sprintf("Último factor estimado: %s\n", ultimo_factor$fecha))
cat(sprintf("β₀=%.4f  β₁=%.4f  β₂=%.4f\n\n",
            ultimo_factor$beta0, ultimo_factor$beta1,
            ultimo_factor$beta2))

panel_2026 <- panel_diario[
  panel_diario$fecha >= as.Date("2026-01-01") &
    panel_diario$fecha <= as.Date("2026-04-30") &
    !panel_diario$outlier, ]
panel_2026 <- panel_2026[order(panel_2026$fecha), ]

cat(sprintf("Días hábiles en 2026: %d\n\n", nrow(panel_2026)))

cat(sprintf("%-12s %8s %8s %8s %8s\n",
            "Fecha","h(días)","MAPE_RW","MAPE_AR","MAPE_VAR"))
cat(strrep("-",52),"\n")

fechas_2026     <- as.Date(c("2026-01-30","2026-02-27",
                             "2026-03-31","2026-04-21"))
resultados_2026 <- list()

for (fecha_pred in fechas_2026) {
  fecha_pred <- as.Date(fecha_pred)
  
  h_dias <- nrow(panel_diario[
    panel_diario$fecha > as.Date("2025-12-31") &
      panel_diario$fecha <= fecha_pred, ])
  if (h_dias <= 0) next
  
  idx        <- which.min(abs(panel_2026$fecha - fecha_pred))
  obs        <- as.numeric(panel_2026[idx, BONOS_GD])
  fecha_real <- panel_2026$fecha[idx]
  
  # RW: último factor conocido
  betas_RW <- c(ultimo_factor$beta0, ultimo_factor$beta1,
                ultimo_factor$beta2, TAU_FIJO)
  
  # AR: predicción dinámica desde fin de muestra
  betas_AR_v <- setNames(numeric(3), c("beta0","beta1","beta2"))
  for (nm in c("beta0","beta1","beta2")) {
    nd  <- diag_din$adf[[nm]]$necesita_diff
    ul  <- tail(factores_B[[nm]], 1)
    serie_full <- if(nd) diff(factores_B[[nm]]) else factores_B[[nm]]
    ord <- ar_modelos[[nm]]$modelo$arma[1]
    betas_AR_v[nm] <- predecir_AR_h(serie_full, ord,
                                    min(h_dias,21), nd, ul)
  }
  betas_AR <- c(betas_AR_v, TAU_FIJO)
  
  # VAR: especificación mixta
  mat_full <- cbind(
    beta0 = if(diag_din$adf$beta0$necesita_diff)
      diff(factores_B$beta0) else factores_B$beta0,
    beta1 = if(diag_din$adf$beta1$necesita_diff)
      diff(factores_B$beta1) else factores_B$beta1,
    beta2 = if(diag_din$adf$beta2$necesita_diff)
      diff(factores_B$beta2) else factores_B$beta2)
  
  var_full <- tryCatch(
    VAR(mat_full, p=diag_din$orden_var, type="const"),
    error=function(e) NULL)
  
  if (!is.null(var_full)) {
    fc_var <- predict(var_full, n.ahead=min(h_dias,21))
    h_use  <- min(h_dias, 21)
    betas_VAR <- c(
      beta0 = if(diag_din$adf$beta0$necesita_diff)
        tail(factores_B$beta0,1) +
        fc_var$fcst$beta0[h_use,"fcst"]
      else fc_var$fcst$beta0[h_use,"fcst"],
      beta1 = fc_var$fcst$beta1[h_use,"fcst"],
      beta2 = if(diag_din$adf$beta2$necesita_diff)
        tail(factores_B$beta2,1) +
        fc_var$fcst$beta2[h_use,"fcst"]
      else fc_var$fcst$beta2[h_use,"fcst"],
      TAU_FIJO)
  } else {
    betas_VAR <- betas_RW
  }
  
  calc_mape <- function(betas, fecha, obs) {
    p_teo <- sapply(BONOS_GD, function(b)
      precio_teo(betas, fecha, cf_list[[b]]))
    mean(abs((p_teo - obs) / obs * 100), na.rm=TRUE)
  }
  
  mape_RW  <- calc_mape(betas_RW,  fecha_real, obs)
  mape_AR  <- calc_mape(betas_AR,  fecha_real, obs)
  mape_VAR <- calc_mape(betas_VAR, fecha_real, obs)
  
  cat(sprintf("%-12s %8d %8.3f %8.3f %8.3f\n",
              as.character(fecha_real), h_dias,
              mape_RW, mape_AR, mape_VAR))
  
  resultados_2026[[as.character(fecha_pred)]] <- list(
    fecha    = fecha_real, h = h_dias,
    mape_RW  = mape_RW, mape_AR = mape_AR, mape_VAR = mape_VAR,
    betas_RW = betas_RW, betas_AR = betas_AR, betas_VAR = betas_VAR,
    obs      = obs)
}

# ── ETAPA 2b: Valuación de cartera hipotética ─────────────────────────────────

cat("\n=== ETAPA 2b: Cartera hipotética — NS vs mercado ===\n")
cat("Cartera inicial: USD 1.000.000 al 3 de abril de 2024\n\n")

cal_GD_sp  <- readRDS("output/calibracion_GD_sp.rds")
fecha_inicio <- as.Date("2024-04-03")

p_inicio <- as.numeric(panel_diario[
  panel_diario$fecha == fecha_inicio, BONOS_GD])
names(p_inicio) <- BONOS_GD

VN_TOTAL <- 1000000
VN_bono  <- sapply(BONOS_GD, function(b)
  (VN_TOTAL / 6) / (p_inicio[b] / 100))

cat(sprintf("Verificación valor inicial: USD %.2f\n\n",
            sum(p_inicio / 100 * VN_bono)))

cat(sprintf("%-12s %12s %12s %12s %8s\n",
            "Fecha","V_NS","V_Obs","Error_NS","Error%"))
cat(strrep("-",58),"\n")

fechas_check <- as.Date(c(
  "2024-04-03","2024-06-28","2024-09-27","2024-12-27",
  "2025-03-28","2025-06-27","2025-09-30","2025-12-29"))

error_vals     <- c()
error_pct_vals <- c()
resultados_cartera <- list()

for (fecha_check in fechas_check) {
  fecha_check <- as.Date(fecha_check)
  
  p_obs_raw <- as.numeric(panel_diario[
    panel_diario$fecha == fecha_check, BONOS_GD])
  if (all(is.na(p_obs_raw))) next
  
  res_cal <- Filter(function(r)
    !is.null(r) && !is.na(r$fecha) &&
      abs(as.numeric(r$fecha - fecha_check)) <= 3,
    cal_GD_sp)
  if (length(res_cal) == 0) next
  
  r      <- res_cal[[1]]
  fecha  <- r$fecha
  params <- c(r$beta0, r$beta1, r$beta2, TAU_FIJO)
  
  p_obs <- as.numeric(panel_diario[
    panel_diario$fecha == fecha, BONOS_GD])
  p_ns  <- sapply(BONOS_GD, function(b)
    precio_teo(params, fecha, cf_list[[b]]))
  
  v_obs   <- sum(p_obs / 100 * VN_bono, na.rm=TRUE)
  v_ns    <- sum(p_ns  / 100 * VN_bono, na.rm=TRUE)
  err     <- v_ns - v_obs
  err_pct <- err / v_obs * 100
  
  error_vals     <- c(error_vals,     err)
  error_pct_vals <- c(error_pct_vals, err_pct)
  
  resultados_cartera[[as.character(fecha)]] <- list(
    fecha=fecha, v_ns=v_ns, v_obs=v_obs,
    error=err, error_pct=err_pct)
  
  cat(sprintf("%-12s %12.0f %12.0f %12.0f %7.3f%%\n",
              as.character(fecha), v_ns, v_obs, err, err_pct))
}

cat(sprintf("\nPromedio |Error_NS|: USD %.0f (%.4f%%)\n",
            mean(abs(error_vals)),
            mean(abs(error_pct_vals))))
cat(sprintf("Sesgo promedio:      USD %+.0f\n", mean(error_vals)))

# ── Guardar ───────────────────────────────────────────────────────────────────

saveRDS(list(
  etapa1     = rmse_resultados,
  etapa2     = resultados_2026,
  etapa2b    = resultados_cartera,
  horizontes = horizontes,
  n_train    = n_train,
  n_test     = n_test,
  fecha_corte= factores_B$fecha[n_train]
), "output/prediccion_resultados.rds")

cat("\n✓ output/prediccion_resultados.rds\n")