# ==============================================================================
# 05_dinamica_var.R
# Análisis de la dinámica temporal de los factores NS — Subperíodo B
#
# Contenido:
#   1. Tests ADF de raíz unitaria por factor
#   2. Modelos AR univariados por factor (orden por AIC)
#   3. Modelo VAR conjunto (orden por BIC)
#   4. Tests de causalidad de Granger — par a par (lmtest::grangertest)
#      NOTA: se usa grangertest() y NO vars::causality(), que testea
#      causalidad conjunta (X causa al sistema {Y1,Y2}) y no par a par,
#   5. Funciones impulso-respuesta (IRF)
#   6. Descomposición de varianza (FEVD)
#
# Requiere: output/factores_GD.rds
# Output:   output/ar_modelos.rds
#           output/var_modelo.rds
#           output/diagnostico_dinamica.rds
# ==============================================================================

set.seed(42)

library(dplyr)
library(tseries)
library(vars)
library(forecast)
library(lmtest)

source("scripts/00_funciones.R")

factores_full <- readRDS("output/factores_GD.rds")

factores <- factores_full %>%
  filter(fecha >= SUBPERIODOS$B$inicio,
         fecha <= SUBPERIODOS$B$fin,
         !is.na(beta0), beta0 > 0.01) %>%
  arrange(fecha)

cat(sprintf("=== Dinámica de factores — Subperíodo B ===\n"))
cat(sprintf("Período: %s a %s\n", min(factores$fecha), max(factores$fecha)))
cat(sprintf("Observaciones: %d días hábiles\n\n", nrow(factores)))

b0_ts <- ts(factores$beta0, frequency = 252)
b1_ts <- ts(factores$beta1, frequency = 252)
b2_ts <- ts(factores$beta2, frequency = 252)

# ── 1. Tests ADF ──────────────────────────────────────────────────────────────
# H₀: tiene raíz unitaria (no estacionaria)
# Si p < 0.05 → estacionaria → modelar en niveles
# Si p > 0.05 → I(1) → diferenciar antes de modelar

cat("=== 1. Tests ADF de raíz unitaria ===\n")
cat("H₀: la serie tiene raíz unitaria\n")
cat("p < 0.05 → estacionaria → usar niveles\n")
cat("p > 0.05 → I(1) → diferenciar\n\n")
cat(sprintf("%-8s %12s %10s %10s %10s\n",
            "Factor","Estadístico","p-valor","Orden","Acción"))
cat(strrep("-",55),"\n")

adf_resultados <- list()
series_modelo  <- list()

for (nm in c("beta0","beta1","beta2")) {
  x   <- factores[[nm]]
  res <- test_raiz_unitaria(x)
  adf_resultados[[nm]] <- res
  accion <- ifelse(res$necesita_diff, "Diferenciar", "Niveles")
  cat(sprintf("%-8s %12.4f %10.4f %10s %10s\n",
              nm, res$estadistico, res$p_valor_niveles,
              res$orden, accion))
  series_modelo[[nm]] <- preparar_serie(x, res$necesita_diff)
}

# ── 2. Modelos AR univariados ─────────────────────────────────────────────────

cat("\n=== 2. Modelos AR univariados (orden por AIC) ===\n")
cat("AR(p): el valor de hoy depende de los p valores anteriores\n\n")

ar_modelos <- list()

for (nm in c("beta0","beta1","beta2")) {
  serie <- series_modelo[[nm]]$serie
  mod   <- auto.arima(serie, d=0, max.p=10, max.q=0,
                      ic="aic", stepwise=FALSE, approximation=FALSE)
  lb    <- Box.test(residuals(mod), lag=10, type="Ljung-Box")
  
  ar_modelos[[nm]] <- list(
    modelo        = mod,
    necesita_diff = series_modelo[[nm]]$necesita_diff,
    ultimo_nivel  = series_modelo[[nm]]$ultimo_nivel
  )
  
  cat(sprintf("Factor %s: AR(%d)  AIC=%.2f  Ljung-Box p=%.4f %s\n",
              nm, mod$arma[1], mod$aic, lb$p.value,
              ifelse(lb$p.value > 0.05, "✓ residuos OK",
                     "⚠ autocorrelación en residuos")))
  cat(sprintf("  Coeficientes φ: "))
  cat(round(coef(mod), 4))
  cat("\n")
}

# ── 3. Modelo VAR ─────────────────────────────────────────────────────────────
# Especificación mixta: Δβ₀, β₁ en niveles, Δβ₂ — según ADF
# Orden por BIC (más parsimonioso que AIC)

cat("\n=== 3. Modelo VAR conjunto (orden por BIC) ===\n")
cat("Especificación mixta: Δβ₀, β₁ en niveles, Δβ₂\n")
cat("VAR(p): los 3 factores modelados simultáneamente\n\n")

# CORRECCIÓN DE ALINEAMIENTO TEMPORAL:
# diff(x) tiene n-1 observaciones vs x que tiene n.
# Si beta1 es I(0) (en niveles) y beta0/beta2 son I(1) (en diferencias),
# hay que recortar beta1 al mismo largo que diff(beta0) y diff(beta2).
# Sin este recorte R hace recycling silencioso: beta1[1] (fecha t=1)
# queda emparejado con Δbeta0[1] (que es el cambio de t=1 a t=2),
# introduciendo un desalineamiento de un período.
# Solución: usar tail() para tomar las últimas n_min observaciones.

n_min <- min(
  length(if(adf_resultados$beta0$necesita_diff)
    diff(factores$beta0) else factores$beta0),
  length(if(adf_resultados$beta1$necesita_diff)
    diff(factores$beta1) else factores$beta1),
  length(if(adf_resultados$beta2$necesita_diff)
    diff(factores$beta2) else factores$beta2))

mat_var <- cbind(
  beta0 = tail(
    if(adf_resultados$beta0$necesita_diff)
      diff(factores$beta0) else factores$beta0,
    n_min),
  beta1 = tail(
    if(adf_resultados$beta1$necesita_diff)
      diff(factores$beta1) else factores$beta1,
    n_min),
  beta2 = tail(
    if(adf_resultados$beta2$necesita_diff)
      diff(factores$beta2) else factores$beta2,
    n_min))

cat(sprintf("n_min (observaciones alineadas): %d\n", n_min))

seleccion <- VARselect(mat_var, lag.max=10, type="const")
orden_bic  <- seleccion$selection["SC(n)"]
cat(sprintf("Orden seleccionado por BIC: VAR(%d)\n\n", orden_bic))

var_mod <- VAR(mat_var, p=orden_bic, type="const")

cat("R² por ecuación:\n")
for (nm in c("beta0","beta1","beta2")) {
  r2 <- summary(var_mod$varresult[[nm]])$r.squared
  cat(sprintf("  %-8s R² = %.4f\n", nm, r2))
}

# ── 4. Causalidad de Granger — par a par ─────────────────────────────────────
# MÉTODO: lmtest::grangertest() — testea cada par (X→Y) de forma independiente
#
#
# grangertest(Y ~ X, order=p) estima dos modelos:
#   Restringido: Y_t = f(Y_{t-1},...,Y_{t-p})
#   No restringido: Y_t = f(Y_{t-1},...,Y_{t-p}, X_{t-1},...,X_{t-p})
#   F testea si los rezagos de X agregan poder predictivo sobre Y
#
# H₀: X no causa Granger a Y
# Si p < 0.05 → X tiene poder predictivo sobre Y en el margen

cat("\n=== 4. Causalidad de Granger — par a par ===\n")
cat("Método: lmtest::grangertest() | Orden:", orden_bic, "\n")
cat("Series: Δβ₀, β₁ en niveles, Δβ₂ (según ADF)\n")
cat("H₀: X no causa Granger a Y\n\n")
cat(sprintf("%-20s %10s %10s %12s\n",
            "Par (X→Y)","F-stat","p-valor","Conclusión"))
cat(strrep("-",55),"\n")

# Alinear longitudes — β₁ en niveles tiene n obs, Δβ₀ y Δβ₂ tienen n-1
n_comun <- min(sapply(c("beta0","beta1","beta2"), function(nm)
  length(series_modelo[[nm]]$serie)))

series_gt <- lapply(c("beta0","beta1","beta2"), function(nm)
  tail(as.numeric(series_modelo[[nm]]$serie), n_comun))
names(series_gt) <- c("beta0","beta1","beta2")

cat(sprintf("(n alineado = %d observaciones)\n\n", n_comun))

pares <- list(
  c("beta1","beta0"), c("beta2","beta0"),
  c("beta0","beta1"), c("beta2","beta1"),
  c("beta0","beta2"), c("beta1","beta2"))

granger_resultados <- list()
n_signif <- 0

for (par in pares) {
  causa  <- par[1]; efecto <- par[2]
  label  <- sprintf("%s → %s", causa, efecto)
  
  res <- tryCatch(
    grangertest(series_gt[[efecto]] ~ series_gt[[causa]],
                order = orden_bic),
    error = function(e) NULL)
  
  if (!is.null(res)) {
    fstat  <- res$F[2]
    pval   <- res$`Pr(>F)`[2]
    signif <- pval < 0.05
    if (signif) n_signif <- n_signif + 1
    
    cat(sprintf("%-20s %10.4f %10.4f %12s\n",
                label, fstat, pval,
                ifelse(signif, "Sí causa ✓", "No causa")))
    
    granger_resultados[[label]] <- list(
      causa         = causa,
      efecto        = efecto,
      f_stat        = fstat,
      p_valor       = pval,
      significativo = signif)
  }
}

cat(sprintf("\nPares significativos: %d de 6\n", n_signif))

if (n_signif == 0) {
  cat("→ Factores independientes entre sí\n")
  cat("→ AR univariados son suficientes\n")
} else if (n_signif <= 3) {
  cat("→ Interdependencias selectivas entre factores\n")
  cat("→ VAR captura relaciones que los AR univariados no pueden\n")
  cat("→ Pero R² bajos sugieren magnitud económica limitada\n")
  cat("→ AR univariados siguen siendo una alternativa parsimoniosa válida\n")
} else {
  cat("→ Interdependencia generalizada entre factores\n")
  cat("→ VAR claramente preferible a AR univariados\n")
}

# ── 5. IRF ────────────────────────────────────────────────────────────────────

irf_resultados <- NULL
if (n_signif > 0) {
  cat("\n=== 5. Funciones impulso-respuesta (21 días) ===\n")
  irf_resultados <- irf(var_mod, n.ahead=21, boot=TRUE,
                        runs=200, ci=0.90)
  print(irf_resultados)
}

# ── 6. FEVD ───────────────────────────────────────────────────────────────────

cat("\n=== 6. Descomposición de varianza (FEVD a 21 días) ===\n")
fevd_resultados <- fevd(var_mod, n.ahead=21)
print(fevd_resultados)

# ── 7. Guardar ────────────────────────────────────────────────────────────────

if (!dir.exists("output")) dir.create("output")
saveRDS(ar_modelos, "output/ar_modelos.rds")
saveRDS(var_mod,    "output/var_modelo.rds")
saveRDS(list(
  adf              = adf_resultados,
  granger          = granger_resultados,
  irf              = irf_resultados,
  fevd             = fevd_resultados,
  n_obs            = nrow(factores),
  periodo          = c(min(factores$fecha), max(factores$fecha)),
  orden_var        = orden_bic,
  n_granger_signif = n_signif
), "output/diagnostico_dinamica.rds")

cat("\n✓ output/ar_modelos.rds\n")
cat("✓ output/var_modelo.rds\n")
cat("✓ output/diagnostico_dinamica.rds\n")