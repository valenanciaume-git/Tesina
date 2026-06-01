# ==============================================================================
# 00_funciones.R
# Funciones compartidas por todos los scripts del pipeline
#
# ==============================================================================

library(tseries)

# ── 1. Funciones Nelson-Siegel ────────────────────────────────────────────────

ns_spot <- function(beta0, beta1, beta2, tau, t) {
  lt     <- t / tau
  carga1 <- ifelse(lt < 1e-6, 1.0, (1 - exp(-lt)) / lt)
  carga2 <- ifelse(lt < 1e-6, 0.0, carga1 - exp(-lt))
  beta0 + beta1 * carga1 + beta2 * carga2
}

ns_forward <- function(beta0, beta1, beta2, tau, t) {
  lt <- t / tau
  beta0 + beta1 * exp(-lt) + beta2 * (lt / tau) * exp(-lt)
}

# Fecha a partir de la cual cambia el régimen de liquidación
FECHA_CAMBIO_LIQUIDACION <- as.Date("2024-06-03")

#' Precio teórico descontando desde la fecha de LIQUIDACIÓN
#' Antes del 3-jun-2024: T+2 (liquida 2 días hábiles después)
#' Desde el 3-jun-2024:  T+1 (liquida 1 día hábil después)
#' Se calcula el cupón corrido hasta la fecha de liquidación,
#' por lo que el precio observado corresponde a esa fecha.
precio_teo <- function(params, val_date, cf) {
  beta0 <- params[1]; beta1 <- params[2]
  beta2 <- params[3]; tau   <- params[4]
  
  val_date <- as.Date(val_date)
  
  # Determinar días hábiles a agregar según régimen de liquidación
  dias_liq <- ifelse(val_date < FECHA_CAMBIO_LIQUIDACION, 2, 1)
  
  # Fecha de liquidación: avanzar días hábiles
  # Aproximación: sumar días calendario (fines de semana raramente caen
  # exactamente en la liquidación, el error es < 1 día sobre 365)
  liq_date <- val_date + dias_liq
  
  # Si cae en fin de semana, correr al lunes
  dow <- as.integer(format(liq_date, "%u"))  # 6=sábado, 7=domingo
  if (dow == 6) liq_date <- liq_date + 2
  if (dow == 7) liq_date <- liq_date + 1
  
  cf_fut <- cf[cf$fecha > liq_date, ]
  if (nrow(cf_fut) == 0) return(NA_real_)
  
  t <- as.numeric(cf_fut$fecha - liq_date) / 365
  s <- ns_spot(beta0, beta1, beta2, tau, t)
  sum(cf_fut$ff * exp(-s * t))
}

construir_curva <- function(res, plazos = seq(0.25, 30, by = 0.25)) {
  if (is.null(res)) return(NULL)
  tau <- if (!is.null(res$tau)) res$tau else 1/res$lambda
  data.frame(
    fecha   = res$fecha,
    plazo   = plazos,
    spot    = ns_spot(res$beta0, res$beta1, res$beta2, tau, plazos),
    forward = ns_forward(res$beta0, res$beta1, res$beta2, tau, plazos)
  )
}

# ── 2. Funciones de cash flows ────────────────────────────────────────────────

flujos_futuros <- function(bono, val_date, cf_list) {
  val_date <- as.Date(val_date)
  cf <- cf_list[[bono]]
  if (is.null(cf)) return(NULL)
  cf_fut <- cf[cf$fecha > val_date, ]
  if (nrow(cf_fut) == 0) return(NULL)
  cf_fut$t <- as.numeric(cf_fut$fecha - val_date) / 365
  cf_fut
}

get_ancla <- function(fecha, ancla_df) {
  if (is.null(ancla_df)) return(NULL)
  fecha <- as.Date(fecha)
  idx   <- which(ancla_df$fecha <= fecha & ancla_df$tasa > 0.0001)
  if (length(idx) == 0) return(NULL)
  ancla_df$tasa[max(idx)]
}

# ── 3. Funciones de calibración ───────────────────────────────────────────────

duration_mod <- function(bono, val_date, cf_list, ytm_bono) {
  cf <- flujos_futuros(bono, val_date, cf_list)
  if (is.null(cf)) return(NA_real_)
  pv  <- cf$ff * exp(-ytm_bono * cf$t)
  p   <- sum(pv)
  if (p <= 0) return(NA_real_)
  sum(cf$t * pv) / p
}

calcular_ponderadores <- function(bonos, val_date, cf_list,
                                  tirs_obs) {
  durs <- sapply(bonos, function(b) {
    tir_col  <- paste0("tir_", b)
    ytm_comp <- tirs_obs[[tir_col]]
    
    # Verificar que la TIR sea válida
    if (is.null(ytm_comp) || is.na(ytm_comp) || ytm_comp <= 0.001)
      return(NA_real_)
    
    # Convertir TIR compuesta a continua: r_cont = ln(1 + r_comp)
    ytm_cont <- log(1 + ytm_comp)
    
    duration_mod(b, val_date, cf_list, ytm_bono = ytm_cont)
  })
  
  inv <- ifelse(is.na(durs) | durs <= 0, NA, 1 / durs)
  tot <- sum(inv, na.rm = TRUE)
  if (tot <= 0) return(rep(1/length(bonos), length(bonos)))
  inv / tot
}

loss_precio_ancla <- function(params, val_date, bonos, precios_obs,
                              cf_list, ancla = NULL,
                              t_ancla = 30/365,
                              ponderadores = NULL) {
  beta0 <- params[1]; beta1 <- params[2]
  beta2 <- params[3]; tau   <- params[4]
  
  if (tau <= 0 || tau > 50) return(1e10)
  
  sse <- 0; n <- 0
  
  for (i in seq_along(bonos)) {
    if (is.na(precios_obs[i])) next
    p <- precio_teo(params, val_date, cf_list[[bonos[i]]])
    if (!is.na(p)) {
      w   <- if (!is.null(ponderadores) && !is.na(ponderadores[i]))
        ponderadores[i] else 1
      sse <- sse + w * (precios_obs[i] - p)^2
      n   <- n + 1
    }
  }
  if (n == 0) return(1e10)
  
  if (!is.null(ancla) && !is.na(ancla) && ancla >= 0.001) {
    s_corto <- ns_spot(beta0, beta1, beta2, tau, t_ancla)
    sse     <- sse + ((s_corto - ancla) / ancla)^2
  }
  sse
}

# ── 4. Series de tiempo ───────────────────────────────────────────────────────

test_raiz_unitaria <- function(x, alpha = 0.05) {
  x <- na.omit(x)
  res <- tryCatch(adf.test(x, alternative="stationary"), error=function(e) NULL)
  if (is.null(res))
    return(list(p_valor_niveles=NA, estadistico=NA,
                necesita_diff=FALSE, orden="I(0)", p_valor_diff=NA))
  nd <- res$p.value >= alpha
  pd <- NA
  if (nd) {
    rd <- tryCatch(adf.test(diff(x), alternative="stationary"),
                   error=function(e) NULL)
    if (!is.null(rd)) pd <- rd$p.value
  }
  list(p_valor_niveles=res$p.value, estadistico=as.numeric(res$statistic),
       p_valor_diff=pd, necesita_diff=nd,
       orden=ifelse(nd,"I(1)","I(0)"))
}

preparar_serie <- function(x, necesita_diff) {
  x <- na.omit(x)
  list(serie=ts(if(necesita_diff) diff(x) else x, frequency=252),
       ultimo_nivel=tail(x,1), necesita_diff=necesita_diff)
}

# ── 5. Constantes ─────────────────────────────────────────────────────────────

# τ calibrado 
# mediana de calibración libre sobre período estable may-dic 2025
# Resultado: mediana=1.13, media=1.14, SD=0.28 (n=33 obs diarias)
TAU_FIJO    <- 1.13
LAMBDA_FIJO <- 1 / TAU_FIJO

# Serie GD: bonos ley Nueva York
BONOS_GD <- c("GD29","GD30","GD35","GD38","GD41","GD46")

# Serie AL: bonos ley Argentina

BONOS_AL <- c("AL29","AL30","AE38","AL35","AL41") 

# Grilla de puntos de partida — Alberti y Zincenko (2019)
GRILLA_INICIO <- expand.grid(
  beta0 = c(0.07, 0.09, 0.10, 0.11),
  beta1 = c(-0.065, -0.055, -0.045, -0.025),
  beta2 = c(0.10, 0.20, 0.30),
  tau   = TAU_FIJO
)
# 4 × 4 × 3 × 1 = 48 puntos de partida

BOUNDS_LOWER <- c(0.00, -1.50, -1.00, TAU_FIJO * 0.999)
BOUNDS_UPPER <- c(1.50,  0.50,  1.00, TAU_FIJO * 1.001)

# Definición de subperíodos
SUBPERIODOS <- list(
  A = list(
    nombre = "Distressed debt",
    inicio = as.Date("2021-01-01"),
    fin    = as.Date("2024-03-31"),
    desc   = "Yields 30-50%, dos regímenes de pricing"
  ),
  B = list(
    nombre = "Normalización",
    inicio = as.Date("2024-04-01"),
    fin    = as.Date("2025-12-31"),
    desc   = "Curva suave invertida, yields 7-18%"
  )
)

# Fechas ex-cupón conocidas (T+1 desde jun-2024)
# Días donde el modelo falla por caída diferencial de cupones
FECHAS_EXCUPON <- as.Date(c(
  "2021-01-07","2021-07-07",
  "2022-01-07","2022-07-07",
  "2023-01-06","2023-07-07",
  "2024-01-05","2024-07-05",  # T+2
  "2025-01-08","2025-07-08"   # T+1 desde jun-2024
))