# SCRIPT HEADER ---------------------------------------------------------------------------------------------
# This script contains helper functions required for model analysis

# Lydia Braunack-Mayer
#l.brauanckmayer@gmail.com
#lydia.braunack-mayer@swisstph.ch

# 2024


# SCRIPT ----------------------------------------------------------------------------------------------------

## Define interventions -------------------------------------------------------------------------------------

# Define effect decay for medical interventions
weibull <- function(t, init.pe, L, k) init.pe * exp(-(t/L)^k * log(2))

# Define a generic preventative intervention function with a decaying effect
intervention.pe <- function(t.start, t.end, FUN.decay, id, ...) {
  
  # Inputs:
  #  t.start: numeric vector, times (in days) when the intervention is deployed
  #  t.end: integer, end time of monitoring period
  #  FUN.decay: function, protective efficacy decay profile
  #  id: string, ID for intervention
  #  ...: additional arguments passed to FUN.decay
  
  # Outputs: data frame, protective efficacy over time
  
  # Create time index
  t <- 0:t.end
  
  # Calculate protective efficacy over time
  pe <- FUN.decay(t, ...)
  
  # Initialise outputs
  out <- data.frame("id" = id, "time" = 0:t.end, "pe" = rep(0, t.end + 1))
  
  # Add protective efficacy (assuming additional deployments boost intervention effect to the initial effect)
  for (i in 1:length(t.start)) {
    if (any(out[out$time %in% t.start[i]:t.end, "pe"] > 0)) {
      out[out$time %in% t.start[i]:t.end, "pe"] <- 0
    }
    out[out$time %in% t.start[i]:t.end, "pe"] <- out[out$time %in% t.start[i]:t.end, "pe"] + pe[1:(t.end - t.start[i] + 1)] 
  }
  
  # Correct if protective efficacy greater than 1
  out$pe[out$pe > 1] <- 1
  
  return(out)
}

# # Test function
# df <- intervention.pe(t.start = c(30, 60, 90, 120), t.end = 365, FUN.decay = weibull, id = "ID", init.pe = 0.8, L = 2, k = 2)
# plot(df$time, df$pe, type = "l")

# Define a generic intervention function (preventative or otherwise) with a time-limited, constant effect
intervention.scalar <- function(t.int, t.end, base.effect, int.effect, id) {
  
  # Inputs:
  #  t.int: numeric vector, times (in days) when the intervention is active
  #  t.end: integer, end time of monitoring period
  #  base.effect: scalar, parameter value without intervention
  #  int.effect: scalar, effect of intervention
  #  id: string, ID for intervention
  
  # Outputs: data frame, effect of intervention over time
  
  # Create time index
  t <- 0:t.end
  
  # Generate impact vector
  out <- rep(base.effect, t.end + 1)
  out[t.int] <- int.effect
  
  # Return outputs
  out <- data.frame("id" = id, "time" = 0:t.end, "effect" = out)
  
  return(out)
}

# # Test function
# df <- intervention.scalar(t.int = c(30:365), t.end = 365, id = "ID", effect = 2)
# plot(df$time, df$effect, type = "l")



## Define model equations -----------------------------------------------------------------------------------

# Define baseline SIS model equations
sis.baseline <- function(init, params) {
  
  # Check function inputs
  if (!all(c("I", "N") %in% names(init))) {
    print("init does not contain initial values for either I or N. Check function inputs.")
    stop()
  }
  if (!all(c("tmax", "beta", "delta", "gamma", "threshold") %in% names(params))) {
    print("params does not contain all required parameters. Check function inputs.")
    stop()
  }
  
  # Initialise model outputs
  df <- data.frame(t = 0:params$tmax, S = NA, I = NA, beta = NA)
  
  # Define initial conditions
  df$beta[1] <- params$beta * (1 + params$delta*cos(pi))
  df$I[1] <- init$I
  df$S[1] <- init$N - init$I
  
  # Generate model dynamics
  for (t in 1:params$tmax) {
    df$beta[t + 1] <- params$beta * (1 + params$delta*cos(2*pi*(t + 365/2)/365))
    df$I[t + 1] <- df$beta[t] / init$N * df$I[t] * (init$N - df$I[t]) + params$gamma * df$I[t]
    df$S[t + 1] <- init$N - df$I[t + 1]
    
    # If point prevalence is below the elimination threshold, set infectious population to 0
    if (df$I[t + 1]/init$N < params$threshold) df$I[t + 1] <- 0
  }
  
  # Calculate R0
  df$R0 <- df$beta / (1 - params$gamma)
  
  # Return outputs
  return(df)
  
}

# # Test function
# init.test <- list("I" = 1000, "N" = 10000)
# params.test <- list("tmax" = 365*10,
#                     "beta" = 0.1, 
#                     "delta" = 0,
#                     "gamma" = (1 - 0.30)^(1/14),
#                     "threshold" = 0)
# temp <- sis.baseline(init = init.test, params = params.test)



# Define intervention SIS model equations
sis.intervention <- function(init, params, int) {
  
  # Check function inputs
  if (!all(c("I", "N") %in% names(init))) {
    print("init does not contain initial values for either I or N. Check function inputs.")
    stop()
  }
  if (!all(c("tmax", "beta", "delta", "gamma", "threshold") %in% names(params))) {
    print("params does not contain all required parameters. Check function inputs.")
    stop()
  }
  if (!all(c("drug", "vaccine", "CM", "test") %in% names(int))) {
    print("int does not contain all required interventions Check function inputs.")
    stop()
  }
  
  # Initialise model outputs
  df <- data.frame(t = 0:params$tmax, S = NA, I = NA, beta = NA)
  
  # Define initial conditions
  df$beta[1] <- params$beta * (1 + params$delta*cos(pi))
  df$I[1] <- init$I
  df$S[1] <- init$N - init$I
  
  # Define delta as product of drug- and vaccine-based interventions
  delta <- data.frame("t" = 0:params$tmax, "delta_t" = (1 - int$drug$pe))
  delta$delta_t <- delta$delta_t * (1 - int$vaccine$pe)
  
  # Append intervention dynamics to data frame
  df$pe <- delta$delta_t
  df$drug <- int$drug$pe
  df$vaccine <- int$vaccine$pe
  df$CM <- int$CM$effect
  df$test <- int$test$effect
  
  # Add id for intervention
  df$id <- paste(int$drug$id[1], int$vaccine$id[1], int$CM$id[1], int$test$id[1], sep = ", ")
  
  # Generate model dynamics
  for (t in 1:params$tmax) {
    df$beta[t + 1] <- params$beta * (1 + params$delta*cos(2*pi*(t + 365/2)/365))
    df$I[t + 1] <- delta$delta_t[t] * df$beta[t] / init$N * df$I[t] * (init$N - df$I[t]) + df$test[t] * (max(params$gamma -  df$CM[t], 0)) * df$I[t]
    df$S[t + 1] <- init$N - df$I[t + 1]
    
    # If point prevalence is below the elimination threshold, set infectious population to 0
    if (df$I[t + 1]/init$N < params$threshold) df$I[t + 1] <- 0
  }
  
  # Calculate R0
  df$R0 <- df$beta / (1 - params$gamma)
  
  # Return outputs
  out <- list("df" = df[, c("t", "S", "I", "beta", "R0")],
              "int" = df[, c("id", "t", "pe", "drug", "vaccine", "CM", "test")],
              "init" = init,
              "params" =  params)
  return(out)
  
}

# # Test function
# init.test <- list("I" = 1000, "N" = 10000)
# params.test <- list("tmax" = 365*10,
#                     "beta" = 0.1,
#                     "delta" = 0,
#                     "gamma" = (1 - 0.30)^(1/14),
#                     "threshold" = 0)
# 
# interventions.drug <- intervention.pe(t.start = c(365 + 1),
#                                            t.end = params.test$tmax,
#                                            id = "long acting drug - pulsed",
#                                            FUN.decay = weibull,
#                                            init.pe = 1,
#                                            L = 31.1,
#                                            k = 5.4)
# 
# interventions.vaccine <- intervention.pe(t.start = 365 + 1,
#                                               t.end = params.test$tmax,
#                                               FUN.decay = weibull,
#                                               id = "vaccine - no effect",
#                                               init.pe = 0,
#                                               L = 0,
#                                               k = 0)
# 
# interventions.CM <- intervention.scalar(t.int = 0:params.test$tmax,
#                                              t.end = params.test$tmax,
#                                              base.effect = 1,
#                                              int.effect = 1,
#                                              id = "case management - no effect")
# 
# interventions.test <- intervention.scalar(t.int = 365 + 1,
#                                                t.end = params.test$tmax,
#                                                base.effect = 0,
#                                                int.effect = 0,
#                                                id = "mass test and treat - no effect")
# 
# interventions <- list(interventions.drug, interventions.vaccine, interventions.CM, interventions.test)
# names(interventions) <- c("drug", "vaccine", "CM", "test")
# 
# temp <- sis.intervention(init = init.test, params = params.test, int = interventions)



## Define postprocessing equations --------------------------------------------------------------------------

# Define postprocessing
postprocessing <- function(sis) {
  
  # Inputs:
  #  sis: object generated by sis.intervention(), results from a single model run
  
  # Outputs:
  #  out: 1 x number years tibble, prevalence and elimination outcomes from a single model run
  
  # Load required packages
  require(tidyr)
  require(dplyr)
  require(sfsmisc)
  
  # Define data
  df <- sis$df
  
  # Add year
  df <- df %>%
    mutate(year = floor(t/(365 + 1)) + 1) 
  
  # Calculate mean annual prevalence
  out <- df %>%
    group_by(year) %>%
    summarise(prev = mean(I) / sis$init$N)
  
  # Calculate prevalence reduction
  out <- out %>%
    mutate(prevRed = (prev[year == 1] - prev) / prev[year == 1])
  
  # Calculate total cases and cumulative cases per year
  temp <- df %>%
    group_by(year) %>%
    summarise(cases = sum(I)) %>%
    mutate(casesSum = cumsum(cases))
  out <- merge(out, temp, by = c("year"))
  
  # Identify if elimination was achieved in final year
  out$elimination <- as.logical(out$prev[out$year == max(out$year)] < sis$params$threshold)
  
  # Add intervention id
  out$intervention <- sis$int$id[1]
  
  # Add R0 at t = 0 under no intervention
  out$R0 <- sis$df$R0[1]
  
  # Add intervention number
  out$ID <- sis$int$idnum[1]
  
  # Transform to wide format
  out <- out %>%
    pivot_wider(id_cols = c(ID, intervention, R0, elimination),
                names_from = year,
                values_from = c(prev, prevRed, cases, casesSum)) %>%
    as.data.frame()
  
  # Return outputs
  return(out)
  
}

