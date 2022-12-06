#simulate comparison audit betting 
library(tidyverse)
library(latex2exp)
library(mvtnorm)

######## helper functions #######

shuffle <- function(x){sample(x, size = length(x), replace = FALSE)}

log_expected_value <- function(lambda, a, p_1, p_2){
  log(1 + lambda * (a - 1/2)) * (1 - p_1 - p_2) + log(1 + lambda * (a/2 - 1/2)) * p_1 + log(1 - lambda / 2) * p_2
}

optimal_lambda <- function(a, p_1, p_2){
  temp_log_expected_value <- function(lambda){
    log_expected_value(lambda, a = a, p_1 = p_1, p_2 = p_2)
  }
  derivative <- function(lambda){
    (a - 1/2) * (1 - p_1 - p_2) / (1 + lambda * (a - 1/2)) + (a - 1) * p_1 / (2 - lambda * (1 - a)) + p_2 / (2 - lambda)
  }
  solution <- optimize(temp_log_expected_value, interval = c(0,2), maximum = TRUE)
  solution$maximum
}

#from https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
estimate_beta_params <- function(mu, var){
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  c("alpha" = alpha, "beta" = beta)
}

######### main simulation function ################
simulate_audit <- function(pop, a, risk_limit = 0.05, strategy = "fixed", lambda = NULL, pars = NULL){
  x <- sample(pop, size = length(pop), replace = T)
  
  N <- length(x)
  m <- 1/2
  
  if(strategy == "fixed" & !is.null(lambda)){
    mart <- cumprod(1 + lambda * (x - m))
  } else if(strategy == "grid"){
    lambda <- seq(0, 1/m, length.out = 100)
    mart <- colMeans(t(apply(1 + lambda %o% (x - m), 1, cumprod)))
  } else if(strategy == "adaptive"){
    if(is.null(pars)){stop("Need to specify pars$prior_p_k, pars$d_k, pars$eps_k for k in 1,2")}
    lag_p_hat_1 <- c(0, lag(cummean(x == a/2))[-1])
    lag_p_hat_2 <- c(0, lag(cummean(x == 0))[-1])
    tilde_p_1 <- pmax((pars$d_1 * pars$prior_p_1 + (1:N) * lag_p_hat_1) / (pars$d_1 + 1:N - 1), pars$eps_1)
    tilde_p_2 <- pmax((pars$d_2 * pars$prior_p_2 + (1:N) * lag_p_hat_2) / (pars$d_2 + 1:N - 1), pars$eps_2)
    tilde_lambda <- rep(NA, N)
    for(i in 1:N){
      tilde_lambda[i] <- optimal_lambda(a = a, p_1 = tilde_p_1[i], p_2 = tilde_p_2[i])
    }
    mart <- cumprod(1 + tilde_lambda * (x - m))
  } else if(strategy == "mixture"){
    grid_1 <- seq(1e-6, 2 - 1/a, length.out = 10)
    grid_2 <- seq(1e-6, 1 - 1/(2*a), length.out = 10)
    mu <- c(pars$mean_p_1, pars$mean_p_2)
    Sigma <- rbind(c(pars$sd_p_1^2, pars$cor * pars$sd_p_1 * pars$sd_p_2), 
                      c(pars$cor * pars$sd_p_1 * pars$sd_p_2, pars$sd_p_2^2))
    bivar_grid <- expand.grid(grid_1, grid_2) %>%
      rename("p_1" = "Var1", "p_2" = "Var2") %>% 
      filter(a * p_2 + (a/2) * p_1 < a - 1/2) %>%
      mutate(weights = dmvnorm(cbind(p_1, p_2), mean = mu, sigma = Sigma)) %>%
      mutate(weights = weights/sum(weights))
    B <- nrow(bivar_grid)
    lambda <- rep(NA, B)
    for(b in 1:B){
      lambda[b] <- optimal_lambda(a = a, p_1 = bivar_grid$p_1[b], p_2 = bivar_grid$p_2[b])
    }
    marts <- t(apply(1 + lambda %o% (x - m), 1, cumprod))
    mart <- colSums(bivar_grid$weights * marts)
  } else{
    stop("Input valid strategy for betting")
  }
  
  
  if(any(mart > 1/risk_limit)){
    stopping_time <- min(which(mart > 1/risk_limit))
  } else{
    stopping_time <- length(mart)
  }
  stopping_time
}


############### Table 1: oracle methods ##############
#simulated stopping times for optimal lambda under a range of diluted margins and overstatement rates.

# v <- c(.05, .1, .2)
# a <- 1/(2 - v)
# risk_limit <- 0.05
# #p <- c(.997, .999, 1)
# 
# p <- c(.985, .99, .995, .999, 1)
# n_sims <- 400
# N <- 10000
# 
# stopping_times_optimal <- array(NA, dim = c(n_sims, length(v), length(p)), dimnames = list(1:n_sims, v, p))
# stopping_times_apkelly <- array(NA, dim = c(n_sims, length(v), length(p)), dimnames = list(1:n_sims, v, p))
# stopping_times_gkelly <- array(NA, dim = c(n_sims, length(v), length(p)), dimnames = list(1:n_sims, v, p))
# for(i in 1:length(v)){
#   for(j in 1:length(p)){
#     #reported_votes <- c(rep(0, N * (.5 - v[i]/2)), rep(1, N * (.5 + v[i]/2)))
#     #true_votes <- c(rep(0, N * (.5 - v[i]/2) + N * (1-p[j])), rep(1, N * (.5 + v[i]/2) - N * (1-p[j])))
#     pop <- c(rep(0, N * (1-p[j])), rep(a[i], N*p[j]))
#     
#     lambda_optimal <- max(0, (2 - 4 * a[i] * p[j]) / (1 - 2 * a[i]))
#     eta_apkelly <- (1 - mean(reported_votes - true_votes)) / (2 - (2 * mean(reported_votes) - 1))
#     #using transformation in ALPHA
#     lambda_apkelly <- 4 * eta_apkelly - 2
#     print(lambda_optimal)
#     
#     stopping_times_optimal[,i,j] <- replicate(n_sims, simulate_audit(pop, a =  a[i], lambda = lambda_optimal)) 
#     stopping_times_apkelly[,i,j] <- replicate(n_sims, simulate_audit(pop, a =  a[i], lambda = lambda_apkelly))
#     stopping_times_gkelly[,i,j] <- replicate(n_sims, simulate_audit(pop, a =  a[i], strategy = "grid"))
#   }
# }

# stopping_time_df_optimal <- stopping_times_optimal %>%
#   reshape2::melt(varnames = c("sim","v","p")) %>%
#   mutate(lambda = "optimal")
# stopping_time_df_apkelly <- stopping_times_apkelly %>%
#   reshape2::melt(varnames = c("sim","v","p")) %>%
#   mutate(lambda = "apkelly")
# stopping_time_df_gkelly <- stopping_times_gkelly %>%
#   reshape2::melt(varnames = c("sim","v","p")) %>%
#   mutate(lambda = "gkelly")
# stopping_time_df <- stopping_time_df_optimal %>%
#   bind_rows(stopping_time_df_apkelly) %>%
#   bind_rows(stopping_time_df_gkelly) %>%
#   as_tibble()
# save(stopping_time_df = "table_1_simulations")



##################### Table 2: robust methods #######################
v <- c(.05, .1)
a <- 1/(2 - v)
risk_limit <- 0.05
p_1 <- c(.01, .001) #true 1-vote overstatement rate
p_2 <- c(.01, .001, .0001) #true 2-vote overstatement rate
n_sims <- 500
N <- 20000

prior_p_1 <- c(.01, .001)
prior_p_2 <- c(.001, .0001)
sd_p_1 <- .005
sd_p_2 <- .0025
d_1 <- 100
d_2 <- 1000
eps_1 <- .00001
eps_2 <- .00001

sim_frame <- expand.grid("a" = a, "p_1" = p_1, "p_2" = p_2, "prior_p_1" = prior_p_1, "prior_p_2" = prior_p_2, "sd_p_1" = sd_p_1, "sd_p_2" = sd_p_2, "d_1" = d_1, "d_2" = d_2, "eps_1" = eps_1, "eps_2" = eps_2) %>%
  mutate(pop_mean = a * (1 - p_1/2 - p_2), diluted_margin = 2 - 1/a)
stopping_times_adaptive <- matrix(NA, nrow = nrow(sim_frame), ncol = n_sims)
stopping_times_apriori <- matrix(NA, nrow = nrow(sim_frame), ncol = n_sims)
stopping_times_oracle <- matrix(NA, nrow = nrow(sim_frame), ncol = n_sims)
stopping_times_mixture <- matrix(NA, nrow = nrow(sim_frame), ncol = n_sims)


for(i in 1:nrow(sim_frame)){
  print(i)
  pop <- c(
    rep(0, N * sim_frame$p_2[i]),
    rep(sim_frame$a[i]/2, N * sim_frame$p_1[i]),
    rep(sim_frame$a[i], N * (1 - sim_frame$p_1[i] - sim_frame$p_2[i]))
  )
  pars_adaptive <- list(
    "prior_p_1" = sim_frame$prior_p_1[i], 
    "prior_p_2" = sim_frame$prior_p_2[i], 
    "d_1" = sim_frame$d_1[i], 
    "d_2" = sim_frame$d_2[i], 
    "eps_1" = sim_frame$eps_1[i], 
    "eps_2" = sim_frame$eps_2[i]
    )
  pars_mixture <- list(
    "mean_p_1" = sim_frame$prior_p_1[i],
    "mean_p_2" = sim_frame$prior_p_2[i],
    "sd_p_1" = sim_frame$sd_p_1[i],
    "sd_p_2" = sim_frame$sd_p_2[i],
    "cor" = 0.25
  )
  apriori_lambda <- min(2, max(0, optimal_lambda(a = sim_frame$a[i], p_1 = sim_frame$prior_p_1[i], p_2 = sim_frame$prior_p_2[i])))
  oracle_lambda <- optimal_lambda(a = sim_frame$a[i], p_1 = sim_frame$p_1[i], p_2 = sim_frame$p_2[i])
  
  stopping_times_adaptive[i,] <- replicate(n_sims, simulate_audit(pop, a =  sim_frame$a[i], strategy = "adaptive", pars = pars_adaptive))
  stopping_times_mixture[i,] <- replicate(n_sims, simulate_audit(pop, a =  sim_frame$a[i], strategy = "mixture", pars = pars_mixture))
  stopping_times_apriori[i,] <- replicate(n_sims, simulate_audit(pop, a =  sim_frame$a[i], strategy = "fixed", lambda = apriori_lambda))
  stopping_times_oracle[i,] <- replicate(n_sims, simulate_audit(pop, a =  sim_frame$a[i], strategy = "fixed", lambda = oracle_lambda))
}
#stopping_times_list <- list(stopping_times_adaptive, stopping_times_apriori, stopping_times_oracle)

stopping_time_df_adaptive <- stopping_times_adaptive %>%
  as_tibble() %>%
  bind_cols(sim_frame) %>%
  pivot_longer(cols = starts_with("V", ignore.case = FALSE), names_to = "rep", values_to = "stopping_time", names_prefix = "V") %>%
  mutate(lambda = "adaptive")
stopping_time_df_mixture <- stopping_times_mixture %>%
  as_tibble() %>%
  bind_cols(sim_frame) %>%
  pivot_longer(cols = starts_with("V", ignore.case = FALSE), names_to = "rep", values_to = "stopping_time", names_prefix = "V") %>%
  mutate(lambda = "mixture")
stopping_time_df_apriori <- stopping_times_apriori %>%
  as_tibble() %>%
  bind_cols(sim_frame) %>%
  pivot_longer(cols = starts_with("V", ignore.case = FALSE), names_to = "rep", values_to = "stopping_time", names_prefix = "V") %>%
  mutate(lambda = "apriori")
stopping_time_df_oracle <- stopping_times_oracle %>%
  as_tibble() %>%
  bind_cols(sim_frame) %>%
  pivot_longer(cols = starts_with("V", ignore.case = FALSE), names_to = "rep", values_to = "stopping_time", names_prefix = "V") %>%
  mutate(lambda = "oracle")

stopping_times_df <- stopping_time_df_adaptive %>%
  bind_rows(stopping_time_df_mixture) %>%
  bind_rows(stopping_time_df_apriori) %>%
  bind_rows(stopping_time_df_oracle) 

save(stopping_times_df, file = "table_2_simulations")



