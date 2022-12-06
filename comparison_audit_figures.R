library(tidyverse)
library(latex2exp)
library(gridExtra)

#comparison audit sample sizes
#log scale grid from https://stackoverflow.com/questions/30179442/plotting-minor-breaks-on-a-log-scale-with-ggplot
log10_minor_break <- function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}
geomean <- function(x){
  exp(mean(log(x)))
}

#ballot-polling and comparison audits
N_w <- 500
N_l <- 400
N_u <- 100
N <- N_w + N_l + N_u
v <- (N_w - N_l) / N
assorter_frame <- data.frame(
  assorters = c(rep(1, N_w), rep(0, N_l), rep(1/2, N_u), rep(1 / (2 - v), N)),
  type = c(rep("Polling", N), rep("Comparison", N))
)
ggplot(assorter_frame, aes(assorters)) +
  geom_histogram(binwidth = 0.02) +
  geom_vline(xintercept = 0.5, linetype = 'dashed') +
  facet_grid(~ type) +
  xlab("Population values") +
  ylab("Count") +
  theme_bw() +
  theme(text = element_text(size = 18), axis.text = element_text(size = 14))


#assuming plurality contest, so u = 1
#diluted margin
v <- seq(0.0001, 0.075, length.out = 2000)
a <- 1/(2 - v)
alpha <- c(.001, .01, .05, .1)

n_frame <- expand.grid(a = a, alpha = alpha) %>%
  mutate(v = 2 - 1/a) %>%
  mutate(n = log(1/alpha) / log(2*a)) %>%
  mutate(alpha_pct = factor(paste(100*alpha, "%", sep = ""), levels = c("0.1%", "1%", "5%", "10%")))

ggplot(n_frame, aes(x = v, y = n, color = alpha_pct)) +
  geom_line(size = 1.5) +
  coord_cartesian(ylim = c(50,10000), xlim = c(0,0.05)) +
  scale_y_log10(breaks = 10^(2:4), minor_breaks = log10_minor_break()) +
  scale_x_continuous(labels=scales::percent) +
  xlab("Diluted margin") +
  ylab("Sample size") +
  guides(color=guide_legend(title="Risk limit")) +
  theme_bw() +
  theme(text = element_text(size = 18), axis.text = element_text(size = 14)) +
  annotation_logticks()


# plot lambda for a grid of diluted margins and error rates
p <- seq(.9,1,length.out = 2000)
lambda_frame <- expand.grid(a=a, p=p) %>%
  filter(a * p > 1/2) %>%
  mutate(v = 2 - 1/a) %>%
  mutate(lambda = (2 - 4 * a * p) / (1 - 2 * a)) 

ggplot(lambda_frame, aes(x = v, y = 1-p, fill = lambda)) +
  geom_raster() +
  scale_x_continuous(labels=scales::percent) +
  scale_y_continuous(labels=scales::percent) +
  xlab("Diluted margin") +
  ylab(TeX("Two-vote overstatement rate ($p_2$)")) +
  theme_bw() +
  scale_fill_viridis_c(begin = 0, end = 1, guide_colorbar(title = "Lambda")) +
  theme(text = element_text(size = 18), axis.text = element_text(size = 14)) 


############# create table 1 etc ##########
#table values
load("table_1_simulations")

stopping_time_table <- stopping_time_df %>%
  group_by(lambda, v, p) %>%
  summarize(mean_stopping = round(mean(value)), percentile_stopping = round(quantile(value, .9))) %>%
  mutate(percentile_stopping = paste("(", percentile_stopping, ")", sep = "")) %>%
  unite("stopping", c("mean_stopping", "percentile_stopping"), sep = " ") %>%
  pivot_wider(names_from = "lambda", values_from = "stopping")

print(xtable::xtable(stopping_time_table), include.rownames = FALSE, digits = 4)

ggplot(stopping_time_df %>% filter(v %in% c(0.05)), aes(value, color = lambda)) +
  stat_ecdf(geom = "step") +
  coord_cartesian(xlim = c(0,10000)) +
  facet_grid(. ~ p)

#average across scenarios
average_stopping_times <- stopping_time_df %>%
  group_by(lambda, v, p) %>%
  summarize(mean_stopping = round(mean(value))) %>%
  ungroup() %>%
  pivot_wider(names_from = "lambda", values_from = mean_stopping) %>%
  mutate(ratio_apkelly = optimal / apkelly, ratio_gkelly = optimal / gkelly) %>% 
  summarize(geomean_apkelly = exp(mean(log(ratio_apkelly))), geomean_gkelly = exp(mean(log(ratio_gkelly))))


############# create table 2 etc ##########
load("table_2_simulations")

ggplot(stopping_times_df %>% filter(prior_p_1 == .001, prior_p_2 == .0001, round(diluted_margin, 3) == 0.05), 
       aes(stopping_time, color = lambda)) +
     stat_ecdf(geom = "step") +
     facet_grid(p_1 ~ p_2) +
     coord_cartesian(xlim = c(0,4000))

stopping_time_ratios <- stopping_times_df %>% group_by(lambda, diluted_margin, p_1, p_2, prior_p_1, prior_p_2) %>%
  summarize(mean_stopping = mean(stopping_time)) %>%
  pivot_wider(names_from = "lambda", values_from = "mean_stopping") %>%
  mutate(mixture_oracle = mixture/oracle, apriori_oracle = apriori/oracle, adaptive_oracle = adaptive/oracle) 
  

stopping_times_table <- stopping_times_df %>%
  group_by(lambda, diluted_margin, p_1, p_2, prior_p_1, prior_p_2) %>%
  summarize(mean_stopping = round(mean(stopping_time)), percentile_stopping = round(quantile(stopping_time, .9))) %>%
  mutate(percentile_stopping = paste("(", percentile_stopping, ")", sep = "")) %>%
  unite("stopping", c("mean_stopping", "percentile_stopping"), sep = " ") %>%
  pivot_wider(names_from = "lambda", values_from = "stopping") %>%
  filter(round(diluted_margin, 3) == 0.05) %>%
  ungroup() %>%
  select(p_2, p_1, prior_p_2, prior_p_1, oracle, apriori, adaptive, mixture)  %>%
  arrange(p_2, p_1, prior_p_2, prior_p_1)

print(xtable::xtable(stopping_times_table, digits = 4), include.rownames = FALSE, digits = 4)


average_stopping_times <- stopping_times_df %>%
  group_by(lambda, diluted_margin, p_1, p_2, prior_p_1, prior_p_2)%>%
  summarize(mean_stopping = round(mean(stopping_time))) %>%
  ungroup() %>%
  pivot_wider(names_from = "lambda", values_from = mean_stopping) %>%
  mutate(ratio_adaptive = adaptive/oracle, ratio_apriori = apriori/oracle, ratio_mixture = mixture/oracle) %>% 
  summarize(geomean_adaptive = geomean(ratio_adaptive), 
            geomean_apriori = geomean(ratio_apriori), 
            geomean_mixture = geomean(ratio_mixture))


############ plot induced distribution on lambda from flat prior on overstatements #############

v <- .1
a <- 1 / (2 - v)
grid_1 <- seq(1e-6, 2 - 1/a, length.out = 25)
grid_2 <- seq(1e-6, 1 - 1/(2*a), length.out = 25)
mu <- c(.01, .001)
Sigma <- rbind(c(.02^2, .25 * .02^2 * .01^2), 
               c(.25 * .02^2 * .01^2, .01^2))

bivar_grid_unif <- expand.grid(grid_1, grid_2) %>%
  rename("p_1" = "Var1", "p_2" = "Var2") %>% 
  filter(a * p_2 + (a/2) * p_1 < a - 1/2) %>%
  mutate(weights = 1) %>%
  mutate(weights = weights/sum(weights)) 
lambda_star <- rep(NA, nrow(bivar_grid_unif))
for(i in 1:nrow(bivar_grid_unif)){
  lambda_star[i] <- optimal_lambda(a = a, p_1 = bivar_grid$p_1[i], p_2 = bivar_grid$p_2[i])
}
bivar_grid_unif <- bivar_grid_unif %>%
  mutate(lambda_star = lambda_star) %>%
  mutate(prior = "Uniform")

bivar_grid_normal <- expand.grid(grid_1, grid_2) %>%
  rename("p_1" = "Var1", "p_2" = "Var2") %>% 
  filter(a * p_2 + (a/2) * p_1 < a - 1/2) %>%
  mutate(weights = dmvnorm(cbind(p_1, p_2), mean = mu, sigma = Sigma)) %>%
  mutate(weights = weights/sum(weights)) 
lambda_star <- rep(NA, nrow(bivar_grid_normal))
for(i in 1:nrow(bivar_grid_normal)){
  lambda_star[i] <- optimal_lambda(a = a, p_1 = bivar_grid$p_1[i], p_2 = bivar_grid$p_2[i])
}
bivar_grid_normal <- bivar_grid_normal %>%
  mutate(lambda_star = lambda_star) %>%
  mutate(prior = "Normal")

bivar_grid <- bivar_grid_normal %>%
  bind_rows(bivar_grid_unif)

overstatement_plot_unif <- ggplot(bivar_grid %>% filter(prior == "Uniform"), aes(x = p_1, y = p_2, size = weights)) +
  geom_point() +
  coord_fixed() +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  ylab(TeX("2-vote overstatement rate ($p_2$)")) +
  xlab(TeX("1-vote overstatement rate ($p_1$)")) +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 18), axis.text = element_text(size = 14)) 
overstatement_plot_normal <- ggplot(bivar_grid %>% filter(prior == "Normal"), aes(x = p_1, y = p_2, size = weights)) +
  geom_point() +
  coord_fixed() +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  ylab(TeX("2-vote overstatement rate ($p_2$)")) +
  xlab(TeX("1-vote overstatement rate ($p_1$)")) +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 18), axis.text = element_text(size = 14)) 
lambda_plot_unif <- ggplot(bivar_grid %>% filter(prior == "Uniform"), aes(lambda_star, y = ..density.., weight = weights)) +
  geom_histogram(binwidth = 0.21) +
  xlab(TeX("Bet ($\\lambda$)")) + 
  ylab("Density") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 18), axis.text = element_text(size = 14)) 
lambda_plot_normal <- ggplot(bivar_grid %>% filter(prior == "Normal"), aes(lambda_star, y = ..density.., weight = weights)) +
  geom_histogram(binwidth = 0.25) +
  xlab(TeX("Bet ($\\lambda$)")) + 
  ylab("Density") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 18), axis.text = element_text(size = 14)) 

grid.arrange(
  overstatement_plot_unif, 
  lambda_plot_unif, 
  overstatement_plot_normal, 
  lambda_plot_normal)





