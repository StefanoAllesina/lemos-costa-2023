library(tidyverse)
library(ggpubr)

load("forplots.RData")

#likelihoods <- likelihoods %>% filter(grepl("Cadotte", label))
#observations <- observations %>% filter(grepl("Cadotte", label))

plot_histograms_facet_model <- function(likelihoods){
  # build pvals table
  dtpv <- likelihoods %>% filter(type == "original tree") %>% rename(orLogLik = logLik) %>% select(-type)
  dtpv <- likelihoods %>% filter(type != "original tree", type != "broom tree") %>% left_join(dtpv) 
  dtpv <- dtpv %>% mutate(better = (logLik >= orLogLik) * 1) %>% 
    group_by(type, model, label) %>% summarise(pval = round(mean(better),3), .groups = "drop")
  dtpv <- dtpv %>% pivot_wider(names_from = type, values_from = pval)
  dtpv <- dtpv %>% mutate(textpanel = paste0("p[rnd] ==", `random tree`))
  dtpv <- dtpv %>% mutate(textpanel2 = paste0("p[shuffle] ==", `shuffled tree`))
  dtpv <- dtpv %>% left_join(likelihoods %>% filter(type == "original tree") %>% group_by(model, label) %>% top_n(1, wt = logLik))
  
  pl <- ggplot(likelihoods) + aes(x = logLik) + 
    geom_histogram(data = likelihoods %>% filter(type != "original tree", type != "broom tree"), 
                   aes(fill = type),
                   position = "identity",
                   alpha = 0.65) + 
    geom_vline(data = likelihoods %>% filter(type == "original tree"),
               aes(xintercept = logLik), linetype = 2) +
    geom_vline(data = likelihoods %>% filter(type == "broom tree"),
               aes(xintercept = logLik), linetype = 1) +
    theme_bw() + 
    theme(legend.position = "none") + 
    scale_x_continuous("log-likelihood", limits = c(-355, -330)) + 
    scale_fill_manual(values = c("#2a9d8f", "#264653"))
  pl <- pl + facet_grid(model~.) + geom_text(
    data    = dtpv,
    mapping = aes(x = logLik, y = 250, label = textpanel), parse = TRUE,
    hjust = -0.1, vjust = -1, size = 3)
  pl <- pl + geom_text(
    data    = dtpv,
    mapping = aes(x = logLik, y = 200, label = textpanel2), parse = TRUE,
    hjust = -0.1, vjust = -1, size = 3)
  
  return(pl)
}

plot_histograms_facet_model_label <- function(likelihoods){
  # build pvals table
  likelihoods <- likelihoods %>%  mutate(label = paste(label, model))
  dtpv <- likelihoods %>% filter(type == "original tree") %>% rename(orLogLik = logLik) %>% select(-type)
  dtpv <- likelihoods %>% filter(type != "original tree", type != "broom tree") %>% left_join(dtpv)
  dtpv <- dtpv %>% mutate(better = (logLik >= orLogLik) * 1) %>% 
    group_by(type, model, label) %>% summarise(pval = round(mean(better),3), .groups = "drop")
  dtpv <- dtpv %>% pivot_wider(names_from = type, values_from = pval)
  #dtpv <- dtpv %>% mutate(textpanel = paste0("p[rnd] ==", `random tree`))
  #dtpv <- dtpv %>% mutate(textpanel2 = paste0("p[shuffle] ==", `shuffled tree`))
  #dtpv <- dtpv %>% left_join(likelihoods %>% filter(type == "original tree") %>% group_by(model, label) %>% top_n(1, wt = logLik))
  # pvalues in label
  likelihoods <- likelihoods %>% left_join(dtpv)
  likelihoods <- likelihoods %>%  mutate(label = paste(label, paste0("(", round(`random tree`, 3), ", ", round(`shuffled tree`, 3), ")")))
  
  pl <- ggplot(likelihoods) + aes(x = logLik) + 
    geom_histogram(data = likelihoods %>% filter(type != "original tree", type != "broom tree"), 
                   aes(fill = type),
                   position = "identity",
                   alpha = 0.65) + 
    geom_vline(data = likelihoods %>% filter(type == "original tree"),
               aes(xintercept = logLik), linetype = 2) +
    geom_vline(data = likelihoods %>% filter(type == "broom tree"),
               aes(xintercept = logLik), linetype = 1) +
    theme_bw() + 
    theme(legend.position = "none") + 
    xlab("log-likelihood") + 
    scale_fill_manual(values = c("#2a9d8f", "#264653"))
  pl <- pl + facet_wrap(label~., ncol = 4, scales = "free") 
  #+ geom_text(
  #   data    = dtpv,
  #   mapping = aes(x = logLik, y = 100, label = textpanel), parse = TRUE,
  #   hjust = -0.1, vjust = -1, size = 3)
  # pl <- pl + geom_text(
  #   data    = dtpv,
  #   mapping = aes(x = logLik, y = 80, label = textpanel2), parse = TRUE,
  #   hjust = -0.1, vjust = -1, size = 3)
  # 
  return(pl)
}


plot_predictions_without_faceting <- function(observations){
  pl <- ggplot(observations) + aes(x = observed, y = predicted) + 
    geom_point(alpha = 0.5) + geom_abline(slope = 1, intercept = 0) + 
    scale_x_log10("observed biomass") + scale_y_log10("predicted biomass") + 
    coord_flip() + 
    theme_bw() + 
    theme(legend.position = "bottom")
  labels <- observations %>% group_by(label, model) %>% 
         summarise(corr_log = cor(log(observed), log(predicted)),
                   corr_log_geom = cor(log(mean_observed), log(predicted)), 
                   yval = max(predicted),
                   .groups = "drop")
  marginstibble <- observations %>% group_by(model) %>% summarise(xval = min(observed), .groups = "drop")
  labels <- labels %>% left_join(marginstibble) %>% mutate(cr = paste0("rho == ", round(corr_log, 3)))
  pl <- pl + facet_wrap(~model) + geom_text(
    data    = labels,
    mapping = aes(x = xval, y = yval, label = cr), parse = TRUE,
    hjust = 1, vjust = -1, size = 3)
  return(pl)
}
