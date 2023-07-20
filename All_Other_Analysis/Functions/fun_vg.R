
## plot model assumptions
plot.mod.vg <- function(mod){
  par(mfrow = c(1,2))
  plot1 <- plot(fitted(mod), resid(mod), xlab = "Fitted values", ylab = "Residuals", main = "Tukey-Anscombe plot")
  plot2 <- car::qqPlot(resid(mod), dist = "norm", mean = mean(resid(mod)), sd = sd(resid(mod)),xlab = "Theoretical quantiles", ylab = "Empirical quantiles", main = "Q-Q plot of residuals")
}




# neighborhood analysis fold change
surrounding_fold_change <- function(plot_id, variable, treatment, x_coordinates, y_coordinates, data){
  
  tab.1 <- data %>%
    select({{ plot_id }}, {{ y_coordinates }}, {{ x_coordinates }}, {{ variable }}) %>%
    mutate(pos_x_y = paste0({{ x_coordinates }}, ";",{{ y_coordinates }}),
           pos_x_plus1_y = paste0({{ x_coordinates }}+1, ";",{{ y_coordinates }}),
           pos_x_minus1_y = paste0({{ x_coordinates }}-1, ";",{{ y_coordinates }}),
           pos_x_y_plus1 = paste0({{ x_coordinates }}, ";",{{ y_coordinates }}+1),
           pos_x_y_minus1 = paste0({{ x_coordinates }}, ";",{{ y_coordinates }}-1))
  
  tab.2 <- tibble(plot_id = character(), 
                  variable_focal_plot = numeric(), 
                  surrounding_mean = numeric(),
                  fold_change = numeric())
  tab.2[1:nrow(tab.1), 1:ncol(tab.2)] <- NA
  
  for (i in 1: nrow(tab.1)) {
    results <- tab.1 %>%
      filter(pos_x_y %in% c(pos_x_plus1_y[i], pos_x_minus1_y[i],
                            pos_x_y_plus1[i], pos_x_y_minus1[i])) %>% 
      summarise(
        plot_id = pull(tab.1, {{ plot_id }})[i],
        variable_focal_plot = pull(tab.1, {{ variable }})[i],
        surrounding_mean = mean({{ variable }}, na.rm = TRUE),
        fold_change = case_when(str_detect({{ plot_id }}, "bx") ~ surrounding_mean/variable_focal_plot,
                                str_detect({{ plot_id }}, "WT") ~ variable_focal_plot/surrounding_mean))
    
    tab.2$plot_id[i] <- results$plot_id[]
    tab.2$variable_focal_plot[i] <- results$variable_focal_plot
    tab.2$surrounding_mean[i] <- results$surrounding_mean
    tab.2$fold_change[i] <- results$fold_change
  }
  tab.3 <- tab.2 %>% 
    rename("fold_change_{{variable}}" := fold_change) %>% 
    left_join(data, ., by = c("plot_id" = "plot_id")) 
  
  output <- if_else(identical(paste0(tab.3$variable_focal_plot), paste0(pull(tab.1, {{ variable }}))), "Matching Successful!", "Error!!! Not matching")
  
  tab.4 <- tab.3 %>% 
    select( - variable_focal_plot, -surrounding_mean )%>%
    as_tibble()
  
  message(output)
  
  return(tab.4)
  
}







# neighborhood analysis log-response ratio (LRR)
surrounding_lrr <- function(plot_id = plot_id, variable, treatment, x_coordinates = length,
                            y_coordinates = width, data){
  
  # make table with coordinates of surrounding plots
  tab.1 <- data %>%
    select({{ plot_id }}, {{ y_coordinates }}, {{ x_coordinates }}, {{ variable }}) %>%
    mutate(pos_x_y = paste0({{ x_coordinates }}, ";",{{ y_coordinates }}),
           pos_x_plus1_y = paste0({{ x_coordinates }}+1, ";",{{ y_coordinates }}),
           pos_x_minus1_y = paste0({{ x_coordinates }}-1, ";",{{ y_coordinates }}),
           pos_x_y_plus1 = paste0({{ x_coordinates }}, ";",{{ y_coordinates }}+1),
           pos_x_y_minus1 = paste0({{ x_coordinates }}, ";",{{ y_coordinates }}-1))
  
  # make output table
  tab.2 <- tibble(plot_id = character(), 
                  variable_focal_plot = numeric(), 
                  surrounding_mean = numeric(),
                  fold_change = numeric())
  tab.2[1:nrow(tab.1), 1:ncol(tab.2)] <- NA
  
  # calculate fold change relative to surrounding mean of other treatment (WT/bx)
  for (i in 1: nrow(tab.1)) {
    results <- tab.1 %>%
      filter(pos_x_y %in% c(pos_x_plus1_y[i], pos_x_minus1_y[i],
                            pos_x_y_plus1[i], pos_x_y_minus1[i])) %>% 
      summarise(
        plot_id = pull(tab.1, {{ plot_id }})[i],
        variable_focal_plot = pull(tab.1, {{ variable }})[i],
        surrounding_mean = mean({{ variable }}, na.rm = TRUE),
        fold_change = case_when(str_detect({{ plot_id }}, "bx") ~ 
                                  log(surrounding_mean/variable_focal_plot),
                                str_detect({{ plot_id }}, "WT") ~ 
                                  log(variable_focal_plot/surrounding_mean)))
    
    tab.2$plot_id[i] <- results$plot_id[]
    tab.2$variable_focal_plot[i] <- results$variable_focal_plot
    tab.2$surrounding_mean[i] <- results$surrounding_mean
    tab.2$fold_change[i] <- results$fold_change
  }
  
  # change name of variable to fold_change_ + variable name and append to original data frame
  tab.3 <- tab.2 %>% 
    rename("lrr_{{variable}}" := fold_change) %>% 
    left_join(data, ., by = c("plot_id" = "plot_id")) 
  
  # remove unneeded columns
  tab.4 <- tab.3 %>% 
    select( - variable_focal_plot, -surrounding_mean )%>%
    as_tibble()

  return(tab.4)
  
}






# Plot LRR
# Base function 
ggplot_lrr <- function(data, x = PCA_soil_chem_1, y = !!LRR, fill, ...) {
  
  data %>% 
    mutate(psf = case_when(
      {{y}} < 0  ~ "neg",
      {{y}} >= 0  ~ "pos")) %>% 
    ggplot(aes(x = {{x}}, y = {{y}})) +
    geom_rect(aes(xmin = {{x}}-0.06, 
                  xmax = {{x}}+0.06,
                  ymax = {{y}}, 
                  ymin = 0,
                  fill = psf),
              alpha = 0.3) +
    geom_point(size = 2, alpha = 0.9) 
}

# Make visually more appealing
impr_graph_lrr <- function(data, x = PCA_soil_chem_1, y = !!LRR) {
  
  # Calculate position for legend 
  # pos_legend <- data %>% 
  #   summarise(
  #     max = max({{y}}),
  #     min = min({{y}}),
  #     range = (max - min),
  #     y = (0-min)/range) %>% 
  #   pull(y) %>% 
  #   c(1.23, .)
  
  # Calculate linear regression for label
  mod <- data %>%
    rename("LRR" = {{y}},
           x = {{x}}) %>% 
    lm(LRR ~ x, data = .)
  
  r.squared <- summary(mod)$r.squared %>% signif(digits = 2)
  p.val <- anova(mod)$'Pr(>F)'[1] %>% signif(digits = 2)
  # p.val <- if_else(p.val < 0.001, paste0("p < 0.001"),
  #           paste0("p = ", sprintf("%.3f", p.val)))
  # 
  
  # Add geoms
  list(
    geom_hline(aes(yintercept = 0), color = "black", alpha = 0.3, size = 0.7),
    scale_fill_manual(labels = c(pos = "Positive feedback", 
                                 neg = "Negative feedback"),
                      values = c(pos = "gold2", 
                                 neg = "darkgreen")),
    geom_smooth(method = lm, se = FALSE, formula = 'y ~ x',
                color = "black"),
    xlab("Soil chemistry PC1"),
    theme(legend.title = element_blank(),
#          legend.position = pos_legend,
#          plot.margin = unit(c(0.2, 5, 0.2, 0.2),"cm"),
          legend.position = "top"),
    annotate('text', x = -Inf, y = Inf, 
             label = paste("italic(R)^2==", r.squared, "~','~",
                           "~italic(p)==", p.val), 
             parse = TRUE, hjust = -0.1, vjust = 1.5, size = 3.5)
  )
}

