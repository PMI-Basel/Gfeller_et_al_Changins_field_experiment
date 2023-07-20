## Wheat phenotyping

# Output for anova
output_anova_phen <- function(data) {
  t.anova <- data %>% 
    rownames_to_column(var = "Variable") %>%
    tibble() %>%
    mutate(Variable = str_replace_all(Variable, "trt:PCA_soil_chem_1", "C x C"), 
           Variable = str_replace_all(Variable, "trt", "Cond"),
           Variable = str_replace_all(Variable, "PCA_soil_chem_1", "Chem"),
           p_new = case_when(`Pr(>F)` < 0.001   ~  paste0("<sup>_***_</sup>"),
                             `Pr(>F)` < 0.01    ~  paste0("<sup>_**_</sup>"),
                             `Pr(>F)` < 0.05    ~  paste0("<sup>_*_</sup>"),
                             `Pr(>F)` >= 0.05   ~  paste0("<sup>ns</sup>")),
           across(where(is.numeric), ~ as.character(signif(., 2))),
           across(everything(), ~ replace_na(., ""))) 
  
  # tidy anova table for plot
  lab_anova <- paste0(t.anova$Variable[1], t.anova$p_new[1], ", ", 
                      t.anova$Variable[2], t.anova$p_new[2], ", ",
                      t.anova$Variable[3], t.anova$p_new[3])
  return(lab_anova)
}





# Add anova table to plot
add_anova_tab_phen <- function(tab, ...) {
  geom_richtext(data = data.frame(), mapping =  aes(x = -Inf, y = Inf, label = {{tab}}),
                fill = alpha(c("white"), 0.7), color = "black",
                label.padding = grid::unit(rep(0.1, 3), "lines"),
                hjust = -0.05, label.color = alpha(c("white"), 0.7),
                vjust = 1.5, size = 3.5,
                label.r = unit(0, "lines"))
}

# Boxplot wheat phenotyping
plot_box_wheat_feed <- function(my_dat, var_measured, var_col_grad = PCA_soil_chem_1, 
                                y_lab = NULL, tab_anova = NULL){
  p <- my_dat %>% 
    ggplot(aes(x = trt, y ={{ var_measured }})) +
    geom_boxplot(aes(fill = trt), outlier.shape = NA, alpha = 0.3, 
                 show.legend = FALSE) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 6, 
                 color = "black") +
    stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.8, 
                 color = "black", width = 0.25, alpha = 0.65) + 
    geom_quasirandom(aes(color = {{ var_col_grad }}), size = 3, shape = 16, 
                     alpha = 0.9, width = 0.2) +
    scale_x_discrete(name = "Conditioning", 
                     labels = c(W = "WT", b = expression(italic(bx1)))) +
    scale_fill_manual(values = c("W" = "gold", "b" = "darkgreen")) +
    scale_colour_continuous(name = "Soil chemistry PC1")  +
    # scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
    {if(!is.null(y_lab))ylab(y_lab)} +
    {if(!is.null(tab_anova))add_anova_tab_phen(tab_anova)} +
    theme(legend.position = "top") 
  
  return(p)
}




# plot regression plot position
plot_pos_reg <- function(my_dat, var_measured, y_lab = NULL){
  p <- my_dat %>% 
    ggplot(aes(x = pos_lat, y ={{ var_measured }}, color = trt)) + 
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("W" = "gold3", "b" = "darkgreen"),
                       name = "Conditioning", 
                       labels = c("W" = "W22", "b" = expression(italic(bx1)))) +
    xlab("Latitudinal position")  +
    # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    {if(!is.null(y_lab))ylab(y_lab)} +
    theme(legend.position = "top") +
    geom_smooth(method = lm, se = FALSE, fullrange = TRUE) +
    ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
                     label.x.npc = "middle")
  
  
  return(p)
}

# plot regression PCA soil chem
plot_pc1_reg <- function(my_dat, var_measured, y_lab = NULL,
                         tab_anova = NULL){
  p <- my_dat %>% 
    ggplot(aes(x = PCA_soil_chem_1, y ={{ var_measured }}, color = trt)) + 
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("W" = "gold3", "b" = "darkgreen"),
                       name = "Conditioning", 
                       labels = c("W" = "W22", "b" = expression(italic(bx1)))) +
    xlab("Soil chemistry PC1")  +
    # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    {if(!is.null(y_lab))ylab(y_lab)} +
    {if(!is.null(tab_anova))add_anova_tab_phen(tab_anova)} +
    theme(legend.position = "top") +
    geom_smooth(method = lm, se = FALSE, fullrange = TRUE, formula = 'y ~ x') +
    ggpubr::stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
                     label.x.npc = "left", label.y.npc = "top")
  
  return(p)
}


# combine plots long
plot_comb_long <- function(p.cond, p.soil, p.psf){
  p <- plot_grid(p.cond, p.soil, p.psf, labels=LETTERS[1:3], align = "h", 
                 ncol = 3, rel_widths = c(1, 1.6, 1.6))
  return(p)
}

# combine plots
plot_comb <- function(p.cond, p.psf){
  p <- plot_grid(p.cond, p.psf, 
                 align = "h", rel_widths = c(1, 1.6))
  return(p)
}


# Plot BXs at maize harvest along soil chemistry PC1
plot_lin_reg_bx_cond <- function() {
  list(
    geom_point(size = 2, alpha = 3/4), 
    stat_smooth(method = "lm", formula = 'y ~ x',
                color = "black", se = FALSE, fullrange = TRUE),
    ylab("Concentration (ng/ml)"),
    xlab("Soil chemistry PC1"),
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
             vjust = 2),
    # scale_fill_manual(labels = c(W = "WT", b = expression(italic(bx1))),
    #                   values = c(W = "gold3", b = "darkgreen")),
    facet_wrap(vars(Compound_name), nrow = 2, scales = "free"))
}

