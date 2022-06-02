#################
### Libraries ###
#################
# install.packages(c("ggrepel","viridis","survival","survminer","tidyverse")) # install needed packages if not installed already
library(ggrepel)
library(viridis)
library(survival)
library(survminer)
library(tidyverse)

####################
### Import files ###
####################
df <- read_csv("./Data.csv") %>% filter(Censored == 0)
meta <- read_csv("./Metainfo.csv")

########################
#### Diagonal plots ####
########################
df_plot <- df %>%
  group_by(Trial,Strain,RNAi) %>%
  summarise(RNAi_ls = median(Time, na.rm = TRUE)) %>%
  ungroup() %>%
  spread(key = RNAi,value = RNAi_ls) %>%
  mutate(delta_lifespan = `daf-2` - L4440,
         control_lifespan = L4440,
         extended_lifespan = `daf-2`) %>%
  replace_na(list(delta_lifespan = 0,control_lifespan = 0,extended_lifespan = 0)) %>%
  gather(key = "RNAi",value = "RNAi_ls",-Trial, -Strain, -delta_lifespan,-control_lifespan,-extended_lifespan) %>%
  left_join(x = df,y = .,by = c("Trial", "Strain", "RNAi")) %>%
  arrange(desc(control_lifespan)) %>% # decide how to order data for plots
  mutate(Strain = factor(Strain,unique(Strain)))

df_plot_label <- df_plot %>%
  group_by(Strain,RNAi) %>%
  summarise(Time = mean(Time,na.rm = T))

stat_ls <- df %>%
  group_by(Strain,RNAi) %>%
  summarise(RNAi_mean = mean(Time),
            ls_ste = sd(Time) / sqrt(length(Time))) %>%
  ungroup()

df_mean <- stat_ls %>%
  dplyr::select(Strain,RNAi,RNAi_mean) %>%
  spread(key = RNAi,value = RNAi_mean) %>%
  rename(`daf-2 Mean` = `daf-2`, `L4440 Mean` = L4440) %>%
  mutate(delta_lifespan = `daf-2 Mean` - `L4440 Mean`,
         control_lifespan = `L4440 Mean`,
         extended_lifespan = `daf-2 Mean`)

df_ste <- stat_ls %>%
  dplyr::select(Strain,RNAi,ls_ste) %>%
  spread(key = RNAi,value = ls_ste) %>%
  rename(`daf-2 STE` = `daf-2`, `L4440 STE` = L4440)

df_range <- full_join(x = df_mean,y = df_ste) %>%
  mutate(`daf-2 min` = `daf-2 Mean` - `daf-2 STE`,
         `daf-2 max` = `daf-2 Mean` + `daf-2 STE`,
         `L4440 min` = `L4440 Mean` - `L4440 STE`,
         `L4440 max` = `L4440 Mean` + `L4440 STE`,
         noeffect_y = `L4440 Mean`,
         deviation_from_no_effect = `daf-2 Mean` - `L4440 Mean`
  ) %>%
  drop_na()

N2_slope <- df_range %>%
  filter(Strain == "N2") %>%
  transmute(N2_slope = `daf-2 Mean` / `L4440 Mean`) %>%
  pull(N2_slope)

# Prepare diagonal plot
quantile_cutoff <- 0.4
plot_range <- df_range %>%
  left_join(x = .,y = meta, by = "Strain") %>%
  mutate(label_all = case_when(is.na(combined_genes) ~ Strain,TRUE ~ paste0(Strain, " (",combined_genes,")"))) %>%
  mutate(label_show = ifelse(deviation_from_no_effect < quantile(deviation_from_no_effect,quantile_cutoff) | Strain == "N2", label_all, NA)) %>%
  mutate(combined_category = case_when(combined_category == "" ~ "Adhesion Receptors and Downstream Signaling", TRUE ~ combined_category),
         combined_division = case_when(combined_division == "" ~ "Adhesion Receptors and Downstream Signaling", TRUE ~ combined_division)
  )
range_ls <- plot_range %>% dplyr::select(`L4440 Mean`,`daf-2 Mean`) %>% unlist() %>% range(.,na.rm = TRUE)

# Make diagonal plots
coord_range <- c(range_ls[1],range_ls[2])
show_raw <- plot_range %>%
  filter(combined_division != "Not in matrisome") %>%
  filter(Strain != "WS3403") #(is an overexpressor)

# Save supplementary table
show_raw %>% select(-contains("list"),-Mutation) %>% write_csv(path = paste0("./diag_plot_extended_data.csv"))

# Add N2 as overlay to each panel.
show_category <- show_raw %>%
  filter(Strain != "N2") %>%
  select(combined_category) %>%
  distinct() %>%
  mutate(Strain = "N2") %>%
  arrange(combined_category) %>%
  arrange(combined_category %in% c("Adhesion Receptors and Downstream Signaling"))

show_division <- show_raw %>%
  filter(Strain != "N2") %>%
  select(combined_division) %>%
  distinct() %>%
  mutate(Strain = "N2") %>%
  arrange(combined_division) %>%
  arrange(combined_division %in% c("Adhesion Receptors and Downstream Signaling"))

show_N2_category <- show_raw %>%
  filter(Strain == "N2") %>%
  select(-combined_category) %>%
  full_join(show_category, by = "Strain")

show_N2_division <- show_raw %>%
  filter(Strain == "N2") %>%
  select(-combined_division) %>%
  full_join(show_division, by = "Strain")

show <- show_raw %>%
  filter(Strain != "N2") %>%
  mutate(combined_category_fct = factor(combined_category,levels = show_category %>% pull(combined_category)),
         combined_division_fct = factor(combined_division,levels = show_division %>% pull(combined_division)))

#############################
### Main plots and export ###
#############################
p_division <- ggplot(data = show,mapping = aes(x = `L4440 Mean`, y = `daf-2 Mean`, color = combined_division, alpha = deviation_from_no_effect)) +
  # N2 data
  geom_segment(data = show_N2_division, aes(x = `L4440 Mean`, y = noeffect_y, xend = `L4440 Mean`, yend = `daf-2 Mean`), size = 0.7, linetype = 1, color = "black", alpha = 0.2) +
  geom_point(data = show_N2_division,size = 3, shape = 18, color = "black", alpha = 0.7) +
  geom_errorbarh(data = show_N2_division,aes(xmin = `L4440 Mean` - `L4440 STE`, xmax = `L4440 Mean` + `L4440 STE`, y = `daf-2 Mean`), size = 0.6, height = 0.6, color = "black", alpha = 0.7) +
  geom_errorbar(data = show_N2_division,aes(ymin = `daf-2 Mean` - `daf-2 STE`, ymax = `daf-2 Mean` + `daf-2 STE`, x = `L4440 Mean`), size = 0.6, width = 0.6, color = "black", alpha = 0.7) +
  # diagonal lines
  geom_abline(aes(intercept = 0, slope = N2_slope), linetype =2) +
  geom_abline(aes(intercept = 0, slope = 1, color = combined_division), linetype = 1) +
  # Remaining strains
  geom_point(size = 3, shape = 18, alpha = 1) +
  coord_cartesian(xlim = coord_range, ylim = coord_range) +
  geom_segment(aes(x = `L4440 Mean`, y = noeffect_y, xend = `L4440 Mean`, yend = `daf-2 Mean`), size = 1, linetype = 1) +
  geom_errorbarh(data = show,aes(xmin = `L4440 Mean` - `L4440 STE`, xmax = `L4440 Mean` + `L4440 STE`, y = `daf-2 Mean`), size = 0.6, height = 0.6, alpha = 1) +
  geom_errorbar(aes(ymin = `daf-2 Mean` - `daf-2 STE`, ymax = `daf-2 Mean` + `daf-2 STE`, x = `L4440 Mean`), size = 0.6, width = 0.6, alpha = 1) +
  geom_label_repel(aes(label = label_all),show.legend = FALSE, box.padding = 1,  size = 2, segment.size = 0.3, segment.alpha = 0.5, force = 5, nudge_x = range_ls[2] * 2, fontface = "italic", direction    = "y", alpha = 0.7) +
  scale_color_viridis_d(option = "D",direction = -1,end = 0.9) +
  scale_alpha_continuous(range = c(1,0),guide = F) +
  facet_wrap(. ~ combined_division) +
  guides(color = guide_legend(title = "Deviation from the zero-effect line")) +
  labs(title = "Survival by matrisome division",subtitle = expression(paste(italic("solid line: "), "no-effect trajectory / ", italic("dashed line: "), "N2 extrapolation"))) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(face="bold.italic"))
ggsave(plot = p_division,filename = paste0("./p_diag_plot_divisions.pdf"),width = 10,height = 10)

p_category <- ggplot(data = show,mapping = aes(x = `L4440 Mean`, y = `daf-2 Mean`, color = combined_category, alpha = deviation_from_no_effect)) +
  # N2 data
  geom_segment(data = show_N2_category, aes(x = `L4440 Mean`, y = noeffect_y, xend = `L4440 Mean`, yend = `daf-2 Mean`), size = 0.7, linetype = 1, color = "black", alpha = 0.2) +
  geom_point(data = show_N2_category,size = 3, shape = 18, color = "black", alpha = 0.7) +
  geom_errorbarh(data = show_N2_category,aes(xmin = `L4440 Mean` - `L4440 STE`, xmax = `L4440 Mean` + `L4440 STE`, y = `daf-2 Mean`), size = 0.6, height = 0.6, color = "black", alpha = 0.7) +
  geom_errorbar(data = show_N2_category,aes(ymin = `daf-2 Mean` - `daf-2 STE`, ymax = `daf-2 Mean` + `daf-2 STE`, x = `L4440 Mean`), size = 0.6, width = 0.6, color = "black", alpha = 0.7) +
  # diagonal lines
  geom_abline(aes(intercept = 0, slope = N2_slope), linetype =2) +
  geom_abline(aes(intercept = 0, slope = 1, color = combined_category), linetype = 1) +
  # Remaining strains
  geom_point(size = 3, shape = 18, alpha = 1) +
  coord_cartesian(xlim = coord_range, ylim = coord_range) +
  geom_segment(aes(x = `L4440 Mean`, y = noeffect_y, xend = `L4440 Mean`, yend = `daf-2 Mean`), size = 1, linetype = 1) +
  geom_errorbarh(data = show,aes(xmin = `L4440 Mean` - `L4440 STE`, xmax = `L4440 Mean` + `L4440 STE`, y = `daf-2 Mean`), size = 0.6, height = 0.6, alpha = 1) +
  geom_errorbar(aes(ymin = `daf-2 Mean` - `daf-2 STE`, ymax = `daf-2 Mean` + `daf-2 STE`, x = `L4440 Mean`), size = 0.6, width = 0.6, alpha = 1) +
  geom_label_repel(aes(label = label_all),show.legend = FALSE, box.padding = 1,  size = 2, segment.size = 0.3, segment.alpha = 0.5, force = 5, nudge_x = range_ls[2] * 2, fontface = "italic", direction    = "y", alpha = 0.7) +
  scale_color_viridis_d(option = "D",direction = -1,end = 0.9) +
  scale_alpha_continuous(range = c(1,0),guide = F) +
  facet_wrap(. ~ combined_category_fct) +
  guides(color = guide_legend(title = "Deviation from the zero-effect line")) +
  labs(title = "Survival by matrisome category",subtitle = expression(paste(italic("solid line: "), "no-effect trajectory / ", italic("dashed line: "), "N2 extrapolation"))) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(face="bold.italic"))
ggsave(plot = p_category,filename = paste0("./p_diag_plot_categories.pdf"),width = 10,height = 10)


