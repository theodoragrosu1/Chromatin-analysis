library(dplyr)
library(ggplot2)


#AP-1  #change this based on which footprint needs to be plotted
input_file <- "./tobias/tf_specific_footprints/Jun-AP1_Jun-AP1_footprint.txt"

plot_aggregate(input_file, title = "Jun-AP-1", labels_facets = c("WT", "KO")) +
  scale_colour_manual(values = c("#000000", "#E69F00"), guide = "none")


plot_aggregate <- function(input_file, title = "", labels_facets){
  # read aggregate data into a data frame and make it into a "long" format for easier plotting
  df <- read.table(input_file, sep = "\t", col.names = c("signal", "regions", "aggregate"))
  n_pos <- stringr::str_count(df$aggregate[1], pattern = ",") + 1
  df_long <- tidyr::separate(df, aggregate, sep = ",", convert = TRUE, into = paste0("pos", 1:n_pos)) %>% 
    tidyr::pivot_longer(cols = starts_with("pos"), 
                        names_to = "position", names_prefix = "pos", names_transform = list(position = as.numeric), 
                        values_to = "score") %>% 
    tidyr::extract(regions, into = c("output_prefix", "regions_group"), 
                   regex = "(.*)_(bound|unbound|all)", remove = FALSE) 
  
  # remove extra text from "signal" column to just get cluster ID
  df_long <- df_long %>% 
    mutate(signal = gsub("_.mRp.clN.sorted_corrected", "", signal)) %>%
    mutate(signal = gsub("Normal_integrated_with_ambRNAremoval_cluster_", "", signal)) %>%
    mutate(signal = factor(signal, levels = labels_facets))
  # you should here use mutate(signal = factor(signal, levels = cluster_order))
  # to force the clusters to be plotted in the right developmental order
  # you can also use left_join() to add your cluster annotations if you want to
  # have cluster names instead of just numbers
  
  # get mean signal in flank and centre regions
  df_mean_signal <- df_long %>% 
    mutate(position_group = case_when(position <= 50 ~ "flank",
                                      position > 50 & position <= 70 ~ "centre",
                                      position > 70 ~ "flank")) %>% 
    group_by(signal, position_group) %>% 
    summarise(mean_score = mean(score)) %>% 
    tidyr::pivot_wider(names_from = position_group, values_from = mean_score) %>% 
    mutate(footprint_depth = flank - centre) 
  
  
  # construct plot
  p <- df_long %>% 
    left_join(df_mean_signal) %>% 
    ggplot(aes(x = position, y = score, colour = signal)) +
    geom_line() +
    theme_bw() +
    geom_vline(xintercept = c(50, 70), lty = 3) +
    geom_segment(aes(y = centre, yend = centre, x = 50, xend = 70), lty = 3, colour = "#C44536") +
    geom_segment(aes(y = flank, yend = flank, x = 0, xend = 50), lty = 3, colour = "#C44536") +
    facet_wrap(~signal, nrow = 1) +
    labs(x = "Position", y = "Score") +
    ggtitle(title) + 
    xlim(0, 125) +
    guides(colour = "none") +
    theme_bw(base_size = 22)
  return(p)
}


