
# Experimental Treatment Outcome Summary

### Purpose: To visualize the experimental design and phenotypic outcome of the sctld transmission experiment as a foundation for 16s and rnaseq work.

![image](https://github.com/user-attachments/assets/2dc84e3b-7ec0-4080-a2a6-9b0f688c34dd)


## Summary Experimental Outcome Table
> Using the dataset 'metadata'
```
colnames(metadata)

### Load necessary packages
library(dplyr)
library(ggplot2)
library(ggpubr)

### Assuming your dataset is named 'metadata'
### Create the summary table
summary_table <- metadata %>%
  group_by(Coral.Species) %>%
  summarise(
    Control = sum(Outcome == "Control", na.rm = TRUE),
    Healthy = sum(Outcome == "Healthy", na.rm = TRUE),
    Diseased = sum(Outcome == "Disease", na.rm = TRUE)
  )

### Display the summary table
print(summary_table)

### Visualize summary table using ggplot
summary_table_plot <- ggtexttable(summary_table, rows = NULL, theme = ttheme("minimal"))
print(summary_table_plot)
```

## Calculate percent infected
```
infection_data <- metadata %>%
  group_by(Coral.Species) %>%
  summarise(
    Exposed = sum(Treatment == "Disease", na.rm = TRUE),
    Diseased = sum(Outcome == "Disease", na.rm = TRUE)
  ) %>%
  mutate(Percent.Infected = (Diseased / Exposed) * 100) %>%
  arrange(desc(Percent.Infected)) ### Sort in descending order

### Create a new label including sample count
infection_data <- infection_data %>%
  mutate(Label = paste0(Coral.Species, "\n(n=", Exposed, ")"))  ### Format: "Species (n=Samples)"

### Bar plot of percent infected with sample count labels
infection_plot <- ggplot(infection_data, aes(x = reorder(Label, -Percent.Infected), y = Percent.Infected)) +
  geom_bar(stat = "identity", fill = "#D8BFD8") +  ### Pastel purple fill
  geom_text(aes(label = round(Percent.Infected, 1)), vjust = -0.5) +  ### Add percent labels
  theme_minimal() +
  labs(
    title = "Percent Infected by Coral Species",
    x = "Coral Species",
    y = "Percent Infected"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylim(0, 100)  ### Set y-axis to 0-100

### Print updated infection plot
print(infection_plot)
```


```
### Load necessary libraries
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(patchwork)
library(ggtext)  ### Required for markdown formatting

### Map abbreviated species names to full italicized names
species_labels <- c(
  "Acer" = "A. cervicornis",
  "Mcav" = "M. cavernosa",
  "Ofav" = "O. faveolata",
  "Past" = "P. astreoides"
)

### Update species names in summary table
summary_table_renamed <- summary_table %>%
  mutate(Species = species_labels[Coral.Species]) %>%
  select(Species, Control, Healthy, Diseased) %>%
  rename("Control" = Control, "Disease-Exposed" = Healthy, "Disease-Infected" = Diseased)

### Convert summary table to a ggplot table (Panel A)
summary_table_plot <- ggtexttable(summary_table_renamed, rows = NULL, theme = ttheme("minimal"))

### Add a title to the table
summary_table_plot <- annotate_figure(summary_table_plot, 
                                      top = text_grob("Experimental Design", face = "bold", size = 14)
)

### Print table with title
print(summary_table_plot)
### Update species names in infection data
infection_data <- infection_data %>%
  mutate(Label = paste0(species_labels[Coral.Species], "\n(n=", Exposed, ")"))

### Create Panel B (infection plot with sample counts)
infection_plot <- ggplot(infection_data, aes(x = reorder(Label, -Percent.Infected), y = Percent.Infected)) +
  geom_bar(stat = "identity", fill = "#D8BFD8") +  ### Pastel purple fill
  geom_text(aes(label = round(Percent.Infected, 1)), vjust = -0.5) +  ### Add percent labels
  theme_minimal() +
  labs(
    title = "**Disease Prevalence**",
    x = "Coral Species",
    y = "Percent Infected"
  ) +
  theme(
    plot.title = element_markdown(size = 14),  ### Make "B" bold
    axis.text.x = element_markdown(size = 10, angle = 45, hjust = 1, vjust = 1)  ### Italicize species names
  ) +
  ylim(0, 100)  ### Set y-axis to 0-100

### Combine the two panels using patchwork with bold A and B
final_plot <- (summary_table_plot / infection_plot) + 
  plot_annotation(tag_levels = 'A', title = "Experimental Design")  ### Make "A" bold

### Print the final multi-panel figure
print(final_plot)
```

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

##  Relative Risk
```
### Load necessary libraries
library(dplyr)
library(rjags)
library(coda)
library(tidyr)
library(ggplot2)

### Prepare data for relative risk calculation
rr_data <- metadata %>%
  group_by(Coral.Species) %>%
  summarise(
    Exposed_Diseased = sum(Treatment == "Disease" & Outcome == "Disease", na.rm = TRUE),
    Exposed_Total = sum(Treatment == "Disease", na.rm = TRUE),
    NonExposed_Diseased = sum(Treatment != "Disease" & Outcome == "Disease", na.rm = TRUE),
    NonExposed_Total = sum(Treatment != "Disease", na.rm = TRUE)
  ) %>%
  mutate(
    Risk_Exposed = Exposed_Diseased / Exposed_Total,
    Risk_NonExposed = NonExposed_Diseased / NonExposed_Total,
    RR = Risk_Exposed / Risk_NonExposed
  )


### Prepare data for JAGS
jags_data <- list(
  N = nrow(rr_data),
  Exposed_Diseased = rr_data$Exposed_Diseased,
  Exposed_Total = rr_data$Exposed_Total,
  NonExposed_Diseased = rr_data$NonExposed_Diseased,
  NonExposed_Total = rr_data$NonExposed_Total
)

### Define JAGS model using log transformation for numerical stability
jags_model <- "
model {
  for (i in 1:N) {
    Exposed_Diseased[i] ~ dbin(p_exposed[i], Exposed_Total[i])
    NonExposed_Diseased[i] ~ dbin(p_nonexposed[i], NonExposed_Total[i])
    
    p_exposed[i] ~ dbeta(1,1)  ### Non-informative priors
    p_nonexposed[i] ~ dbeta(1,1)

    logRR[i] <- log(p_exposed[i]) - log(p_nonexposed[i])  ### Log Relative Risk
  }
}
"

### Run JAGS model
jags_fit <- jags.model(textConnection(jags_model), data = jags_data, n.chains = 3, n.adapt = 1000)
update(jags_fit, 2000)  ### Burn-in period
samples <- coda.samples(jags_fit, variable.names = c("logRR"), n.iter = 10000)

### Convert logRR to RR
rr_results <- as.data.frame(do.call(rbind, lapply(samples, as.matrix)))

### Reshape from wide to long format and exponentiate to get RR
rr_long <- rr_results %>%
  pivot_longer(cols = everything(), names_to = "Species_Index", values_to = "logRR") %>%
  mutate(
    Species_Index = as.numeric(gsub("logRR\\[(\\d+)\\]", "\\1", Species_Index)),
    RR = exp(logRR)  ### Convert logRR back to RR
  )


### Compute 95% Credible Intervals per species
### > Because our nonexposed (control) samples never got disease, this introduces a zero which is being divided by which makes the calculations go haywire because its infinite. So people use relative risk but the point of it now becomes whether or not the lower credible interval limit overlaps the 1.0 line or not as a signficant or non signficant determiniation of the relative risk of disease relative to control samples. So to consdier all that, here is how the CrI was calculated:
rr_summary <- rr_long %>%
  group_by(Species_Index) %>%
  summarise(
    Median_RR = median(RR),
    Lower_CrI = quantile(RR, 0.025),
    Upper_CrI = median(RR)+(median(RR)-Lower_CrI),
    .groups = "drop"
  )

### Merge species names back
rr_final <- rr_summary %>%
  left_join(rr_data %>% mutate(Species_Index = row_number()), by = "Species_Index") %>%
  select(Coral.Species, Median_RR, Lower_CrI, Upper_CrI) %>%
  mutate(Significant = ifelse(Lower_CrI > 1 | Upper_CrI < 1, "Yes", "No"))

### Print Bayesian RR estimates
print(rr_final)

### Manual Validation: Check if Bayesian RR matches expected RR
rr_validation <- rr_data %>%
  mutate(
    Risk_Exposed = Exposed_Diseased / Exposed_Total,
    Risk_NonExposed = NonExposed_Diseased / NonExposed_Total,
    RR_check = Risk_Exposed / 0.1  ### Manually computed RR
  )

### Print validation check
print(rr_validation)

median_rr <- median(rr_long$RR)
print(median_rr)

### Visualization: Bayesian RR with Credible Intervals
rr_plot <- ggplot(rr_final, aes(y = reorder(Coral.Species, Median_RR), x = Median_RR, color = Significant)) +
  geom_point(size = 3) +  ### Dots for RR values
  geom_errorbarh(aes(xmin = Lower_CrI, xmax = Upper_CrI), height = 0.2) +  ### Horizontal credible interval bars
  geom_vline(xintercept = 1, linetype = "solid", color = "black") +  ### Significance reference line
  geom_vline(xintercept = median_rr, linetype = "dashed", color = "grey") +  ### Significance reference line
  scale_color_manual(values = c("Yes" = "#BA55D3", "No" = "gray")) +  ### Color by significance
  theme_minimal() +
  labs(
    title = "Relative Risk of Disease Incidence",
    y = "Coral Species",
    x = "Relative Risk (RR)",
    color = "Significant?"
  ) +
  #xlim(0,15)+
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        legend.position = "none"  ### Hides the legend
        )

### Print RR plot
print(rr_plot)



median_rr <- median(rr_long$RR)
print(median_rr)




final_plot <- (summary_table_plot / (infection_plot + rr_plot)) + 
  plot_annotation(tag_levels = 'A')  ### Bold title

### Print the final multi-panel figure
print(final_plot)
```

### Save rr_final to a CSV file
```
write.csv(rr_final, "rr_final_output.csv", row.names = FALSE)
```
