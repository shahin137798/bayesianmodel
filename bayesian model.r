library(readr)
library(brms)
library(reshape2)
library(bayesplot)
library(ggplot2)
library(igraph)
library(ggraph)


data <- read_csv("./long_format_data.csv")

colnames(data)[colnames(data) == "abundance"] <- "Abundance"
colnames(data)[colnames(data) == "taxa"] <- "Taxa"

data$Abundance <- data$Abundance / max(data$Abundance, na.rm = TRUE)


data <- na.omit(data)



data$Abundance <- data$Abundance / max(data$Abundance, na.rm = TRUE)


data$Abundance <- ifelse(data$Abundance >= 1, 0.999, data$Abundance)
data$Abundance <- ifelse(data$Abundance <= 0, 0.001, data$Abundance)


print(range(data$Abundance, na.rm = TRUE))

library(brms)

model <- brm(
  bf(Abundance ~ Birthmode + Breastmilk + Antibiotic + (1|Taxa)),
  data = data,
  family = zero_inflated_beta(),
  chains = 4,          
  iter = 8000,         
  warmup = 2000,       
  cores = 2            
)


# Save RDS file
saveRDS(model, file = "./fitted_model.rds")

# Load RDS file
model <- readRDS("./fitted_model.rds")




# summary  text
capture.output(summary(model), file = "./model_summary.txt")




#  plots to PDF
pdf("./model_plots.pdf")
plot(model)
dev.off()




#  marginal effects
marginal_effects_data <- marginal_effects(model)

plot_marginal_effects <- plot(marginal_effects_data, ask = FALSE)

ggsave("./marginal_effects_plot.pdf", plot = plot_marginal_effects[[1]], width = 10, height = 8, dpi = 300)




interaction_graph <- graph_from_data_frame(edges_data, directed = FALSE)


E(interaction_graph)$weight <- abs(E(interaction_graph)$Weight)

# Plot the network
network_plot <- ggraph(interaction_graph, layout = "fr") +
  geom_edge_link(aes(edge_alpha = weight, edge_width = weight), color = "blue") +
  geom_node_point(color = "skyblue", size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, max.overlaps = 50) +
  ggtitle("Bacterial Interaction Network") +
  theme_void()


ggsave("./network_analysis_plot.pdf", plot = network_plot, width = 12, height = 8, dpi = 300)
