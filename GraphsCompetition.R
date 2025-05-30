

# Install if necessary
install.packages("paletteer")
library(paletteer)
library(ggplot2)

gge <- ggplot(Opt_competitor_per_species_potential, 
              aes(x = Corr_Act_Biserial, 
                  y = Correlation_Richness, 
                  color = Num_Competitors, 
                  label = Species)) +
  geom_point(size = 1.5, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradientn(
    colors = as.character(rev(paletteer::paletteer_c("ggthemes::Classic Red-Blue", 30))),
    name = "Num. Competitors"
  ) +
  ggrepel::geom_text_repel(size = 1.5, color = "black") +
  labs(
    title = "Effect of Competitive Pressure vs. Richness and Actual Distribution",
    x = "Biserial Correlation (Actual)",
    y = "Body Mass (g)"
  ) +
  facet_wrap(~ Class,scales = "free") +
  theme_linedraw() +
  theme(legend.position = "right")

print(gge)

ggsave("Figures/Correlation_potential_facet_bm.jpg", plot = gge, width = 10, height = 6, dpi = 600)


ggm <- ggplot(Opt_competitor_per_species_potential, 
              aes(x = TrophicPosition, 
                  y = BodyMass_g, 
                  color = Correlation_Richness, 
                  label = Species)) +
  geom_point(size = 1.5, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradientn(
    colors = as.character(rev(paletteer::paletteer_c("ggthemes::Classic Red-Blue", 30))),
    name = "Correlation Richness"
  ) +
  ggrepel::geom_text_repel(size = 1.5, color = "black") +
  labs(
    title = "Effect of Competitive Pressure vs. Richness and Actual Distribution",
    x = "Trophic Position",
    y = "Body Mass (g)"
  ) +
  facet_wrap(~ Class,scales = "free") +
  theme_linedraw() +
  theme(legend.position = "right")

print(ggm)
ggsave("Figures/Correlation_potential_tp_bm_facet.jpg", plot = ggm, width = 10, height = 6, dpi = 600)



gge <- ggplot(Opt_competitor_per_species_potential, 
              aes(x = Corr_Act_Biserial, 
                  y = Corr_hfp, 
                  color = TrophicPosition, 
                  label = Species)) +
  geom_point(size = 1.5, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_gradientn(
    colors = as.character(rev(paletteer::paletteer_c("ggthemes::Classic Red-Blue", 30))),
    name = "Num. Competitors"
  ) +
  ggrepel::geom_text_repel(size = 1.5, color = "black") +
  labs(
    title = "Effect of Competitive Pressure vs. Richness and Actual Distribution",
    x = "Biserial Correlation (Actual)",
    y = "Human Footprint"
  ) +
  facet_wrap(~ Class,scales = "free") +
  theme_linedraw() +
  theme(legend.position = "right")

print(gge)

ggsave("Figures/Correlation_potential_facet_bm.jpg", plot = gge, width = 10, height = 6, dpi = 600)


######
library(ggplot2)
library(dplyr)
library(FactoMineR)
library(factoextra)

# Select and scale only numeric columns for PCA
pca_data_scaled <- Opt_competitor_per_species_potential[,c(2:6)]%>%scale()

# Proceed with PCA
pca_result <- prcomp(pca_data_scaled, center = TRUE, scale. = TRUE)

# Combine with metadata for plotting
pca_df <- as.data.frame(pca_result$x)
pca_df$Species <- Opt_competitor_per_species_potential$Species
pca_df$Class <- Opt_competitor_per_species_potential$Class
# Extract loadings (rotation matrix) and rescale for arrows
loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$Variable <- rownames(loadings)

# Optionally scale arrows for visualization
arrow_scale <- 2  # Increase if arrows are too short
loadings <- loadings %>%
  mutate(PC1 = PC1 * arrow_scale,
         PC2 = PC2 * arrow_scale)

# Base PCA plot with points and labels
ggplot(pca_df, aes(x = PC1, y = PC2, color = Class, label = Species)) +
  geom_point(size = 2, alpha = 0.8) +
  ggrepel::geom_text_repel(size = 2, max.overlaps = 10) +
  
  # Add arrows
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", inherit.aes = FALSE) +
  geom_text(data = loadings, aes(x = PC1, y = PC2, label = Variable),
            color = "black", vjust = -1, inherit.aes = FALSE, size = 3) +
  
  # Formatting
  theme_minimal() +
  labs(
    title = "PCA of Competitive and Trophic Traits",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "% variance)")
  )
