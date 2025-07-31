install.packages("factoextra")
library(factoextra)

# Read in soil data 
soil_data <- read.csv("/.../soil_data_normalised.csv", header=TRUE)
row.names(soil_data) <- soil_data$Population
#remove populations D32 and H12
soil_data <- soil_data[!soil_data$Population %in% c("D32", "H12"), ]

# PCA 
soil_pca <- prcomp(soil_data[,3:41])
pca_summary <- summary(soil_pca)

# Create a data frame with the key info: Proportion of Variance and Cumulative Proportion
pc_variance_df <- data.frame(
  PC = paste0("PC", seq_along(pca_summary$importance[2, ])),
  Proportion_Variance = pca_summary$importance[2, ],
  Cumulative_Variance = pca_summary$importance[3, ]
)

# Write to CSV file
write.csv(pc_variance_df, file = "/.../soil_pca_variance_summary.csv", row.names = FALSE)

#extract the PCA co-ord for each variable for first three PCs 
PCA_out <- data.frame(PC1 = soil_pca$x[,1], PC2 = soil_pca$x[,2], PC3 = soil_pca$x[,3], 
                      Population = soil_data$Population, Ecotype = soil_data$Ecotype)

write.csv(PCA_out, "/.../pop_PCA_coords.csv")

#transpose 
PCA_data <- t(PCA_out)
write.csv(PCA_data, "/.../PCA_soil_data_names.csv")

#extract the loadings 
PCA_loadings <- soil_pca$rotation
write.csv(PCA_loadings, "/.../soil_PCA/PCA_loadings.csv")

# plot the % variation explained by each eigenvector (PC axis)
fviz_eig(soil_pca)

# plot loadings for PCA 
fviz_pca_var(soil_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# Create the biplot and save it as a PNG file
png(filename = "/.../soil_pca_biplot_PC2PC3.png", width = 1200, height = 1000, res = 150)

fviz_pca_biplot(soil_pca,
                repel = TRUE,
                axes = c(1, 2),
                col.var = "contrib",  # Color variables by contributions
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                col.ind = "black"     # Color individuals in black
)

fviz_pca_biplot(soil_pca,
                repel = TRUE,
                axes = c(1, 2),
                col.ind = "black",          # Color individuals in black
                geom.ind = "point",         # Show individuals as points
                pointsize = 4,              # Larger dots
                col.var = "contrib",     # Color variables by contribution
) +
  scale_color_viridis_c(option = "D", end = 0.9)


fviz_pca_biplot(soil_pca,
                repel = TRUE,
                axes = c(1, 2),
                col.ind = "black",                 # Color individuals in black
                geom.ind = c("point", "text"),     # Show both points and labels
                pointsize = 4,                     # Larger dots
                labelsize = 4,                     # Larger labels
                label.fontface = "bold",           # Bold text
                col.var = "contrib"                # Color variables by contribution
) +
  scale_color_viridis_c(option = "D", end = 0.9) +
  theme(panel.grid = element_blank())

dev.off()
