install.packages("factoextra")
library(factoextra)

# Read in soil data 
soil_data <- read.csv("/Users/kathleenmclay/google_drive/PhD/Chapter_2_still_inversions/2.1_datasets/soil_data_normalised.csv", header=TRUE)
row.names(soil_data) <- soil_data$Population
# PCA 
soil_pca <- prcomp(soil_data[,3:41])

# the first three PCs explain 80% of variation 
PCA_out <- summary(soil_pca)
write.csv(PCA_out, "/Users/kathleenmclay/google_drive/PhD/Chapter_2_still_inversions/2.1_datasets/GEA/PCA_soil_data/PCA_out.txt", sep = "/t")

#extract the loadings 
PCA_loadings <- soil_pca$rotation
write.csv(PCA_loadings, "/Users/kathleenmclay/google_drive/PhD/Chapter_2_still_inversions/2.1_datasets/GEA/PCA_soil_data/PCA_loadings.csv")

#extract the PCA co-ord for each variable for first three PCs 
PCA_data <- soil_pca$x[,1:3]
write.csv(PCA_data, "/Users/kathleenmclay/google_drive/PhD/Chapter_2_still_inversions/2.1_datasets/GEA/PCA_soil_data/PCA_soil_data.csv")

# plot the % variation explained by each eigenvector (PC axis)
fviz_eig(soil_pca)

# plot samples on PCA 
fviz_pca_ind(soil_pca,
             repel = TRUE
)

# plot loadings for PCA 
fviz_pca_var(soil_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# combine them in a biplot 
fviz_pca_biplot(soil_pca, repel = TRUE,
                col.var = "contrib", # Color by contributions to the PC
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                col.ind = "black"  # Individuals color
)