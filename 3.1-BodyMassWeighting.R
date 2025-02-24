


#Body sizes from https://github.com/fernandacaron/body_size_evol
bm.amp<-read.csv2("metaweb/body_size/BodySizeAmphibia_RMA_17jan24.csv",sep = ",")
bm.mam<-read.csv2("metaweb/body_size/BodySizeMammalia_09set21.csv",sep = ",")
bm.aves<-read.csv2("metaweb/body_size/BodySizeAves_10set22.csv",sep = ",")
bm.rep<-read.csv2("metaweb/body_size/BodySizeReptilia_15set21.csv",sep = ",")

#Subsample interesting columns
bm.amp<-bm.amp[,c(1,5,8)]
bm.mam<-bm.mam[,c(1,5,10)]
bm.aves<-bm.aves[,c(2,6,11)]
bm.rep<-bm.rep[,c(1,6,10)]

# Rename columns to ensure uniformity
colnames(bm.amp) <- c("Class", "Species", "BodyMass_g")
colnames(bm.mam) <- c("Class", "Species", "BodyMass_g")
colnames(bm.aves) <- c("Class", "Species", "BodyMass_g")
colnames(bm.rep) <- c("Class", "Species", "BodyMass_g")

# Combine all datasets into one
bm_combined <- bind_rows(bm.amp, bm.mam, bm.aves, bm.rep)

# Remove underscores in the Species column
bm_combined <- bm_combined %>%
  mutate(Species = gsub("_", " ", Species))  # Replace "_" with space


# Join bodymass_df to assign mass to Predator_A
Competition_weighted <- Competition_spain_df %>%
  left_join(bm_combined, by = c("Predator_A" = "Species")) %>%
  rename(BodyMass_A = BodyMass_g)%>%  # Rename to avoid conflict
  rename(Class_A = Class)

# Join again for Predator_B
Competition_weighted <- Competition_weighted %>%
  left_join(bm_combined, by = c("Predator_B" = "Species")) %>%
  rename(BodyMass_B = BodyMass_g) %>%  # Rename to avoid conflict
  rename(Class_B = Class)

# Convert new columns to numeric (if not already)
Competition_weighted <- Competition_weighted %>%
  mutate(BodyMass_A = as.numeric(BodyMass_A),
         BodyMass_B = as.numeric(BodyMass_B))

#sigmoid <- function(x) {
#  1 / (1 + exp(-x))  # Sigmoid function
#}

#Competition_weighted <- Competition_weighted %>%
#  mutate(
#    SigmoidRatioAtoB = sigmoid(BodyMass_A / BodyMass_B - 1),  
#    SigmoidRatioBtoA = sigmoid(BodyMass_B / BodyMass_A - 1)
 # )

#Calculate the amount of pressure competitors exert over the others
#Competition_weighted <- Competition_weighted %>%
 # mutate(
  #  PressureBtoA = round(
   #   (alphaA * Proportion_A_Shares + (1 - alphaA) * RealizedDietOverlapA) * SigmoidRatioBtoA, 3
    #),
    #PressureAtoB = round(
    #  (alphaB * Proportion_B_Shares + (1 - alphaB) * RealizedDietOverlapB) * SigmoidRatioAtoB, 3
  #)
    #)


# Define Gaussian function
gaussian_similarity <- function(mass_A, mass_B, sigma = 1.8) {
  exp(-((log(mass_A) - log(mass_B))^2 / (2 * sigma^2)))  # Gaussian decay
}

# Apply Gaussian weighting to favor similar-sized competitors
Competition_weighted <- Competition_weighted %>%
  mutate(
    MassSimilarity = gaussian_similarity(BodyMass_A, BodyMass_B, sigma = 1.8),  # Adjust sigma as needed
    PressureBtoA = round(
      (alphaA * Proportion_A_Shares + (1 - alphaA) * RealizedDietOverlapA) * MassSimilarity, 3
    ),
    PressureAtoB = round(
      (alphaB * Proportion_B_Shares + (1 - alphaB) * RealizedDietOverlapB) * MassSimilarity, 3
    )
  )


