

setwd("~/Desktop/Github/BioticSDMs")


#Load the list from the atlas, these are the ones with atlas data as well.
save(names_sps, file = "Atlas_data/names_atlas_sps.RData")
load("Atlas_data/names_atlas_sps.RData")

# List of vertebrates according to the ministry, it includes many poorly represented sp
vert<-read.csv2("metaweb/Listado_vertebrados_terrestres_completo.csv",)

# Extract first two words from Sinónimo (filtering out empty strings)
sinonimos <- unique(trimws(vert$Sinónimo))
sinonimos <- sinonimos[sinonimos != ""]  # Remove empty strings
sinonimos_first2 <- unique(sapply(strsplit(sinonimos, "\\s+"), function(x) paste(head(x, 2), collapse = " ")))

# Get unique values from WithoutAutorship
without_authorship <- unique(trimws(vert$WithoutAutorship))

# Combine and get final unique list
all_names <- unique(c(sinonimos_first2, without_authorship))

# View result
all_names


# Keep only the first two words of each element
all_names <- sapply(strsplit(all_names, "\\s+"), function(x) paste(head(x, 2), collapse = " "))

# Get unique values
all_names <- unique(all_names)

library(rgbif)

# Initialize an empty data frame to store results
species_observation_count <- data.frame(species = character(), count = integer(), stringsAsFactors = FALSE)

# Loop over each species
for (species in all_names) {
  
  # Get count
  count_es <- occ_count(scientificName = species, country = "ES;PT", year = "2000,2020")

  # Total count
  total_count <- count_es
  
  # Append to the results data frame
  species_observation_count <- rbind(
    species_observation_count,
    data.frame(species = species, count = total_count, stringsAsFactors = FALSE)
  )
}

# Print the result
print(species_observation_count)

# Filter out species with fewer than 100 records
filtered_species_df <- species_observation_count[species_observation_count$count >= 200, ]

all_names_filtered<-unique(filtered_species_df$species)


diff_1 <- setdiff(all_names_filtered, names_sps)

# Filter filtered_species_df to only include species in diff_1
filtered_diff_df <- filtered_species_df[filtered_species_df$species %in% diff_1, ]

#This is a more complete list of species that are abundant in iberia.
all_sp<-c(names_sps,diff_1)
save(all_sp, file = "Atlas_data/Complete_sp_list_iberia.RData")


# View the result
print(filtered_diff_df)

