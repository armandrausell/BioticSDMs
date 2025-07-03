


library(rentrez)
library(stringr)
library(dplyr)
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(ape)
#This script retrieves genetic data for a set of species, aligns the sequences and
#computes a genetic distance matrix, which would be later used to check if competitive
#exclusion is related to genetic distance

# Function to fetch and clean a COI sequence for one species
get_coi_sequence <- function(species_name) {
  query <- paste0(species_name, "[ORGN] AND (COI[All Fields] OR CO1[All Fields] OR COX1[All Fields])")
  
  search_results <- entrez_search(db = "nucleotide", term = query, retmax = 1)
  
  if (length(search_results$ids) == 0) {
    return(data.frame(species = species_name, sequence = NA, length = 0))
  }
  
  raw_fasta <- entrez_fetch(db = "nucleotide", id = search_results$ids[1], rettype = "fasta", retmode = "text")
  
  # Remove FASTA header and collapse sequence lines
  lines <- strsplit(raw_fasta, "\n")[[1]]
  seq_lines <- lines[!grepl("^>", lines)]
  sequence <- paste(seq_lines, collapse = "")
  
  return(data.frame(species = species_name, sequence = sequence, length = nchar(sequence)))
}

# Apply to all species
results <- lapply(names_gbif_unf, get_coi_sequence)

# Combine into a dataframe
coi_df <- bind_rows(results)

# View the result
print(coi_df)

library(Biostrings)
library(DECIPHER)
library(ape)
# Remove rows with missing sequences
coi_clean <- coi_df[!is.na(coi_df$sequence) & coi_df$length > 0, ]

# Convert to DNAStringSet with species names as names
dna_set <- DNAStringSet(coi_clean$sequence)
names(dna_set) <- coi_clean$species

alignment <- AlignSeqs(dna_set, anchor = NA)
alignment_trimmed <- subseq(alignment, start = 9990, end = 10603)

# View alignment interactively in browser
BrowseSeqs(alignment)
BrowseSeqs(alignment_trimmed)

# Save and re-load as DNAbin for compatibility
writeXStringSet(alignment_trimmed, filepath = "coi_alignment_trimmed.fasta")
aligned <- read.nexus.data("sequences/Armand_Cox1.nex")
dist_matrix <- dist.dna(aligned, model = "raw", pairwise.deletion = TRUE)

# View the matrix
c<-as.matrix(dist_matrix)

saveRDS(names_gbif_unf, "names_species.RData")
write.table(names_gbif_unf,"names_species.txt", sep=";",col.names = FALSE, row.names = FALSE)

########
library(ape)
#Read nexus from annie
aligned <- read.nexus.data("sequences/Armand_Cox1.nex")

# Assuming your Nexus data is stored in a variable named 'aligned'

# Extract species names
species_names <- names(aligned)

species_names_clean <- gsub("_", " ", species_names)
missing_names <- setdiff(names_gbif_unf, species_names_clean)
print(missing_names)


# Create a data frame
df <- data.frame(
  species = species_names,
  sequence = sapply(aligned, function(seq) paste(seq, collapse = "")),
  stringsAsFactors = FALSE
)

# Add a column counting non-gap characters
df$num_non_gaps <- sapply(df$sequence, function(seq) sum(strsplit(seq, "")[[1]] != "-"))

# View the data frame
print(df)

update_sequence <- function(df, species_name, raw_sequence, add_if_missing = FALSE) {
  # Remove numbers, spaces, and any non-ATCGN characters (e.g., dashes, line breaks)
  cleaned_sequence <- gsub("[^ATCGNatcgn-]", "", raw_sequence)
  cleaned_sequence <- toupper(cleaned_sequence)
  non_gap_count <- sum(strsplit(cleaned_sequence, "")[[1]] != "-")
  
  # Check if species exists
  if (species_name %in% df$species) {
    df[df$species == species_name, "sequence"] <- cleaned_sequence
    df[df$species == species_name, "num_non_gaps"] <- non_gap_count
  } else if (add_if_missing) {
    # Add new row
    df <- rbind(df, data.frame(
      species = species_name,
      sequence = cleaned_sequence,
      num_non_gaps = non_gap_count,
      stringsAsFactors = FALSE
    ))
  } else {
    warning(paste("Species", species_name, "not found in DataFrame and not added."))
  }
  
  return(df)
}

#Fill the sequence to be included
raw_seq<-"gtgtca cttattcgtt gattattttc aacaaaccac aaagacatcg gcaccctata
     5281 cttattattt ggcgcctgag ccggtatagc cggcacagct cttagccttc tgatccgaac
     5341 agagctaagt caacccggca cccttcttgg ggacgaccag gtatataacg tggttgtcac
     5401 agctcatgct ttcgttataa ttttcttctt agttatacct gtaataattg gcgggtttgg
     5461 gaactgactt gtcccattaa taattggtgc ccctgacata gcatttccac gaataaataa
     5521 cataagcttt tgactccttc ccccatctct tcttctgctt ctatcttctt ctggaattga
     5581 agctggtgcc gggaccggtt gaactgtcta ccccccccta gccggaaatc ttgcccacgc
     5641 aggggcatca gtcgatctaa ctattttttc acttcactta gctggagttt cttcgatttt
     5701 aggcgcaatt aattttatca ccacctgcat caacataaaa ccccccaaca taacacaata
     5761 tcaaacccct ttatttgttt gatccgtctt gattacagcc gtattactat tgctctctct
     5821 ccccgtacta gccgcaggca ttacaatgct actaacagac cgcaatctaa acacatcatt
     5881 ttttgacccc gcgggagggg gagacccgat cctctaccaa cacttatttt gattctttgg
     5941 gcaccctgaa gtctacattc ttatcctccc aggttttggc ataatttctc atattgttac
     6001 atactatgca ggaaaaaaag aacccttcgg ctacatagga atagtctgag ccataatgtc
     6061 aattgggttt ttaggcttca tcgtatgagc tcatcatata tttaccgtag gaatggatgt
     6121 tgacacccgg gcctacttta catcagctac aataattatt gctattccca caggggtaaa
     6181 agtctttagc tgacttgcaa ctcttcatgg cggaactatt aaatgagacg cagctatact
     6241 ttgggctcta ggctttatct tcctgtttac tgttgggggt ctaacaggca ttattctagc
     6301 caactcctca ttagatattg tccttcatga tacatattac gtagttgccc acttccacta
     6361 tgttttgtcc ataggagctg tctttgccat tataggcgga tttgttcatt gattcccact
     6421 ctttacgggt ttcaccctac atagctcatg aacaaaagct caatttggtg ttatattcac
     6481 tggcgtcaat ataacattct ttcctcaaca cttcctgggt ttagctggca taccccgacg
     6541 ctactctgac tacccagatg catatactct ttgaaactcc atttcatcga ttggctccct
     6601 aatctcatta acagctgtaa ttataataat atttattatt tgagaagccc tagcagctaa
     6661 acgcgaagta cttaccctcg aacttactag cactaatcta gagtgacttc acggctgccc
     6721 acctccatac cacacctacg aagaagcaac ccatgtacaa acctcaaggg"

df <- update_sequence(df, "Lacerta_bilineata", raw_seq, add_if_missing = F)

write.csv(df, "FilledDNASequences.csv", row.names = FALSE)
df <- read.csv("FilledDNASequences.csv", stringsAsFactors = FALSE)

df1<-df
df1$sequence_clean <- gsub("-", "", df1$sequence)

dna_set <- DNAStringSet(df1$sequence_clean)
names(dna_set) <- df1$species
alignment <- AlignSeqs(dna_set, anchor = NA)

# Save and re-load as DNAbin for compatibility
writeXStringSet(alignment, filepath = "sequences/coi_alignment_filled.fasta")
aligned <- read.FASTA("sequences/coi_alignment_filled.fasta")
dist_matrix <- dist.dna(aligned, model = "raw", pairwise.deletion = TRUE)

c<-as.matrix(dist_matrix)

aligned <- readDNAStringSet("sequences/coi_alignment_filled.fasta")
tree_nj <- nj(dist_matrix)
plot(tree_nj, main = "Neighbor-Joining Tree (ape)", cex = 0.5)
