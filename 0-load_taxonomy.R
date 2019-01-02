require(dplyr)

# Taxonomy from he file
taxonomy_file <- read.csv("OTU taxa information.txt", sep='\t', header=T, row.names=1)

# Get the taxonomy
taxonomy <- observation_metadata(otus_table_raw)

# Samples with no nulls 
tax_sums <- apply(otus_initial, 1, sum)

# Filter out zeroes from taxonomies
taxonomy <- taxonomy[tax_sums > 0]

# Get the names
names_of_samples <- names(taxonomy)

# Metod to process the names
process_name <- function(name) {
  if (name == "Unassigned") {
    return(NA)
  } else {
    # Delete not needed values
    name <- gsub("\"", "", name)
    # Find the start
    index <- regexpr("__", name)
    # Cut the substring
    return(substring(name, index + 2))
  }
}

# Get the taxes
new_taxes <- lapply(taxonomy, sapply, process_name)
new_taxes <- lapply(new_taxes, unname)
raw_taxes <- lapply(taxonomy, unname)

# Splitted taxes
splitted_new_taxes <- lapply(taxonomy, sapply, process_name)
splitted_new_taxes <- unname(lapply(splitted_new_taxes, unname))
splitted_new_taxes <- lapply(splitted_new_taxes, function(x) { return(x[x != ""])})

# Get maxnumber of taxes
splitted_new_taxes_max <- max(sapply(splitted_new_taxes, length))

# Method that gets the taxonomy
get_taxonomies <- function(vector_of_taxes, number_of_taxes) {
  
  # Initialize vector
  trial_vector <- vector(mode = "character", length = number_of_taxes)
  # Vector length
  taxes_number <- length(vector_of_taxes)
  # Assign
  vector_of_taxes[is.na(vector_of_taxes)] <- "Unassigned"
  trial_vector[1:taxes_number] <- vector_of_taxes
  if (taxes_number != 14) {
    trial_vector[(taxes_number + 1) : number_of_taxes] <- "Unassigned"
  }
  
  return(trial_vector)
}

# Method gets the last taxa
get_last_taxonomy <- function(vector_of_taxas, number_of_taxes) {
  # Check if the first taxonomy is unassigned
  if (vector_of_taxas[1] == "Unassigned") {
    return("Unassigned")
  } else {
    # If the vector has the last element
    if (vector_of_taxas[number_of_taxes] != "Unassigned") {
      return(vector_of_taxas[number_of_taxes])
    }
    # Find the index of the first unassigned
    ind_of_first_na <- min(which(vector_of_taxas == "Unassigned"))
    # Return the previous to unassigned
    return(vector_of_taxas[ind_of_first_na - 1])
  }
}

# Standartize
splitted_new_taxes <- lapply(splitted_new_taxes, get_taxonomies, splitted_new_taxes_max)

# Colnames for splitted df
splitted_colnames <- sapply(seq(splitted_new_taxes_max), function(x) {return(paste0("tax_", x))})
# Make a matrix from a list of vectors
splitted_matrix <- matrix(unlist(splitted_new_taxes), byrow = TRUE, ncol = splitted_new_taxes_max)
# Last non-NA taxes
splitted_last_nonna <- unlist(lapply(splitted_new_taxes, get_last_taxonomy, splitted_new_taxes_max))
# Make a dataframe
splitted_df <- data.frame(samples = names_of_samples)

# Add splitted colnames
for (i in seq(splitted_colnames)) {
  splitted_df[[splitted_colnames[i]]] <- splitted_matrix[, i]
}
# Add the last tax colname
splitted_df[["last_tax"]] <- unname(unlist(splitted_last_nonna))

# Get taxes as a string
new_taxes_merged <- unname(lapply(new_taxes, paste, collapse = ", "))
raw_taxes_merged <- unname(lapply(raw_taxes, paste, collapse = ", "))

# Processed taxes dataframe
processed_df <- data.frame(samples = names_of_samples, tax = unlist(new_taxes_merged))
raw_taxes_df <- data.frame(samples = names_of_samples, tax = unlist(raw_taxes_merged))