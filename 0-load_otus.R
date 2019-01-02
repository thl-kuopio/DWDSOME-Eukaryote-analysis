# Import libraries
library("biomformat")
library("vegan")

# Load 
otus_table_raw <- read_biom("otu.table.json.fmt.biom")
mapping_raw <- read.csv("map_new_JI mod.txt", sep='\t', header=T, row.names=1)

# OTUs initial
otus_initial <- as.matrix(biom_data(otus_table_raw))
# Manipulate
otus_preprocessed <- t(as.matrix(biom_data(otus_table_raw)))

# Find same ids
same_ids <- intersect(rownames(otus_preprocessed), rownames(mapping_raw))

# Get the required samples
otus <- otus_preprocessed[same_ids,]
mapping <- mapping_raw[same_ids,]

# Mapping colnames
column_names <- c(
  "coding",
  "template",
  "waterworks",
  "Cu",
  "Al",
  "Fe",
  "Mn",
  "Turbidity",
  "Absorbance",
  "Absorbance.420nm",
  "AOC_unknwn",
  "AOC",
  "MAP",
  "HPC",
  "SybrGreen",
  "DAPI",
  "Temperature",
  "FreeChl",
  "Chlorine",
  "EC",
  "pH"
)
# Change colnames
colnames(mapping) <- column_names

# Initial OTUs rownames
otus_rownames <- rownames(otus_preprocessed)

# Method gets samples by mapping parameters
get_samples_by <- function(mapping, template = NULL, waterworks = NULL) {
  # Reinitialize the mapping
  mapping_result <- mapping
  # Filter by the template
  if (!is.null(template)) {
    mapping_result <- mapping_result[mapping_result[["template"]] == template, ]
  }
  # Filter by the field
  if (!is.null(waterworks)) {
    mapping_result <- mapping_result[mapping_result[["waterworks"]] == waterworks, ]
  }
  
  return(mapping_result)
}