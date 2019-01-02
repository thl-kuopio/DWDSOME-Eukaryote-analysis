require(ggplot2)
require(purrr)

# Get the samples
current_mapping <- get_samples_by(mapping = mapping, template = "DNA")

# Get all otus
current_otus <- otus[rownames(current_mapping), ]
# Make cca
current_cca <- cca(current_otus ~ Fe + Turbidity + Absorbance + AOC + MAP + HPC + DAPI + Temperature + Chlorine + EC + pH, data=current_mapping, na.action=na.exclude)

# Species
species_all <- scores(current_cca, scaling = 2)$species
# Sites
sites_all <- scores(current_cca, scaling = 2)$sites

# Get summary for explained variances
cca_summary <- summary(current_cca)$concont$importance

# Get biplot
biplot_all <- current_cca$CCA$biplot

# CCA species
cca_species <- rownames(species_all)

# Make a dataframe for sites
data_species <- data.frame(species_all)
data_species[["Genera"]] <- taxonomy_file[rownames(taxonomy_file) %in% rownames(species_all),]$Genera
data_species[["Kingdom"]] <- taxonomy_file[rownames(taxonomy_file) %in% rownames(species_all),]$Kingdom_class

# COloring for sites
group_colors <- c("orange2", "blue3", "red3", "green3", "purple3")
group_names <- c("A","B","C","D","E")

# Data sites
data_sites <- data.frame(sites_all)
data_sites[["waterworks"]] <- current_mapping$waterworks

# ----------------------------------- FILTERING OUT FORBIDDEN TAXES ----------------------------------------

# To filter out
if (TRUE) {
  # The list of genera's to exclude
  forbidden_generas <- c("Unassigned","Metazoa")
  data_species <- data_species %>% filter(! Kingdom %in% forbidden_generas)
}

species_to_plot <- data_species %>% group_by(Kingdom) %>% summarise_all(funs(mean(., na.rm=TRUE)))

# ----------------------------------------------------------------------------------------------------------

# Scaling coefficient (visualization of the species)
multi_coef <- 1.5

# Initialize the machine
pdf(paste0("dna", ".pdf"))

if (TRUE) {
  # Add all the points to plots
  p <- ggplot(data = data_sites) +
    geom_point(mapping = aes(x = CCA1, y = CCA2, color = waterworks)) +
    labs(x = paste0("CCA1 (", round(cca_summary[2,1]*100, digits = 2), "%)"), y = paste0("CCA2 (", round(cca_summary[2,2]*100, digits = 2), "%)")) +
    scale_color_manual(values = setNames(group_colors, group_names), name = "Waterworks:")
} else {
  p <- ggplot() +
    labs(x = paste0("CCA1 (", round(cca_summary[2,1]*100, digits = 2), "%)"), y = paste0("CCA2 (", round(cca_summary[2,2]*100, digits = 2), "%)"))
}



# Remove legend?
if (FALSE) {
  p <- p + theme(legend.position = "none")
}

# Add species?
if (FALSE) {
  p <- p + geom_point(data = data_species, mapping = aes(x = CCA1, y = CCA2, color = Genera))
}

# Find point with max length
max_point_length <- max(apply(data_sites[c("CCA1", "CCA2")], 1, function(x) {return(sqrt(x[1]*x[1] + x[2]*x[2]))}))
max_arrow_length <- max(apply(biplot_all[,c("CCA1", "CCA2")], 1, function(x) {return(sqrt(x[1]*x[1] + x[2]*x[2]))}))

# Arrow multiplier
arrow_multi <- (max_point_length / max_arrow_length) * 0.8
# Change arrows length
arrows_all <- data.frame(biplot_all)

# Improve the arrors length
arrows_all$CCA1 <- arrows_all$CCA1 * arrow_multi
arrows_all$CCA2 <- arrows_all$CCA2 * arrow_multi

# Add arrows
for (i in seq(nrow(biplot_all))) {
  p <- p + geom_segment(data=arrows_all, mapping=aes(x=0, y=0, xend=CCA1, yend=CCA2), arrow=arrow(type = "closed", angle = 15, length = unit(2, "mm")), size=0.1, color="black")
}

# Varialbes coordinates
labels_coordinates <- data.frame(name = as.character(rownames(biplot_all)), x = unname(arrows_all[,1]), y = unname(arrows_all[,2]))
rownames(labels_coordinates) <- rownames(biplot_all)


# Manually setting the coordinates
labels_coordinates[labels_coordinates$name == "Fe", "x"] <- -0.0 + labels_coordinates[labels_coordinates$name == "Fe", "x"]
labels_coordinates[labels_coordinates$name == "Fe", "y"] <- 0.1 + labels_coordinates[labels_coordinates$name == "Fe", "y"]

labels_coordinates[labels_coordinates$name == "Turbidity", "x"] <- -0.5 + labels_coordinates[labels_coordinates$name == "Turbidity", "x"]
labels_coordinates[labels_coordinates$name == "Turbidity", "y"] <- labels_coordinates[labels_coordinates$name == "Turbidity", "y"]

labels_coordinates[labels_coordinates$name == "Absorbance", "x"] <- 0.0 + labels_coordinates[labels_coordinates$name == "Absorbance", "x"]
labels_coordinates[labels_coordinates$name == "Absorbance", "y"] <- -0.1 + labels_coordinates[labels_coordinates$name == "Absorbance", "y"]

labels_coordinates[labels_coordinates$name == "AOC", "x"] <- 0.3 + labels_coordinates[labels_coordinates$name == "AOC", "x"]
labels_coordinates[labels_coordinates$name == "AOC", "y"] <- labels_coordinates[labels_coordinates$name == "AOC", "y"]

labels_coordinates[labels_coordinates$name == "MAP", "x"] <- 0.1 + labels_coordinates[labels_coordinates$name == "MAP", "x"]
labels_coordinates[labels_coordinates$name == "MAP", "y"] <- -0.2 + labels_coordinates[labels_coordinates$name == "MAP", "y"]

labels_coordinates[labels_coordinates$name == "HPC", "x"] <- 0.5 + labels_coordinates[labels_coordinates$name == "HPC", "x"]
labels_coordinates[labels_coordinates$name == "HPC", "y"] <- 0.1 + labels_coordinates[labels_coordinates$name == "HPC", "y"]

labels_coordinates[labels_coordinates$name == "DAPI", "x"] <- 0.5 + labels_coordinates[labels_coordinates$name == "DAPI", "x"]
labels_coordinates[labels_coordinates$name == "DAPI", "y"] <- labels_coordinates[labels_coordinates$name == "DAPI", "y"]

labels_coordinates[labels_coordinates$name == "Temperature", "x"] <- -1.2 + labels_coordinates[labels_coordinates$name == "Temperature", "x"]
labels_coordinates[labels_coordinates$name == "Temperature", "y"] <- 1 + labels_coordinates[labels_coordinates$name == "Temperature", "y"]

# Add segment for label
if (TRUE) {
  p <- p + geom_segment(data=arrows_all, mapping=aes(x=arrows_all["Temperature", "CCA1"], y=arrows_all["Temperature", "CCA2"], xend=labels_coordinates["Temperature", "x"], yend=(labels_coordinates["Temperature", "y"])- 0.05), size=0.1, color="darkgrey", linetype = "dashed")
}


labels_coordinates[labels_coordinates$name == "Chlorine", "x"] <- 0.1 + labels_coordinates[labels_coordinates$name == "Chlorine", "x"]
labels_coordinates[labels_coordinates$name == "Chlorine", "y"] <- 0.2 + labels_coordinates[labels_coordinates$name == "Chlorine", "y"]

labels_coordinates[labels_coordinates$name == "EC", "x"] <- -0.2 + labels_coordinates[labels_coordinates$name == "EC", "x"]
labels_coordinates[labels_coordinates$name == "EC", "y"] <- labels_coordinates[labels_coordinates$name == "EC", "y"]

labels_coordinates[labels_coordinates$name == "pH", "x"] <- labels_coordinates[labels_coordinates$name == "pH", "x"]
labels_coordinates[labels_coordinates$name == "pH", "y"] <- -0.1 + labels_coordinates[labels_coordinates$name == "pH", "y"]

# Add annotations
for (i in seq(nrow(labels_coordinates))) {
  p <- p + annotate("text", x = labels_coordinates$x[i], y = labels_coordinates$y[i], label = labels_coordinates$name[i])
}

# Plot species ----------------------------------------------------------------------------------------

if (TRUE) {
  
  # Find point with max length
  max_point_length <- max(apply(data_sites[c("CCA1", "CCA2")], 1, function(x) {return(sqrt(x[1]*x[1] + x[2]*x[2]))}))
  max_species_length <- max(apply(species_to_plot[,c("CCA1", "CCA2")], 1, function(x) {return(sqrt(x[1]*x[1] + x[2]*x[2]))}))
  
  # Arrow multiplier
  species_multi <- (max_point_length / max_species_length) * multi_coef
  
  # Make a species DF
  species_annotation <- data.frame(x = species_to_plot$CCA1 * species_multi,
                                   y = species_to_plot$CCA2 * species_multi,
                                   text = species_to_plot$Kingdom,
                                   to_plot = !vector(mode = "logical", length = nrow(species_to_plot)),
                                   to_highlight = vector(mode = "logical", length = nrow(species_to_plot)))
  species_annotation$text <- as.character(species_annotation$text)
  
  # Overlapping names
  overlapping_labels <- c("Discosea")
  
  if (TRUE) {
    # Mutate overlapping
    species_annotation[species_annotation$text %in% overlapping_labels, "to_plot"] <- FALSE
    
    overlapped <- species_annotation %>% filter(text %in% overlapping_labels)
  }
  
  need_background <- c("Copepoda", "Embryophyta")
  # Need backgrounds?
  if (TRUE) {
    species_annotation[species_annotation$text %in% need_background, "to_highlight"] <- TRUE
  }
  
  # Way to improve
  # ways: "move", "line", "point"
  way <- "point"
  
  if (way == "move") {
    
    overlapped[overlapped$text == overlapping_labels[1], "x"] <- -0.0 + species_annotation[species_annotation$text == overlapping_labels[1], "x"]
    overlapped[overlapped$text == overlapping_labels[1], "y"] <- 0.05 + species_annotation[species_annotation$text == overlapping_labels[1], "y"]
    
  } else if (way == "line") {
    
    overlapped[overlapped$text == overlapping_labels[1], "x"] <- -0.5 + species_annotation[species_annotation$text == overlapping_labels[1], "x"]
    overlapped[overlapped$text == overlapping_labels[1], "y"] <- 0.1 + species_annotation[species_annotation$text == overlapping_labels[1], "y"]
    
  } else if (way == "point") {
    
    overlapped[overlapped$text == overlapping_labels[1], "x"] <- -0.5 + species_annotation[species_annotation$text == overlapping_labels[1], "x"]
    overlapped[overlapped$text == overlapping_labels[1], "y"] <- 0.1 + species_annotation[species_annotation$text == overlapping_labels[1], "y"]
    
  }
  
  # Add taxes
  for (i in seq(nrow(species_to_plot))) {
    # If it is needed to be plotted
    if (species_annotation[i, "to_plot"]) {
      # If it is needed to be displayed as it is
      if (!species_annotation[i, "to_highlight"]) {
        p <- p + geom_text(data = species_annotation[i,], 
                           mapping = aes(x = x, y = y, label = text), 
                           color = "slategrey", 
                           size = 3)
      } else {
        p <- p + geom_label(data = species_annotation[i,], 
                            mapping = aes(x = x, y = y, label = text), 
                            color = "slategrey",
                            fill = "white",
                            alpha = 0.7,
                            size = 3)
      }
    } else if (way == "move") {
      this_data <- overlapped[overlapped$text == species_annotation[i,]$text, ]
      
      p <- p + geom_text(data = this_data, 
                         mapping = aes(x = x, y = y, label = text), 
                         color = "slategrey", 
                         size = 3)
    } else if (way == "line") {
      this_data <- overlapped[overlapped$text == species_annotation[i,]$text, ]
      # Draw the label
      p <- p + geom_text(data = this_data, 
                         mapping = aes(x = x, y = y, label = text), 
                         color = "slategrey", 
                         size = 3)
      # Line data
      line_data <- data.frame(x_s = species_annotation[i,]$x, y_s = species_annotation[i,]$y, x_e = this_data$x, y_e = this_data$y)
      # Draw the line
      p <- p + geom_segment(data = line_data, mapping=aes(x=x_s, y=y_s, xend=x_e, yend=y_e), size=0.1, color="darkgrey", linetype = "dashed")
      
    } else if (way == "point") {
      this_data <- overlapped[overlapped$text == species_annotation[i,]$text, ]
      # Draw the label
      p <- p + geom_text(data = this_data, 
                         mapping = aes(x = x, y = y, label = text), 
                         color = "slategrey", 
                         size = 3)
      # Line data
      line_data <- data.frame(x_s = species_annotation[i,]$x, y_s = species_annotation[i,]$y, x_e = this_data$x, y_e = this_data$y)
      # Draw the line
      p <- p + geom_segment(data = line_data, mapping=aes(x=x_s, y=y_s, xend=x_e, yend=y_e), size=0.1, color="darkgrey", linetype = "dashed")
      # Draw the point
      p <- p + geom_point(data = line_data, mapping = aes(x = x_s, y = y_s), shape = 21, color = "darkgrey", fill = "white")
      
    }
  }
  
}

if (FALSE) {
  data_to_visiolize <- data_species %>% filter(Genera %in% c("Cryptomycota", "Basidiomycota", "Copepoda"))
  
  p <- p + stat_ellipse(data = data_to_visiolize, aes(x=CCA1, y=CCA2, fill=Genera),
                        geom="polygon", level=0.95, alpha=0.2)
  p <- p + geom_point(data = data_to_visiolize, mapping = aes(x = CCA1, y = CCA2, fill = Genera), shape = 21) 
}

p <- p + xlim(-5.2, 3.8)

print(p)

dev.off()


