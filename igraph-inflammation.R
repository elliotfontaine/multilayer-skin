# devtools::install_github("npjc/biogridr")
setwd("/home/elliot/RWork")
getwd()

### ____ Requetes REST APIs pour retrouver les pathways ____###
library(httr)
library(jsonlite)
library(biogridr)
source("rest_functions.R")


# Lecture du contenu du fichier des proteines impliques inflammation peau
uniprot_ids <- read.table("prots_skin_infl.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
uniprot_ids <- as.vector(uniprot_ids$UniProtKB)

## KEGG Pathways ##
# Récupérer les pathways KEGG pour les IDs UniProt
kegg_ids <- getKeggIdFromUniprot(uniprot_ids)
kegg_ids <- as.vector(kegg_ids$To)
pathway_ids <- unique(getKeggPathways(kegg_ids))
# Recuperer tous les genes associes au pathways
genes_ids <- unique(getGenesFromPathways(pathway_ids))
# Retour aux IDs UniProt pour les genes
enz_ids <- getUniprotIdFromKegg(genes_ids)
enz_ids <- subset(enz_ids, Reviewed == "reviewed")


## PPI reactome.org ##
# Enter you personal BioGrid API access key
bg_get_key("Elliot", "Fontaine", "elliot.fontaine@gmail.com", "CHUL")
interactions_df <- bg("interactions") %>%
  bg_constrain(geneList = paste(uniprot_ids, collapse = "|")) %>%
  bg_constrain(taxId = taxId("Homo sapiens")) %>%
  bg_get_results()



### ____ Graphes ____###

# Creation sous-reseau inflammation
library(tidyverse)
library(igraph)

# https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.221/BIOGRID-ALL-4.4.221.tab3.zip
biogrid <- read_tsv("./BIOGRID-ALL-4.4.221.tab3.txt") %>%
  dplyr::select(
    "SWISS-PROT Accessions Interactor A",
    "SWISS-PROT Accessions Interactor B",
    "Organism Name Interactor A"
  ) %>%
  set_names("A", "B", "Species") %>%
  filter(Species == "Homo sapiens")

# Interactions PPI avec au moins l'une des 2 enzymes provenant de la liste excel
biogrid.small <- biogrid %>% filter(A %in% uniprot_ids | B %in% uniprot_ids)

# Interactions PPI avec au moins l'une des 2 enzymes provenant de la liste excel,
# ou bien les 2 enzymes dans l'environnement metabolique proche
biogrid.big <- biogrid %>% filter(
  (A %in% uniprot_ids | B %in% uniprot_ids) | (A %in% enz_ids$Entry & B %in% enz_ids$Entry)
)

# Reseau igraph small
graph.small <- igraph::graph_from_data_frame(biogrid.small, directed = FALSE)
graph.small <- igraph::simplify(graph.small, remove.multiple = TRUE, remove.loops = TRUE)
summary(graph.small)
# Voir ligne 145
tkplot(graph.small)
plot(
  graph.small,
  vertex.label = NA,
  edge.width = 0.5,
  layout = layout_with_kk
)


# Reseau igraph big
graph.big <- igraph::graph_from_data_frame(biogrid.big, directed = FALSE)
graph.big <- igraph::simplify(graph.big, remove.multiple = TRUE, remove.loops = TRUE)
summary(graph.big)

layout <- layout_in_circle(graph.big)
plot(graph.big, layout = layout)
tkplot(graph.big)












### ____ Graphes OLD ____###

# Creation sous-reseau inflammation
library(tidyverse)
library(igraph)

# https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.221/BIOGRID-ALL-4.4.221.tab3.zip
biogrid <- read_tsv("./BIOGRID-ALL-4.4.221.tab3.txt") %>%
  dplyr::select(
    "SWISS-PROT Accessions Interactor A",
    "SWISS-PROT Accessions Interactor B",
    "Organism Name Interactor A"
  ) %>%
  set_names("A", "B", "Species")
# Extraction du reseau PPI humain
biogrid <- biogrid %>% filter(Species == "Homo sapiens")


# https://reactome.org/PathwayBrowser/#/R-HSA-168256&DTAB=MT
reactome.immunity <- read_tsv("./reactomeDB_immunity_homo.tsv") %>%
  dplyr::select("Identifier", "MoleculeName") %>%
  mutate(MoleculeName = str_remove(MoleculeName, "^[^\\s]*\\s"))

# https://reactome.org/PathwayBrowser/#/R-HSA-168249&DTAB=MT
reactome.innate <- read_tsv("./reactomeDB_innate_homo.tsv") %>%
  dplyr::select("Identifier", "MoleculeName") %>%
  mutate(MoleculeName = str_remove(MoleculeName, "^[^\\s]*\\s"))

# Sous-tableau contenant seulement les interactions entre proteines de l'immunite
biogrid.immunity <- biogrid %>% filter(A %in% reactome.immunity$Identifier & B %in% reactome.immunity$Identifier)
# Sous-tableau contenant seulement les interactions entre proteines de l'immunite innee
biogrid.innate <- biogrid %>% filter(A %in% reactome.innate$Identifier & B %in% reactome.immunity$Identifier)


# Reseau igraph immunite
graph.immunity <- igraph::graph_from_data_frame(biogrid.immunity, directed = FALSE)
graph.immunity <- igraph::simplify(graph.immunity, remove.multiple = TRUE, remove.loops = TRUE)
summary(graph.immunity)
plot(graph.immunity)

# Reseau igraph immunite innee
graph.innate <- igraph::graph_from_data_frame(biogrid.innate, directed = FALSE)
graph.innate <- igraph::simplify(graph.innate, remove.multiple = TRUE, remove.loops = TRUE)
summary(graph.innate)

layout <- layout_in_circle(graph.innate)
plot(graph.innate, layout = layout)



# Calculate eigen centrality and check the distribution We're attaching the
# result of eigen_centrality() straight onto the vertices as verticy-attributes
V(graph.small)$ec <- eigen_centrality(graph.small, directed = T, weights = NA)$vector
hist(V(graph.small)$ec)

# You could use the scales package, or define this normalisation function:
normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
(V(graph.small)$ec_index <- round(normalize(V(graph.small)$ec) * 2) + 1)
# ec_index should now be a category between 1 and 10 for your centralities

# Build a color-mapping with 10 categories and set the color of each
# node to the corresponding color of each centrality-measure category
V(graph.small)$color <- colorRampPalette(c("turquoise", "yellow", "red"))(3)[V(graph.small)$ec_index]
V(graph.small)$size <- V(graph.small)$ec_index
# Look at what we did
table(V(graph.small)$color)
plot(graph.small, vertex.label = NA, vertex.size = 5)
