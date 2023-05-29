# Requetes REST APIs pour retrouver les pathways
library(httr)
library(jsonlite)

# Fonction pour rechercher les IDs KEGG à partir des IDs UniProt
getKeggIdFromUniprot <- function(uniprot_ids) {
  # URL de l'API UniProt
  api_url <- "https://rest.uniprot.org/idmapping/run"

  uniprot_ids <- vectorToString(uniprot_ids, sep = ",")
  # Corps de la requête
  request_body <- list(
    ids = uniprot_ids,
    from = "UniProtKB_AC-ID",
    to = "KEGG"
  )

  # Faire la requête HTTP POST à l'API UniProt
  response <- POST(url = api_url, body = request_body, encode = "multipart", accept_json())
  submission <- content(response, as = "parsed")

  if (!isJobReady(submission[["jobId"]])) {
    stop("Le job n'a pas pu etre retrouve sur UniProt.")
  } else {
    url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
    r <- GET(url = url, accept_json())
    details <- content(r, as = "parsed")
    url <- getResultsURL(details[["redirectURL"]])
    # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
    url <- paste(url, "?format=tsv", sep = "")
    r <- GET(url = url, accept_json())
    resultsTable <- read.table(text = content(r), sep = "\t", header = TRUE, encoding = "UTF-8")
  }
}


# Fonction pour vérifier et transformer un vecteur de strings en une seule string avec des éléments séparés par une virgule
vectorToString <- function(obj, sep = ",") {
  if (is.vector(obj) && all(sapply(obj, is.character))) {
    return(paste(obj, collapse = sep))
  } else {
    stop("L'objet doit être un vecteur de strings.")
  }
}


isJobReady <- function(jobId) {
  pollingInterval <- 5
  nTries <- 20
  for (i in 1:nTries) {
    url <- paste("https://rest.uniprot.org/idmapping/status/", jobId, sep = "")
    r <- GET(url = url, accept_json())
    status <- content(r, as = "parsed")
    if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
      return(TRUE)
    }
    if (!is.null(status[["messages"]])) {
      print(status[["messages"]])
      return(FALSE)
    }
    Sys.sleep(pollingInterval)
  }
  return(FALSE)
}


getResultsURL <- function(redirectURL) {
  if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
    url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
  } else {
    url <- gsub("/results/", "/results/stream/", redirectURL)
  }
}


# Récupérer les ID de pathways dans un vecteur à partir de l'identifiant de gene KEGG
getKeggPathways <- function(kegg_ids) {
  # Construire l'URL de l'API KEGG pour l'identifiant KEGG
  api_url <- paste0("http://rest.kegg.jp/link/pathway/", vectorToString(kegg_ids, sep = "+"))

  # Faire la requête HTTP à l'API KEGG
  response <- GET(api_url)

  # Vérifier si la requête est réussie
  if (http_status(response)$category == "Success") {
    # Extraire les pathways de la réponse
    pathways <- content(response, "text")
    # Nettoyer la réponse et extraire les identifiants de pathways
    pathway_ids <- strsplit(pathways, "\n")[[1]]
    pathway_ids <- gsub("^.*:(.*)", "\\1", pathway_ids)
    return(pathway_ids)
  } else {
    # En cas d'erreur, afficher un message d'erreur
    stop("Une erreur s'est produite lors de la récupération des pathways.")
  }
}

getGenesFromPathways <- function(pathway_ids, species = "hsa") {
  # Construire l'URL de l'API KEGG pour l'identifiant KEGG
  api_url <- paste0("http://rest.kegg.jp/link/", species, "/", vectorToString(pathway_ids, sep = "+"))

  # Faire la requête HTTP à l'API KEGG
  response <- GET(api_url)

  # Vérifier si la requête est réussie
  if (http_status(response)$category == "Success") {
    # Extraire les genes de la réponse
    genes <- content(response, "text")
    # Nettoyer la réponse et extraire les identifiants de genes
    genes_ids <- strsplit(genes, "\n")[[1]]
    genes_ids <- gsub("^.*\t(\\w+:\\d+)$", "\\1", genes_ids)
    return(genes_ids)
  } else {
    # En cas d'erreur, afficher un message d'erreur
    stop("Une erreur s'est produite lors de la récupération des genes associes aux pathways KEGG.")
  }
}


# Fonction pour rechercher les IDs UniProt à partir des IDs KEGG
getUniprotIdFromKegg <- function(kegg_ids) {
  # URL de l'API UniProt
  api_url <- "https://rest.uniprot.org/idmapping/run"

  kegg_ids <- vectorToString(kegg_ids, sep = ",")
  # Corps de la requête
  request_body <- list(
    ids = kegg_ids,
    to = "UniProtKB",
    from = "KEGG"
  )

  # Faire la requête HTTP POST à l'API UniProt
  response <- POST(url = api_url, body = request_body, encode = "multipart", accept_json())
  submission <- content(response, as = "parsed")

  if (!isJobReady(submission[["jobId"]])) {
    stop("Le job n'a pas pu etre retrouve sur UniProt.")
  } else {
    url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
    r <- GET(url = url, accept_json())
    details <- content(r, as = "parsed")
    url <- getResultsURL(details[["redirectURL"]])
    # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
    url <- paste(url, "?format=tsv", sep = "")
    r <- GET(url = url, accept_json())
    resultsTable <- read.table(text = content(r), sep = "\t", header = TRUE, encoding = "UTF-8")
  }
}
