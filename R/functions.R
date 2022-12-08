#' get all known O.sative gene ids - the universe
#' genefile csv file with GeneID column of all O. sativa ids
#' @param which subspecies with either "indica" or "japonica". Default is "indica"
all_genes <- function(which = "indica"){
  
  col_types <- readr::cols(
    GeneID = readr::col_character()
  )
  
  if (which == "indica"){
    filename = "indica_group_uniq_ids.txt"
  }
  else if (which == "japonica") {
    filename ="japonica_group_uniq_ids.txt"
  }
  
  path_to_file <- system.file("extdata", 
    filename, 
    package = "osgo", 
    mustWork = TRUE)
  
  readr::read_csv(path_to_file, col_types = col_types)$GeneID
}

#' read the terms file
#' @param termsfile TSV file from export of GO terms for indica group at plants.ensembl.org Oryza sativa Japonica Group genes (IRGSP-1.0)
read_terms <- function(termsfile = NULL) {
readr::read_tsv(termsfile, col_types = readr::cols(
    `Gene stable ID` = readr::col_character(),
    `Transcript stable ID` = readr::col_character(),
    `GO term accession` = readr::col_character(),
    `GO term name` = readr::col_character(),
    `GO term definition` = readr::col_character(),
    `GO term evidence code` = readr::col_character(),
    `GO domain` = readr::col_character()
  ))
  
}

#' get mapping between Os Gene IDs and GO Terms
#' @param termsfile TSV file from export of GO terms for indica group at plants.ensembl.org Oryza sativa Japonica Group genes (IRGSP-1.0)
#' @param which subspecies with either "indica" or "japonica". Default is "indica"
mapping <- function(termsfile = NULL, which = "indica") {

  terms <- read_terms(termsfile)
  
  term2gene <- data.frame(
    term = terms$`GO term accession`,
    gene = terms$`Transcript stable ID`
  )
  
  term2go <- data.frame(
    term = terms$`GO term accession`,
    go = terms$`GO term accession`
  )
  
  term2name <- data.frame(
    term = terms$`GO term accession`,
    name = terms$`GO term name`
  )
  
  term2description <- data.frame(
    term = terms$`GO term accession`,
    description = terms$`GO term definition`
  )

  all_genes <- all_genes(which = which)

  return(list(
    term2go = term2go,
    term2gene = term2gene,
    term2name = term2name,
    term2description = term2description,
    all_genes = all_genes
  ))
}

#' run clusterProfiler::enricher on vector of gene ids,
#' @importFrom clusterProfiler enricher
#' @param genes character vector of gene IDs of interest
#' @param which subspecies with either "indica" or "japonica". Default is "indica"
#' @param termsfile TSV file from biomart export of GO terms
#' @param label_type string for the type of labels wanted; default is "go" for the go term (GO:number), "name" for the GO name, and "description" for the GO long description.
#' @param ... more options described in clusterProfiler::enricher()
#' @return enricher object
#' @export
do_enrich <- function(genes, termsfile = NULL,
                      which = "indica", label_type = "go",
                      ...) {
  
  if (is.null(termsfile) & which == "indica"){
    filename = "indica_group_mart_export.txt"
    termsfile <- system.file("extdata", filename, package = "osgo", mustWork = TRUE)
  } else if (is.null(termsfile) & which == "japonica"){
    filename = "japonica_group_mart_export.txt"
    termsfile <- system.file("extdata", filename, package = "osgo", mustWork = TRUE)
  }

  info <- mapping(termsfile, which = which)
  
  if(label_type == "go"){
    labels = info$term2go
  } else if(label_type == "name") {
    labels = info$term2name
  } else if(label_type == "description") {
    labels = info$term2description
  }
  
  clusterProfiler::enricher( genes,
    universe=info$all_genes,
    TERM2GENE=info$term2gene,
    TERM2NAME=labels,
    ...
  )
}

#' convert enricher object to DAVID format
#' @param enrich enricher object from do_enrich
#' @param which subspecies with either "indica" or "japonica". Default is "indica"
#' @return DAVID format dataframe
#' @export
enricher_to_david <- function(enrich, which = "indica"){
  
  
  comma_sep_genes = gsub("/", ", ", enrich@result$geneID)

  terms <- NULL
  if (which == "indica" ){
    terms <- read_terms(termsfile = here::here("inst", "extdata", "indica_group_mart_export.txt") )
  }
  else if (which == "japonica" ){
    terms <- read_terms(here::here("inst", "extdata", "japonica_group_mart_export.txt") )
  }
  
  terms <- dplyr::mutate(terms,
    short_category = dplyr::if_else(`GO domain` == 'biological_process', "BP",
                                    dplyr::if_else(`GO domain` == "molecular function", "MF", "CC"))
  )

  id_to_category <- terms$short_category
  names(id_to_category) <- terms$`GO term accession`
  fixed_category <- id_to_category[enrich@result$ID]

  data.frame(
    category = fixed_category,
    ID = enrich@result$ID,
    term = enrich@result$Description,
    genes = comma_sep_genes,
    adj_pval = enrich@result$p.adjust
  )

}
