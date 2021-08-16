
#' get all known O.sative gene ids - the universe
#' @param genefile csv file with GeneID column of all Mo ids
all_genes <- function(which="indica"){
  
  col_types <- readr::cols(
    GeneID = readr::col_character()
  )
  
  if (which == "indica"){
    readr::read_csv(here::here("inst", "extdata", "indica_group_uniq_ids.txt"), col_types = col_types)$GeneID
    }
  else if (which == "japonica") {
    readr::read_csv(here::here("inst", "extdata", "japonica_group_uniq_ids.txt"), col_types = col_types)$GeneID 
  }
}


read_terms <- function(termsfile=NULL) {
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
mapping <- function(termsfile=NULL, which = "indica") {

  terms <- read_terms(termsfile)
  
  term2gene <- data.frame(
    term = terms$`GO term accession`,
    gene = terms$`Transcript stable ID`
  )

  term2name <- data.frame(
    term = terms$`GO term accession`,
    name = terms$`GO term definition`
  )

  all_genes <- all_genes(which = which)

  return(list(
    term2gene = term2gene,
    term2name = term2name,
    all_genes = all_genes
  ))
}

#' run clusterProfiler::enricher on vector of gene ids,
#' @param genes character vector of gene IDs of interest
#' @param termsfile TSV file from biomart export of GO terms
#' @return enricher object
#' @export
do_enrich <- function(genes, termsfile=NULL,
                      which = "indica", 
                      ...) {
  
  if (is.null(termsfile) & which == "indica"){
    termsfile <- here::here("inst", "extdata", "indica_group_mart_export.txt")
  }
  else if (is.null(termsfile) & which == "japonica"){
    termsfile <- here::here("inst", "extdata", "japonica_group_mart_export.txt")
  }
  
  info <- mapping(termsfile, which = which)

  clusterProfiler::enricher( genes,
                             universe=info$all_genes,
                             TERM2GENE=info$term2gene,
                             TERM2NAME=info$term2name,
                             ...
                             )
}

#' conver enricher object to DAVID format
#' @param enrich enricher object from do_enrich
#' @param termsfile TSV file from biomart export of GO terms
#' @return DAVID format dataframe
#' @export
enricher_to_david <- function(enrich, which="indica"){
  
  
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
