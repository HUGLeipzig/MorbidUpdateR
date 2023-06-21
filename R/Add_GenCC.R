#' Add Genes from the GenCC Database
#'
#' @param .gencc_tsv The GenCC tsv file.
#' Should be in your global environment after running \code{\link{StartNewVersion}}
#' @param classification The confidence of GenCC-Genes. One of "Definitive" (default),
#' "Limited", "Moderate", "Refuted" or "All"
#'
#' @return A dataframe with the GenCC genes
#' @export
#'
#' @import stringr
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' gencc = Add_GenCC()
#' }

Add_GenCC = function(.gencc_tsv = gencc_tsv, classification = "Definitive"){

  gencc = gencc_tsv %>%
    mutate(hgnc_id = str_remove(gene_curie, "HGNC:")) %>%
    filter(classification_title == classification) %>%
    select(gene_symbol, classification_title, hgnc_id) %>%
    distinct() %>%
    mutate(isGenCCGene = T)

  return(gencc)

}
