#' Add the HGNC Data and genomic positions (GRCh38) to the Varvis Gene Management file
#'
#' @param .VarvisGeneManagement The Varvis Gene Management file.
#' Should be in your global environment after running \code{\link{StartNewVersion}}
#'
#' @return
#' @export
#'
#' @import dplyr
#' @import biomaRt
#'
#' @examples
#' \dontrun{
#' VarvisGeneManagement_HGNC = Add_HGNC(.VarvisGeneManagement = VarvisGeneManagement)
#' }

Add_HGNC = function(.VarvisGeneManagement = VarvisGeneManagement){

  # get genomic position via Ensembl/BiomaRt
  mart = useMart("ensembl")
  mart = useDataset("hsapiens_gene_ensembl", mart)
  attributes = c("ensembl_gene_id", "start_position","end_position","hgnc_symbol",
                 "chromosome_name")
  all.genes = getBM(attributes=attributes, mart=mart)

  all.genes = all.genes %>%
    filter(chromosome_name %in% c(seq(1:22), "X", "Y", "MT")) %>%
    dplyr::select(-chromosome_name) %>%
    distinct()

  VarvisGeneManagement_HGNC <- VarvisGeneManagement %>%
    left_join(hgnc_complete_set, by = c("NCBIID" = "entrez_id")) %>%
    select("SYMBOL","NAME","CHROMOSOME","CHROMOSOMELOCATION","TRANSCRIPT",
           "NCBIID","OMIMID","LRGID","ENSEMBLID","HGNCID","symbol") %>%
    mutate(HGNC_symbol_corrected =
             case_when(
               is.na(symbol) ~ SYMBOL,
               !is.na(symbol) ~ symbol,
             )
    ) %>%
    left_join(all.genes, by = c("ENSEMBLID" = "ensembl_gene_id"),
              na.matches = "never")

  return(VarvisGeneManagement_HGNC)
}

