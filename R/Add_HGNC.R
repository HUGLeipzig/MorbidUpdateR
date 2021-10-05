#' Add the HGNC Data to the Varvis Gene Management file
#'
#' @param .VarvisGeneManagement The Varvis Gene Management file.
#' Should be in your global environment after running \code{\link{StartNewVersion}}
#'
#' @return
#' @export
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' VarvisGeneManagement_HGNC = Add_HGNC(.VarvisGeneManagement = VarvisGeneManagement)
#' }

Add_HGNC = function(.VarvisGeneManagement = VarvisGeneManagement){
  VarvisGeneManagement_HGNC <- VarvisGeneManagement %>%
    left_join(hgnc_complete_set, by = c("NCBIID" = "entrez_id")) %>%
    select("SYMBOL","NAME","CHROMOSOME","CHROMOSOMELOCATION","TRANSCRIPT",
           "NCBIID","OMIMID","LRGID","ENSEMBLID","HGNCID","symbol") %>%
    mutate(HGNC_symbol_corrected =
             case_when(
               is.na(symbol) ~ SYMBOL,
               !is.na(symbol) ~ symbol,
             )
    )
  return(VarvisGeneManagement_HGNC)
}
