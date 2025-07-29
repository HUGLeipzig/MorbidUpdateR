#' Add genes which have to be added manually from the Excel spreadsheet
#'
#' @param ManualFile The path to the file where the Manual Genes are stored
#'
#' @return A dataframe with the manually added genes
#' @export
#'
#' @import dplyr
#' @import readxl
#'
#' @examples
#' \dontrun{
#' Manual = Add_ManualGenes()
#' }

Add_ManualGenes = function(ManualFile = "W:/HUG/04 Klinische Genomik/10 Panels/MorbidGenes-Panel/GenesToBeAddedManually.xlsx"){

  Manual = read_excel(ManualFile,
                      sheet = 1, col_names = c("Gene", "PMID", "comment"), trim_ws = F)

  Manual = Manual %>%
    mutate(Gene = str_remove_all(string = Gene, pattern = "\\s")) %>%
    mutate(addedManually = TRUE)

  return(Manual)
}
