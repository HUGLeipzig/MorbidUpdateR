#' Add the PanelApp Genes
#'
#' @param PanelAppFile The path to the file where the PanelApp Genes are stored
#'
#' @return
#' @export
#'
#' @import dplyr
#' @import readr
#'
#' @examples
#' \dontrun{
#' PanelAppGenes = Add_PanelApp()
#' }

Add_PanelApp = function(PanelAppFile = "W:/HUG/04 Klinische Genomik/10 Panels/MorbidGenes-Panel/PanelAppGenes/2021_08_25_PA_all_genes.csv"){

  PanelAppGenes = read_tsv(PanelAppFile,
                           col_names = T)

  PanelAppGenes = PanelAppGenes %>%
    group_by(hgnc) %>%
    summarise(sum = sum(confidence)) %>%
    ungroup() %>%
    filter(sum > 1) %>%
    mutate(isPanelAppGene = TRUE)

  return(PanelAppGenes)
}
