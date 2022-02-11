#' Add Genes from the SysNDD Database
#'
#' @param category The confidence of SysNDD-Genes. One of "Definitive" (default),
#' "Limited", "Moderate", "Refuted" or "All"
#' @param inheritance Filter for specific inheritance patterns. One of "All" (default),
#' "Dominant", "Other", "Recessive" or "X-linked"
#' @param download Should the latest SysNDD be downloaded? If not, previous version will be used
#'
#' @return
#' @export
#'
#' @import httr
#' @import jsonlite
#' @import readr
#'
#' @examples
#' \dontrun{
#' SysNDD = Add_SysNDD()
#' }

Add_SysNDD = function(category = "Definitive", inheritance = "All",
                      download = T){

  if(download == T){
    httr::set_config(config(ssl_verifypeer = 0L))
    sysndd_res = GET(paste0("http://sysndd.org/alb/api/panels/browse?",
                            "category_input=", category,
                            "&inheritance_input=", inheritance))
    sysndd_res_text = suppressMessages(content(sysndd_res, as = "text", type = "text/csv"))
    sysnddTest = fromJSON(sysndd_res_text)$data

    genes = unique(sysnddTest$symbol)

    SysNDD = data.frame(SysNDDGene = genes, isSysNDDGene = T)

    write_tsv(SysNDD, "W:/HUG/04 Klinische Genomik/10 Panels/MorbidGenes-Panel/SysNDD_Genes.tsv",
              col_names = T)
  } else {

    SysNDDfile = read_tsv("W:/HUG/04 Klinische Genomik/10 Panels/MorbidGenes-Panel/SysNDD_Genes.tsv",
                      col_names = T)

    genes = unique(SysNDDfile$symbol)

    SysNDD = data.frame(SysNDDGene = genes, isSysNDDGene = T)
  }

  return(SysNDD)
}
