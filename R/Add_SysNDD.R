#' Add Genes from the SysNDD Database
#'
#' @param category The confidence of SysNDD-Genes. One of "Definitive" (default),
#' "Limited", "Moderate", "Refuted" or "All"
#' @param inheritance Filter for specific inheritance patterns. One of "All" (default),
#' "Dominant", "Other", "Recessive" or "X-linked"
#' @param download Should the latest SysNDD be downloaded? If not, previous version will be used
#'
#' @return A dataframe with the SysNDD genes
#' @export
#'
#' @import httr
#' @import jsonlite
#' @import readr
#' @import config
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

    write_tsv(SysNDD, config::get("sysndd_tsv_path"),
              col_names = T)
  } else {

    SysNDD = read_tsv(config::get("sysndd_tsv_path"),
                      col_names = T)

  }

  return(SysNDD)
}
