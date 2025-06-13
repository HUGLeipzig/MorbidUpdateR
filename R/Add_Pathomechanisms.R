#' Add the pathogenic variant types from Clinvar as well as STRs
#'
#' @param .clinvar_tsv_filtered ClinVar's filtered .tsv file.
#' Should be in your global environment after running \code{\link{StartNewVersion}}
#'
#' @return A dataframe with pathogenic variant types per gene
#' @export
#'
#' @import dplyr
#' @import stringr
#' @import readr
#'
#' @examples
#' \dontrun{
#' PathoClinVarVars = Add_Pathomechanisms(.clinvar_tsv_filtered = clinvar_tsv_filtered)
#' }
Add_Pathomechanisms = function(.clinvar_tsv_filtered = clinvar_tsv_filtered){

  # Read ClinVar tsv and select pathogenic variants
  # Subset to GeneSymbol and Variant Type
  PathoClinVarVars = .clinvar_tsv_filtered %>%
    filter(GeneID != -1) %>%
    filter(ClinSigSimple == 1) %>%
    select(GeneSymbol, Type) %>%
    unique()

  # Read STR database and subset to GeneSymbol
  strs_path = config::get("strs_tsv_path")
  strs = read_tsv(strs_path) %>%
    select("gene") %>%
    distinct() %>%
    rename("GeneSymbol" = "gene") %>%
    mutate("Type" = "STR")

  # bind strs to PathoClinVarVars and collapse by variant type
  PathoClinVarVars_STRs = PathoClinVarVars %>%
    rbind(strs) %>%
    mutate(Type = str_replace_all(Type, "single nucleotide variant", "SNV")) %>%
    mutate(Type = str_replace_all(Type, "copy number gain", "CNV gain")) %>%
    mutate(Type = str_replace_all(Type, "copy number loss", "CNV loss")) %>%
    group_by(GeneSymbol) %>%
    arrange(Type) %>%
    summarise(PathoVars = paste(Type, collapse = ", "))

  return(PathoClinVarVars_STRs)

}


