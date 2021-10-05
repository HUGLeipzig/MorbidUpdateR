#' Add HGMD Data and determine pathogenic cutoff
#'
#' @param cutoff minimum number of pathogenic variants per gene to be included into the Morbid Genes panel
#'
#' @return
#' @export
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' HGMD_count_HGNC = Add_HGMD(cutoff = 4)
#' }

Add_HGMD = function(cutoff = 4){
  HGMD = read_csv("W:/HUG/04 Klinische Genomik/09 HGMD Download and Analyses/HGMD_Advanced_Substitutions.csv.gz",
                  col_select = c("gene", "entrezid", "Variant_class"))

  HGMD_count = HGMD %>%
    group_by(gene, entrezid) %>%
    summarise(HGMD_pathogenic_variant_count = sum(Variant_class == "DM")) %>%
    ungroup() %>%
    mutate(HGMD_pathogenic_variant_count_cutoff = (HGMD_pathogenic_variant_count >= cutoff)) %>%
    transmute(gene = gene, HGMD_pathogenic_variant_count = HGMD_pathogenic_variant_count,
              HGMD_pathogenic_variant_count_cutoff = HGMD_pathogenic_variant_count_cutoff,
              entrezid = as.integer(entrezid))


  HGMD_count_HGNC <- HGMD_count %>%
    left_join(hgnc_complete_set, by = c("entrezid" = "entrez_id")) %>%
    select("gene", "HGMD_pathogenic_variant_count",
           "HGMD_pathogenic_variant_count_cutoff", "entrezid", "symbol") %>%
    mutate(HGNC_symbol_corrected =
             case_when(
               is.na(symbol) ~ gene,
               !is.na(symbol) ~ symbol,
             )
    ) %>%
    select("gene","HGMD_pathogenic_variant_count","HGMD_pathogenic_variant_count_cutoff","entrezid","HGNC_symbol_corrected")

  return(HGMD_count_HGNC)
}
