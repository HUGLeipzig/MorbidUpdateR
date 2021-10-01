#' Title
#'
#' @param directory The same you provided in \code{\link{StartNewVersion}}
#' @param version The same you provided in \code{\link{StartNewVersion}}
#' @param .VarvisGeneManagement_HGNC The result of the function \code{\link{Add_HGNC}}
#' or in your Global Env after running \code{\link{Add_all}}
#' @param .clinvar_tsv_filtered_patho_HGNC The result of the function \code{\link{Add_ClinVar}}
#' or in your Global Env after running \code{\link{Add_all}}
#' @param .morbidmap_reshape The result of the function \code{\link{Add_OMIM}}
#' or in your Global Env after running \code{\link{Add_all}}
#' @param .Manual The result of the function \code{\link{Add_ManualGenes}}
#' or in your Global Env after running \code{\link{Add_all}}
#' @param .PanelAppGenes The result of the function \code{\link{Add_PanelApp}}
#' or in your Global Env after running \code{\link{Add_all}}
#' @param save should the Panel be saved as a .csv? Defaults to Yes in the corresponding directory
#'
#' @return
#' @export
#'
#' @import readr
#' @import dplyr
#'
#' @examples
#' Build_MorbidGenesPanel()

Build_MorbidGenesPanel = function(directory = "W:/HUG/04 Klinische Genomik/10 Panels/MorbidGenes-Panel/",
                                  version = format(Sys.time(), "%Y-%m"),
                                  .VarvisGeneManagement_HGNC = VarvisGeneManagement_HGNC,
                                  .clinvar_tsv_filtered_patho_HGNC = clinvar_tsv_filtered_patho_HGNC,
                                  .morbidmap_reshape = morbidmap_reshape,
                                  .Manual = Manual,
                                  .PanelAppGenes = PanelAppGenes,
                                  save = T){

  MorbidGenes_Panel = VarvisGeneManagement_HGNC %>%
    left_join(HGMD_count_HGNC,
              by = c("HGNC_symbol_corrected" = "HGNC_symbol_corrected"),
              na_matches = "never") %>%
    left_join(clinvar_tsv_filtered_patho_HGNC,
              by = c("HGNC_symbol_corrected" = "Approved_Gene_Symbol_HGNC"),
              na_matches = "never") %>%
    left_join(morbidmap_reshape,
              by = c("HGNC_symbol_corrected" = "Approved_Gene_Symbol_HGNC"),
              na_matches = "never") %>%
    left_join(Manual,
              by = c("HGNC_symbol_corrected" = "Gene"),
              na_matches = "never") %>%
    left_join(PanelAppGenes,
              by = c("HGNCID" = "hgnc"),
              na_matches = "never") %>%
    select("SYMBOL", "NAME", "CHROMOSOME", "CHROMOSOMELOCATION", "TRANSCRIPT",
           "NCBIID", "OMIMID", "LRGID", "ENSEMBLID", "HGNCID", "symbol",
           "HGNC_symbol_corrected", "HGMD_pathogenic_variant_count",
           "HGMD_pathogenic_variant_count_cutoff", "ClinVarPathogenicCount",
           "ClinVarPathogenicCount_cutoff","Phenotype","Phenotype_MIM_Numbers",
           "MIM_Numbers", "addedManually", "isPanelAppGene") %>%
    mutate(has_Phenotype_MIM_Number = !is.na(Phenotype_MIM_Numbers)) %>%
    mutate(keep = case_when(
      # Genes with HGMD cutoff
      HGMD_pathogenic_variant_count_cutoff ~ TRUE,
      # Genes with ClinVar cutoff
      ClinVarPathogenicCount_cutoff ~ TRUE,
      # PanelAppGenes
      isPanelAppGene ~ TRUE,
      # Genes added manually
      addedManually ~ TRUE,
      # has OMIM Phenotype
      has_Phenotype_MIM_Number ~ TRUE
    )) %>%
    distinct()

  # check if the directory string ends with a "\" or "/" and append if not
  directory = ifelse(grepl("/$|\\$", directory),
                     directory,
                     paste0(directory, "/"))

  # check if a subversion is already present
  # if no, define subversion = 1
  # if yes, add new subversion
  files = list.files(path = directory)

  filename = tail(sort(files[grep(version, files)]), n = 1)

  if(save == T){
    write_csv2(MorbidGenes_Panel,
               paste0(directory, filename, "/", filename, ".csv"),
               na = "NA", append = FALSE, col_names = T, escape = "double")
  }

}
