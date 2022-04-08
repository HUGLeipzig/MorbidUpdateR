#' Build the final Morbid Genes Panel
#'
#' @param directory The same you provided in \code{\link{StartNewVersion}}
#' @param version The same you provided in \code{\link{StartNewVersion}}
#' @param .VarvisGeneManagement_HGNC The result of the function \code{\link{Add_HGNC}}
#' or in your Global Env after running \code{\link{Add_all}}
#' @param .clinvar_tsv_filtered_patho_HGNC The result of the function \code{\link{Add_ClinVar}}
#' or in your Global Env after running \code{\link{Add_all}}
#' @param .HGMD_count_HGNC The result of the function \code{\link{Add_HGMD}}
#' or in your Global Env after running \code{\link{Add_all}}
#' @param .morbidmap_reshape The result of the function \code{\link{Add_OMIM}}
#' or in your Global Env after running \code{\link{Add_all}}
#' @param .Manual The result of the function \code{\link{Add_ManualGenes}}
#' or in your Global Env after running \code{\link{Add_all}}
#' @param .SysNDD The result of the function \code{\link{Add_SysNDD}}
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
#' \dontrun{
#' Build_MorbidGenesPanel()
#' }

Build_MorbidGenesPanel = function(directory = "W:/HUG/04 Klinische Genomik/10 Panels/MorbidGenes-Panel/",
                                  version = format(Sys.time(), "%Y-%m"),
                                  .VarvisGeneManagement_HGNC = VarvisGeneManagement_HGNC,
                                  .HGMD_count_HGNC = HGMD_count_HGNC,
                                  .clinvar_tsv_filtered_patho_HGNC = clinvar_tsv_filtered_patho_HGNC,
                                  .morbidmap_reshape = morbidmap_reshape,
                                  .Manual = Manual,
                                  .PanelAppGenes = PanelAppGenes,
                                  .SysNDD = SysNDD,
                                  save = T){

  MorbidGenes_Panel = VarvisGeneManagement_HGNC %>%
    left_join(.HGMD_count_HGNC,
              by = c("HGNC_symbol_corrected" = "HGNC_symbol_corrected"),
              na_matches = "never") %>%
    left_join(.clinvar_tsv_filtered_patho_HGNC,
              by = c("HGNC_symbol_corrected" = "Approved_Gene_Symbol_HGNC"),
              na_matches = "never") %>%
    left_join(.morbidmap_reshape,
              by = c("HGNC_symbol_corrected" = "Approved_Gene_Symbol_HGNC"),
              na_matches = "never") %>%
    left_join(.Manual,
              by = c("HGNC_symbol_corrected" = "Gene"),
              na_matches = "never") %>%
    left_join(.PanelAppGenes,
              by = c("HGNC_symbol_corrected" = "entity_name"),
              na_matches = "never") %>%
    left_join(.SysNDD,
              by = c("HGNC_symbol_corrected" = "SysNDDGene"),
              na_matches = "never") %>%
    select("HGNC_symbol_corrected", "NAME", "CHROMOSOME", "bed_hg19",
           "bed_hg38", "CHROMOSOMELOCATION", "TRANSCRIPT",
           "NCBIID", "OMIMID", "LRGID", "ENSEMBLID", "HGNCID",
           "HGMD_pathogenic_variant_count",
           "HGMD_pathogenic_variant_count_cutoff", "ClinVarPathogenicCount",
           "ClinVarPathogenicCount_cutoff","Phenotype","Phenotype_MIM_Numbers",
           "MIM_Numbers", "addedManually", "isPanelAppGene",
           "isUKPanelAppGene", "isAustraliaPanelAppGene", "isSysNDDGene") %>%
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
      # SysID Genes
      isSysNDDGene ~ TRUE,
      # has OMIM Phenotype
      has_Phenotype_MIM_Number ~ TRUE
    )) %>%
    # add a score counting the "TRUE" values to determine which gene has most evidence
    mutate(MorbidScore = rowSums(.[c("HGMD_pathogenic_variant_count_cutoff",
                                     "ClinVarPathogenicCount_cutoff",
                                     "isPanelAppGene",
                                     "isSysNDDGene",
                                     "has_Phenotype_MIM_Number")],
                                 na.rm = T)) %>%
    distinct()

  colnames(MorbidGenes_Panel) = c("symbol", "name", "chromosome", "bed_hg19",
                                  "bed_hg38", "chromosome_location", "transcript",
                                  "id_ncbi", "id_omim", "id_lrg", "id_ensembl", "id_hgnc",
                                  "hgmd_pathogenic_count",
                                  "hgmd_pathogenic_cutoff", "clinvar_pathogenic_count",
                                  "clinvar_pathogenic_cutoff","phenotype","phenotype_mim_numbers",
                                  "mim_numbers", "manually_added", "panelapp",
                                  "panelapp_UK", "panelapp_australia", "sysndd",
                                  "omim_phenotype", "morbidgene", "morbidscore")

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

  assign("MorbidGenes_Panel", MorbidGenes_Panel, envir = .GlobalEnv)

  # Write the gene list for Varvis
  Varvis_GeneList = MorbidGenes_Panel %>%
    filter(morbidgene == T) %>%
    select(symbol)

  write_tsv(Varvis_GeneList,
            paste0(directory, filename, "/Varvis_GeneList.tsv"),
            na = "NA", append = FALSE, col_names = F, escape = "double")

  # Write the Genes into a Varvis-readable Filter structure
  MG_Genes = Varvis_GeneList$symbol
  MG_vector = vector()
  c = 1
  max = length(MG_Genes)
  for(Gene in MG_Genes){
    if(c == 1){
      Genename = paste0("'Gene' #'(^", Gene, "$|")
      MG_vector = c(MG_vector, Genename)
      c = c+1
    }else if(c==max){
      Genename = paste0("^", Gene, "$)'")
      MG_vector = c(MG_vector, Genename)
    }else{
      Genename = paste0("^", Gene, "$|")
      MG_vector = c(MG_vector, Genename)
      c = c+1
    }
  }
  ## save list
  MG = paste(MG_vector,collapse="")
  cat(MG, file = paste0(directory, filename, "/Varvis_Filter.txt"))

}
