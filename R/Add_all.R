#' Run all relevant scripts before you build the Morbid Genes Panel
#'
#' @description This is a convenience wrapper script utilizing \code{\link{Add_ManualGenes}},
#' \code{\link{Add_PanelApp}}, \code{\link{Add_ClinVar}}, \code{\link{Add_OMIM}},
#' \code{\link{Add_HGMD}}, \code{\link{Add_HGNC}} and \code{\link{Add_Pathomechanisms}}.
#'
#' @section Warning:
#' This script will add all the relevant variables to your global environment for convenience reasons
#'
#' @param add_coordinates Add genomic coordinates from Ensembl, see \code{\link{Add_HGNC}} for details
#' @param ClinVarCutoff ClinVar's pathogenic cutoff value. See \code{\link{Add_ClinVar}} for details
#' @param HGMDCutoff HGMD's pathogenic cutoff value. See \code{\link{Add_HGMD}} for details
#' @param .morbidmap See \code{\link{Add_HGNC}} for details
#' @param .mim2gene See \code{\link{Add_OMIM}} for details
#' @param downloadPanelApp See \code{\link{Add_PanelApp}} for details
#' @param PanelApp_Panel See \code{\link{Add_PanelApp}} for details
#' @param PanelApp_confidence See \code{\link{Add_PanelApp}} for details
#' @param SysNDD_category See \code{\link{Add_SysNDD}} for details
#' @param SysNDD_inheritance See \code{\link{Add_SysNDD}} for details
#' @param download_SysNDD See \code{\link{Add_SysNDD}} for details
#' @param manualFile See \code{\link{Add_ManualGenes}} for details
#' @param .gencc_tsv See \code{\link{Add_GenCC}} for details
#' @param gencc_classification See \code{\link{Add_GenCC}} for details
#' @param .clinvar_tsv_filtered2 See \code{\link{Add_Pathomechanisms}} for details
#' @param ... Additional arguments passed on to \code{\link{Add_ClinVar}}
#'
#' @return Assigns all values to the global environment
#' @export
#'
#' @examples
#' \dontrun{
#' Add_all(ClinVarCutoff = 4, mim2gene = mim2gene)
#' }

Add_all = function(add_coordinates = T,
                   ClinVarCutoff = 4,
                   HGMDCutoff = 4,
                   .morbidmap = morbidmap,
                   .mim2gene = mim2gene,
                   PanelApp_Panel = "both",
                   PanelApp_confidence = 3,
                   downloadPanelApp = T,
                   SysNDD_category = "Definitive",
                   SysNDD_inheritance = "All",
                   download_SysNDD = T,
                   .gencc_tsv = gencc_tsv,
                   gencc_classification = "Definitive",
                   .clinvar_tsv_filtered2 = clinvar_tsv_filtered,
                   manualFile = "W:/HUG/04 Klinische Genomik/10 Panels/MorbidGenes-Panel/GenesToBeAddedManually.xlsx",
                   ...){

  variables = c("VarvisGeneManagement_HGNC",
                "HGMD_count_HGNC",
                "morbidmap_reshape",
                "clinvar_tsv_filtered_patho_HGNC",
                "PanelAppGenes",
                "Manual", "SysNDD")

  variablePresent = TRUE %in% (variables %in% ls(envir = .GlobalEnv))

  if(variablePresent == TRUE){
    answer = askYesNo(msg = paste0("Beware! One or more of your global variables will be overwritten. The variables are: ",
                                   paste(variables, collapse = ", "),
                                   ". Do you wish to continue and overwrite your current variables? "),
                      prompts = c("yes", "no", "hell no!"))

    if(answer != TRUE){
      stop("Clear your global environment from the abovementioned variables or consent to overwriting")
    }
  }

  message("\nImporting Varvis Gene Management File\n")
  VarvisGeneManagement_HGNC = Add_HGNC(add_coordinates = add_coordinates)
  assign("VarvisGeneManagement_HGNC", VarvisGeneManagement_HGNC, envir = .GlobalEnv)

  message("\nAdding HGMD Data\n")
  HGMD_count_HGNC = Add_HGMD(cutoff = HGMDCutoff)
  assign("HGMD_count_HGNC", HGMD_count_HGNC, envir = .GlobalEnv)

  message("\nAdding OMIM Data\n")
  morbidmap_reshape = Add_OMIM(.morbidmap = morbidmap, .mim2gene = mim2gene)
  assign("morbidmap_reshape", morbidmap_reshape, envir = .GlobalEnv)

  message("\nAdding ClinVar Data\n")
  clinvar_tsv_filtered_patho_HGNC = Add_ClinVar(cutoff = ClinVarCutoff, .mim2gene = mim2gene)
  assign("clinvar_tsv_filtered_patho_HGNC", clinvar_tsv_filtered_patho_HGNC, envir = .GlobalEnv)

  message("\nAdding PanelApp Genes\n")
  PanelAppGenes = Add_PanelApp(Panel = PanelApp_Panel, confidence = PanelApp_confidence,
                               download = downloadPanelApp)
  assign("PanelAppGenes", PanelAppGenes, envir = .GlobalEnv)

  message("\nAdding Manual Genes\n")
  Manual = Add_ManualGenes(ManualFile = manualFile)
  assign("Manual", Manual, envir = .GlobalEnv)

  message("\nAdding SysNDD Genes\n")
  SysNDD = Add_SysNDD(category = SysNDD_category, inheritance = SysNDD_inheritance,
                      download = download_SysNDD)
  assign("SysNDD", SysNDD, envir = .GlobalEnv)

  message("\nAdding GenCC Genes\n")
  gencc = Add_GenCC(.gencc_tsv = gencc_tsv, classification = gencc_classification)
  assign("gencc", gencc, envir = .GlobalEnv)

  message("\nAdding Pathomechanisms\n")
  PathoClinVarVars = Add_Pathomechanisms(.clinvar_tsv_filtered = .clinvar_tsv_filtered2)
  assign("PathoClinVarVars", PathoClinVarVars, envir = .GlobalEnv)

}
