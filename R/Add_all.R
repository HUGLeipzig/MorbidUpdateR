#' Run all relevant scripts before you build the Morbid Genes Panel
#'
#' @description This is a convenience wrapper script utilizing \code{\link{Add_ManualGenes}},
#' \code{\link{Add_PanelApp}}, \code{\link{Add_ClinVar}}, \code{\link{Add_OMIM}},
#' \code{\link{Add_HGMD}}, \code{\link{Add_HGNC}}.
#'
#' @section Warning:
#' This script will add all the relevant variables to your global environment for convenience reasons
#'
#' @param ClinVarCutoff ClinVar's pathogenic cutoff value. See \code{\link{Add_ClinVar}} for details
#' @param HGMDCutoff HGMD's pathogenic cutoff value. See \code{\link{Add_HGMD}} for details
#' @param .VarvisGeneManagement See \code{\link{Add_HGNC}} for details
#' @param .morbidmap See \code{\link{Add_HGNC}} for details
#' @param .mim2gene See \code{\link{Add_OMIM}} for details
#' @param PanelAppFile See \code{\link{Add_PanelApp}} for details
#' @param ManualFile See \code{\link{Add_ManualGenes}} for details
#' @param ... Additional arguments passed on to \code{\link{Add_ClinVar}}
#'
#' @return
#' @export
#'
#' @examples
#' Add_all(ClinVarCutoff = 4, mim2gene = mim2gene)
Add_all = function(ClinVarCutoff = 4,
                   HGMDCutoff = 4,
                   .VarvisGeneManagement = VarvisGeneManagement,
                   .morbidmap = morbidmap,
                   .mim2gene = mim2gene,
                   panelAppFile = "W:/HUG/04 Klinische Genomik/10 Panels/MorbidGenes-Panel/PanelAppGenes/2021_08_25_PA_all_genes.csv",
                   manualFile = "W:/HUG/04 Klinische Genomik/10 Panels/MorbidGenes-Panel/GenesToBeAddedManually.xlsx",
                   ...){

  variables = c("VarvisGeneManagement_HGNC",
                "HGMD_count_HGNC",
                "morbidmap_reshape",
                "clinvar_tsv_filtered_patho_HGNC",
                "PanelAppGenes",
                "Manual")

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

  VarvisGeneManagement_HGNC = Add_HGNC()
  assign("VarvisGeneManagement_HGNC", VarvisGeneManagement_HGNC, envir = .GlobalEnv)

  HGMD_count_HGNC = Add_HGMD(cutoff = HGMDCutoff)
  assign("HGMD_count_HGNC", HGMD_count_HGNC, envir = .GlobalEnv)

  morbidmap_reshape = Add_OMIM(.morbidmap = morbidmap, .mim2gene = mim2gene)
  assign("morbidmap_reshape", morbidmap_reshape, envir = .GlobalEnv)

  clinvar_tsv_filtered_patho_HGNC = Add_ClinVar(cutoff = ClinVarCutoff, .mim2gene = mim2gene)
  assign("clinvar_tsv_filtered_patho_HGNC", clinvar_tsv_filtered_patho_HGNC, envir = .GlobalEnv)

  PanelAppGenes = Add_PanelApp(PanelAppFile = panelAppFile)
  assign("PanelAppGenes", PanelAppGenes, envir = .GlobalEnv)

  Manual = Add_ManualGenes(ManualFile = manualFile)
  assign("Manual", Manual, envir = .GlobalEnv)

}
