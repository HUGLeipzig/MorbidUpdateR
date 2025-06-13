#' Load all relevant files after building the Morbid Genes Panel in
#' case they are not in your environment anymore
#'
#' @description This is a convenience wrapper script utilizing \code{\link{Add_ManualGenes}},
#' \code{\link{Add_PanelApp}}, \code{\link{Add_ClinVar}}, \code{\link{Add_OMIM}},
#' \code{\link{Add_HGMD}}, \code{\link{Add_HGNC}}.
#'
#' @section Warning:
#' This script will add all the relevant variables to your global environment for convenience reasons
#'
#' @param directory The directory where all the files are saved, e.g. "MorbidGenes-Panel-v2021-12.1
#' @param .ClinVarCutoff ClinVar's pathogenic cutoff value. See \code{\link{Add_ClinVar}} for details
#' @param .HGMDCutoff HGMD's pathogenic cutoff value. See \code{\link{Add_HGMD}} for details
#'
#' @return Loads all files to the global environment
#' @export
#'
#' @import readr
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' Load_Morbidfiles(directory = "path/to/morbidgenesfiles")
#' }

Load_Morbidfiles = function(directory,
                            .ClinVarCutoff = 4, .HGMDCutoff = 4){

  variables = c("VarvisGeneManagement_HGNC", "VarvisGeneManagement_HGNC",
                "HGMD_count_HGNC",
                "morbidmap_reshape",
                "clinvar_tsv_filtered_patho_HGNC",
                "PanelAppGenes",
                "Manual", "SysNDD", "MorbidGenes_Panel")

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

  ############################ DIRECTORY STUFF ############################

  # check if the directory string ends with a "\" or "/" and append if not
  directory = ifelse(grepl("/$|\\$", directory),
                     directory,
                     paste0(directory, "/"))

  downloadDir = paste0(directory, "downloads/")

  message("\nReading Varvis Gene Management File\n")
  VarvisGeneManagement = read_csv(file = paste0(downloadDir, "VarvisGeneManagement.csv"),
                                    col_names = T)
  assign("VarvisGeneManagement", VarvisGeneManagement, envir = .GlobalEnv)

  message("\nreading HGNC\n")
  hgnc_complete_set = read_delim(paste0(downloadDir, "hgnc_complete_set.txt"), "\t",
                                 escape_double = FALSE, trim_ws = TRUE,
                                 show_col_types = FALSE)
  assign("hgnc_complete_set", hgnc_complete_set, envir = .GlobalEnv)

  message("\nreading OMIM Data")
  mim2gene = read_delim(paste0(downloadDir, "mim2gene.txt"),
                        "\t", escape_double = FALSE, trim_ws = TRUE, comment = "#",
                        col_names = c("MIM_Number", "MIM_Entry_Type",
                                      "Entrez_Gene_ID_NCBI",
                                      "Approved_Gene_Symbol_HGNC",
                                      "Ensembl_Gene_ID_Ensembl"),
                        show_col_types = FALSE)

  mimTitles = read_delim(paste0(downloadDir, "mimTitles.txt"), "\t",
                         escape_double = FALSE, trim_ws = TRUE, comment = "#",
                         col_names = c("Prefix","MIM_Number","Preferred_Title_symbol",
                                       "Alternative_Titles_symbols",
                                       "Included_Titles_ symbols"),
                         show_col_types = FALSE)

  genemap2 = read_delim(paste0(downloadDir, "genemap2.txt"), "\t",
                        escape_double = FALSE, trim_ws = TRUE, comment = "#",
                        col_names = c("Chromosome", "Genomic_Position_Start",
                                      "Genomic_Position_End",
                                      "Cyto_Location", "Computed_Cyto_Location",
                                      "MIM_Number", "Gene_Symbols",
                                      "Gene_Name", "Approved_Symbol", "Entrez_Gene_ID",
                                      "Ensembl_Gene_ID",
                                      "Comments", "Phenotypes", "Mouse_Gene_Symbol_ID"),
                        show_col_types = FALSE)

  morbidmap = read_delim(paste0(downloadDir, "morbidmap.txt"), "\t",
                         escape_double = FALSE, trim_ws = TRUE, comment = "#",
                         col_names = c("Phenotype", "Gene_Symbols",
                                       "MIM_Number", "Cyto_Location"),
                         show_col_types = FALSE)

  assign("mim2gene", mim2gene, envir = .GlobalEnv)
  assign("mimTitles", mimTitles, envir = .GlobalEnv)
  assign("genemap2", genemap2, envir = .GlobalEnv)
  assign("morbidmap", morbidmap, envir = .GlobalEnv)

  message("\nReading ClinVar vcf (This might take a while...)\n")

  clinvar_vcf = read.vcfR(paste0(downloadDir, "clinvar.vcf.gz"), verbose = T)

  tidy_clinvar_vcf_meta = extract_info_tidy(clinvar_vcf, info_fields = NULL,
                                            info_types = TRUE, info_sep = ";")

  assign("tidy_clinvar_vcf_meta", tidy_clinvar_vcf_meta, envir = .GlobalEnv)

  message("\nReading ClinVar tsv (This might take a while...)\n")

  clinvar_tsv = read_tsv(paste0(downloadDir, "variant_summary.txt.gz"),
                         col_select = c("#AlleleID", "Name", "GeneSymbol",
                                        "ClinicalSignificance", "Assembly",
                                        "GeneID", "HGNC_ID"),
                         show_col_types = FALSE)

  colnames(clinvar_tsv)[1] = "AlleleID"

  # subset clinvar table to GrCH37
  clinvar_tsv_filtered = clinvar_tsv %>%
    filter(Assembly == "GRCh37")

  assign("clinvar_tsv_filtered", clinvar_tsv_filtered, envir = .GlobalEnv)

  Add_all(add_coordinates = F,
          ClinVarCutoff = .ClinVarCutoff, HGMDCutoff = .HGMDCutoff,
          downloadPanelApp = F)


  message("\nReading MorbidGenesPanel\n")
  files = list.files(directory)
  Morbidfilename = files[grep("MorbidGenes-Panel", files)]
  MorbidGenes_Panel = read_csv2(file = paste0(directory, Morbidfilename),
                                col_names = T)
  assign("MorbidGenes_Panel", MorbidGenes_Panel, envir = .GlobalEnv)
}
