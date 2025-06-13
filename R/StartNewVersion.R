#' Start New Morbid Genes Panel Version
#'
#' @param directory Where the files will be saved. Defaults to the Working Group directory where the other panel versions are saved
#' @param version The panel version. Defaults to "Year-Month"
#' @param download_varvis Should the Varvis Gene Management table be downloaded?
#' @param download_HGNC Should the current HGNC nomenclature be downloaded?
#' @param download_OMIM Should the current OMIM files be downloaded?
#' @param download_ClinVar_tsv Should the current ClinVar tsv file be downloaded?
#' @param download_ClinVar_vcf Should the current ClinVar vcf file be downloaded?
#' @param download_GenCC Should the current ClinVar vcf file be downloaded?
#'
#' @return Generates a new folder and downloads data
#' @export
#'
#' @import jsonlite
#' @import stringr
#' @import utils
#' @import vcfR
#' @import readr
#' @import httr
#' @import dplyr
#' @import stringi
#' @import config
#'
#' @details The files needed for the next step will also be added to the global environment.
#' If you do not download the latest datasets, the ones from the previous version
#' will be used.
#'
#' It's best to stick with the default values, otherwise things will get
#' confusing. This script will also automatically add a subversion in case there
#' is a previous version from the same month, to avoid overwriting.
#'
#' @examples
#' \dontrun{
#' StartNewVersion(download_varvis = F, download_HGNC = F,
#' download_OMIM = F, download_ClinVar_tsv = F, download_ClinVar_vcf = F)
#' # This will create a panel directory in the default HUG panel folder,
#' # the version will be "Year-Month" and a subversion
#' }
#'
StartNewVersion = function(directory = "W:/HUG/04 Klinische Genomik/10 Panels/MorbidGenes-Panel/",
                           version = format(Sys.time(), "%Y-%m"),
                           download_varvis = T,
                           download_HGNC = T,
                           download_OMIM = T,
                           download_ClinVar_tsv = T,
                           download_ClinVar_vcf = T,
                           download_GenCC = T){

  ############################ PREFACE ####################################

  # check if a variable is already present in the global environment
  # these will be assigned within this funtion (with "assigned")
  # if yes, confirm to overwrite it

  options(download.file.method="curl",
          download.file.extra = "-k -L")

  variables = c("VarvisGeneManagement",
                "hgnc_complete_set",
                "mim2gene", "mimTitles",
                "genemap2", "morbidmap",
                "tidy_clinvar_vcf_meta",
                "clinvar_tsv_filtered", "genecc_tsv")

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


  # check if a subversion is already present
  # if no, define subversion = 1
  # if yes, add new subversion
  files = list.files(path = directory)

  prevVersion = tail(sort(files[grep("MorbidGenes-Panel", files)]), n = 1)
  prevVersion = paste0(directory, prevVersion)

  subversion = ifelse(identical(files[grep(version, files)], character(0)),
                      1, # subverion = 1 if no subversion present yet
                      # if more than one version is present:
                      ifelse(length(files[grep(version, files)]) == 1,
                             as.numeric(str_sub(files[grep(version, files)], start = -1)) + 1,
                             as.numeric(str_sub(sort(files[grep(version, files)])[-1], start = -1)) + 1))

  # create directory where the files will be saved
  directoryNew = paste0(directory, "MorbidGenes-Panel-v", version, ".", subversion)
  dir.create(directoryNew)
  downloadDir = paste0(directoryNew, "/downloads/")
  dir.create(downloadDir)


  ########################### VARVIS GENE MANAGEMENT ########################

  # download Varvis with the API
  if(download_varvis == T){
    message("\nDownloading Varvis Data\n")

    target=config::get("varvis_target")
    user_name=config::get("varvis_user")

    # parse the password as unicode because of the special characters
    password = stri_unescape_unicode(config::get("varvis_password"))

      # 1) Get CSRF token and session ID to log in
    res = GET(paste0("https://", target, ".varvis.com/authenticate"))
    token = res$headers$`x-csrf-token`
    sessionID = res$cookies$value

      # 2) build json body and update token in another step
    json_body = toJSON(list("_csrf" = token, "username" = user_name,
                            "password" = password), auto_unbox = TRUE)
    res = POST(paste0("https://", target, ".varvis.com/login"),
               body = list("_csrf" = token, "username" = user_name,
                           "password" = password))

      ## get Varvis gene IDs
    res = GET(paste0("https://", target, ".varvis.com/virtual-panel-genes"))
    res_text = content(res, as = "text", type = "text/csv")

    VarvisGeneManagement = fromJSON(res_text)$response
    colnames(VarvisGeneManagement) = toupper(colnames(VarvisGeneManagement))

    write_csv(VarvisGeneManagement, file = paste0(downloadDir, "VarvisGeneManagement.csv"),
              col_names = T)

    assign("VarvisGeneManagement", VarvisGeneManagement, envir = .GlobalEnv)

  } else { # get the VarvisGeneManagement from the previous version
    message(paste0("\nReading Varvis Gene Management file from ", prevVersion, "\n"))

    VarvisGeneManagement = read.csv(file = paste0(prevVersion, "/downloads/VarvisGeneManagement.csv"),
                                    header = T, sep = ",")

    file.copy(from = paste0(prevVersion, "/downloads/VarvisGeneManagement.csv"),
              to = paste0(downloadDir, "VarvisGeneManagement.csv"),
              copy.date = T)

    assign("VarvisGeneManagement", VarvisGeneManagement, envir = .GlobalEnv)
  }


  ############################# HGNC DATA #################################

  # download HGNC files
  if(download_HGNC == T){
    message("\nDownloading HGNC Data\n")

    download.file("https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
                  paste0(downloadDir, "hgnc_complete_set.txt"), quiet=F)
    hgnc_path = downloadDir

  } else {
    message(paste0("\nReading HGNC Data from ", prevVersion, "\n"))

    file.copy(from = paste0(prevVersion, "/downloads/hgnc_complete_set.txt"),
              to = paste0(downloadDir, "hgnc_complete_set.txt"),
              copy.date = T)

    hgnc_path = paste0(prevVersion, "/downloads/")

  }

  hgnc_complete_set = read_delim(paste0(hgnc_path, "hgnc_complete_set.txt"), "\t",
                                  escape_double = FALSE, trim_ws = TRUE,
                                 show_col_types = FALSE)
  assign("hgnc_complete_set", hgnc_complete_set, envir = .GlobalEnv)


  ############################ OMIM DATA ###################################

  # download OMIM files
  if(download_OMIM == T){
    message("\nDownloading OMIM Data\n")

    OMIM_ID = config::get("omim_id")

    download.file("https://omim.org/static/omim/data/mim2gene.txt",
                  paste0(downloadDir, "mim2gene.txt"), quiet=F)
    download.file(paste0("https://data.omim.org/downloads/", OMIM_ID, "/mimTitles.txt"),
                  paste0(downloadDir, "mimTitles.txt"), quiet=F)
    download.file(paste0("https://data.omim.org/downloads/", OMIM_ID, "/genemap2.txt"),
                  paste0(downloadDir, "genemap2.txt"), quiet=F)
    download.file(paste0("https://data.omim.org/downloads/", OMIM_ID, "/morbidmap.txt"),
                  paste0(downloadDir, "morbidmap.txt"), quiet=F)

    Omim_path = downloadDir

  } else {
    message(paste0("\nReading OMIM data from ", prevVersion, "\n"))

    file.copy(from = paste0(prevVersion, "/downloads/mim2gene.txt"),
              to = paste0(downloadDir, "mim2gene.txt"),
              copy.date = T)
    file.copy(from = paste0(prevVersion, "/downloads/mimTitles.txt"),
              to = paste0(downloadDir, "mimTitles.txt"),
              copy.date = T)
    file.copy(from = paste0(prevVersion, "/downloads/genemap2.txt"),
              to = paste0(downloadDir, "genemap2.txt"),
              copy.date = T)
    file.copy(from = paste0(prevVersion, "/downloads/morbidmap.txt"),
              to = paste0(downloadDir, "morbidmap.txt"),
              copy.date = T)

    Omim_path = paste0(prevVersion, "/downloads/")

  }


  mim2gene = read_delim(paste0(Omim_path, "mim2gene.txt"),
                         "\t", escape_double = FALSE, trim_ws = TRUE, comment = "#",
                         col_names = c("MIM_Number", "MIM_Entry_Type",
                                       "Entrez_Gene_ID_NCBI",
                                       "Approved_Gene_Symbol_HGNC",
                                       "Ensembl_Gene_ID_Ensembl"),
                        show_col_types = FALSE)

  mimTitles = read_delim(paste0(Omim_path, "mimTitles.txt"), "\t",
                          escape_double = FALSE, trim_ws = TRUE, comment = "#",
                          col_names = c("Prefix","MIM_Number","Preferred_Title_symbol",
                                        "Alternative_Titles_symbols",
                                        "Included_Titles_ symbols"),
                         show_col_types = FALSE)

  genemap2 = read_delim(paste0(Omim_path, "genemap2.txt"), "\t",
                         escape_double = FALSE, trim_ws = TRUE, comment = "#",
                         col_names = c("Chromosome", "Genomic_Position_Start",
                                       "Genomic_Position_End",
                                       "Cyto_Location", "Computed_Cyto_Location",
                                       "MIM_Number", "Gene_Symbols",
                                       "Gene_Name", "Approved_Symbol", "Entrez_Gene_ID",
                                       "Ensembl_Gene_ID",
                                       "Comments", "Phenotypes", "Mouse_Gene_Symbol_ID"),
                        show_col_types = FALSE)

  morbidmap = read_delim(paste0(Omim_path, "morbidmap.txt"), "\t",
                          escape_double = FALSE, trim_ws = TRUE, comment = "#",
                          col_names = c("Phenotype", "Gene_Symbols",
                                        "MIM_Number", "Cyto_Location"),
                         show_col_types = FALSE)

  assign("mim2gene", mim2gene, envir = .GlobalEnv)
  assign("mimTitles", mimTitles, envir = .GlobalEnv)
  assign("genemap2", genemap2, envir = .GlobalEnv)
  assign("morbidmap", morbidmap, envir = .GlobalEnv)



  ############################### CLINVAR DATA #################################

  # download ClinVar vcf
  if(download_ClinVar_vcf == T){
    message("\nDownloading ClinVar vcf (This might take a while...)\n")

    download.file("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
                  paste0(downloadDir, "/clinvar.vcf.gz"), quiet=F)

    clinvarvcf_path = downloadDir

  } else {
    message(paste0("\nReading ClinVar vcf from ", prevVersion, "\n"))

    file.copy(from = paste0(prevVersion, "/downloads/clinvar.vcf.gz"),
              to = paste0(downloadDir, "clinvar.vcf.gz"),
              copy.date = T)

    clinvarvcf_path = paste0(prevVersion, "/downloads/")

  }

  message("\nReading ClinVar vcf (This might take a while...)\n")

  clinvar_vcf = read.vcfR(paste0(clinvarvcf_path, "clinvar.vcf.gz"), verbose = T,
                          limit = 6e+09)

  tidy_clinvar_vcf_meta = extract_info_tidy(clinvar_vcf, info_fields = NULL,
                                             info_types = TRUE, info_sep = ";")

  assign("tidy_clinvar_vcf_meta", tidy_clinvar_vcf_meta, envir = .GlobalEnv)


  # download ClinVar tsv
  if(download_ClinVar_tsv == T){
    message("\nDownloading ClinVar tsv (This might take a while...)\n")

    download.file("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
                  paste0(downloadDir, "variant_summary.txt.gz"), quiet = F)

    clinvartsv_path = downloadDir

  } else {
    message(paste0("\nReading ClinVar tsv from ", prevVersion, "\n"))

    file.copy(from = paste0(prevVersion, "/downloads/variant_summary.txt.gz"),
              to = paste0(downloadDir, "variant_summary.txt.gz"),
              copy.date = T)

    clinvartsv_path = paste0(prevVersion, "/downloads/")

  }


  message("\nReading ClinVar tsv (This might take a while...)\n")

  clinvar_tsv = read_tsv(paste0(clinvartsv_path, "variant_summary.txt.gz"),
                          col_select = c("#AlleleID", "Name", "GeneSymbol",
                                         "ClinicalSignificance", "Assembly",
                                         "GeneID", "HGNC_ID", "ClinSigSimple",
                                         "Type"),
                         show_col_types = FALSE)

  colnames(clinvar_tsv)[1] = "AlleleID"

  # subset clinvar table to GrCH37
  clinvar_tsv_filtered = clinvar_tsv %>%
    filter(Assembly == "GRCh37")

  assign("clinvar_tsv_filtered", clinvar_tsv_filtered, envir = .GlobalEnv)


  ############################# GENCC DATA #####################################

  if(download_GenCC == T){
    message("\nDownloading GenCC Data\n")

    download.file("https://search.thegencc.org/download/action/submissions-export-tsv",
                  paste0(downloadDir, "GenCC_Data.tsv"), quiet = F)

    gencc_path = downloadDir

  } else {
    message(paste0("\nReading GenCC tsv from ", prevVersion, "\n"))

    file.copy(from = paste0(prevVersion, "/downloads/GenCC_Data.tsv"),
              to = paste0(downloadDir, "GenCC_Data.tsv"),
              copy.date = T)

    gencc_path = paste0(prevVersion, "/downloads/")

  }

  message("\nReading GenCC Data\n")

  gencc = read_tsv(paste0(gencc_path, "GenCC_Data.tsv"),
                         col_select = c("uuid", "gene_curie", "gene_symbol",
                                        "classification_title"),
                         show_col_types = FALSE)

  gencc = gencc %>%
    mutate(hgnc_id = str_remove(gene_curie, "HGNC:"))


  assign("gencc_tsv", gencc, envir = .GlobalEnv)

  ############################# LookUp Table ###############################
  lookup = read_csv2(paste0(prevVersion, "/LookUpTable_FirstOccurrence.csv"))
  assign("lookup", lookup, envir = .GlobalEnv)

  ############################# Pseudogenes ###############################
  pseudogenes_path = config::get("pseudogenes_tsv_path")
  pseudogenes = read_tsv(pseudogenes_path) %>%
    select("Parent name") %>%
    distinct() %>%
    rename("symbol" = "Parent name") %>%
    mutate("has_pseudogene" = T)

  assign("pseudogenes", pseudogenes, envir = .GlobalEnv)


  ############################# FINALE ########################################

  message(paste0("\nSuccess! All files have been successfully downloaded to ",
                 downloadDir, "\n"))
}
