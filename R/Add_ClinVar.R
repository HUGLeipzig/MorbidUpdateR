#' Add the data from both Clinvar files
#'
#' @param .clinvar_tsv_filtered ClinVar's filtered .tsv file.
#' Should be in your global environment after running \code{\link{StartNewVersion}}
#' @param .tidy_clinvar_vcf_meta ClinVar's filtered .vcf file.
#' Should be in your global environment after running \code{\link{StartNewVersion}}
#' @param .mim2gene mim2gene file from OMIM.
#' Should be in your global environment after running \code{\link{StartNewVersion}}
#' @param cutoff minimum number of pathogenic variants per gene to be included into the Morbid Genes panel
#' @param print_NonMatching Should the Transcriptgenes be printed which do not match the HGNC nomenclature?
#'
#' @return A dataframe with the ClinVar genes
#' @export
#'
#' @import dplyr
#' @import stringr
#' @importFrom purrr discard
#'
#' @examples
#' \dontrun{
#' clinvar_tsv_filtered_patho_HGNC = Add_ClinVar(cutoff = 4)
#' }

Add_ClinVar = function(.clinvar_tsv_filtered = clinvar_tsv_filtered,
                       .tidy_clinvar_vcf_meta = tidy_clinvar_vcf_meta,
                       .mim2gene = mim2gene,
                       cutoff = 4,
                       print_NonMatching = T){

  # identify multiple genes
  tidy_clinvar_vcf_meta_multiGene =
    tidy_clinvar_vcf_meta[str_count(pattern = ":", tidy_clinvar_vcf_meta$GENEINFO) > 1, ]

  tidy_clinvar_vcf_meta_multiGene = tidy_clinvar_vcf_meta_multiGene[!is.na(tidy_clinvar_vcf_meta_multiGene$ALLELEID), ]

  clinvar_tsv_filtered$TranscriptGene = str_extract(pattern = "\\([A-Za-z][A-Za-z0-9-_]+\\)",
                                                    clinvar_tsv_filtered$Name) %>%
    str_extract(pattern = "[A-Za-z0-9-_]+")

  # filter tsv
  clinvar_tsv_filtered_noGeneID = clinvar_tsv_filtered %>%
    filter(grepl(";", GeneSymbol)) %>%
    filter(str_detect(ClinicalSignificance,
                      "((?!>(known_|non-))[Pp]athogenic(?![A-z]))")) %>%
    filter(!is.na(TranscriptGene))

  clinvar_tsv_filtered_GeneID = clinvar_tsv_filtered %>%
    filter(GeneID != -1) %>%
    filter(str_detect(ClinicalSignificance,
                      "((?!>(known_|non-))[Pp]athogenic(?![A-z]))")) %>%
    filter(!is.na(TranscriptGene)) %>%
    mutate(TranscriptGene = if_else(TranscriptGene == GeneSymbol,
                                    TranscriptGene,
                                    GeneSymbol))

  # look for corresponding gene in vcf
  get_GeneID = function(Gene){
    Gene = paste0(Gene, ":")
    if(identical(tidy_clinvar_vcf_meta$GENEINFO[grepl(Gene, tidy_clinvar_vcf_meta$GENEINFO)],
                 character(0))){
      as.numeric(-1)
    } else if(identical(tidy_clinvar_vcf_meta$GENEINFO[grepl(Gene, tidy_clinvar_vcf_meta$GENEINFO)],
                        numeric(0))){
      as.numeric(-1)
    } else{as.numeric(unique(tidy_clinvar_vcf_meta$GENEINFO[grepl(Gene, tidy_clinvar_vcf_meta$GENEINFO)]) %>%
                        str_split("\\|", simplify = T) %>%
                        str_extract(paste0("^", Gene, "[0-9]+")) %>%
                        discard(is.na) %>%
                        unique() %>%
                        str_extract("[0-9]+$") %>%
                        discard(is.na))}

  }

  # run the GetGeneID function on all Genes which are missing the Gene ID
  # this might take a while

  clinvar_tsv_filtered_noGeneID$GeneID = sapply(clinvar_tsv_filtered_noGeneID$TranscriptGene, get_GeneID, simplify = "vector", USE.NAMES = F)

  #clinvar_tsv_filtered_noGeneID = clinvar_tsv_filtered_noGeneID %>%
  #  filter(GeneID == -1)

  # combine clinvar tables
  clinvar_tsv_filtered_patho = rbind(clinvar_tsv_filtered_GeneID,
                                     clinvar_tsv_filtered_noGeneID)


  clinvar_tsv_filtered_patho <- clinvar_tsv_filtered_patho %>%
    select(AlleleID, Name, ClinicalSignificance, GeneID, TranscriptGene) %>%
    add_count(GeneID) %>%
    group_by(GeneID) %>%
    summarise(CinSigs = paste(ClinicalSignificance, collapse=","),
                     ALLELEIDs = paste(AlleleID, collapse=","),
                     NcbiGeneID = unique(GeneID), ClinVarPathogenicCount = unique(n)) %>%
    ungroup() %>%
    #filter(!is.na(TranscriptGene)) %>%
    mutate(ClinVarPathogenicCount_cutoff = (ClinVarPathogenicCount >= cutoff))

  clinvar_tsv_filtered_patho$NcbiGeneID = as.double(clinvar_tsv_filtered_patho$NcbiGeneID)

  # add hgnc and omim
  clinvar_tsv_filtered_patho_HGNC <- clinvar_tsv_filtered_patho %>%
    left_join((mim2gene %>%
                 select(Entrez_Gene_ID_NCBI, Approved_Gene_Symbol_HGNC)),
              by = c("NcbiGeneID" = "Entrez_Gene_ID_NCBI"))


  return(clinvar_tsv_filtered_patho_HGNC)
}
