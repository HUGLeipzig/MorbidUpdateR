#' Reshape the morbidmap to be added to the Panel
#'
#' @param .morbidmap The morbidmap file from OMIM.
#' Should be in your global environment after running \code{\link{StartNewVersion}}
#' @param .mim2gene The mim2gene file from OMIM.
#' Should be in your global environment after running \code{\link{StartNewVersion}}
#'
#' @return A dataframe with the OMIM data
#' @export
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' \dontrun{
#' morbidmap_reshape = Add_OMIM(.morbidmap = morbidmap, .mim2gene = mim2gene)
#' }


Add_OMIM = function(.morbidmap = morbidmap, .mim2gene = mim2gene){


  morbidmap_reshape <- morbidmap %>%
    left_join(mim2gene, by = "MIM_Number") %>%
    separate(Phenotype, c("Phenotype", "phenotype_mapping_key"),
                    sep = " \\((?=[0-9]\\)$)") %>%
    mutate(phenotype_mapping_key = str_replace(phenotype_mapping_key, "\\)",
                                                      "")) %>%
    separate(Phenotype, c("Phenotype", "Phenotype_MIM_Number"),
                    sep = "\\, (?=[0-9]{6}$)") %>%
    mutate(OMIM_association_doubtful = grepl("(\\?|\\[|\\{)", Phenotype)) %>%
    select(Phenotype, Phenotype_MIM_Number, phenotype_mapping_key,
                  Approved_Gene_Symbol_HGNC, MIM_Number, MIM_Entry_Type,
                  OMIM_association_doubtful) %>%
    filter(phenotype_mapping_key == 3, OMIM_association_doubtful == FALSE,
                  !is.na(Approved_Gene_Symbol_HGNC)) %>%
    select(Phenotype, Phenotype_MIM_Number,
                  Approved_Gene_Symbol_HGNC, MIM_Number) %>%
    group_by(Approved_Gene_Symbol_HGNC) %>%
    summarise(Phenotype = paste(Phenotype, collapse=","),
                     Phenotype_MIM_Numbers = paste(Phenotype_MIM_Number, collapse=","),
                     MIM_Numbers = paste(unique(MIM_Number), collapse=",")) %>%
    ungroup()

  return(morbidmap_reshape)
}
