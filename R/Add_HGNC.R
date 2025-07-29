#' Add the HGNC Data and genomic positions (GRCh38) to the Varvis Gene Management file
#'
#' @param .VarvisGeneManagement The Varvis Gene Management file.
#' Should be in your global environment after running \code{\link{StartNewVersion}}
#' @param directory The same you provided in \code{\link{StartNewVersion}}
#' @param version The same you provided in \code{\link{StartNewVersion}}
#' @param add_coordinates Add genomic coordinates from Ensembl. Can be set to F if Ensembl is down
#'
#' @return A dataframe with the HGNC Data
#' @export
#'
#' @import dplyr
#' @import biomaRt
#'
#' @examples
#' \dontrun{
#' VarvisGeneManagement_HGNC = Add_HGNC(.VarvisGeneManagement = VarvisGeneManagement)
#' }

Add_HGNC = function(.VarvisGeneManagement = VarvisGeneManagement,
                    directory = "W:/HUG/04 Klinische Genomik/10 Panels/MorbidGenes-Panel/",
                    version = format(Sys.time(), "%Y-%m"),
                    add_coordinates = T){

  ##### Directory Stuff #####
  # check if the directory string ends with a "\" or "/" and append if not
  directory = ifelse(grepl("/$|\\$", directory),
                     directory,
                     paste0(directory, "/"))

  # check if a subversion is already present
  # if no, define subversion = 1
  # if yes, add new subversion
  files = list.files(path = directory)

  filename = tail(sort(files[grep(version, files)]), n = 1)

  ##### get genomic position via Ensembl/BiomaRt #####
  ##### define functions #####

  httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))

  if(add_coordinates == T){
    #mart_hg19 <- useMart("ensembl", host="https://grch37.ensembl.org")
    #mart_hg19 <- useDataset("hsapiens_gene_ensembl", mart_hg19)

    mart_hg38 <- useMart("ensembl", host="https://www.ensembl.org")
    mart_hg38 <- useDataset("hsapiens_gene_ensembl", mart_hg38)

    # function to retrive bed format style gene coordinates
    gene_coordinates_from_symbol <- function(gene_symbols, reference = "hg38") {
      gene_symbol_list <- as_tibble(gene_symbols) %>%
        dplyr::select(hgnc_symbol = value)

      if (reference == "hg19") {
        mart <- mart_hg19
      } else {
        mart <- mart_hg38
      }

      attributes <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position")
      filters <- c("hgnc_symbol")

      values <- list(hgnc_symbol = gene_symbol_list$hgnc_symbol)

      gene_coordinates_hg19 <- getBM(attributes=attributes, filters=filters, values=values, mart=mart) %>%
        group_by(hgnc_symbol) %>%
        summarise(hgnc_symbol = max(hgnc_symbol), chromosome_name = max(chromosome_name), start_position = max(start_position), end_position = max(end_position)) %>%
        mutate(bed_format = paste0("chr", chromosome_name, ":", start_position, "-", end_position)) %>%
        dplyr::select(hgnc_symbol, bed_format)

      gene_symbol_list_return <- gene_symbol_list %>%
        left_join(gene_coordinates_hg19, by = ("hgnc_symbol"))

      return(gene_symbol_list_return)
    }

    #
    gene_coordinates_from_ensembl <- function(ensembl_id, reference = "hg38") {
      ensembl_id_list <- as_tibble(ensembl_id) %>%
        dplyr::select(ensembl_gene_id = value)

      if (reference == "hg19") {
        mart <- mart_hg19
      } else {
        mart <- mart_hg38
      }

      attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position")
      filters <- c("ensembl_gene_id")

      values <- list(ensembl_gene_id = ensembl_id_list$ensembl_gene_id)

      gene_coordinates_hg19 <- getBM(attributes=attributes, filters=filters, values=values, mart=mart) %>%
        group_by(ensembl_gene_id) %>%
        summarise(ensembl_gene_id = max(ensembl_gene_id), chromosome_name = max(chromosome_name), start_position = max(start_position), end_position = max(end_position)) %>%
        mutate(bed_format = paste0("chr", chromosome_name, ":", start_position, "-", end_position)) %>%
        dplyr::select(ensembl_gene_id, bed_format)

      ensembl_id_list_return <- ensembl_id_list %>%
        left_join(gene_coordinates_hg19, by = ("ensembl_gene_id"))

      return(ensembl_id_list_return)
    }

    non_alt_loci_set_coordinates <- hgnc_complete_set %>%
      #mutate(hg19_coordinates_from_ensembl = gene_coordinates_from_ensembl(ensembl_gene_id)) %>%
      #mutate(hg19_coordinates_from_symbol = gene_coordinates_from_symbol(symbol)) %>%
      mutate(hg19_coordinates_from_ensembl = NA) %>%
      mutate(hg19_coordinates_from_symbol = NA) %>%
      mutate(hg38_coordinates_from_ensembl = gene_coordinates_from_ensembl(ensembl_gene_id, reference = "hg38")) %>%
      mutate(hg38_coordinates_from_symbol = gene_coordinates_from_symbol(symbol, reference = "hg38")) %>%
      #mutate(bed_hg19 =
      #         case_when(
      #           !is.na(hg19_coordinates_from_ensembl$bed_format) ~ hg19_coordinates_from_ensembl$bed_format,
      #           is.na(hg19_coordinates_from_ensembl$bed_format) ~ hg19_coordinates_from_symbol$bed_format,
      #         )
      #) %>%
      mutate(bed_hg19 = NA) %>%
      mutate(bed_hg38 =
               case_when(
                 !is.na(hg38_coordinates_from_ensembl$bed_format) ~ hg38_coordinates_from_ensembl$bed_format,
                 is.na(hg38_coordinates_from_ensembl$bed_format) ~ hg38_coordinates_from_symbol$bed_format,
               )
      ) %>%
      dplyr::select(-hg19_coordinates_from_ensembl, -hg19_coordinates_from_symbol, -hg38_coordinates_from_ensembl, -hg38_coordinates_from_symbol)

  } else {

    non_alt_loci_set_coordinates <- hgnc_complete_set %>%
      mutate(bed_hg19 = NA) %>%
      mutate(bed_hg38 = NA)
  }


  mb_genes_hgnc_connect <- non_alt_loci_set_coordinates %>%
    dplyr::select(hgnc_id) %>%
    mutate(is_active = TRUE)

  ############################################



  ############################################
  ## export table as csv with date of creation
  creation_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

  write_csv(non_alt_loci_set_coordinates, file = paste0(directory, filename, "/non_alt_loci_set_coordinates.",creation_date,".csv"))

  write_csv(mb_genes_hgnc_connect, file = paste0(directory, filename,"/mb_genes_hgnc_connect.",creation_date,".csv"))

  ############################################

  VarvisGeneManagement_HGNC <- VarvisGeneManagement %>%
    left_join(non_alt_loci_set_coordinates, by = c("NCBIID" = "entrez_id")) %>%
    dplyr::select("SYMBOL","NAME","CHROMOSOME","CHROMOSOMELOCATION","TRANSCRIPT",
           "NCBIID","OMIMID","LRGID","ENSEMBLID","HGNCID","symbol", "bed_hg19", "bed_hg38") %>%
    mutate(HGNC_symbol_corrected =
             case_when(
               is.na(symbol) ~ SYMBOL,
               !is.na(symbol) ~ symbol,
             )
    )

  return(VarvisGeneManagement_HGNC)
}

