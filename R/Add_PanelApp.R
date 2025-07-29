#' Add the PanelApp Genes
#'
#' @param Panel One of "UK", "AUS" or "both" (default). This is to specify which PanelApp to use
#' @param confidence The minimum confidence of genes to keep. Can be 1 (minimum red genes),
#' 2 (minimum amber genes) or the default 3 (minimum green genes)
#' @param download Should the latest data be downloaded? If not, the previous versions will be read
#'
#' @return A dataframe with the PanelApp genes
#' @export
#'
#' @import dplyr
#' @import readr
#' @import httr
#' @import jsonlite
#' @import config
#'
#' @examples
#' \dontrun{
#' PanelAppGenes = Add_PanelApp()
#' }

Add_PanelApp = function(Panel = "both", confidence = 3,
                        download = T){

  httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))

  if(download == T){
    # Determine which PanelApp to use (UK or Australia)
    baseURL_UK = "panelapp.genomicsengland.co.uk"
    baseURL_AUS = "panelapp-aus.org"

    # confidence will be the gene levels, i.e. 3=green, 2=amber, 1=red
    .confidence = confidence

    # this function will get the Panels and then the Genes from the corresponding PanelApp
    PanelAppAPI = function(baseURL, confidence = .confidence){
      # The results come in pages, so start with page 1
      print("Reading page 1")
      panelres = GET(paste0(baseURL, "/api/v1/panels/?page=1"))
      panelres_text = suppressMessages(content(panelres, as = "text", type = "text/csv"))

      panelTest = fromJSON(panelres_text)

      # create a vector with the panel IDs which will be appended later and used to get the genes
      IDs = panelTest$results$id

      page = 2

      # This while loop checks if there will be a next page.
      # If not, the last page is reached and the loop stops
      while(!is.null(panelTest$`next`)){
        print(paste0("Reading page ", page))
        panelres = GET(paste0(baseURL, "/api/v1/panels/?page=", page))
        panelres_text = suppressMessages(content(panelres, as = "text", type = "text/csv"))
        panelTest = fromJSON(panelres_text)

        IDs = c(IDs, panelTest$results$id)

        page = page + 1
      }

      # add a counter for convenience printing as a status bar
      i = 0

      # iterate over the IDs and get the corresponding genes
      for(ID in IDs){
        panelres = GET(paste0(baseURL, "/api/v1/panels/", ID))
        panelres_text = suppressMessages(content(panelres, as = "text", type = "text/csv"))
        panelTest = fromJSON(panelres_text)
        i = i+1

        # determine which PanelApp is used
        message_Country = ifelse(grepl(".uk", baseURL),
                                 " UK Panels", " Australia Panels")
        print(paste0("Reading Panel ", i, " of ", length(IDs), message_Country))

        # some panels do not contain genes but CNVs, skip these here
        if(length(panelTest$genes) == 0){
          next
        } else {
          # in the first loop, create a dataframe where the genes will be stored
          if(i == 1){
            PanelAppGenes = panelTest$genes %>%
              select(entity_name, confidence_level) %>%
              # filter for the confidence level (minimum)
              filter(confidence_level >= confidence)
          } else {
            # after the first loop, append the genes to the df created in the first loop
            a = panelTest$genes %>%
              select(entity_name, confidence_level) %>%
              filter(confidence_level >= confidence)

            PanelAppGenes = rbind(PanelAppGenes, a)
          }
        }
      }

      # genes can be in multiple panesl, so keep only distinct genes
      PanelAppGenes = PanelAppGenes %>% distinct()
      return(PanelAppGenes)
    }

    # depending on the specified country, download the data from the corresponding PanelApp
    if(Panel == "UK"){
      PanelAppResult = PanelAppAPI(baseURL = baseURL_UK) %>%
        mutate(isUKPanelAppGene = TRUE)
    } else if (Panel == "AUS"){
      PanelAppResult = PanelAppAPI(baseURL = baseURL_AUS) %>%
        mutate(isAustraliaPanelAppGene = TRUE)

      # if "both" is specified, download both and combine them here
    } else if (Panel == "both"){
      PanelAppResult_UK = PanelAppAPI(baseURL = baseURL_UK) %>%
        mutate(isUKPanelAppGene = TRUE)
      PanelAppResult_AUS = PanelAppAPI(baseURL = baseURL_AUS) %>%
        mutate(isAustraliaPanelAppGene = TRUE)
      PanelAppResult = full_join(PanelAppResult_UK, PanelAppResult_AUS)
    }

    PanelAppResult = PanelAppResult %>%
      mutate(isPanelAppGene = TRUE) %>%
      select(-confidence_level)

    write_tsv(PanelAppResult, config::get("panelapp_tsv_path"),
              col_names = T)
  } else {

    PanelAppResult = read_tsv(config::get("panelapp_tsv_path"),
                              col_names = T)

  }

  return(PanelAppResult)

}
