StartNewVersion = function(directory = "W:/HUG/04 Klinische Genomik/02 Dokumentation/SNV_CNV_Doku/01 Panelzusammenstellung/MorbidGenes-Panel/",
                           version = format(Sys.time(), "%Y-%m")){
  directory = ifelse(grepl("/$", directory),
                     directory,
                     paste0(directory, "/"))
  dir.create(paste0(directory, "MorbidGenes-Panel-v", version))
  directoryNew = paste0(directory, "MorbidGenes-Panel-v", version)
  #setwd(paste0(directory, "MorbidGenes-Panel-v", version))
  dir.create(paste0(directoryNew, "/downloads"))
  dir.create(paste0(directoryNew, "/_OtherFiles"))
  print(paste0("Success! Please download the Varvis Gene Management file and copy it into ",
              directoryNew, "/_OtherFiles"))
}
