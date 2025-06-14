% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Add_Pathomechanisms.R
\name{Add_Pathomechanisms}
\alias{Add_Pathomechanisms}
\title{Add the pathogenic variant types from Clinvar as well as STRs}
\usage{
Add_Pathomechanisms(.clinvar_tsv_filtered = clinvar_tsv_filtered)
}
\arguments{
\item{.clinvar_tsv_filtered}{ClinVar's filtered .tsv file.
Should be in your global environment after running \code{\link{StartNewVersion}}}
}
\value{
A dataframe with pathogenic variant types per gene
}
\description{
Add the pathogenic variant types from Clinvar as well as STRs
}
\examples{
\dontrun{
PathoClinVarVars = Add_Pathomechanisms(.clinvar_tsv_filtered = clinvar_tsv_filtered)
}
}
