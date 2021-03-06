% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Enrichments.R
\name{get_hm_enrichments}
\alias{get_hm_enrichments}
\title{Get dataframe for plotting heatmap}
\usage{
get_hm_enrichments(exprs_data, use_markers, groups, summarize_mode = "median")
}
\arguments{
\item{exprs_data}{matrix with enrichment data}

\item{use_markers}{list of used markers}

\item{groups}{vector with grouped info for each cell (clusters, samples or other)}

\item{summarize_mode}{method to calculate enrichment of cells groups(median or mean)}
}
\value{

}
\description{
Get dataframe for plotting heatmap
}
