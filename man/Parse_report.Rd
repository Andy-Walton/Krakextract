% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{parse_report}
\alias{parse_report}
\title{parse_report}
\usage{
parse_report(
  file,
  sci_names,
  rank_method = "reads",
  rank_no = 5,
  min_reads = 0,
  out_type = "names",
  outdir
)
}
\arguments{
\item{file}{This is the file path to your Kraken report file}

\item{sci_names}{This is a vector of the form c("Lactobacillaceae", "Streptococcaceae") ONLY FAMILY NAMES WORK}

\item{rank_method}{Sort the kraken hits by percentage or number of reads args = "reads" or "percent"}

\item{rank_no}{The top x number of species (default is 5)}

\item{min_reads}{The minimum number of reads for a species to be included (default is 0)}

\item{out_type}{Specifies what you want to output args: IDs or names}

\item{outdir}{The filepath you want to save your text files to}
}
\value{
Text files with NCBI IDs for each species of interest above thresholds in your kraken output
}
\description{
Extracts species of interest from Kraken reports.

License: MIT + file LICENSE
}
\examples{
\dontrun{
families = c("Lactobacillaceae", "Aerococcaceae")
path_out = "man/Example_Data/"

parse_report("man/Example_Data/Example.kreport",
sci_names=families,
rank_method="reads",
rank_no=10,
out_type="IDs",
outdir=path_out)}

}
