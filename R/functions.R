#' @title Parse_report
#'
#' @description Extracts species of interest from Kraken reports
#'
#' License: MIT + file LICENSE
#'
#' @param file This is the file path to your Kraken report file
#' @param filetype This is the type of file, only "Kraken" works at this point
#' @param sci_names This is a vector of the form c("Lactobacillaceae", "Streptococcaceae") ONLY FAMILY NAMES WORK
#' @param rank_method Sort the kraken hits by percentage or number of reads args = "reads" or "percent"
#' @param rank_no The top x number of species (default is 5)
#' @param min_reads The minimum number of reads for a species to be included (default is 0)
#' @param outdir The filepath you want to save your text files to
#'
#' @return Text files with NCBI IDs for each species of interest above thresholds in your kraken output
#'
#' @examples
#' \dontrun{
#' families = c("Lactobacillaceae", "Aerococcaceae")
#' path_out = "man/Example_Data/"
#'
#' Parse_report("man/Example_Data/Example.kreport",
#' filetype="Kraken",
#' sci_names=families,
#' rank_method="reads",
#' rank_no=10,
#' outdir=path_out)}
#'
#' @importFrom utils read.table
#' @export
Parse_report <- function(file, filetype, sci_names, rank_method, rank_no = 5, min_reads = 0, outdir){
  df <- read.table(file, sep = '\t', header = FALSE)

  print("report loaded")

  if (filetype == "Kraken"){
    names(df)[1] <- "percent"
    names(df)[2] <- "no_reads_clade"
    names(df)[3] <- "no_reads_taxon"
    names(df)[4] <- "rank_code"
    names(df)[5] <- "NCBI_id"
    names(df)[6] <- "sci_name"

  for (s in 1:length(sci_names)){

    Family_name <- paste0("              ", sci_names[s])

    start_taxa <- which(Family_name == df$sci_name)

    df_trim<- df[-(1:(start_taxa)),]

    Interest_IDs <- y <- data.frame(percent=numeric(),
                                    no_reads_clade=numeric(),
                                    no_reads_taxon=numeric(),
                                    rank_code=character(),
                                    NCBI_id=numeric(),
                                    sci_name=character())

    # Define valid rank codes
    valid_codes <- c("G", "G1", "G2", "S", "S1", "S2")

    for (i in 1:nrow(df_trim)) {
      current_code <- df_trim$rank_code[i]

      if (current_code %in% valid_codes) {
        Interest_IDs[i,] <- df_trim[i,]
      } else {
        break
      }
    }

    Interest_IDs_Species <- subset(Interest_IDs, rank_code == "S")

    print(paste0(nrow(Interest_IDs_Species), " taxa of ", sci_names[s]," found"))

    if (min_reads != 0 ){

      Interest_IDs_Species <- subset(Interest_IDs_Species, no_reads_taxon > min_reads)

      print("min reads applied")

    }



    if (rank_method == "percent" || rank_method == "Percent"){
      Interest_IDs_Species$rank <- rank(-Interest_IDs_Species$percent, ties.method = "first")
      Interest_IDs_Species <- subset(Interest_IDs_Species, rank >= rank_no)

      print("Species ranked and subsetted")
    }

    else if (rank_method == "reads" || rank_method == "Reads"){
      Interest_IDs_Species$rank <- rank(-Interest_IDs_Species$no_reads_taxon, ties.method = "first")
      Interest_IDs_Species <- subset(Interest_IDs_Species, rank <= rank_no)

      print("Species ranked and subsetted")

    }

    print(paste0(nrow(Interest_IDs_Species), " species of ", sci_names[s]," remain"))

    Interest_IDs_Species$sci_name <- gsub("^\\s+", "", Interest_IDs_Species$sci_name)

    writeLines(as.character(Interest_IDs_Species$sci_name), file.path(outdir, paste0(sci_names[s], ".txt")))
  }

  }


}
