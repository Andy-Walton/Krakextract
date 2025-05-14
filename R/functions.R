#' @title parse_report
#'
#' @description Extracts species of interest from Kraken reports.
#'
#' License: MIT + file LICENSE
#'
#' @param file This is the file path to your Kraken report file
#' @param sci_names This is a vector of the form c("Lactobacillaceae", "Streptococcaceae") ONLY FAMILY NAMES WORK
#' @param rank_method Sort the kraken hits by percentage or number of reads args = "reads" or "percent"
#' @param rank_no The top x number of species (default is 5)
#' @param min_reads The minimum number of reads for a species to be included (default is 0)
#' @param out_type Specifies what you want to output args: IDs or names
#' @param outdir The filepath you want to save your text files to
#'
#' @return Text files with NCBI IDs for each species of interest above thresholds in your kraken output
#'
#' @examples
#' \dontrun{
#' families = c("Lactobacillaceae", "Aerococcaceae")
#' path_out = "man/Example_Data/"
#'
#' parse_report("man/Example_Data/Example.kreport",
#' sci_names=families,
#' rank_method="reads",
#' rank_no=10,
#' out_type="IDs",
#' outdir=path_out)}
#'
#' @importFrom utils read.table
#' @export
parse_report <- function(file, sci_names, rank_method = "reads", rank_no = 5, min_reads = 0, out_type = "names", outdir){
  df <- read.table(file, sep = '\t', header = FALSE)

  print("report loaded")


    names(df)[1] <- "percent"
    names(df)[2] <- "no_reads_clade"
    names(df)[3] <- "no_reads_taxon"
    names(df)[4] <- "rank_code"
    names(df)[5] <- "NCBI_id"
    names(df)[6] <- "sci_name"

  for (s in seq_along(sci_names)){

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

    if (out_type == "names"){

      writeLines(as.character(Interest_IDs_Species$sci_name), file.path(outdir, paste0(sci_names[s], ".txt")))

      print("File with names produced")

    }

    if (out_type == "IDs") {
      writeLines(as.character(Interest_IDs_Species$NCBI_id), file.path(outdir, paste0(sci_names[s], ".txt")))

      print("File with IDs produced")
    }


  }


}

#' @title get_metadata
#'
#' @description downloads refseq metadata - you need to run this once after downloading this package for every kingdom
#' you are interested in (e.g. bacteria, fungi).
#'
#' License: MIT + file LICENSE
#'
#' @param kingdom the kingdoms that the refseq genomes you want to obtain for downstream analysis
#' come from. These can be any of these options:  c("archaea",
#' "fungi",
#' "bacteria",
#' "invertebrate",
#' "plant",
#' "protozoa",
#' "vertebrate_mammalian",
#' "vertebrate_other",
#' "viral")
#' @param outdir The filepath you want to save your metadata files to
#'
#' @return saves refseq metadata to your machine - necessary to download genomes
#'
#' @examples
#' \dontrun{
#' kings <- "fungi"
#' file = "~/Full_test"
#' get_metadata(kings, file)}
#'
#' @importFrom utils read.table
#' @importFrom utils read.delim
#' @importFrom utils download.file
#' @export
get_metadata <- function(kingdom, outdir){

  options(timeout = 100000)

  for (k in 1:length(kingdom)){

    site <- paste("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/", kingdom[k], "/assembly_summary.txt", sep = "")
    out <- paste(outdir,"/assembly_summary_refseq_", kingdom[k], ".txt", sep = "")

    if (file.exists(out)) {print("Summary already downloaded in that location")
    } else {
      download.file(url = site, destfile = out)

      lines <- readLines(out)

      # Find where the header starts (line that starts with "#assembly_accession")
      header_line_index <- which(grepl("^#assembly_accession", lines))

      if (length(header_line_index) == 0) {
        stop("Header line not found in file: ", original_file)
      }

      # Keep from header line to end
      data_lines <- lines[header_line_index:length(lines)]

      # Remove the leading "#" from the header line
      data_lines[1] <- sub("^#", "", data_lines[1])

      # Write to a new file
      writeLines(data_lines, out)

    }
  }
}

#' @title download_refs
#'
#' @description downloads reference genomes using species lists generated by the parse_report() function. It will only download full, complete genomes from the refseq database.
#'
#' License: MIT + file LICENSE
#'
#' @param ID_file the location of the file with list of organisms generated by parse_report()
#' @param assembly_location the location of the refseq metadata files you downloaded using get_metadata()
#' @param outdir The filepath you want to save your reference genomes to
#'
#' @return downloades reference genomes using lists of organism names from parse_report() files
#'
#' @examples
#' \dontrun{
#' path = "~/Test_Data/Lactobacillaceae.txt"
#' location = "~/Full_test"
#' Download_Refs(path, assembly_location = location, outdir = location)}
#'
#' @importFrom utils read.table
#' @importFrom utils read.delim
#' @importFrom utils download.file
#' @export
download_refs <- function(ID_file, assembly_location, outdir){

  options(timeout = 100000)

  df <- read.table(ID_file, sep = '\t', header = FALSE)

  names(df)[1] <- "ID"

  for (i in 1:nrow(df)){

    kingdoms <- c("archaea",
                  "fungi",
                  "bacteria",
                  "invertebrate",
                  "plant",
                  "protozoa",
                  "vertebrate_mammalian",
                  "vertebrate_other",
                  "viral")

    # This checks to see what assembly file to use
    for (k in 1:length(kingdoms)){

      file_name <- paste(assembly_location,"/assembly_summary_refseq_", kingdoms[k], ".txt", sep = "")

      if (file.exists(file_name)) {

        chk <- read.delim(file_name, sep = '\t', header = TRUE, quote = "", comment.char = "#", fill = TRUE)

        status <- sum(grepl(df$ID[i], chk$organism_name, fixed = TRUE))

        if(status>0) {
          break}
      }
    }

    chk <- subset(chk, organism_name == df$ID[i])

    print(paste("Found", nrow(chk), df$ID[i], "genomes", sep = " "))

    chk <- subset(chk, assembly_level == "Complete Genome")
    chk <- subset(chk, genome_rep == "Full")

    print(paste(nrow(chk),"are complete and full", sep = " "))

    chk$seq_rel_date <- gsub("-", "", chk$seq_rel_date)

    chk$seq_rel_date <- as.numeric(chk$seq_rel_date)

    sequence <- chk[which.max(chk$seq_rel_date), ]

    print(paste("Downloading most recent from", sequence$seq_rel_date, sep = " "))

    file_prefix <- gsub(" ", "_", df$ID[i])

    genome_destination = paste0(outdir, "/",  file_prefix, "_genomic_refseq.fna.gz")

    download.file(url = sequence$ftp_path, destfile = genome_destination)
  }

}
