#' Generate Primer Binding Report
#'
#' @param plasmid_path Path to the plasmid text file.
#' @param primers_path Path to the primers Excel file.
#' @param output_name Base name for the output HTML and PDF files.
#' @export
generate_primer_report <- function(plasmid_path, primers_path, output_name = "primer_report") {
  # 1. Read inputs
  primer_df <- readxl::read_excel(primers_path, col_names = TRUE)
  primer_df$Sequence <- toupper(gsub("[[:space:]]", "", primer_df$Sequence))

  plasmid_dna_set <- Biostrings::readDNAStringSet(plasmid_path)
  plasmid_seq <- Biostrings::DNAString(toupper(gsub("[[:space:]]", "", as.character(plasmid_dna_set[[1]]))))

  # [Insert all the logic we built for matching, color assignment, and HTML generation here]
  # ...

  # 7. Build HTML report and save
  html_file <- paste0(output_name, ".html")
  pdf_file <- paste0(output_name, ".pdf")

  htmltools::save_html(html_content, file = html_file)
  webshot2::webshot(html_file, file = pdf_file, delay = 2)

  message(sprintf("Success! Reports saved as %s and %s", html_file, pdf_file))
}
