#' Check Primer Bindings Against Plasmid
#'
#' Internal helper function to calculate matching coordinates.
#' @param plasmid DNAString object of the plasmid.
#' @param primers_df Dataframe containing 'Name' and 'Sequence' columns.
#' @return Dataframe of match results.
#' @import Biostrings
#' @noRd
check_primers <- function(plasmid, primers_df) {
  results <- list()

  for (i in 1:nrow(primers_df)) {
    p_name <- primers_df$Name[i]
    p_seq_str <- primers_df$Sequence[i]

    if (is.na(p_seq_str) || nchar(p_seq_str) == 0) next

    p_dna <- Biostrings::DNAString(p_seq_str)
    fwd_matches <- Biostrings::matchPattern(p_dna, plasmid)
    rc_dna <- Biostrings::reverseComplement(p_dna)
    rc_matches <- Biostrings::matchPattern(rc_dna, plasmid)

    if (length(fwd_matches) > 0) {
      results[[length(results) + 1]] <- data.frame(
        Name = p_name, Sequence = p_seq_str, Strand = "+ (Forward)",
        Match_Count = length(fwd_matches),
        Positions = paste(BiocGenerics::start(fwd_matches), collapse=", "),
        stringsAsFactors = FALSE
      )
    }
    if (length(rc_matches) > 0) {
      results[[length(results) + 1]] <- data.frame(
        Name = p_name, Sequence = p_seq_str, Strand = "- (Reverse Comp)",
        Match_Count = length(rc_matches),
        Positions = paste(BiocGenerics::start(rc_matches), collapse=", "),
        stringsAsFactors = FALSE
      )
    }
    if (length(fwd_matches) == 0 && length(rc_matches) == 0) {
      results[[length(results) + 1]] <- data.frame(
        Name = p_name, Sequence = p_seq_str, Strand = "None",
        Match_Count = 0, Positions = "No Match",
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, results)
}

#' Generate Primer Binding Report
#'
#' Reads a plasmid sequence and a list of primers, calculates binding sites,
#' and exports a color-coded HTML and PDF report.
#'
#' @param plasmid_path Path to the plasmid text file (e.g., "plasmid.txt").
#' @param primers_path Path to the primers Excel file (e.g., "primers.xlsx").
#' @param output_name Base name for the output files (default: "primer_report").
#' @return Saves an HTML and PDF report to the current working directory.
#' @import Biostrings
#' @import gt
#' @import htmltools
#' @importFrom dplyr select %>%
#' @importFrom readxl read_excel
#' @importFrom webshot2 webshot
#' @importFrom scales col_numeric
#' @export
generate_primer_report <- function(plasmid_path, primers_path, output_name = "primer_report") {

  # 1. Load and clean inputs
  message("Reading inputs...")
  primer_df <- readxl::read_excel(primers_path, col_names = TRUE)
  primer_df$Sequence <- toupper(gsub("[[:space:]]", "", primer_df$Sequence))

  plasmid_dna_set <- Biostrings::readDNAStringSet(plasmid_path)
  plasmid_seq <- Biostrings::DNAString(toupper(gsub("[[:space:]]", "", as.character(plasmid_dna_set[[1]]))))
  plasmid_char <- as.character(plasmid_seq)
  plasmid_len <- nchar(plasmid_char)

  # 2. Run the matching algorithm
  message("Matching primers...")
  match_results <- check_primers(plasmid_seq, primer_df)

  # 3. Strictly map colors by UNIQUE PRIMER NAME
  human_palette <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
                     "#b15928", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f")

  matched_names <- unique(match_results$Name[match_results$Match_Count > 0])
  primer_color_map <- character(length(unique(match_results$Name)))
  names(primer_color_map) <- unique(match_results$Name)
  primer_color_map[] <- "#ffffff" # Initialize all with white/neutral

  if(length(matched_names) > 0) {
    palette_assigned <- rep(human_palette, ceiling(length(matched_names) / length(human_palette)))[1:length(matched_names)]
    primer_color_map[matched_names] <- palette_assigned
  }

  match_results$Row_Color <- "transparent"
  for (i in 1:nrow(match_results)) {
    if (match_results$Match_Count[i] > 0) {
      match_results$Row_Color[i] <- primer_color_map[match_results$Name[i]]
    }
  }

  # 4. Create gt table
  message("Generating table...")
  gt_table <- match_results %>%
    dplyr::select(-Row_Color) %>%
    gt::gt() %>%
    gt::tab_header(
      title = "Primer Binding Verification & Sequences",
      subtitle = paste("Vector:", basename(plasmid_path))
    ) %>%
    gt::data_color(
      columns = Name,
      colors = primer_color_map
    ) %>%
    gt::data_color(
      columns = Match_Count,
      colors = scales::col_numeric(
        palette = c("#ffcccc", "#d4edda"),
        domain = c(0, max(match_results$Match_Count, na.rm = TRUE))
      )
    )

  # 5. Build the color-coded & strand-coded plasmid sequence string
  message("Building sequence map...")
  base_colors <- rep("transparent", plasmid_len)
  base_strands <- rep("none", plasmid_len)

  for (i in 1:nrow(match_results)) {
    if (match_results$Match_Count[i] > 0 && match_results$Positions[i] != "No Match") {
      p_color <- match_results$Row_Color[i]
      is_reverse <- grepl("-", match_results$Strand[i])
      p_len <- nchar(match_results$Sequence[i])
      pos_list <- as.numeric(strsplit(match_results$Positions[i], ", ")[[1]])

      for (p in pos_list) {
        end_pos <- min(p + p_len - 1, plasmid_len)
        base_colors[p:end_pos] <- p_color
        base_strands[p:end_pos] <- ifelse(is_reverse, "reverse", "forward")
      }
    }
  }

  # Helper for HTML span generation
  generate_span <- function(chunk, color, strand) {
    if (color == "transparent") {
      return(chunk)
    } else {
      style_str <- sprintf("background-color: %s; color: #fff; padding: 1px 0; font-weight: bold;", color)
      if (strand == "reverse") {
        style_str <- paste0(style_str, " border-bottom: 3px dashed #2c3e50; padding-bottom: 1px;")
      }
      return(sprintf('<span style="%s" title="%s match">%s</span>', style_str, strand, chunk))
    }
  }

  seq_chars <- strsplit(plasmid_char, "")[[1]]
  html_seq_parts <- c()
  current_color <- base_colors[1]
  current_strand <- base_strands[1]
  current_chunk <- seq_chars[1]

  if (plasmid_len > 1) {
    for (i in 2:plasmid_len) {
      if (base_colors[i] == current_color && base_strands[i] == current_strand) {
        current_chunk <- paste0(current_chunk, seq_chars[i])
      } else {
        html_seq_parts <- c(html_seq_parts, generate_span(current_chunk, current_color, current_strand))
        current_color <- base_colors[i]
        current_strand <- base_strands[i]
        current_chunk <- seq_chars[i]
      }
    }
  }
  html_seq_parts <- c(html_seq_parts, generate_span(current_chunk, current_color, current_strand))
  highlighted_sequence_html <- paste0(html_seq_parts, collapse = "")

  # 6. Build HTML content
  html_content <- htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$style(htmltools::HTML("
        body { font-family: Arial, sans-serif; margin: 20px; color: #333; -webkit-print-color-adjust: exact; print-color-adjust: exact; }
        h2, h3 { color: #2c3e50; }
        .sequence-box { font-family: monospace; white-space: pre-wrap; word-break: break-all; background: #fdfdfd; padding: 15px; border: 1px solid #ddd; max-height: 400px; overflow-y: scroll; line-height: 1.6; font-size: 14px; }
        .legend { font-size: 14px; margin-bottom: 10px; padding: 10px; background: #f0f0f0; border-radius: 5px; display: inline-block; }
        .legend-rev { border-bottom: 3px dashed #2c3e50; padding-bottom: 1px; font-weight: bold; }
      "))
    ),
    htmltools::tags$h2("Primer Binding Verification & Sequences"),
    htmltools::as_raw_html(gt_table),
    htmltools::tags$hr(),
    htmltools::tags$h3("Color-Coded Plasmid Sequence Reference"),
    htmltools::tags$div(class = "legend",
                        htmltools::HTML("<strong>Legend:</strong> Highlighted regions denote matched primer binding sites. <span class='legend-rev'>Dashed underline</span> indicates a Reverse Complement match.")),
    htmltools::tags$p(paste("Total Plasmid Length:", plasmid_len, "bp")),
    htmltools::tags$div(class = "sequence-box", htmltools::HTML(highlighted_sequence_html))
  )

  # 7. Save Outputs
  message("Saving reports...")
  out_html <- paste0(output_name, ".html")
  out_pdf <- paste0(output_name, ".pdf")

  htmltools::save_html(html_content, file = out_html)
  webshot2::webshot(out_html, file = out_pdf, delay = 2)

  message(sprintf("Success! Reports saved as:\n  - %s\n  - %s", out_html, out_pdf))
  return(invisible(TRUE))
}
