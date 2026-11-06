############################################################################################################
#' Copy de-go-analysis HTML reports into Shiny app www/, add a DE/GO iframe tab, and add an
#' interactive volcano tab (DE TSVs) under Annotations like other ShinyCell tabs (optional).
############################################################################################################

esc_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

is_dego_shiny_yes <- function(x) {
  if (is.null(x) || length(x) != 1L || is.na(x)) {
    return(FALSE)
  }
  if (isTRUE(x)) {
    return(TRUE)
  }
  identical(toupper(trimws(as.character(x))), "YES")
}

resolve_dego_shiny_plots_dir <- function(root_dir, dego_mode) {
  dego_analysis_dir <- file.path(root_dir, "analyses", "de-go-analysis")
  config_util <- file.path(dego_analysis_dir, "util", "function-de-go-config.R")
  if (file.exists(config_util)) {
    sys.source(config_util, envir = environment(), chdir = FALSE)
    dego_mode <- normalize_dego_mode_label(dego_mode)
    plots_dir <- dego_module_plots_dir(dego_analysis_dir, dego_mode)
  } else {
    plots_dir <- file.path(dego_analysis_dir, "plots", paste0("plots_", dego_mode))
  }
  if (!dir.exists(plots_dir)) {
    legacy <- file.path(dego_analysis_dir, "plots")
    if (dir.exists(legacy)) {
      return(list(dir = legacy, mode = dego_mode, legacy = TRUE))
    }
    return(NULL)
  }
  list(dir = plots_dir, mode = dego_mode, legacy = FALSE)
}

resolve_dego_results_de_dir_for_shiny <- function(root_dir, dego_mode) {
  dego_analysis_dir <- file.path(root_dir, "analyses", "de-go-analysis")
  config_util <- file.path(dego_analysis_dir, "util", "function-de-go-config.R")
  if (file.exists(config_util)) {
    sys.source(config_util, envir = environment(), chdir = FALSE)
    mode <- normalize_dego_mode_label(dego_mode)
    return(dego_results_de_dir(dego_analysis_dir, mode))
  }
  mode <- tolower(gsub("-", "_", as.character(dego_mode)))
  file.path(dego_analysis_dir, "results", paste0("results_", mode), "01-de-analysis")
}

is_placeholder_markers_tsv <- function(df) {
  if (is.null(df) || nrow(df) == 0L) {
    return(TRUE)
  }
  if (identical(names(df), "note")) {
    return(TRUE)
  }
  FALSE
}

prepare_volcano_frame_from_de_tsv <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  df <- tryCatch(
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || is_placeholder_markers_tsv(df)) {
    return(NULL)
  }
  req <- c("avg_log2FC", "p_val", "gene")
  if (!all(req %in% names(df))) {
    return(NULL)
  }
  df <- df[!is.na(df$gene) & nzchar(as.character(df$gene)), , drop = FALSE]
  if (nrow(df) == 0L) {
    return(NULL)
  }
  eps <- 1e-300
  df$p_val <- pmax(as.numeric(df$p_val), eps)
  df$y_neg_log10_p <- -log10(df$p_val)
  if ("p_val_adj" %in% names(df)) {
    df$p_val_adj <- suppressWarnings(as.numeric(df$p_val_adj))
    df$p_val_adj[which(is.na(df$p_val_adj))] <- 1
    df$p_val_adj <- pmax(df$p_val_adj, eps)
    df$y_neg_log10_padj <- -log10(df$p_val_adj)
  }
  df
}

human_label_de_marker_file <- function(fname) {
  base <- basename(fname)
  if (grepl("_Markers_within_celltype_", base, fixed = TRUE)) {
    suf <- sub(paste0("^Res_.*?_Markers_within_celltype_(.*)\\.tsv$"), "\\1", base)
    return(paste0("Within cell type (", suf, ")"))
  }
  if (grepl("_Markers_all_", base, fixed = TRUE) && !grepl("_Markers_all\\.tsv$", base)) {
    suf <- sub(paste0("^Res_.*?_Markers_all_(.*)\\.tsv$"), "\\1", base)
    return(paste0("Global interaction (", suf, ")"))
  }
  if (grepl("_Markers_all\\.tsv$", base)) {
    return("All cell types (FindAllMarkers)")
  }
  sub("\\.tsv$", "", base, ignore.case = TRUE)
}

de_marker_tsv_paths_for_shiny <- function(results_dir, resolution_list) {
  if (!dir.exists(results_dir)) {
    return(character(0L))
  }
  esc_res <- gsub("([.|()\\^{}+$*?[]\\\\])", "\\\\\\1", as.character(resolution_list))
  pat <- paste0("^Res_", esc_res, "_(Markers_all.*\\.tsv|Markers_within_celltype_.*\\.tsv)$")
  sort(list.files(results_dir, pattern = pat, full.names = TRUE, ignore.case = FALSE))
}

build_de_go_volcano_index_list <- function(tsv_paths) {
  out <- list()
  for (p in tsv_paths) {
    df <- prepare_volcano_frame_from_de_tsv(p)
    if (is.null(df)) {
      next
    }
    lab <- human_label_de_marker_file(p)
    key <- paste0(lab, " — ", basename(p))
    out[[key]] <- df
    attr(out[[key]], "source_tsv") <- normalizePath(p, winslash = "/", mustWork = FALSE)
  }
  out
}

remove_static_volcano_iframe_tab <- function(ui_lines) {
  if (!any(grepl("src = \"volcano_reports/index.html\"", ui_lines, fixed = TRUE))) {
    return(ui_lines)
  }
  idx <- grep("src = \"volcano_reports/index.html\"", ui_lines, fixed = TRUE)
  if (length(idx) != 1L) {
    return(ui_lines)
  }
  i <- idx[1L]
  start <- i
  while (start > 1L && !grepl("tabPanel\\(", ui_lines[start], perl = TRUE)) {
    start <- start - 1L
  }
  end <- i
  n <- length(ui_lines)
  while (end < n && !grepl("^\\s*\\),\\s*$", ui_lines[end], perl = TRUE)) {
    end <- end + 1L
  }
  if (end > n) {
    return(ui_lines)
  }
  ui_lines[-(start:end)]
}

detect_shiny_dataset_prefix <- function(server_lines) {
  pat <- "^([a-zA-Z0-9]+)meta\\s*=\\s*readRDS"
  hits <- grep(pat, server_lines, perl = TRUE, value = TRUE)
  if (length(hits) == 0L) {
    return("sc1")
  }
  m <- regmatches(hits[length(hits)], regexec(pat, hits[length(hits)], perl = TRUE))
  if (length(m[[1L]]) >= 2L) {
    return(m[[1L]][2L])
  }
  "sc1"
}

inject_volcano_readRDS_after_meta_load <- function(server_path) {
  lines <- readLines(server_path, encoding = "UTF-8", warn = FALSE)
  if (any(grepl("de_go_volcano_index", lines, fixed = TRUE))) {
    return(invisible(TRUE))
  }
  idx <- grep("^([a-zA-Z0-9]+)meta\\s*=\\s*readRDS", lines, perl = TRUE)
  if (length(idx) == 0L) {
    message("Could not find *meta = readRDS in server.R; skipping de_go_volcano_index load")
    return(invisible(FALSE))
  }
  insert_after <- idx[length(idx)]
  block <- c(
    paste0(
      "de_go_volcano_index <- tryCatch(readRDS(\"de_go_volcano_index.rds\"), ",
      "error = function(e) NULL)"
    )
  )
  lines <- append(lines, block, after = insert_after)
  writeLines(lines, server_path, useBytes = TRUE)
  message("server.R: load de_go_volcano_index.rds after Seurat meta load")
  invisible(TRUE)
}

inject_volcano_server_observers <- function(server_path, prefix) {
  lines <- readLines(server_path, encoding = "UTF-8", warn = FALSE)
  if (any(grepl("### Interactive DE volcano", lines, fixed = TRUE))) {
    return(invisible(TRUE))
  }
  close_candidates <- grep("^\\s*\\}\\)\\s*$", lines, perl = TRUE)
  if (length(close_candidates) == 0L) {
    message("Could not find closing )} for shinyServer; skipping volcano server code")
    return(invisible(FALSE))
  }
  close_idx <- close_candidates[length(close_candidates)]
  if (close_idx < 3L) {
    message("Could not find closing )} for shinyServer; skipping volcano server code")
    return(invisible(FALSE))
  }
  v <- prefix
  block <- c(
    "",
    paste0("  ### Interactive DE volcano (sc-rna-seq-snap; prefix ", v, ")"),
    paste0("  observe({"),
    paste0("    if (is.null(de_go_volcano_index) || length(de_go_volcano_index) == 0L) return(invisible(NULL))"),
    paste0("    updateSelectInput(session, \"", v, "v1contrast\", choices = names(de_go_volcano_index), "),
    paste0("                      selected = names(de_go_volcano_index)[[1]])"),
    paste0("  })"),
    paste0("  output$", v, "v1plot <- renderPlot({"),
    paste0("    validate(need(!is.null(de_go_volcano_index) && length(de_go_volcano_index) > 0, "),
    paste0("            \"No DE tables for volcano (run de-go precompute / run-de-computation).\"))"),
    paste0("    req(input$", v, "v1contrast)"),
    paste0("    df0 <- de_go_volcano_index[[input$", v, "v1contrast]]"),
    paste0("    req(!is.null(df0), nrow(df0) > 0)"),
    paste0("    df <- df0"),
    paste0("    if (!is.null(input$", v, "v1ct) && input$", v, "v1ct != \"All\" && \"cluster\" %in% names(df)) {"),
    paste0("      df <- df[as.character(df$cluster) == input$", v, "v1ct, , drop = FALSE]"),
    paste0("    }"),
    paste0("    validate(need(nrow(df) > 0, \"No rows for this selection.\"))"),
    paste0("    use_adj <- identical(input$", v, "v1y, \"FDR (-log10 adj. p)\") && \"y_neg_log10_padj\" %in% names(df)"),
    paste0("    df$yplot <- if (use_adj) df$y_neg_log10_padj else df$y_neg_log10_p"),
    paste0("    fc_cut <- req(input$", v, "v1fc)"),
    paste0("    p_cut <- req(input$", v, "v1p)"),
    paste0("    y_metric <- if (use_adj) as.numeric(df$p_val_adj) else as.numeric(df$p_val)"),
    paste0("    sig <- abs(as.numeric(df$avg_log2FC)) >= fc_cut & y_metric <= p_cut"),
    paste0("    sig[is.na(sig)] <- FALSE"),
    paste0("    ggplot2::ggplot(df, ggplot2::aes(avg_log2FC, yplot)) +"),
    paste0("      ggplot2::geom_point(ggplot2::aes(color = sig), alpha = 0.35, size = ggplot2::rel(1.2), na.rm = TRUE) +"),
    paste0("      ggplot2::scale_color_manual(values = c(`FALSE` = \"grey70\", `TRUE` = \"firebrick\"), guide = \"none\") +"),
    paste0("      ggplot2::labs(x = \"avg_log2FC\", y = if (use_adj) \"-log10 FDR\" else \"-log10 p\", "),
    paste0("                    title = input$", v, "v1contrast) +"),
    paste0("      ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = 14))"),
    paste0("  })"),
    paste0("  output$", v, "v1ct_ui <- renderUI({"),
    paste0("    if (is.null(de_go_volcano_index) || length(de_go_volcano_index) == 0) return(NULL)"),
    paste0("    req(input$", v, "v1contrast)"),
    paste0("    df0 <- de_go_volcano_index[[input$", v, "v1contrast]]"),
    paste0("    if (!is.null(df0) && \"cluster\" %in% names(df0)) {"),
    paste0("      lev <- sort(unique(as.character(df0$cluster[!is.na(df0$cluster)])))"),
    paste0("      if (length(lev) > 0L) {"),
    paste0("        return(selectInput(\"", v, "v1ct\", \"Cell type filter:\", "),
    paste0("                           choices = c(\"All\", lev), selected = \"All\"))"),
    paste0("      }"),
    paste0("    }"),
    paste0("    NULL"),
    paste0("  })"),
    paste0("  output$", v, "v1table <- DT::renderDT({"),
    paste0("    validate(need(!is.null(de_go_volcano_index) && length(de_go_volcano_index) > 0, \"\"))"),
    paste0("    req(input$", v, "v1contrast)"),
    paste0("    df0 <- de_go_volcano_index[[input$", v, "v1contrast]]"),
    paste0("    req(!is.null(df0), nrow(df0) > 0)"),
    paste0("    df <- df0"),
    paste0("    if (!is.null(input$", v, "v1ct) && input$", v, "v1ct != \"All\" && \"cluster\" %in% names(df)) {"),
    paste0("      df <- df[as.character(df$cluster) == input$", v, "v1ct, , drop = FALSE]"),
    paste0("    }"),
    paste0("    use_adj <- identical(input$", v, "v1y, \"FDR (-log10 adj. p)\") && \"y_neg_log10_padj\" %in% names(df)"),
    paste0("    df$yplot <- if (use_adj) df$y_neg_log10_padj else df$y_neg_log10_p"),
    paste0("    br <- input$", v, "v1brush"),
    paste0("    if (!is.null(br)) {"),
    paste0("      keep <- df$avg_log2FC >= br$xmin & df$avg_log2FC <= br$xmax &"),
    paste0("        df$yplot >= br$ymin & df$yplot <= br$ymax"),
    paste0("      df <- df[keep, , drop = FALSE]"),
    paste0("    }"),
    paste0("    show_cols <- intersect(c(\"gene\", \"cluster\", \"avg_log2FC\", \"p_val\", \"p_val_adj\", "),
    paste0("                          \"pct.1\", \"pct.2\"), names(df))"),
    paste0("    DT::datatable(df[, show_cols, drop = FALSE], rownames = FALSE, "),
    paste0("                  options = list(scrollX = TRUE, pageLength = 15))"),
    paste0("  })"),
    ""
  )
  lines <- append(lines, block, after = close_idx - 1L)
  writeLines(lines, server_path, useBytes = TRUE)
  message("server.R: added interactive DE volcano outputs")
  invisible(TRUE)
}

volcano_dynamic_tab_lines <- function(prefix) {
  v <- prefix
  c(
    "  tabPanel(",
    paste0("    HTML(\"Volcano (DE)\"),"),
    paste0("    h4(\"Interactive DE volcano\"),"),
    "    \"Select a DE table, brush the plot to list genes, adjust FC and p-value cutoffs.\",",
    "    br(), br(),",
    "    fluidRow(",
    "      column(",
    "        3, style = \"border-right: 2px solid black\",",
    paste0("        selectInput(\"", v, "v1contrast\", \"DE result:\", choices = NULL),"),
    paste0("        uiOutput(\"", v, "v1ct_ui\"),"),
    paste0("        sliderInput(\"", v, "v1fc\", \"|log2 FC| threshold (highlight):\", min = 0, max = 4, value = 0.25, step = 0.05),"),
    paste0("        sliderInput(\"", v, "v1p\", \"P-value cutoff (highlight):\", min = 1e-20, max = 0.2, value = 0.05),"),
    paste0("        selectInput(\"", v, "v1y\", \"Y-axis:\", "),
    paste0("                    choices = c(\"P-value (-log10)\", \"FDR (-log10 adj. p)\"),"),
    paste0("                    selected = \"P-value (-log10)\"),"),
    paste0("        br(), helpText(\"Data: de-go-analysis precomputed TSVs in results/results_<mode>/01-de-analysis/\")"),
    "      ),",
    "      column(",
    "        9,",
    paste0("        plotOutput(\"", v, "v1plot\", height = \"600px\", "),
    paste0("                   brush = brushOpts(\"", v, "v1brush\", resetOnNew = TRUE)),"),
    paste0("        DTOutput(\"", v, "v1table\")"),
    "      )",
    "    )",
    "  ),"
  )
}

write_dego_reports_index_html <- function(report_files, dest_dir, dego_mode) {
  items <- vapply(
    report_files,
    function(f) {
      label <- sub("\\.html$", "", f, ignore.case = TRUE)
      sprintf(
        "<li><a href=\"%s\" target=\"_blank\" rel=\"noopener\">%s</a></li>",
        esc_html(f),
        esc_html(label)
      )
    },
    character(1L)
  )
  body <- paste0(
    "<!DOCTYPE html><html><head><meta charset=\"utf-8\">",
    "<title>de-go-analysis reports (", esc_html(dego_mode), ")</title>",
    "<style>body{font-family:Arial,sans-serif;margin:1.5rem;} ul{line-height:1.8;}</style>",
    "</head><body>",
    "<h1>de-go-analysis reports</h1>",
    "<p>Mode: <strong>", esc_html(dego_mode), "</strong></p>",
    "<ul>", paste(items, collapse = "\n"), "</ul>",
    "<p><em>Generated by sc-rna-seq-snap rshiny-app module.</em></p>",
    "</body></html>"
  )
  writeLines(body, file.path(dest_dir, "index.html"), useBytes = TRUE)
}

write_volcano_reports_index_html <- function(report_files, dest_dir, dego_mode) {
  items <- vapply(
    report_files,
    function(f) {
      label <- sub("\\.html$", "", f, ignore.case = TRUE)
      sprintf(
        "<li><a href=\"%s\" target=\"_blank\" rel=\"noopener\">%s</a></li>",
        esc_html(f),
        esc_html(label)
      )
    },
    character(1L)
  )
  body <- paste0(
    "<!DOCTYPE html><html><head><meta charset=\"utf-8\">",
    "<title>Volcano plots (", esc_html(dego_mode), ")</title>",
    "<style>body{font-family:Arial,sans-serif;margin:1.5rem;} ul{line-height:1.8;}</style>",
    "</head><body>",
    "<h1>Volcano plots (DE module)</h1>",
    "<p>Mode: <strong>", esc_html(dego_mode), "</strong></p>",
    "<p>Knitted from <code>02-volcano-plots.Rmd</code>. Open a report or use the links below.</p>",
    "<ul>", paste(items, collapse = "\n"), "</ul>",
    "<p><em>Generated by sc-rna-seq-snap rshiny-app module.</em></p>",
    "</body></html>"
  )
  writeLines(body, file.path(dest_dir, "index.html"), useBytes = TRUE)
}

inject_dego_tab_in_ui <- function(ui_path) {
  if (!file.exists(ui_path)) {
    message("ui.R not found; skipping DE/GO tab injection: ", ui_path)
    return(invisible(FALSE))
  }
  ui <- readLines(ui_path, encoding = "UTF-8", warn = FALSE)
  if (any(grepl("DE & GO reports", ui, fixed = TRUE))) {
    return(invisible(TRUE))
  }

  tab_lines <- c(
    "  tabPanel(",
    "    title = \"DE & GO reports\",",
    "    icon = icon(\"file-alt\"),",
    "    tags$iframe(",
    "      src = \"de_go_reports/index.html\",",
    "      width = \"100%\", height = \"900px\",",
    "      style = \"border: none;\"",
    "    )",
    "  ),"
  )

  insert_at <- NA_integer_
  nav_start <- grep("^navbarPage\\(", ui)[1]
  if (!is.na(nav_start)) {
    footer_line <- grep('p\\("Dyer Lab"', ui)[1]
    if (is.na(footer_line)) {
      footer_line <- grep("ShinyCell", ui)[1]
    }
    if (!is.na(footer_line)) {
      menu_close <- grep("^   \\),\\s*$", ui)
      menu_close <- menu_close[menu_close > nav_start & menu_close < footer_line]
      if (length(menu_close) > 0L) {
        insert_at <- menu_close[length(menu_close)]
      }
    }
  }

  if (is.na(insert_at)) {
    shiny_app_line <- grep("^shinyApp\\(", ui)[1]
    if (!is.na(shiny_app_line)) {
      insert_at <- shiny_app_line - 1L
      while (insert_at > 1L && grepl("^\\s*$", ui[insert_at])) {
        insert_at <- insert_at - 1L
      }
      tab_lines[1] <- "    , tabPanel("
      tab_lines[length(tab_lines)] <- "    )"
    }
  }

  if (is.na(insert_at)) {
    message("Could not find navbarPage / shinyApp anchor in ui.R; skipping DE/GO tab")
    return(invisible(FALSE))
  }

  ui <- append(ui, tab_lines, after = insert_at)
  writeLines(ui, ui_path, useBytes = TRUE)
  message("Added \"DE & GO reports\" tab to ui.R")
  invisible(TRUE)
}

#' Insert anchor line after which a companion tab should be added (after DE & GO tab).
find_insert_after_de_go_tab <- function(ui) {
  de_title <- grep("title = \"DE & GO reports\"", ui, fixed = TRUE)[1]
  if (is.na(de_title)) {
    return(NA_integer_)
  }
  n <- length(ui)
  close_lines <- grep("^  \\),\\s*$", ui, perl = TRUE)
  close_lines <- close_lines[close_lines > de_title & close_lines <= n]
  if (length(close_lines) < 1L) {
    return(NA_integer_)
  }
  close_lines[1L]
}

#' ShinyCell (makeShinyCodesMulti) places dataset tabs inside \code{navbarMenu("<header>", ...)}.
#' Insert the Volcano \code{tabPanel} inside the Annotations menu next to Proportion plot and
#' Bubbleplot / Heatmap — before the \code{### Tab1.d1} / Bubbleplot block when possible.
find_insert_volcano_in_annotations_menu <- function(ui) {
  if (any(grepl("Interactive DE volcano", ui, fixed = TRUE)) ||
      any(grepl("Volcano (DE)", ui, fixed = TRUE))) {
    return(NA_integer_)
  }

  d1_comment <- grep("### Tab1.d1", ui, fixed = TRUE)
  bubble_html <- grep("HTML(\"Bubbleplot / Heatmap\")", ui, fixed = TRUE)

  target <- NA_integer_
  if (length(d1_comment) > 0L) {
    target <- min(d1_comment)
  } else if (length(bubble_html) > 0L) {
    b <- min(bubble_html)
    prev <- grep("^\\s*###\\s+Tab", ui, perl = TRUE)
    prev <- prev[prev < b]
    if (length(prev) > 0L) {
      target <- max(prev)
    } else {
      ## tabPanel line immediately above Bubbleplot title
      tp <- grep("^\\s*tabPanel\\(\\s*$", ui, perl = TRUE)
      tp <- tp[tp < b]
      if (length(tp) > 0L) {
        target <- max(tp)
      }
    }
  }

  if (is.na(target)) {
    return(NA_integer_)
  }
  if (target <= 1L) {
    return(NA_integer_)
  }
  target - 1L
}

inject_volcano_tab_in_ui <- function(ui_path, prefix = "sc1") {
  if (!file.exists(ui_path)) {
    message("ui.R not found; skipping Volcano tab injection: ", ui_path)
    return(invisible(FALSE))
  }
  ui <- readLines(ui_path, encoding = "UTF-8", warn = FALSE)
  if (any(grepl("Interactive DE volcano", ui, fixed = TRUE))) {
    return(invisible(TRUE))
  }

  ui <- remove_static_volcano_iframe_tab(ui)
  tab_lines <- volcano_dynamic_tab_lines(prefix)

  ## Prefer Annotations navbarMenu (ShinyCell): same menu as Proportion plot / Heatmap.
  insert_at <- find_insert_volcano_in_annotations_menu(ui)

  if (is.na(insert_at)) {
    insert_at <- find_insert_after_de_go_tab(ui)
  }

  if (is.na(insert_at)) {
    nav_start <- grep("^navbarPage\\(", ui)[1]
    if (!is.na(nav_start)) {
      footer_line <- grep('p\\("Dyer Lab"', ui)[1]
      if (is.na(footer_line)) {
        footer_line <- grep("ShinyCell", ui)[1]
      }
      if (!is.na(footer_line)) {
        menu_close <- grep("^   \\),\\s*$", ui)
        menu_close <- menu_close[menu_close > nav_start & menu_close < footer_line]
        if (length(menu_close) > 0L) {
          insert_at <- menu_close[length(menu_close)]
        }
      }
    }
  }

  if (is.na(insert_at)) {
    shiny_app_line <- grep("^shinyApp\\(", ui)[1]
    if (!is.na(shiny_app_line)) {
      insert_at <- shiny_app_line - 1L
      while (insert_at > 1L && grepl("^\\s*$", ui[insert_at])) {
        insert_at <- insert_at - 1L
      }
      tab_lines[1] <- "    , tabPanel("
      tab_lines[length(tab_lines)] <- "    )"
    }
  }

  if (is.na(insert_at)) {
    message("Could not find ui.R anchor for Volcano tab; skipping")
    return(invisible(FALSE))
  }

  ui <- append(ui, tab_lines, after = insert_at)
  writeLines(ui, ui_path, useBytes = TRUE)
  message(
    "Added interactive \"Volcano (DE)\" tab to ui.R (under Annotations when ShinyCell layout detected)"
  )
  invisible(TRUE)
}

attach_interactive_volcano_to_shiny_app <- function(shiny_dir,
                                                  root_dir,
                                                  dego_mode = "global_interaction",
                                                  resolution_list = NULL) {
  server_path <- file.path(shiny_dir, "server.R")
  ui_path <- file.path(shiny_dir, "ui.R")
  if (!file.exists(server_path) || !file.exists(ui_path)) {
    message("server.R or ui.R missing; skipping interactive volcano tab")
    return(invisible(FALSE))
  }

  prefix <- detect_shiny_dataset_prefix(readLines(server_path, encoding = "UTF-8", warn = FALSE))

  results_dir <- resolve_dego_results_de_dir_for_shiny(root_dir, dego_mode)
  if (!dir.exists(results_dir)) {
    message("DE results dir not found: ", results_dir, "; skipping interactive volcano")
    return(invisible(FALSE))
  }

  res_l <- resolution_list
  if (is.null(res_l) || length(res_l) == 0L || !nzchar(as.character(res_l[[1L]]))) {
    cfg <- file.path(root_dir, "project_parameters.Config.yaml")
    if (file.exists(cfg)) {
      y <- yaml::read_yaml(cfg)
      if (!is.null(y$resolution_list_find_markers)) {
        res_l <- y$resolution_list_find_markers
      }
    }
  }
  if (is.null(res_l) || length(res_l) == 0L || !nzchar(as.character(res_l[[1L]]))) {
    message("resolution_list_find_markers not set; skipping interactive volcano (no TSV glob)")
    return(invisible(FALSE))
  }

  tsv_paths <- de_marker_tsv_paths_for_shiny(results_dir, res_l)
  if (length(tsv_paths) == 0L) {
    message("No Res_*Markers*.tsv in ", results_dir, "; run run-de-computation first")
    return(invisible(FALSE))
  }

  volcano_index <- build_de_go_volcano_index_list(tsv_paths)
  if (length(volcano_index) == 0L) {
    message("DE TSVs found but none parsed for volcano plots; skipping")
    return(invisible(FALSE))
  }

  saveRDS(volcano_index, file.path(shiny_dir, "de_go_volcano_index.rds"), compress = "xz")

  inject_volcano_readRDS_after_meta_load(server_path)
  inject_volcano_server_observers(server_path, prefix)
  inject_volcano_tab_in_ui(ui_path, prefix = prefix)

  manifest <- data.frame(
    tsv = unlist(lapply(volcano_index, function(z) attr(z, "source_tsv")), use.names = FALSE),
    dego_mode = dego_mode,
    resolution_list = as.character(res_l[[1L]]),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(dirname(shiny_dir), "de_go_volcano_shiny_manifest.tsv")
  utils::write.table(
    manifest,
    manifest_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  message(
    "Interactive volcano: ", length(volcano_index), " DE table(s) → de_go_volcano_index.rds; ",
    "Annotations → \"Volcano (DE)\""
  )
  invisible(TRUE)
}

attach_volcano_html_reports_optional <- function(shiny_dir,
                                               root_dir,
                                               dego_mode = "global_interaction") {
  resolved <- resolve_dego_shiny_plots_dir(root_dir, dego_mode)
  if (is.null(resolved)) {
    return(invisible(FALSE))
  }
  plots_dir <- resolved$dir
  dego_mode <- resolved$mode

  all_reports <- sort(list.files(
    plots_dir,
    pattern = "^Report-.*\\.html$",
    ignore.case = TRUE,
    full.names = FALSE
  ))
  volcano_files <- all_reports[grepl("volcano", all_reports, ignore.case = TRUE)]
  if (length(volcano_files) == 0L) {
    message("No Report-*volcano*.html in ", plots_dir, " (optional www copy skipped)")
    return(invisible(FALSE))
  }

  www_dir <- file.path(shiny_dir, "www", "volcano_reports")
  if (!dir.exists(www_dir)) {
    dir.create(www_dir, recursive = TRUE)
  }

  for (f in volcano_files) {
    file.copy(
      file.path(plots_dir, f),
      file.path(www_dir, f),
      overwrite = TRUE
    )
  }

  write_volcano_reports_index_html(volcano_files, www_dir, dego_mode)
  message(
    "Copied ", length(volcano_files), " static volcano HTML report(s) to www/volcano_reports/ (optional)"
  )
  invisible(TRUE)
}

#' @param shiny_dir path to results/shinyApp (contains ui.R)
#' @export
attach_de_go_reports_to_shiny_app <- function(shiny_dir,
                                             root_dir,
                                             dego_mode = "global_interaction",
                                             resolution_list = NULL) {
  resolved <- resolve_dego_shiny_plots_dir(root_dir, dego_mode)
  if (is.null(resolved)) {
    message(
      "de-go plots directory not found for mode=", dego_mode,
      " (skipping DE/GO HTML copy; interactive volcano may still work if DE TSVs exist)."
    )
  }

  plots_ok <- FALSE
  report_files <- character(0L)
  dego_mode_use <- dego_mode

  if (!is.null(resolved)) {
    plots_dir <- resolved$dir
    dego_mode_use <- resolved$mode
    report_files <- sort(list.files(
      plots_dir,
      pattern = "^Report-.*\\.html$",
      ignore.case = TRUE,
      full.names = FALSE
    ))
    plots_ok <- length(report_files) > 0L

    if (plots_ok) {
      www_dir <- file.path(shiny_dir, "www", "de_go_reports")
      if (!dir.exists(www_dir)) {
        dir.create(www_dir, recursive = TRUE)
      }
      for (f in report_files) {
        file.copy(
          file.path(plots_dir, f),
          file.path(www_dir, f),
          overwrite = TRUE
        )
      }
      write_dego_reports_index_html(report_files, www_dir, dego_mode_use)

      manifest <- data.frame(
        report = report_files,
        source_path = file.path(plots_dir, report_files),
        shiny_www_path = file.path("www/de_go_reports", report_files),
        dego_mode = dego_mode_use,
        stringsAsFactors = FALSE
      )
      if (!is.null(resolution_list) && nzchar(as.character(resolution_list))) {
        manifest$resolution_list <- as.character(resolution_list)
      }
      manifest_path <- file.path(dirname(shiny_dir), "de_go_shiny_manifest.tsv")
      utils::write.table(
        manifest,
        manifest_path,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )

      inject_dego_tab_in_ui(file.path(shiny_dir, "ui.R"))
      message(
        "Attached ", length(report_files), " de-go HTML report(s) to Shiny app (mode=",
        dego_mode_use, "). Open tab \"DE & GO reports\" or www/de_go_reports/index.html"
      )
    } else {
      message("No Report-*.html files in ", plots_dir, " (DE/GO iframe tab not added)")
    }
  }

  volcano_ok <- attach_interactive_volcano_to_shiny_app(
    shiny_dir = shiny_dir,
    root_dir = root_dir,
    dego_mode = dego_mode_use,
    resolution_list = resolution_list
  )

  attach_volcano_html_reports_optional(
    shiny_dir = shiny_dir,
    root_dir = root_dir,
    dego_mode = dego_mode_use
  )

  if (!plots_ok && !volcano_ok) {
    message(
      "include_de_go_analysis_rshiny_app_module is YES but neither de-go HTML reports nor ",
      "DE TSVs for interactive volcano were found. Run de-go-analysis / run-de-computation."
    )
    return(invisible(FALSE))
  }

  invisible(TRUE)
}

############################################################################################################
