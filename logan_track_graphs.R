library(shiny)
library(dplyr)
library(ggplot2)
library(rlang)
library(mgcv)

# Allow uploads up to 8 GB
options(shiny.maxRequestSize = 8 * 1024^3)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

safe_num <- function(x, default) {
  val <- suppressWarnings(as.numeric(x))
  if (length(val) == 0 || is.na(val)) default else val
}

sample_tracks_per_strata <- function(df, track_var = "track_id", strata_vars = character(), max_n = Inf) {
  if (!track_var %in% names(df) || is.infinite(max_n) || is.null(max_n) || max_n <= 0) {
    return(df)
  }
  
  strata_vars <- strata_vars[strata_vars %in% names(df)]
  distinct_tracks <- df %>% distinct(across(all_of(c(strata_vars, track_var))))
  
  sampled_tracks <- if (length(strata_vars) == 0) {
    if (nrow(distinct_tracks) <= max_n) distinct_tracks else slice_sample(distinct_tracks, n = max_n)
  } else {
    distinct_tracks %>%
      group_by(across(all_of(strata_vars))) %>%
      group_modify(~ if (nrow(.x) <= max_n) .x else slice_sample(.x, n = max_n)) %>%
      ungroup()
  }
  
  df %>%
    inner_join(sampled_tracks, by = c(strata_vars, track_var))
}

ui <- fluidPage(
  titlePanel("LoganTrack Interactive Grapher"),
  fluidRow(
    column(
      width = 12,
      div(
        style = "margin-bottom: 10px;",
        actionButton("draw_plot", "Generate plot", class = "btn-primary"),
        downloadButton("download_pdf", "Export PDF")
      ),
      uiOutput("plot_ui"),
      verbatimTextOutput("data_note")
    )
  ),
  hr(),
  fluidRow(
    column(
      width = 3,
      wellPanel(
        h4("Data"),
        fileInput("data_file", "Upload tracks data (CSV/RDS)", accept = c(".csv", ".rds")),
        helpText("If no file is uploaded, the app tries to use an existing object named 'tracks_filt' in the R environment."),
        h4("Mapping"),
        selectInput("x_var", "X variable", choices = NULL),
        selectInput("y_var", "Y variable", choices = NULL),
        checkboxInput("use_multi_measure", "Show multiple Y measurements side-by-side", value = FALSE),
        selectizeInput("multi_y_vars", "Y variables for side-by-side plots", choices = NULL, multiple = TRUE),
        numericInput("multi_cols", "Number of side-by-side columns", value = 2, min = 1, max = 4, step = 1),
        radioButtons(
          "plot_mode",
          "Plot mode",
          choices = c(
            "Individual tracks (one line per track_id)" = "individual",
            "Grouped average + smoothing" = "summary"
          ),
          selected = "individual"
        ),
        selectInput("track_var", "Track ID column", choices = NULL, selected = "track_id"),
        checkboxInput("use_track_select", "Select specific track IDs", value = FALSE),
        selectizeInput("selected_track_ids", "Track IDs to include", choices = NULL, multiple = TRUE),
        selectInput("group_var", "Group column (for summary/downsampling strata)", choices = NULL),
        selectInput("color_var", "Color by", choices = NULL)
      )
    ),
    column(
      width = 3,
      wellPanel(
        h4("Summary & Uncertainty"),
        selectInput(
          "summary_stat",
          "Summary line statistic (summary mode)",
          choices = c("Mean" = "mean", "Median" = "median"),
          selected = "mean"
        ),
        selectInput(
          "smooth_method",
          "Smoothing method (summary mode)",
          choices = c("None" = "none", "loess" = "loess", "gam" = "gam"),
          selected = "loess"
        ),
        selectInput(
          "uncertainty_type",
          "Uncertainty overlay (summary mode)",
          choices = c("None" = "none", "Error bars (SE)" = "error_bars", "Error bars (95% CI)" = "error_bars_ci95", "CV cloud" = "cv_cloud"),
          selected = "none"
        ),
        sliderInput("uncertainty_alpha", "Uncertainty alpha", min = 0.05, max = 1, value = 0.2, step = 0.05),
        numericInput("error_mult", "Error bar multiplier", value = 1, min = 0.1, step = 0.1),
        sliderInput("smooth_span", "Loess span", min = 0.05, max = 1, value = 0.3, step = 0.05),
        sliderInput("smooth_se", "Smoothing CI alpha", min = 0, max = 0.5, value = 0.2, step = 0.05),
        h4("Filtering"),
        selectizeInput("filter_cols", "Columns to filter", choices = NULL, multiple = TRUE),
        uiOutput("filter_ui"),
        h4("Downsampling"),
        numericInput("max_lines", "Max lines per group/facet", value = 250, min = 1, step = 1)
      )
    ),
    column(
      width = 3,
      wellPanel(
        h4("Facets"),
        checkboxInput("use_facet", "Use faceting", value = FALSE),
        selectInput("facet_rows", "Facet rows", choices = NULL),
        selectInput("facet_cols", "Facet columns", choices = NULL),
        h4("Axes & Titles"),
        fluidRow(
          column(6, numericInput("x_min", "X min", value = NA)),
          column(6, numericInput("x_max", "X max", value = NA))
        ),
        fluidRow(
          column(6, numericInput("y_min", "Y min", value = NA)),
          column(6, numericInput("y_max", "Y max", value = NA))
        ),
        textInput("plot_title_main", "Main title", value = ""),
        textInput("plot_title_x", "X-axis title", value = ""),
        textInput("plot_title_y", "Y-axis title", value = ""),
        h4("Plot display"),
        sliderInput("plot_width", "Plot width (px)", min = 100, max = 2000, value = 1100, step = 10),
        sliderInput("plot_height", "Plot height (px)", min = 100, max = 2000, value = 780, step = 10),
        sliderInput("plot_scale", "Plot scale factor", min = 0.5, max = 3, value = 1, step = 0.1),
        numericInput("plot_res", "Render resolution (DPI)", value = 96, min = 72, max = 600, step = 12)
      )
    ),
    column(
      width = 3,
      wellPanel(
        h4("Theme controls"),
        selectInput("base_theme", "Base theme", choices = c("classic", "bw", "minimal", "light", "gray"), selected = "classic"),
        sliderInput("text_size", "Text size", min = 6, max = 48, value = 12, step = 1),
        sliderInput("axis_text_size", "Axis text size", min = 6, max = 48, value = 12, step = 1),
        sliderInput("legend_title_size", "Legend title size", min = 6, max = 48, value = 12, step = 1),
        checkboxInput("show_color_legend", "Show color legend", value = TRUE),
        sliderInput("line_width", "Line width", min = 0.05, max = 3, value = 0.2, step = 0.05),
        sliderInput("line_alpha", "Line alpha", min = 0.05, max = 1, value = 0.3, step = 0.05),
        sliderInput("aspect_ratio", "Aspect ratio", min = 0.2, max = 2, value = 1, step = 0.1),
        sliderInput("title_hjust", "Title hjust", min = 0, max = 1, value = 0.5, step = 0.05),
        checkboxInput("show_axis_line", "Show axis line", value = FALSE),
        checkboxInput("show_strip_bg", "Show strip background", value = FALSE),
        textInput("panel_fill", "Panel background fill", value = "white"),
        textInput("panel_border", "Panel border color", value = "black"),
        numericInput("panel_border_size", "Panel border linewidth", value = 1, min = 0, step = 0.1),
        checkboxInput("use_custom_palette", "Use custom line colors", value = FALSE),
        textInput("custom_colors", "Comma-separated colors", value = "deepskyblue,deepskyblue4,green1,green4")
      )
    )
  )
)

server <- function(input, output, session) {
  raw_data <- reactive({
    if (!is.null(input$data_file)) {
      ext <- tools::file_ext(input$data_file$name)
      if (tolower(ext) == "csv") {
        return(read.csv(input$data_file$datapath, check.names = FALSE, stringsAsFactors = FALSE))
      }
      if (tolower(ext) == "rds") {
        return(readRDS(input$data_file$datapath))
      }
      return(NULL)
    }
    
    default_df <- get0("tracks_filt", envir = .GlobalEnv, inherits = TRUE, ifnotfound = NULL)
    if (!is.null(default_df)) {
      return(default_df)
    }
    
    NULL
  })
  
  observeEvent(raw_data(), {
    df <- raw_data()
    req(df)
    cols <- names(df)
    factor_like <- cols[vapply(df, function(x) is.factor(x) || is.character(x), logical(1))]
    
    updateSelectInput(session, "x_var", choices = cols, selected = if ("h_rel_fm" %in% cols) "h_rel_fm" else cols[1])
    updateSelectInput(session, "y_var", choices = cols, selected = if ("mean_foci" %in% cols) "mean_foci" else cols[min(2, length(cols))])
    updateSelectizeInput(session, "multi_y_vars", choices = cols, selected = if ("mean_foci" %in% cols) "mean_foci" else cols[min(2, length(cols))], server = TRUE)
    updateSelectInput(session, "track_var", choices = cols, selected = if ("track_id" %in% cols) "track_id" else cols[1])
    updateSelectInput(session, "group_var", choices = c("None" = "", cols), selected = if ("treat" %in% cols) "treat" else "")
    updateSelectInput(session, "color_var", choices = c("None" = "", cols), selected = if ("treat" %in% cols) "treat" else "")
    updateSelectizeInput(session, "filter_cols", choices = cols, selected = NULL, server = TRUE)
    updateSelectInput(session, "facet_rows", choices = c("None" = "", cols), selected = "")
    updateSelectInput(session, "facet_cols", choices = c("None" = "", cols), selected = "")
    
    if (length(factor_like) == 0) {
      updateSelectInput(session, "color_var", selected = "")
    }
  }, ignoreInit = FALSE)
  
  observeEvent(list(filtered_data(), input$track_var), {
    df <- filtered_data()
    req(df, input$track_var)
    req(input$track_var %in% names(df))
    
    track_choices <- sort(unique(as.character(df[[input$track_var]])))
    
    current <- isolate(input$selected_track_ids %||% character())
    current <- current[current %in% track_choices]
    
    updateSelectizeInput(
      session,
      "selected_track_ids",
      choices = track_choices,
      selected = current,
      server = TRUE
    )
  }, ignoreInit = FALSE)
  
  output$filter_ui <- renderUI({
    df <- raw_data()
    req(df)
    
    if (length(input$filter_cols %||% character()) == 0) {
      return(helpText("Select one or more columns to add interactive filters."))
    }
    
    controls <- lapply(input$filter_cols, function(col) {
      v <- df[[col]]
      id <- paste0("flt_", col)
      
      if (is.numeric(v)) {
        rng <- range(v, na.rm = TRUE)
        sliderInput(id, paste("Range:", col), min = floor(rng[1]), max = ceiling(rng[2]), value = rng)
      } else {
        choices <- sort(unique(as.character(v)))
        selectizeInput(id, paste("Values:", col), choices = choices, selected = choices, multiple = TRUE)
      }
    })
    
    do.call(tagList, controls)
  })
  
  filtered_data <- reactive({
    df <- raw_data()
    req(df)
    
    for (col in input$filter_cols %||% character()) {
      id <- paste0("flt_", col)
      filter_value <- input[[id]]
      
      if (is.null(filter_value)) next
      
      if (is.numeric(df[[col]])) {
        df <- df %>% filter(.data[[col]] >= filter_value[1], .data[[col]] <= filter_value[2])
      } else {
        df <- df %>% filter(as.character(.data[[col]]) %in% as.character(filter_value))
      }
    }
    
    df
  })
  
  plot_data <- reactive({
    df <- filtered_data()
    req(df)
    req(input$x_var, input$track_var)
    
    group_var <- input$group_var
    color_var <- input$color_var
    facet_vars <- c(input$facet_rows, input$facet_cols)
    strata <- unique(c(group_var[group_var != ""], color_var[color_var != ""], facet_vars[facet_vars != ""]))
    
    if (isTRUE(input$use_track_select)) {
      chosen_ids <- input$selected_track_ids %||% character()
      if (length(chosen_ids) > 0) {
        df <- df %>% filter(as.character(.data[[input$track_var]]) %in% as.character(chosen_ids))
      }
    }
    
    sample_tracks_per_strata(
      df = df,
      track_var = input$track_var,
      strata_vars = strata,
      max_n = input$max_lines
    )
  })
  
  plot_state <- reactive({
    req(input$draw_plot > 0)
    
    df <- plot_data()
    choices <- names(df)
    
    y_vars <- if (isTRUE(input$use_multi_measure)) {
      vars <- input$multi_y_vars %||% character()
      vars <- vars[vars %in% choices]
      if (length(vars) == 0 && input$y_var %in% choices) vars <- input$y_var
      vars
    } else {
      if (input$y_var %in% choices) input$y_var else character()
    }
    
    list(df = df, y_vars = y_vars)
  })
  
  build_plot <- function(df, y_var) {
    req(nrow(df) > 1)
    
    x_var <- input$x_var
    track_var <- input$track_var
    group_var <- input$group_var
    color_var <- input$color_var
    
    p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]]))
    
    if (input$plot_mode == "individual") {
      if (color_var != "") {
        p <- p + geom_line(
          aes(group = .data[[track_var]], color = .data[[color_var]]),
          linewidth = safe_num(input$line_width, 0.2),
          alpha = safe_num(input$line_alpha, 0.3)
        )
      } else {
        p <- p + geom_line(
          aes(group = .data[[track_var]]),
          linewidth = safe_num(input$line_width, 0.2),
          alpha = safe_num(input$line_alpha, 0.3)
        )
      }
    } else {
      facet_active <- character()
      if (isTRUE(input$use_facet)) {
        facet_active <- c(input$facet_rows, input$facet_cols)
        facet_active <- facet_active[facet_active != "" & facet_active %in% names(df)]
      }
      
      summary_fun <- if (identical(input$summary_stat, "median")) stats::median else mean
      base_group_cols <- unique(c(group_var[group_var != ""], color_var[color_var != ""], facet_active))
      grouping_cols <- unique(c(base_group_cols, x_var))
      
      sum_df <- df %>%
        group_by(across(all_of(grouping_cols))) %>%
        summarise(
          .y = summary_fun(.data[[y_var]], na.rm = TRUE),
          .sd = stats::sd(.data[[y_var]], na.rm = TRUE),
          .n = sum(!is.na(.data[[y_var]])),
          .groups = "drop"
        ) %>%
        mutate(
          .se = ifelse(.n > 0, .sd / sqrt(.n), NA_real_),
          .t_crit = ifelse(.n > 1, stats::qt(0.975, df = .n - 1), NA_real_),
          .ci95 = .se * .t_crit,
          .cv = ifelse(is.na(.y) | .y == 0, NA_real_, .sd / abs(.y)),
          .cv = ifelse(is.finite(.cv), .cv, NA_real_)
        )
      
      if (length(base_group_cols) == 0) {
        p <- ggplot(sum_df, aes(x = .data[[x_var]], y = .y))
      } else {
        p <- ggplot(
          sum_df,
          aes(
            x = .data[[x_var]],
            y = .y,
            group = interaction(!!!syms(base_group_cols), drop = TRUE)
          )
        )
      }
      
      if (color_var != "") {
        p <- p + aes(color = .data[[color_var]])
      }
      
      if (identical(input$uncertainty_type, "error_bars")) {
        p <- p + geom_errorbar(
          aes(
            ymin = .y - (.se * safe_num(input$error_mult, 1)),
            ymax = .y + (.se * safe_num(input$error_mult, 1))
          ),
          alpha = safe_num(input$uncertainty_alpha, 0.2),
          width = 0
        )
      } else if (identical(input$uncertainty_type, "error_bars_ci95")) {
        p <- p + geom_errorbar(
          aes(
            ymin = .y - .ci95,
            ymax = .y + .ci95
          ),
          alpha = safe_num(input$uncertainty_alpha, 0.2),
          width = 0
        )
      } else if (identical(input$uncertainty_type, "cv_cloud")) {
        if (color_var != "") {
          p <- p + geom_ribbon(
            aes(
              ymin = .y * (1 - .cv),
              ymax = .y * (1 + .cv),
              fill = .data[[color_var]]
            ),
            alpha = safe_num(input$uncertainty_alpha, 0.2),
            colour = NA
          )
        } else {
          p <- p + geom_ribbon(
            aes(ymin = .y * (1 - .cv), ymax = .y * (1 + .cv)),
            alpha = safe_num(input$uncertainty_alpha, 0.2),
            fill = "grey70",
            colour = NA
          )
        }
      }
      
      if (input$smooth_method == "none") {
        p <- p + geom_line(linewidth = safe_num(input$line_width, 0.2), alpha = safe_num(input$line_alpha, 0.3))
      } else if (input$smooth_method == "loess") {
        p <- p + geom_smooth(
          method = "loess",
          span = safe_num(input$smooth_span, 0.3),
          se = TRUE,
          level = 0.95,
          alpha = safe_num(input$smooth_se, 0.2),
          linewidth = safe_num(input$line_width, 0.2)
        )
      } else {
        p <- p + geom_smooth(
          method = "gam",
          formula = y ~ s(x, bs = "cs"),
          se = TRUE,
          level = 0.95,
          alpha = safe_num(input$smooth_se, 0.2),
          linewidth = safe_num(input$line_width, 0.2)
        )
      }
    }
    
    if (isTRUE(input$use_custom_palette) && color_var != "") {
      colors <- trimws(strsplit(input$custom_colors, ",")[[1]])
      if (length(colors) > 0) {
        p <- p + scale_color_manual(values = colors)
        if (identical(input$plot_mode, "summary") && identical(input$uncertainty_type, "cv_cloud")) {
          p <- p + scale_fill_manual(values = colors)
        }
      }
    }
    
    if (isTRUE(input$use_facet)) {
      row_var <- if (input$facet_rows == "") "." else input$facet_rows
      col_var <- if (input$facet_cols == "") "." else input$facet_cols
      facet_formula <- as.formula(paste(row_var, "~", col_var))
      p <- p + facet_grid(facet_formula)
    }
    
    x_limits <- c(suppressWarnings(as.numeric(input$x_min)), suppressWarnings(as.numeric(input$x_max)))
    y_limits <- c(suppressWarnings(as.numeric(input$y_min)), suppressWarnings(as.numeric(input$y_max)))
    
    x_title <- if (nzchar(input$plot_title_x %||% "")) input$plot_title_x else x_var
    y_title <- if (nzchar(input$plot_title_y %||% "")) input$plot_title_y else y_var
    main_title <- if (nzchar(input$plot_title_main %||% "")) input$plot_title_main else y_var
    
    p <- p + labs(
      x = x_title,
      y = y_title,
      title = main_title,
      color = ifelse(color_var == "", "", color_var)
    )
    
    x_lim_use <- if (!anyNA(x_limits)) x_limits else NULL
    y_lim_use <- if (!anyNA(y_limits)) y_limits else NULL
    if (!is.null(x_lim_use) || !is.null(y_lim_use)) {
      p <- p + coord_cartesian(xlim = x_lim_use, ylim = y_lim_use)
    }
    
    base_theme <- switch(
      input$base_theme,
      classic = theme_classic(),
      bw = theme_bw(),
      minimal = theme_minimal(),
      light = theme_light(),
      gray = theme_gray(),
      theme_classic()
    )
    
    p +
      base_theme +
      theme(
        text = element_text(size = safe_num(input$text_size, 12)),
        axis.text = element_text(size = safe_num(input$axis_text_size, 12)),
        legend.title = element_text(size = safe_num(input$legend_title_size, 12)),
        legend.position = if (isTRUE(input$show_color_legend)) "right" else "none",
        aspect.ratio = safe_num(input$aspect_ratio, 1),
        plot.title = element_text(hjust = safe_num(input$title_hjust, 0.5)),
        axis.line = if (isTRUE(input$show_axis_line)) element_line() else element_blank(),
        strip.background = if (isTRUE(input$show_strip_bg)) element_rect(fill = "grey90", color = "grey40") else element_blank(),
        panel.background = element_rect(fill = input$panel_fill, colour = input$panel_border, linewidth = safe_num(input$panel_border_size, 1)),
        legend.key = element_rect(color = NA)
      )
  }
  
  plot_width_fn <- function() {
    round(
      safe_num(input$plot_width, 1100) *
        safe_num(input$plot_scale, 1) *
        (safe_num(input$plot_res, 96) / 96)
    )
  }
  
  plot_height_fn <- function() {
    round(
      safe_num(input$plot_height, 780) *
        safe_num(input$plot_scale, 1) *
        (safe_num(input$plot_res, 96) / 96)
    )
  }
  
  panel_plot_width_fn <- function() {
    n_cols <- if (isTRUE(input$use_multi_measure)) max(1, min(4, round(safe_num(input$multi_cols, 2)))) else 1
    max(240, floor(plot_width_fn() / n_cols))
  }
  
  output$plot_ui <- renderUI({
    req(input$draw_plot > 0)
    st <- plot_state()
    y_vars <- st$y_vars
    req(length(y_vars) > 0)
    
    display_width_px <- round(safe_num(input$plot_width, 1100) * safe_num(input$plot_scale, 1))
    display_height_px <- round(safe_num(input$plot_height, 780) * safe_num(input$plot_scale, 1))
    
    if (!isTRUE(input$use_multi_measure) || length(y_vars) == 1) {
      return(plotOutput("line_plot", width = paste0(display_width_px, "px"), height = paste0(display_height_px, "px")))
    }
    
    n_cols <- max(1, min(4, round(safe_num(input$multi_cols, 2))))
    panel_width_px <- max(240, floor(display_width_px / n_cols))
    
    plot_outputs <- lapply(seq_along(y_vars), function(i) {
      tags$div(
        style = paste0("flex: 0 0 ", panel_width_px, "px;"),
        plotOutput(
          outputId = paste0("line_plot_", i),
          width = paste0(panel_width_px, "px"),
          height = paste0(display_height_px, "px")
        )
      )
    })
    
    tags$div(
      style = "display:flex; flex-wrap:wrap; gap:12px; align-items:flex-start;",
      do.call(tagList, plot_outputs)
    )
  })
  
  output$line_plot <- renderPlot({
    req(input$draw_plot > 0)
    st <- plot_state()
    y_vars <- st$y_vars
    req(!isTRUE(input$use_multi_measure) || length(y_vars) == 1)
    req(length(y_vars) >= 1)
    build_plot(st$df, y_vars[[1]])
  }, width = plot_width_fn, height = plot_height_fn)
  
  observe({
    req(input$draw_plot > 0)
    st <- plot_state()
    req(isTRUE(input$use_multi_measure))
    y_vars <- st$y_vars
    req(length(y_vars) > 0)
    
    for (i in seq_along(y_vars)) {
      local({
        idx <- i
        y_var_i <- y_vars[[idx]]
        output[[paste0("line_plot_", idx)]] <- renderPlot({
          build_plot(st$df, y_var_i)
        }, width = panel_plot_width_fn, height = plot_height_fn)
      })
    }
  })
  
  output$download_pdf <- downloadHandler(
    filename = function() {
      paste0("logantrack_plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(input$draw_plot > 0)
      st <- plot_state()
      req(length(st$y_vars) > 0)
      
      grDevices::pdf(
        file = file,
        width = safe_num(input$plot_width, 1100) / 96,
        height = safe_num(input$plot_height, 780) / 96
      )
      on.exit(grDevices::dev.off(), add = TRUE)
      
      for (y_var_i in st$y_vars) {
        print(build_plot(st$df, y_var_i))
      }
    }
  )
  
  output$data_note <- renderText({
    df <- raw_data()
    if (is.null(df)) {
      return("No data available. Upload a CSV/RDS or preload an object named 'tracks_filt'.")
    }
    if (input$draw_plot <= 0) {
      return(paste0("Rows: ", nrow(df), " | Columns: ", ncol(df), " | Click 'Generate plot' to render."))
    }
    st <- plot_state()
    y_info <- paste(st$y_vars, collapse = ", ")
    paste0("Rows: ", nrow(df), " | Columns: ", ncol(df), " | Y displayed: ", y_info)
  })
}

shinyApp(ui, server)
