##############################################################################
# ELISA UNIVERSAL ANALYSIS
##############################################################################

suppressPackageStartupMessages({
  
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(drc)
  library(ggplot2)
  library(openxlsx)
  
})

##############################################################################
# HELPER FUNCTIONS
##############################################################################

extract_section <- function(raw,
                            header,
                            nrows = 9,
                            ncols = 13){
  
  idx <- which(raw[[1]] == header)[1]
  
  if(is.na(idx))
    stop(paste("No encontrado:", header))
  
  start <- idx + 2
  
  block <- raw[
    start:(start+nrows-1),
    1:ncols
  ]
  
  colnames(block) <-
    c(
      "Row",
      as.character(1:(ncol(block)-1))
    )
  
  block
}

##############################################################################
# EXTRACT PLATE MAP
##############################################################################

extract_plate_map <- function(raw){
  
  idx <- which(raw[[1]]=="CALCULATION PLATE MAP")[1]
  
  header_row <- idx + 2
  
  plate_rows <- (header_row+1):(header_row+8)
  
  block <- raw[
    plate_rows,
    1:13
  ]
  
  colnames(block) <- c(
    "Row",
    as.character(1:12)
  )
  
  block |>
    tidyr::pivot_longer(
      cols = -Row,
      names_to = "Column",
      values_to = "SampleID"
    ) |>
    dplyr::mutate(
      Well = paste0(Row,Column),
      SampleID = trimws(as.character(SampleID))
    ) |>
    dplyr::select(
      Well,
      SampleID
    )
  
}
##############################################################################
# EXTRACT BLANK CORRECTED OD
##############################################################################

extract_od <- function(raw){
  
  idx <- which(
    raw[[1]] == "CALCULATION OPERATION RESULTS"
  )[1]
  
  block <- raw[
    (idx+1):(idx+20),
    1:13
  ]
  
  header_row <- which(block[[2]] == "1")[1]
  
  od_block <- block[
    (header_row+1):(header_row+8),
    1:13
  ]
  
  colnames(od_block) <- c(
    "Row",
    as.character(1:12)
  )
  
  od_long <- od_block |>
    tidyr::pivot_longer(
      -Row,
      names_to = "Column",
      values_to = "OD"
    ) |>
    dplyr::mutate(
      Well = paste0(Row, Column),
      OD = as.numeric(OD)
    ) |>
    dplyr::select(
      Well,
      OD
    )
  
  return(od_long)
  
}
##############################################################################
# EXTRACT STANDARD CONCENTRATIONS
##############################################################################

extract_standards <- function(raw){
  
  idx <- grep(
    "^Standard [0-9]+ \\(S[0-9]+\\):$",
    raw[[1]]
  )
  
  idx <- idx[1:7]
  
  labels <- as.character(
    raw[[1]][idx]
  )
  
  sample_ids <- stringr::str_match(
    labels,
    "\\((S[0-9]+)\\)"
  )[,2]
  
  concentrations <- as.numeric(
    raw[[2]][idx]
  )
  
  data.frame(
    SampleID = sample_ids,
    Concentration = concentrations,
    stringsAsFactors = FALSE
  )
  
}


##############################################################################
# INVERSE PREDICTION
##############################################################################

predict_concentration <- function(
    fit,
    od,
    min_x,
    max_x,
    min_od,
    max_od
){
  
  if(is.na(od))
    return(NA)
  
  
  root_fun <- function(x){
    
    predict(
      fit,
      newdata =
        data.frame(
          Concentration = x
        )
    ) - od
    
  }
 
  
  tryCatch({
    
    uniroot(
      root_fun,
      interval = c(
        min_x / 1000,
        max_x * 1000
      )
    )$root
    
  },
  error = function(e){
    
    NA
    
  })
  
}

##############################################################################
# MAIN FUNCTION
##############################################################################

analyze_elisa <- function(
    file,
    sheet = "Calculated results",
    dilution_factor = 1,
    output_prefix = "ELISA",
    output_dir = getwd()
){
  if(!dir.exists(output_dir)){
    dir.create(
      output_dir,
      recursive = TRUE
    )
  }
  excel_file <- file.path(
    output_dir,
    paste0(output_prefix, "_results.xlsx")
  )
  
  plot_file <- file.path(
    output_dir,
    paste0(output_prefix, "_curve.png")
  )
  message("Reading Excel...")
  
  raw <- read_excel(
    file,
    sheet = sheet,
    col_names = FALSE
  )
  
  ##########################################################################
  # IMPORT DATA
  ##########################################################################
  
  plate_map <- extract_plate_map(raw)
  
  od <- extract_od(raw)
  
  standard_conc <- extract_standards(raw)
  
  plate <- left_join(
    plate_map,
    od,
    by = "Well"
  )
  
  ##########################################################################
  # STANDARDS
  ##########################################################################
  
  standards <- plate %>%
    filter(
      grepl(
        "^S[0-9]+$",
        SampleID
      )
    ) %>%
    group_by(SampleID) %>%
    summarise(
      Mean_OD = mean(OD),
      SD = sd(OD),
      CV = ifelse(
        Mean_OD == 0,
        NA,
        100 * SD / Mean_OD
      ),
      n = n(),
      .groups = "drop"
    ) %>%
    left_join(
      standard_conc,
      by = "SampleID"
    ) %>%
    arrange(Concentration)
  
  if(any(is.na(standards$Mean_OD))){
    stop(
      "Hay estándares con OD faltantes."
    )
  }
  ##########################################################################
  # STANDARD QC
  ##########################################################################
  
  standards <- standards %>%
    mutate(
      QC = case_when(
        CV > 20 ~ "FAIL",
        CV > 15 ~ "WARNING",
        TRUE ~ "PASS"
      )
    )
  
  ##########################################################################
  # FIT SETTINGS
  ##########################################################################
  
  exclude_failed_standards <- TRUE
  max_removed_standards <- 2
  min_standards_for_fit <- 5
  
  ##########################################################################
  # DATASET FOR QC FIT
  ##########################################################################
  
  standards_fit <- standards
  
  failed_std <- standards %>%
    filter(
      QC == "FAIL" &
        Concentration != min(Concentration) &
        Concentration != max(Concentration)
    )
  
  if(
    exclude_failed_standards &&
    nrow(failed_std) > 0 &&
    nrow(failed_std) <= max_removed_standards &&
    (nrow(standards) - nrow(failed_std)) >= min_standards_for_fit
  ){
    
    standards_fit <- standards %>%
      filter(
        !(SampleID %in% failed_std$SampleID)
      )
    
    message(
      "Excluded ",
      nrow(failed_std),
      " failed standard(s)"
    )
    
  }
  
  ##########################################################################
  # FIT MODELS
  ##########################################################################
  
  fit4_all <- drm(
    Mean_OD ~ Concentration,
    data = standards,
    fct = LL.4()
  )
  fit5_all <- tryCatch(
    
    drm(
      Mean_OD ~ Concentration,
      data = standards,
      fct = LL.5()
    ),
    
    error = function(e) NULL
    
  )
  
  
  fit4_qc <- drm(
    Mean_OD ~ Concentration,
    data = standards_fit,
    fct = LL.4()
  )
  
  fit5_qc <- tryCatch(
    drm(
    Mean_OD ~ Concentration,
    data = standards_fit,
    fct = LL.5()
  ),
  error = function(e) NULL
  )
  
  ##########################################################################
  # MODEL EVALUATION FUNCTION
  ##########################################################################
  
  evaluate_model <- function(fit, data, name){
    if(is.null(fit)){
      
      return(
        data.frame(
          Model = name,
          AIC = Inf,
          R2 = -Inf,
          MAPE = Inf
        )
      )}
    pred <- predict(
      fit,
      newdata = data.frame(
        Concentration = data$Concentration
      )
    )
    rss <- sum(
      (data$Mean_OD - pred)^2
    )
    
    tss <- sum(
      (
        data$Mean_OD -
          mean(data$Mean_OD)
      )^2
    )
    
    r2 <- 1 - rss/tss
    
    residual_pct <- abs(
      100 *
        (data$Mean_OD - pred) /
        data$Mean_OD
    )
    
    mape <- mean(
      residual_pct,
      na.rm = TRUE
    )
    
    data.frame(
      
      Model = name,
      
      AIC = AIC(fit),
      
      R2 = r2,
      
      MAPE = mape
      
    )
    
  }
  
  ##########################################################################
  # COMPARE MODELS
  ##########################################################################
  message("Evaluating fit4_all")
  evaluate_model(
    fit4_all,
    standards,
    "4PL_ALL"
  )
  
  message("Evaluating fit5_all")
  evaluate_model(
    fit5_all,
    standards,
    "5PL_ALL"
  )
  
  message("Evaluating fit4_qc")
  evaluate_model(
    fit4_qc,
    standards_fit,
    "4PL_QC"
  )
  
  message("Evaluating fit5_qc")
  evaluate_model(
    fit5_qc,
    standards_fit,
    "5PL_QC"
  )
  model_comparison <- bind_rows(
    
    evaluate_model(
      fit4_all,
      standards,
      "4PL_ALL"
    ),
    
    evaluate_model(
      fit5_all,
      standards,
      "5PL_ALL"
    ),
    
    evaluate_model(
      fit4_qc,
      standards_fit,
      "4PL_QC"
    ),
    
    evaluate_model(
      fit5_qc,
      standards_fit,
      "5PL_QC"
    )
    
  )
  
  ##########################################################################
  # SELECT BEST MODEL
  ##########################################################################
  
  best_model_row <-
    
    model_comparison %>%
    
    arrange(
      MAPE,
      AIC
    ) %>%
    
    slice(1)
  
  model_used <- best_model_row$Model
  
  fit <- switch(
    
    model_used,
    
    "4PL_ALL" = fit4_all,
    
    "5PL_ALL" = fit5_all,
    
    "4PL_QC" = fit4_qc,
    
    "5PL_QC" = fit5_qc
    
  )
  
  ##########################################################################
  # FINAL R2
  ##########################################################################
  
  if(grepl("_QC$", model_used)){
    
    fit_data <- standards_fit
    
  } else {
    
    fit_data <- standards
    
  }
  
  pred <- predict(
    fit,
    newdata = data.frame(
      Concentration =
        fit_data$Concentration
    )
  )
  
  RSS <- sum(
    (fit_data$Mean_OD - pred)^2
  )
  
  TSS <- sum(
    (
      fit_data$Mean_OD -
        mean(fit_data$Mean_OD)
    )^2
  )
  
  R2 <- 1 - RSS/TSS
  
  ##########################################################################
  # SUMMARY
  ##########################################################################
  
  cat("\n")
  cat("===================================\n")
  cat("MODEL SELECTED:", model_used, "\n")
  cat("R² =", round(R2,5), "\n")
  cat("===================================\n")
  
  
  ##########################################################################
  # SAMPLES
  ##########################################################################
  
  samples <- plate %>%
    filter(
      !grepl("^S[0-9]+$", SampleID)
    ) %>%
    filter(
      SampleID != "B"
    )
  
  sample_results <- samples %>%
    group_by(SampleID) %>%
    summarise(
      N = n(),
      Mean_OD = mean(OD),
      SD = sd(OD),
      CV = ifelse(
        Mean_OD == 0,
        NA,
        100 * SD / Mean_OD
      ),
      .groups = "drop"
    )
  sample_results <- sample_results %>%
    mutate(
      MissingReplicate =
        ifelse(N < 2,"YES","NO")
    )
  
  
  ##########################################################################
  # QUANTIFICATION
  ##########################################################################
  
  min_x <- min(
    standards$Concentration
  )
  
  max_x <- max(
    standards$Concentration
  )
  
  min_od <- min(
    standards$Mean_OD
  )
  
  max_od <- max(
    standards$Mean_OD
  )
  
  sample_results$Calculated_Conc <-
    sapply(
      sample_results$Mean_OD,
      predict_concentration,
      fit = fit,
      min_x = min_x,
      max_x = max_x,
      min_od = min_od,
      max_od = max_od
    )
  
  sample_results$Final_Conc <-
    sample_results$Calculated_Conc *
    dilution_factor
  
  
  
  ##########################################################################
  # QC FLAGS
  ##########################################################################
  standards <- standards %>%
    mutate(
      QC = case_when(
        CV > 20 ~ "FAIL",
        CV > 15 ~ "WARNING",
        TRUE ~ "PASS"
      )
    )
  sample_results <- sample_results %>%
    mutate(
      Duplicate_QC = case_when(
        CV > 20 ~ "FAIL",
        CV > 15 ~ "WARNING",
        TRUE ~ "PASS"
      )
    )
  
  sample_results <- sample_results %>%
    mutate(
      Range_QC = case_when(
        Mean_OD < min_od ~ "BELOW_LLOQ",
        Mean_OD > max_od ~ "ABOVE_ULOQ",
        TRUE ~ "IN_RANGE"
      )
    )
  sample_results <- sample_results %>%
    mutate(
      Extrapolated =
        Range_QC != "IN_RANGE"
    )
  ##########################################################################
  # CURVE DATA
  ##########################################################################
  
  curve_df <- data.frame(
    
    Concentration = exp(
      
      seq(
        log(min_x),
        log(max_x),
        length.out = 1000
      )
      
    )
    
  )
  
  curve_df$Predicted_OD <-
    predict(
      fit,
      newdata =
        curve_df
    )
  
  curve_qc <- standards %>%
    
    mutate(
      
      Predicted_OD = predict(
        fit,
        newdata = data.frame(
          Concentration = Concentration
        )
        ),
      
      Residual =
        Mean_OD -
        Predicted_OD,
      
      PercentResidual =
        100 *
        Residual /
        Mean_OD
      
    )
  curve_qc <- curve_qc %>%
    mutate(
      Curve_QC = case_when(
        abs(PercentResidual) > 20 ~ "FAIL",
        abs(PercentResidual) > 15 ~ "WARNING",
        TRUE ~ "PASS"
      )
    )
  curve_qc <- curve_qc %>%
    mutate(
      UsedInFit =
        SampleID %in%
        standards_fit$SampleID
    )
  
  ##########################################################################
  # PLOT
  ##########################################################################
  
  curve_plot <-
    
    ggplot() +
    
    geom_line(
      data = curve_df,
      aes(
        Concentration,
        Predicted_OD
      ),
      colour = "blue",
      linewidth = 1
    ) +
    
    geom_point(
      data = standards,
      aes(
        Concentration,
        Mean_OD
      ),
      size = 3
    ) +
    
    scale_x_log10() +
    
    theme_bw(base_size = 14) +
    
    labs(
      title =
        paste0(
          model_used,
          " ELISA Curve (R² = ",
          round(R2,4),
          ")"
        )
    )
  
  ##########################################################################
  # EXPORT
  ##########################################################################
  metadata <- data.frame(
    
    ModelUsed = model_used,
    
    R2 = R2,
    
    DilutionFactor = dilution_factor,
    
    TotalStandards = nrow(standards),
    
    StandardsUsedForFit = nrow(fit_data)
    
  )
  write.xlsx(
    list(
      Standards = standards,
      CurveQC = curve_qc,
      Samples = sample_results,
      ModelComparison = model_comparison,
      Metadata = metadata
    ),
     excel_file,
    overwrite = TRUE
  )
  
  ggsave(
    filename=plot_file, plot= curve_plot,
    width = 7,
    height = 5
  )
  
  ##########################################################################
  # REPORT
  ##########################################################################
  
  cat("\n")
  cat("===================================\n")
  cat("ELISA ANALYSIS\n")
  cat("===================================\n")
  cat(
    "R² = ",
    round(R2,5),
    "\n"
  )
  cat(
    "Samples = ",
    nrow(sample_results),
    "\n"
  )
  cat(
    "Standards = ",
    nrow(standards),
    "\n"
  )
  cat("===================================\n")
  
  ##########################################################################
  # RETURN
  ##########################################################################
      return(
    list(
      fit = fit,
      R2 = R2,
      standards = standards,
      curve_qc = curve_qc,
      samples = sample_results,
      plate = plate,
      plot = curve_plot,
      model = model_used
    )
  ) 
}
### select and put the desired output direction, where the files will be saved.
output_dir <- "C:/Users/vmonforte/Documents/Verónica/ELISA/Resultados"

if(!dir.exists(output_dir)){
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

### introduce the name of the file with its file directory, introduce the sample's dilution factor value and the desired prefix for the files that will be generated and saved.
res <- analyze_elisa(
  file = "C:/Users/vmonforte/Documents/Verónica/ELISA/IgM-OD450_20261015_proc.xlsx",
  dilution_factor = 500,
  output_prefix = "IgM_Run01", output_dir = output_dir
    )
  

