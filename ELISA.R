# ============================================================
# PROCESAMIENTO DE RESULTADOS ELISA
# ============================================================

library(tidyverse)
library(nplr)
library(openxlsx)
library(readxl)


# ============================================================
# 1. FUNCIÓN OPCIONAL PARA MODIFICAR EL MAPA DE PLACA
# ============================================================
# Añade "(4)" a las muestras, excepto:
#   - celdas vacías
#   - blancos "B"
#   - estándares que empiezan por "STD"
#
# La función guarda el resultado en un archivo nuevo para no
# sobrescribir accidentalmente el archivo original.
# ============================================================

add_dilution_to_workbook <- function(
    input_file,
    output_file,
    dilution_factor = 4
) {

  if (!file.exists(input_file)) {
    stop("No existe el archivo: ", input_file)
  }

  if (!is.numeric(dilution_factor) ||
      length(dilution_factor) != 1 ||
      is.na(dilution_factor) ||
      dilution_factor <= 0) {
    stop("dilution_factor debe ser un número positivo.")
  }

  wb <- loadWorkbook(input_file)

  for (sheet in names(wb)) {

    plate_data <- read.xlsx(
      input_file,
      sheet = sheet,
      rowNames = TRUE,
      check.names = FALSE
    )

    plate_data[] <- lapply(
      plate_data,
      function(column) {

        vapply(
          column,
          function(value) {

            if (
              is.na(value) ||
              trimws(as.character(value)) == "" ||
              trimws(as.character(value)) == "B" ||
              grepl(
                "^STD",
                trimws(as.character(value)),
                ignore.case = TRUE
              )
            ) {
              return(as.character(value))
            }

            paste0(
              trimws(as.character(value)),
              " (",
              dilution_factor,
              ")"
            )
          },
          FUN.VALUE = character(1)
        )
      }
    )

    writeData(
      wb,
      sheet = sheet,
      x = plate_data,
      rowNames = TRUE,
      colNames = TRUE,
      withFilter = FALSE
    )
  }

  saveWorkbook(
    wb,
    file = output_file,
    overwrite = TRUE
  )

  message("Archivo modificado guardado en: ", output_file)

  invisible(output_file)
}


# ============================================================
# 2. FUNCIÓN PRINCIPAL PARA PROCESAR UNA CARPETA DE ELISA
# ============================================================

process_elisa_folder <- function(
    folder_path,
    sheet_number = 2,
    output_plot_dir = file.path(folder_path, "calibration_plots"),
    min_calibration_points = 4,
    max_cv = 15
) {

  # ----------------------------------------------------------
  # Comprobaciones iniciales
  # ----------------------------------------------------------

  if (!dir.exists(folder_path)) {
    stop("No existe la carpeta: ", folder_path)
  }

  file_list <- list.files(
    path = folder_path,
    pattern = "\\.xlsx$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  # Excluir archivos temporales generados por Excel
  file_list <- file_list[
    !grepl("^~\\$", basename(file_list))
  ]

  if (length(file_list) == 0) {
    stop(
      "No se encontraron archivos .xlsx en: ",
      folder_path
    )
  }

  if (!dir.exists(output_plot_dir)) {
    dir.create(
      output_plot_dir,
      recursive = TRUE
    )
  }


  # ==========================================================
  # FUNCIÓN INTERNA PARA PROCESAR UN ARCHIVO
  # ==========================================================

  process_single_file <- function(file_path) {

    tryCatch({

      message("Procesando: ", basename(file_path))

      raw_data <- read_excel(
        path = file_path,
        sheet = sheet_number,
        col_names = FALSE
      )


      # --------------------------------------------------------
      # Extraer una sección del archivo
      # --------------------------------------------------------

      extract_section <- function(search_pattern) {

        row_text <- apply(
          raw_data,
          1,
          function(row) {
            paste(
              as.character(row),
              collapse = " "
            )
          }
        )

        start_row <- grep(
          search_pattern,
          row_text,
          ignore.case = TRUE
        )

        if (length(start_row) == 0) {
          stop(
            "No se encontró el patrón '",
            search_pattern,
            "'"
          )
        }

        # Si hay varias coincidencias, utilizar la primera
        start_row <- start_row[1]

        first_data_row <- start_row + 2
        last_data_row <- start_row + 10

        if (last_data_row > nrow(raw_data)) {
          stop(
            "La sección '",
            search_pattern,
            "' no contiene suficientes filas."
          )
        }

        if (ncol(raw_data) < 13) {
          stop(
            "La hoja tiene menos de 13 columnas."
          )
        }

        section_data <- raw_data[
          first_data_row:last_data_row,
          1:13
        ]

        column_names <- as.character(
          unlist(section_data[1, ])
        )

        missing_names <- is.na(column_names) |
          trimws(column_names) == ""

        column_names[missing_names] <- paste0(
          "X",
          which(missing_names)
        )

        column_names <- make.unique(column_names)

        names(section_data) <- column_names

        section_data <- section_data[-1, , drop = FALSE]

        section_data %>%
          rename(row_letter = 1) %>%
          mutate(
            row_letter = as.character(row_letter)
          ) %>%
          pivot_longer(
            cols = -row_letter,
            names_to = "column_number",
            values_to = "value"
          ) %>%
          mutate(
            Well = paste0(
              row_letter,
              column_number
            )
          ) %>%
          select(
            Well,
            value
          )
      }


      # --------------------------------------------------------
      # Extraer mapa de placa
      # --------------------------------------------------------

      plate_map <- extract_section(
        "CALCULATION PLATE MAP"
      ) %>%
        rename(
          study_id_original = value
        ) %>%
        mutate(
          study_id_original = as.character(
            study_id_original
          )
        )


      # --------------------------------------------------------
      # Extraer absorbancias
      # --------------------------------------------------------

      absorbance_data <- extract_section(
        "CALCULATION OPERATION RESULT"
      ) %>%
        rename(
          absorbance = value
        ) %>%
        mutate(
          absorbance = suppressWarnings(
            as.numeric(absorbance)
          )
        )


      # --------------------------------------------------------
      # Combinar mapa de placa y absorbancias
      # --------------------------------------------------------

      full_data <- plate_map %>%
        left_join(
          absorbance_data,
          by = "Well"
        ) %>%
        mutate(
          row_id = row_number(),

          study_id_original = str_squish(
            study_id_original
          ),

          # Extrae 4 de "(4)" y 10 de "(10)"
          Dilution_factor = str_extract(
            study_id_original,
            "(?<=\\()[0-9]+(?:[.,][0-9]+)?(?=\\))"
          ),

          Dilution_factor = str_replace(
            Dilution_factor,
            ",",
            "."
          ),

          Dilution_factor = suppressWarnings(
            as.numeric(Dilution_factor)
          ),

          Dilution_factor = replace_na(
            Dilution_factor,
            1
          ),

          # Elimina el indicador de dilución del identificador
          study_id = str_remove(
            study_id_original,
            "\\s*\\([0-9]+(?:[.,][0-9]+)?\\)\\s*$"
          ),

          study_id = str_remove_all(
            study_id,
            "\\s+"
          )
        )


      # --------------------------------------------------------
      # Localizar concentraciones de estándares
      # --------------------------------------------------------

      row_text <- apply(
        raw_data,
        1,
        function(row) {
          paste(
            as.character(row),
            collapse = " "
          )
        }
      )

      standard_start_row <- grep(
        "SAMPLE CONCENTRATIONS",
        row_text,
        ignore.case = TRUE
      )

      if (length(standard_start_row) == 0) {
        stop(
          "No se encontró 'SAMPLE CONCENTRATIONS'."
        )
      }

      standard_start_row <- standard_start_row[1] + 2
      standard_end_row <- standard_start_row + 6

      if (standard_end_row > nrow(raw_data)) {
        stop(
          "La tabla de concentraciones de estándares está incompleta."
        )
      }

      standard_concentrations <- raw_data[
        standard_start_row:standard_end_row,
        1:2
      ]

      names(standard_concentrations) <- c(
        "study_id",
        "Concentration_standard"
      )

      standard_concentrations <-
        standard_concentrations %>%
        mutate(
          study_id = as.character(study_id),

          study_id = str_replace(
            study_id,
            regex(
              "Standard",
              ignore_case = TRUE
            ),
            "STD"
          ),

          study_id = str_remove(
            study_id,
            "\\s*\\([^\\)]+\\)"
          ),

          study_id = str_remove_all(
            study_id,
            ":"
          ),

          study_id = str_remove_all(
            study_id,
            "\\s+"
          ),

          Concentration_standard =
            suppressWarnings(
              as.numeric(
                Concentration_standard
              )
            )
        ) %>%
        filter(
          !is.na(study_id),
          study_id != "",
          !is.na(Concentration_standard)
        ) %>%
        distinct(
          study_id,
          .keep_all = TRUE
        )


      # --------------------------------------------------------
      # Incorporar concentraciones conocidas
      # --------------------------------------------------------

      full_data <- full_data %>%
        left_join(
          standard_concentrations,
          by = "study_id"
        )


      # --------------------------------------------------------
      # Preparar datos de calibración
      # --------------------------------------------------------

      cal_data <- full_data %>%
        filter(
          !is.na(Concentration_standard),
          !is.na(absorbance),
          Concentration_standard > 0,
          is.finite(Concentration_standard),
          is.finite(absorbance)
        )

      if (nrow(cal_data) < min_calibration_points) {
        stop(
          "Solo hay ",
          nrow(cal_data),
          " observaciones válidas de calibración en ",
          basename(file_path),
          ". Se requieren al menos ",
          min_calibration_points,
          "."
        )
      }

      if (
        n_distinct(cal_data$Concentration_standard) <
        min_calibration_points
      ) {
        stop(
          "No hay suficientes concentraciones distintas ",
          "para ajustar la curva en ",
          basename(file_path),
          "."
        )
      }


      # --------------------------------------------------------
      # Normalización consistente
      # --------------------------------------------------------
      # Se utilizan los límites observados en los estándares
      # tanto para estándares como para muestras.
      # No se llama convertToProp() por separado en cada grupo.
      # --------------------------------------------------------

      min_std_abs <- min(
        cal_data$absorbance,
        na.rm = TRUE
      )

      max_std_abs <- max(
        cal_data$absorbance,
        na.rm = TRUE
      )

      absorbance_range <- max_std_abs - min_std_abs

      if (
        !is.finite(absorbance_range) ||
        absorbance_range <= 0
      ) {
        stop(
          "El rango de absorbancia de los estándares es cero ",
          "o no es válido."
        )
      }

      cal_data <- cal_data %>%
        mutate(
          prop = (
            absorbance - min_std_abs
          ) / absorbance_range,

          prop = pmin(
            pmax(prop, 0),
            1
          )
        )


      # --------------------------------------------------------
      # Ajustar modelo logístico de cinco parámetros
      # --------------------------------------------------------

      modelo <- nplr(
        x = cal_data$Concentration_standard,
        y = cal_data$prop,
        npars = 5
      )


      # --------------------------------------------------------
      # Calidad del ajuste
      # --------------------------------------------------------

      fit_quality <- getGoodness(modelo)

      r_squared <- suppressWarnings(
        as.numeric(fit_quality$R2)
      )

      if (
        length(r_squared) == 0 ||
        !is.finite(r_squared)
      ) {
        r_squared <- NA_real_
      }

      message(
        "R² para ",
        basename(file_path),
        ": ",
        ifelse(
          is.na(r_squared),
          "no disponible",
          round(r_squared, 4)
        )
      )


      # --------------------------------------------------------
      # Identificar y cuantificar muestras
      # --------------------------------------------------------

      sample_rows <- full_data %>%
        filter(
          is.na(Concentration_standard),
          !is.na(absorbance),
          is.finite(absorbance)
        ) %>%
        mutate(
          OutOfRange_low =
            absorbance < min_std_abs,

          OutOfRange_high =
            absorbance > max_std_abs,

          OutOfRange =
            OutOfRange_low |
            OutOfRange_high,

          prop = (
            absorbance - min_std_abs
          ) / absorbance_range,

          # Limitar a 0 y 1 para evitar una inversión inválida
          prop_for_prediction = pmin(
            pmax(prop, 0),
            1
          )
        )


      if (nrow(sample_rows) > 0) {

        predicted_values <- getEstimates(
          modelo,
          sample_rows$prop_for_prediction
        )

        predicted_concentration <- suppressWarnings(
          as.numeric(predicted_values$x)
        )

        if (
          length(predicted_concentration) !=
          nrow(sample_rows)
        ) {
          stop(
            "El número de predicciones no coincide con ",
            "el número de muestras."
          )
        }

        sample_rows <- sample_rows %>%
          mutate(
            Concentration_raw =
              predicted_concentration,

            # Concentración corregida por dilución
            Concentration =
              Concentration_raw *
              Dilution_factor
          )

        sample_updates <- sample_rows %>%
          select(
            row_id,
            prop,
            OutOfRange_low,
            OutOfRange_high,
            OutOfRange,
            Concentration_raw,
            Concentration
          )

        full_data <- full_data %>%
          left_join(
            sample_updates,
            by = "row_id"
          )

      } else {

        full_data <- full_data %>%
          mutate(
            prop = NA_real_,
            OutOfRange_low = FALSE,
            OutOfRange_high = FALSE,
            OutOfRange = FALSE,
            Concentration_raw = NA_real_,
            Concentration = NA_real_
          )
      }


      # --------------------------------------------------------
      # Crear curva usando el mismo modelo nplr
      # --------------------------------------------------------

      y_grid <- seq(
        0.001,
        0.999,
        length.out = 300
      )

      curve_estimates <- getEstimates(
        modelo,
        y_grid
      )

      curve_data <- tibble(
        Concentration = suppressWarnings(
          as.numeric(curve_estimates$x)
        ),
        prop = y_grid
      ) %>%
        filter(
          is.finite(Concentration),
          Concentration > 0,
          is.finite(prop)
        ) %>%
        arrange(Concentration)


      calibration_plot <- ggplot() +
        geom_point(
          data = cal_data,
          aes(
            x = Concentration_standard,
            y = prop
          ),
          size = 3,
          color = "black"
        ) +
        geom_line(
          data = curve_data,
          aes(
            x = Concentration,
            y = prop
          ),
          linewidth = 1,
          color = "red"
        ) +
        scale_x_log10() +
        coord_cartesian(
          ylim = c(0, 1)
        ) +
        labs(
          title = basename(file_path),
          subtitle = paste0(
            "Modelo 5PL nplr; R² = ",
            ifelse(
              is.na(r_squared),
              "no disponible",
              round(r_squared, 4)
            )
          ),
          x = "Concentración",
          y = "Absorbancia normalizada"
        ) +
        theme_minimal(base_size = 12)


      # --------------------------------------------------------
      # Guardar gráfico
      # --------------------------------------------------------

      plot_filename <- paste0(
        tools::file_path_sans_ext(
          basename(file_path)
        ),
        "_calibration.png"
      )

      plot_path <- file.path(
        output_plot_dir,
        plot_filename
      )

      ggsave(
        filename = plot_path,
        plot = calibration_plot,
        width = 8,
        height = 6,
        dpi = 300
      )


      # --------------------------------------------------------
      # Resultados por réplica
      # --------------------------------------------------------

      replicate_results <- full_data %>%
        filter(
          !is.na(study_id),
          study_id != "",
          !study_id %in% c(
            "B",
            "BLANK",
            "CT",
            "POOL"
          ),
          !str_detect(
            study_id,
            regex(
              "^STD",
              ignore_case = TRUE
            )
          ),
          !is.na(Concentration)
        ) %>%
        transmute(
          File = basename(file_path),
          Well,
          study_id,
          study_id_original,
          absorbance,
          Dilution_factor,
          Concentration_raw,
          Concentration,
          OutOfRange_low,
          OutOfRange_high,
          OutOfRange
        )


      # --------------------------------------------------------
      # Resumen por muestra
      # --------------------------------------------------------

      sample_results <- replicate_results %>%
        group_by(
          File,
          study_id,
          Dilution_factor
        ) %>%
        summarise(
          n_rep = sum(
            !is.na(Concentration)
          ),

          MeanConc = if (
            all(is.na(Concentration))
          ) {
            NA_real_
          } else {
            mean(
              Concentration,
              na.rm = TRUE
            )
          },

          SDConc = if (
            sum(!is.na(Concentration)) > 1
          ) {
            sd(
              Concentration,
              na.rm = TRUE
            )
          } else {
            NA_real_
          },

          CV = if (
            sum(!is.na(Concentration)) > 1 &&
            is.finite(
              mean(
                Concentration,
                na.rm = TRUE
              )
            ) &&
            mean(
              Concentration,
              na.rm = TRUE
            ) != 0
          ) {
            100 *
              sd(
                Concentration,
                na.rm = TRUE
              ) /
              mean(
                Concentration,
                na.rm = TRUE
              )
          } else {
            NA_real_
          },

          AnyOutOfRange = any(
            OutOfRange,
            na.rm = TRUE
          ),

          .groups = "drop"
        ) %>%
        mutate(
          CV_status = case_when(
            n_rep < 2 ~ "Sin CV: una sola réplica",
            is.na(CV) ~ "CV no calculable",
            CV <= max_cv ~ "Aceptable",
            TRUE ~ paste0(
              "Revisar: CV > ",
              max_cv,
              "%"
            )
          ),

          Range_status = if_else(
            AnyOutOfRange,
            "Fuera del rango de calibración",
            "Dentro del rango"
          ),

          R2 = r_squared
        )


      # --------------------------------------------------------
      # Información del ajuste
      # --------------------------------------------------------

      fit_results <- tibble(
        File = basename(file_path),
        n_calibration_rows = nrow(cal_data),
        n_calibration_levels =
          n_distinct(
            cal_data$Concentration_standard
          ),
        min_standard_concentration =
          min(
            cal_data$Concentration_standard,
            na.rm = TRUE
          ),
        max_standard_concentration =
          max(
            cal_data$Concentration_standard,
            na.rm = TRUE
          ),
        min_standard_absorbance =
          min_std_abs,
        max_standard_absorbance =
          max_std_abs,
        R2 = r_squared,
        Plot_file = plot_path
      )


      list(
        plot = calibration_plot,
        sample_results = sample_results,
        replicate_results = replicate_results,
        fit_results = fit_results,
        model = modelo
      )

    }, error = function(error_condition) {

      message(
        "Error en ",
        basename(file_path),
        ": ",
        conditionMessage(error_condition)
      )

      list(
        plot = NULL,
        sample_results = NULL,
        replicate_results = NULL,
        fit_results = tibble(
          File = basename(file_path),
          Error = conditionMessage(
            error_condition
          )
        ),
        model = NULL
      )
    })
  }


  # ==========================================================
  # PROCESAR TODOS LOS ARCHIVOS
  # ==========================================================

  processed_files <- lapply(
    file_list,
    process_single_file
  )

  names(processed_files) <- basename(file_list)


  # ----------------------------------------------------------
  # Combinar resultados
  # ----------------------------------------------------------

  combined_sample_results <- processed_files %>%
    map("sample_results") %>%
    compact() %>%
    bind_rows()

  combined_replicate_results <- processed_files %>%
    map("replicate_results") %>%
    compact() %>%
    bind_rows()

  combined_fit_results <- processed_files %>%
    map("fit_results") %>%
    compact() %>%
    bind_rows()

  calibration_plots <- processed_files %>%
    map("plot")

  fitted_models <- processed_files %>%
    map("model")


  # ----------------------------------------------------------
  # Devolver todos los objetos
  # ----------------------------------------------------------

  list(
    calibration_plots = calibration_plots,
    results = combined_sample_results,
    replicate_results = combined_replicate_results,
    fit_results = combined_fit_results,
    models = fitted_models
  )
}


# ============================================================
# 3. EJECUCIÓN
# ============================================================
# Sustituye las rutas siguientes por las rutas reales.
# En Windows puedes usar:
# "C:/Users/usuario/Documents/ELISA"
# ============================================================

elisa_folder <- "RUTA/A/LA/CARPETA/ELISA"

resultados <- process_elisa_folder(
  folder_path = elisa_folder,
  sheet_number = 2,
  min_calibration_points = 4,
  max_cv = 15
)


# ============================================================
# 4. CONSULTAR RESULTADOS
# ============================================================

# Resultados resumidos por muestra
results <- resultados$results
print(results)

# Resultados individuales por pocillo/réplica
replicate_results <- resultados$replicate_results
print(replicate_results)

# Calidad del ajuste por archivo
fit_results <- resultados$fit_results
print(fit_results)

# Lista de gráficos
calibration_plots <- resultados$calibration_plots

# Lista de modelos nplr
models <- resultados$models


# ============================================================
# 5. GUARDAR RESULTADOS
# ============================================================

write.csv(
  results,
  file = file.path(
    elisa_folder,
    "ELISA_results_by_sample.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  replicate_results,
  file = file.path(
    elisa_folder,
    "ELISA_results_by_well.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  fit_results,
  file = file.path(
    elisa_folder,
    "ELISA_calibration_quality.csv"
  ),
  row.names = FALSE,
  na = ""
) 
