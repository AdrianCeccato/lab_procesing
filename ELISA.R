
library(tidyverse)
library(nplr)
library(drc)
library(stringr) 
library(openxlsx)
  
  # Cargar el archivo existente
  wb <- loadWorkbook("")
  
  for (sheet in names(wb)) {
    # Leer la hoja como data.frame
    df <- read.xlsx(wb, sheet = sheet, rowNames = TRUE)
    
    # Agregar (2) a cada celda, excepto a "B" y a los que empiezan con "STD"
    df[] <- lapply(df, function(col) {
      sapply(col, function(x) {
        if (is.na(x) || x == "" || x == "B" || grepl("^STD", x)) {
          x
        } else {
          paste0(x, " (4)")
        }
      })
    })
 
    # Escribir de vuelta la hoja
    writeData(wb, sheet = sheet, df, rowNames = TRUE, colNames = TRUE, withFilter = FALSE)
  }
  
  # Guardar el archivo modificado
  saveWorkbook(wb, "", overwrite = TRUE)
  
  df[] <- lapply(df, function(col) {
    gsub("\\.", "-", col)
  })
  
  
  
  

  
  resultados <- process_elisa_folder("")
  resultados[["calibration_plots"]] # Lista de gráficos
  resultados[["results"]] # Lista de gráficos
  
  process_elisa_folder <- function(folder_path) {
    file_list <- list.files(
      path = folder_path,
      pattern = "\\.xlsx$",
      full.names = TRUE
    )
    
    process_single_file <- function(file_path) {
      tryCatch({
        df <- read_excel(file_path, sheet = 2, col_names = FALSE)
        
        # Función para extraer secciones con manejo de errores y nombres únicos
        extract_section <- function(pattern) {
          start_row <- grep(pattern, apply(df, 1, paste, collapse = ""))
          if (length(start_row) == 0) {
            stop("Patrón '", pattern, "' no encontrado")
          }
          start_row <- start_row[1]   # <- para evitar múltiples valores 
          section_data <- df[(start_row + 2):(start_row + 10), 1:13]
          nombres <- as.character(section_data[1, ])
          nombres[is.na(nombres) | nombres == ""] <- 
            paste0("X", seq_along(nombres))[is.na(nombres) | nombres == ""]
          nombres <- make.unique(nombres)
          colnames(section_data) <- nombres
          section_data <- section_data[-1, ]
          section_data %>%
            rename(letra = 1) %>%
            pivot_longer(-letra, names_to = "columna", values_to = "value") %>%
            mutate(celda = paste0(letra, columna)) %>%
            dplyr::select(celda, value)
        }
        
        # Extraer datos
        plate_map <- extract_section("CALCULATION PLATE MAP") %>% 
          rename(study_id = value)
        
        absorbance <- extract_section("CALCULATION OPERATION RESULT") %>% 
          rename(absorbance = value) %>%
          mutate(absorbance = as.numeric(absorbance))
        
        # Combinar datos
        full_data <- plate_map %>%
          left_join(absorbance, by = "celda") %>%
          rename(Well = celda)
        
        # Extraer concentraciones de estándares (bloque robusto)
        fila_inicio <- which(apply(df, 1, function(x) any(grepl("SAMPLE CONCENTRATIONS", x))))
        if (length(fila_inicio) == 0) stop("No se encontró 'SAMPLE CONCENTRATIONS'")
        indice_inicio <- fila_inicio[1] + 2
        concentraciones <- df[indice_inicio:(indice_inicio + 6), 1:2]
        colnames(concentraciones) <- c("study_id", "Concentration")
        concentraciones <- concentraciones %>%
          mutate(
            study_id = as.character(study_id) %>%
              str_replace("Standard", "STD") %>%             # Cambia "Standard" por "STD"
              str_remove("\\s*\\([^\\)]+\\)") %>%            # Elimina lo que está entre paréntesis y el espacio antes
              str_remove(":") %>%                            # Elimina los dos puntos al final
              str_replace_all(" ", ""),                      # Elimina todos los espacios
            Concentration = as.numeric(Concentration)
          )
        
        
        # Procesamiento completo
        full_data <- full_data %>%
          left_join(concentraciones, by = "study_id") %>%
          mutate(
            absorbance = as.numeric(absorbance),
            study_id = str_replace_all(study_id, " ", "")
          )
        
        # Datos de calibración (solo standards)
        cal_data <- full_data %>% filter(!is.na(Concentration))
        if (nrow(cal_data) < 4) stop("No hay suficientes datos de calibración en ", basename(file_path))
        
        # Normalizar standards a proporción
        cal_data <- cal_data %>%
          mutate(prop = convertToProp(absorbance))
        
        # Ajustar modelo 5PL con proporciones
        modelo <- nplr(
          x = cal_data$Concentration,
          y = cal_data$prop,
          npars = 5
        )
        
        # Predecir proporción de muestras y estimar concentración
        muestras <- full_data %>% filter(is.na(Concentration) & !is.na(absorbance))
        if(nrow(muestras) > 0) {
          # Calcular min/max de absorbancia y concentración de los estándares
          min_std <- min(cal_data$absorbance, na.rm = TRUE)
          max_std <- max(cal_data$absorbance, na.rm = TRUE)
          min_conc <- min(cal_data$Concentration, na.rm = TRUE)
          max_conc <- max(cal_data$Concentration, na.rm = TRUE)
          
          # Normalizar muestras
          muestras <- muestras %>%
            mutate(
              prop = (absorbance - min_std) / (max_std - min_std)
            )
          muestras$prop <- pmin(pmax(muestras$prop, 0), 1)
          
          # Calcular concentraciones con el modelo
          preds <- getEstimates(modelo, muestras$prop)$x
          
          # Asignar las concentraciones predichas
          muestras$Concentration <- preds
          
          # Imputar min/max del estándar después del cálculo
          muestras <- muestras %>%
            mutate(
              Concentration = case_when(
                absorbance < min_std ~ min_conc,
                absorbance > max_std ~ max_conc,
                TRUE ~ Concentration
              )
            )
          
          # Actualizar full_data con las concentraciones finales
          full_data$Concentration[which(is.na(full_data$Concentration) & !is.na(full_data$absorbance))] <- muestras$Concentration
        }
        
        
        # Gráfico de calibración
        curva <- ggplot(cal_data, aes(x = Concentration, y = prop)) +
          geom_point(size = 3) +
          geom_smooth(
            method = drm,
            method.args = list(fct = LL.5()),
            se = FALSE,
            color = "red"
          ) +
          scale_x_log10() +
          labs(title = basename(file_path), y = "Absorbancia (proporción)")
        
        # Seleccionar resultados finales
        results <- full_data %>%
          filter(
            !study_id %in% c("B", "STD", "CT", "POOL"),
            !str_detect(study_id, "\\(|\\)")
          ) %>%
          dplyr::select(c(study_id, Concentration))%>%
          group_by(study_id) %>%
          summarise(
            Concentration = mean(Concentration, na.rm = TRUE),
            .groups = "drop"
          )
        
        list(plot = curva, data = results)
        
      }, error = function(e) {
        message(paste("Error en", basename(file_path), ":", e$message))
        return(list(plot = NULL, data = NULL))
      })
    }
    
    resultados <- lapply(file_list, process_single_file)
    
    # Combinar datos excluyendo NULLs
    datos_combinados <- resultados %>%
      map("data") %>%
      compact() %>%
      bind_rows(.id = "Archivo")
    if (nrow(datos_combinados) > 0) {
      datos_combinados <- datos_combinados %>%
        mutate(Archivo = basename(file_list[as.numeric(Archivo)]))
    }
    
    graficos <- map(resultados, "plot")
    
    list(
      calibration_plots = graficos,
      results = datos_combinados
    )
  }
  
  
  
  
  
  
  
  
  
  
  start_row <- grep(pattern, apply(df, 1, paste, collapse = ""))
  if (length(start_row) == 0) {
    stop("Patrón '", pattern, "' no encontrado")
  }
  start_row <- start_row[1]   # <- para evitar múltiples valores 
  
  
  
  results<-resultados[["results"]]
  
  
  
  
  
