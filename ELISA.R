

library(tidyverse)
library(readxl)
library(nplr)
library(drc)
library(stringr)


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

setwd("C:\\Users\\aaceccato\\Documents")
resultados <- process_elisa_folder("C:\\Users\\aaceccato\\Documents\\STREM-1")
resultados[["calibration_plots"]] # Lista de gráficos

results<-resultados[["results"]]

library(dplyr)
library(stringr)

resultados_limpios <- results %>%
  mutate(
    # Extraer solo los dígitos al inicio antes del punto (puede ser 1, 2 o 3 dígitos)
    study_id_num = str_extract(study_id, "^\\d+"),
    # Asignar timepoint según sufijo
    redcap_event_name.factor = case_when(
      grepl("\\.1$", study_id) ~ "0h",
      grepl("\\.2$", study_id) ~ "4h",
      grepl("\\.3$", study_id) | grepl("=3$", study_id) ~ "24h",
      TRUE ~ NA_character_
    )
  ) %>%
  # Filtrar solo IDs numéricos (descarta los que tengan letras o NA)
  filter(!is.na(study_id_num)) %>%
  # Renombrar y seleccionar columnas
  rename(strem1 = Concentration) %>%
  dplyr::select(study_id = study_id_num, redcap_event_name.factor, strem1)
write.csv(resultados_limpios,"strem_septip.csv")




library(readxl)
library(dplyr)
library(nplr)


# Especifica el nombre del archivo
df <- read_excel("170125-STREM1_20250117_154930.xlsx",sheet = 2) 

fila_inicio <- grep("CALCULATION PLATE MAP", apply(df, 1, paste, collapse = ""))

# Calcula el índice de la fila donde comienza la matriz (dos filas después)
indice_matriz <- fila_inicio + 2
indice_matriz_f<- fila_inicio + 10 
# Extrae la matriz desde la fila calculada hasta el final del archivo
matriz <- df[indice_matriz:indice_matriz_f,1:13 ]

# Elimina filas completamente vacías (si es necesario)
matriz <- matriz[rowSums(is.na(matriz)) != ncol(matriz), ]
nombres_columnas <- as.character(matriz[1, ])

# Asigna los nombres a las columnas
colnames(matriz) <- nombres_columnas

# Renombra la primera columna para facilitar el manejo
matriz <- matriz %>% rename(letra = 1)

# La primera fila contiene los números de columna, así que la extraemos
matriz_datos <- matriz[-1, ] 
# Pivotamos a formato largo
df_largo_plate_map <- matriz_datos %>%
  pivot_longer(
    cols = -letra,
    names_to = "columna",
    values_to = "study_id"
  )

df_largo_plate_map <- df_largo_plate_map %>%
  mutate(celda = paste0(letra, columna)) %>%
  dplyr::select(celda, study_id)


#Repetir para absorbancia:

fila_inicio <- grep("CALCULATION OPERATION RESULT", apply(df, 1, paste, collapse = ""))

# Calcula el índice de la fila donde comienza la matriz (dos filas después)
indice_matriz <- fila_inicio + 2
indice_matriz_f<- fila_inicio + 10 
# Extrae la matriz desde la fila calculada hasta el final del archivo
matriz <- df[indice_matriz:indice_matriz_f,1:13 ]

# Elimina filas completamente vacías (si es necesario)
matriz <- matriz[rowSums(is.na(matriz)) != ncol(matriz), ]
nombres_columnas <- as.character(matriz[1, ])

# Asigna los nombres a las columnas
colnames(matriz) <- nombres_columnas

# Renombra la primera columna para facilitar el manejo
matriz <- matriz %>% rename(letra = 1)

# La primera fila contiene los números de columna, así que la extraemos
matriz_datos <- matriz[-1, ]                # elimina la primera fila

# Pivotamos a formato largo
df_largo_abs <- matriz_datos %>%
  pivot_longer(
    cols = -letra,
    names_to = "columna",
    values_to = "absorbance"
  )

df_largo_abs <- df_largo_abs %>%
  mutate(celda = paste0(letra, columna)) %>%
  dplyr::select(celda, absorbance)

absorbance <- df_largo_plate_map %>%
  inner_join(df_largo_abs, by = "celda") %>%
  rename(Well = celda)
  










# 3. Extraer concentraciones de estándares
# Buscar la fila donde está "SAMPLE CONCENTRATIONS"
fila_inicio <- which(apply(df, 1, function(x) any(grepl("SAMPLE CONCENTRATIONS", x))))

# El índice de inicio es dos filas después
indice_inicio <- fila_inicio + 2

# Extraer las 7 filas siguientes y las dos primeras columnas
concentraciones <- df[indice_inicio:(indice_inicio + 6), 1:2]

# Renombrar las columnas
colnames(concentraciones) <- c("study_id", "Concentration")
concentraciones <- concentraciones %>%
  mutate(
    study_id = study_id %>%
      str_replace("Standard", "STD") %>%             # Cambia "Standard" por "STD"
      str_remove("\\s*\\([^\\)]+\\)") %>%            # Elimina lo que está entre paréntesis y el espacio antes
      str_remove(":") %>%                            # Elimina los dos puntos al final
      str_replace_all(" ", "")                       # Elimina todos los espacios
  )


datos_completos <- absorbance %>% 
  left_join(
    concentraciones %>%  dplyr::select(study_id, Concentration),
    by = "study_id"
  )


datos_completos <- datos_completos %>%
  mutate(
    absorbance = as.numeric(absorbance),
    Concentration = as.numeric(Concentration)
  )
datos_completos <- datos_completos %>%
  mutate(absorbance = if_else(absorbance < 0, 0, absorbance))
datos_calibracion <- datos_completos %>%
  filter(!is.na(Concentration))

blanks<-datos_completos%>%
  filter(study_id=="B")





# 2. Datos de calibración
datos_calibracion <- datos_completos %>%
  filter(!is.na(Concentration))
datos_calibracion$absorbance1<-convertToProp(datos_calibracion$absorbance)
# 3. Ajuste 5PL
modelo_5pl <- nplr(x = datos_calibracion$Concentration,
                   y = datos_calibracion$absorbance1,
                   npars = 5)

# 4. Predicción inversa para muestras
muestras <- datos_completos %>%
  filter(is.na(Concentration))

concentraciones_estimadas <- getEstimates(modelo_5pl, muestras$absorbance)

# 5. Asignación de resultados
datos_completos$Concentration[is.na(datos_completos$Concentration)] <- concentraciones_estimadas$x

curva_plot <- ggplot(datos_calibracion, aes(x = Concentration, y = absorbance)) +
  geom_point(color = "blue", size = 3) +
  geom_smooth(
    method = drm,
    method.args = list(fct = LL.5()),
    se = FALSE,
    color = "red"
  ) +
  scale_x_log10() +
  labs(
    title = "Curva de Calibración 4PL",
    x = "Concentración (pg/ml, log10)",
    y = "Absorbancia"
  ) +
  theme_minimal()







# 5. Filtrar estándares válidos
curva_calibracion <- datos_completos %>% 
  filter(str_detect(study_id, "^STD[1-7]$")) %>% 
  mutate(Concentration = as.numeric(Concentration)) %>% 
  drop_na(Concentration)


datos_completos <- datos_completos %>%
  mutate(absorbance = as.numeric(absorbance))

curva_calibracion <- curva_calibracion %>%
  mutate(absorbance = as.numeric(absorbance))

min.absorb <- B
max.absorb <- max(standard$absorbance)




# 6. Ajustar modelo 4PL
modelo <- drm(
  absorbance ~ Concentration,
  data = curva_calibracion,
  fct = LL.4()
)
muestras <- datos_completos %>% 
  filter(is.na(Concentration))%>% 
  mutate(
    Concentration = ifelse(
      absorbance < min(curva_calibracion$absorbance),
      0,
      predict(modelo_4pl, newdata = data.frame(absorbance = absorbance))
    )
  )




datos_completos <- datos_completos%>%
dplyr::mutate(log2concentration = log2(Concentration))
Amin <- median(blanks$absorbance)
Amax = max(datos_completos$absorbance, na.rm=TRUE)
log2EC50 = median(datos_completos$log2concentration, na.rm=TRUE)
n = -1 # o prueba con 1, -0.5, etc.
levm.fit <- nlsLM(
  absorbance ~ Amin + ((Amax - Amin) / (1 + exp(n * (log2concentration - log2EC50)))),
  data = datos_completos,
  start = list(
    Amin = min(datos_completos$absorbance, na.rm=TRUE),
    Amax = max(datos_completos$absorbance, na.rm=TRUE),
    log2EC50 = median(datos_completos$log2concentration, na.rm=TRUE),
    n = -1
  )
)



levm.fit <- nlsLM(
  absorbance ~
    Amin + ((Amax - Amin) / (1 + exp(n * (log2concentration - log2EC50)))),
  data = datos_completos,
  start = list(Amin = 0, Amax = 2, log2EC50 = 0.5, n = -1)
)
# Show summary of the fit
summary(levm.fit)



## Add calibration curve to plot
# Save estimated parameters
Amax.est <- coef(levm.fit)[[1]]
log2EC50.est <- coef(levm.fit)[[2]]
n.est <- coef(levm.fit)[[3]]
# Full range of considered concentrations
log2concentration.new <- log2(seq(0.008, 2, by = 0.001))
# Predict new absorbance values from calibration curve
calibration.curve <- tibble(
  concentration.new = 2^log2concentration.new,
  log2concentration.new = log2concentration.new,
  absorbance.est = predict(levm.fit, tibble(log2concentration = log2concentration.new))
)

## Compute prediction and confidence intervals
# TODO: experimental -- to double-check
fgh <- deriv(
  absorbance ~
    Amin + ((Amax - Amin) / (1 + exp(n * (log2concentration - log2EC50)))),
  c("Amax", "log2EC50", "n"), function(Amax, log2EC50, n, log2concentration) {}
)
f.new <- fgh(Amax.est, log2EC50.est, n.est, log2concentration.new)
g.new <- attr(f.new, "gradient")
cov.fit <- vcov(levm.fit)
gs <- rowSums((g.new %*% cov.fit) * g.new)
alpha <- 0.05
delta.f <- sqrt(gs) * qt(1 - alpha/2, 15)
sigma.est <- summary(levm.fit)$sigma
delta.y <- sqrt(gs + sigma.est^2) * qt(1 - alpha / 2, 15)
calibration.curve <- calibration.curve %>%
  dplyr::mutate(pred.lwr = absorbance.est - delta.y) %>%
  dplyr::mutate(pred.upr = absorbance.est + delta.y) %>%
  dplyr::mutate(conf.lwr = absorbance.est - delta.f) %>%
  dplyr::mutate(conf.upr = absorbance.est + delta.f)

## Plot calibration curve
plot.ccurve <- plot.std(semilog = TRUE) +
  geom_line(
    data = calibration.curve,
    aes(concentration.new, absorbance.est), alpha = 0.8
  ) +
  geom_ribbon(
    data = calibration.curve,
    aes(x = concentration.new, ymin = pred.lwr, ymax = pred.upr),
    alpha = 0.2, fill = "grey"
  ) +
  geom_ribbon(
    data = calibration.curve,
    aes(x = concentration.new, 
        ymin = conf.lwr, 
        ymax = conf.upr),
    alpha = 0.1, fill = NA, color = "grey70", lty = 3
  ) +
  scale_y_continuous(limits = c(-0.3, 3), breaks = seq(-0.3, 3, 0.5)) 


plot.ccurve



blanks<-datos_completos%>%
  filter(study_id=="B")
datos_completos <- datos_completos %>%
  mutate(absorbance = if_else(absorbance < 0, 0, absorbance))
abs_blanco <- mean(blanco)
absorbancias_corregidas <- absorbancias - abs_blanco





library(stringr)

muestras$Dilution <- str_extract(muestras$study_id, "(?<=\\().+?(?=\\))")





# ... después de calcular "muestras" y antes de asignar Concentration ...

# Extraer la dilución del study_id
muestras$Dilution <- str_extract(muestras$study_id, "(?<=\\().+?(?=\\))")
muestras$Dilution_factor <- as.numeric(str_extract(muestras$Dilution, "(?<=:)[0-9]+"))
muestras$Dilution_factor[is.na(muestras$Dilution_factor)] <- 1  # Si no hay dilución, factor 1

# Calcular concentraciones con el modelo
preds <- getEstimates(modelo, muestras$prop)$x

# Ajustar concentración por el factor de dilución
muestras$Concentration <- preds * muestras$Dilution_factor

# Imputar min/max del estándar después del cálculo
muestras <- muestras %>%
  mutate(
    Concentration = case_when(
      absorbance < min_std ~ min_conc * Dilution_factor,
      absorbance > max_std ~ max_conc * Dilution_factor,
      TRUE ~ Concentration
    )
  )

# Actualizar full_data con las concentraciones finales
full_data$Concentration[which(is.na(full_data$Concentration) & !is.na(full_data$absorbance))] <- muestras$Concentration



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
        
        # Extraer la dilución del study_id para todas las muestras
        full_data <- full_data %>%
          mutate(
            Dilution = str_extract(study_id, "(?<=\\().+?(?=\\))"),
            # Extrae TODO después de los dos puntos (si existen) o el valor completo
            Dilution_factor = as.numeric(str_remove(Dilution, "^.*:")),  
            Dilution_factor = ifelse(is.na(Dilution_factor), 1, Dilution_factor)
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
              prop = (absorbance - min_std) / (max_std - min_std),
              prop = pmin(pmax(prop, 0), 1)
            )
          
          # Calcular concentraciones con el modelo
          preds <- getEstimates(modelo, muestras$prop)$x
          
          # Ajustar concentración por el factor de dilución
          muestras$Concentration <- preds * muestras$Dilution_factor
          
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
        dplyr::select(study_id, Concentration, Dilution, Dilution_factor) %>%
          group_by(study_id, Dilution, Dilution_factor) %>%
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
  
  
  
  
  
  
  
  
  process_elisa_folder <- function(folder_path) {
    file_list <- list.files(
      path = folder_path,
      pattern = "\\.xlsx$",
      full.names = TRUE
    )
    
    # Crear directorio para guardar gráficos
    if (!dir.exists("calibration_plots")) {
      dir.create("calibration_plots")
    }
    
    process_single_file <- function(file_path) {
      tryCatch({
        df <- read_excel(file_path, sheet = 2, col_names = FALSE)
        
        # Función para extraer secciones con manejo de errores y nombres únicos
        extract_section <- function(pattern) {
          start_row <- grep(pattern, apply(df, 1, paste, collapse = ""))
          if (length(start_row) == 0) {
            stop("Patrón '", pattern, "' no encontrado")
          }
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
        
        # Extraer la dilución del study_id para todas las muestras
        full_data <- full_data %>%
          mutate(
            Dilution = str_extract(study_id, "(?<=\\().+?(?=\\))"),
            # Extrae TODO después de los dos puntos (si existen) o el valor completo
            Dilution_factor = as.numeric(str_remove(Dilution, "^.*:")),  
            Dilution_factor = ifelse(is.na(Dilution_factor), 1, Dilution_factor)
          )
        full_data <- full_data %>%
          mutate(
            absorbance = as.numeric(absorbance),
            Concentration = as.numeric(Concentration),
            Dilution_factor = as.numeric(Dilution_factor)
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
        
        
        # Obtener valores ajustados y R² del modelo nplr
        fit_values <- getFitValues(modelo)
        fit_quality <- getGoodness(modelo)
        r_cuadrado <- fit_quality$R2
        
        # Texto para el gráfico
        fit_text <- paste("R² =", round(r_cuadrado, 4))
        
        # Gráfico de calibración con texto
        curva <- ggplot(cal_data, aes(x = Concentration, y = prop)) +
          geom_point(size = 3) +
          geom_smooth(
            method = drm,
            method.args = list(fct = LL.5()),
            se = FALSE,
            color = "red"
          ) +
          annotate(
            "text", 
            x = Inf, y = -Inf, 
            label = fit_text,
            hjust = 1.1, vjust = -1.2,
            color = "blue", size = 5
          ) +
          scale_x_log10() +
          labs(
            title = basename(file_path), 
            y = "Absorbancia (proporción)",
            caption = fit_text  # Texto adicional debajo
          ) +
          theme(plot.caption = element_text(size = 12, color = "red"))
        
        # Guardar curva
        ggsave(
          filename = paste0("calibration_plots/", tools::file_path_sans_ext(basename(file_path)), ".png"),
          plot = curva,
          width = 8,
          height = 6,
          dpi = 300
        )
        
        # Imprimir en consola
        message(paste("Calidad del ajuste para", basename(file_path), ":", fit_text))
        
        # Resultados sin R²
        results <- full_data %>%
          dplyr::select(study_id, Concentration, Dilution, Dilution_factor) %>%
          group_by(study_id, Dilution, Dilution_factor) %>%
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
  

  
  
  
  
  library(openxlsx)
  
  # Cargar el archivo existente
  wb <- loadWorkbook("PCT_septip.xlsx")
  
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
  saveWorkbook(wb, "PCT_septip1.xlsx", overwrite = TRUE)
  
  df[] <- lapply(df, function(col) {
    gsub("\\.", "-", col)
  })
  
  
  
  
  mutate(redcap_event_name.factor = case_when(
    grepl("-1$", Sample) ~ "0h",
    grepl("-2$", Sample) ~ "4h",
    grepl("-3$", Sample) | grepl("=3$", Sample) ~ "24h",  # Agregamos la condición para =3
    TRUE ~ NA_character_  # Para otros casos
  ))
  
  
  resultados <- process_elisa_folder("C:\\Users\\aaceccato\\Documents\\PCT")
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
  
  
  resultados_limpios <- results %>%
    mutate(
      # 1. Eliminar cualquier contenido entre paréntesis, incluyendo paréntesis y espacios antes
      study_id_clean = str_remove_all(study_id, "\\s*\\([^\\)]+\\)"),
      
      # 2. Extraer número de paciente o visita
      study_id_num = str_extract(study_id_clean, "^\\d+"),
      
      # 3. Asignar timepoint basado en sufijo (.1, .2, .3 o -1, -2, -3)
      redcap_event_name.factor = case_when(
        grepl("(\\.|-)1$", study_id_clean) ~ "0h",
        grepl("(\\.|-)2$", study_id_clean) ~ "4h",
        grepl("(\\.|-)3$", study_id_clean) ~ "24h",
        TRUE ~ NA_character_
      )
    )%>%
    # Filtrar solo IDs numéricos (descarta los que tengan letras o NA)
    filter(!is.na(study_id_num)) %>%
    # Renombrar y seleccionar columnas
    rename(PCTt = Concentration) %>%
    dplyr::select(study_id = study_id_num, redcap_event_name.factor, PCTt)
  
  
  write.csv(resultados_limpios,"pct_septip.csv")
  