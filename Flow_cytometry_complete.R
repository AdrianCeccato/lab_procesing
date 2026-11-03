# ================================================================
# ANALISIS DE BACTERIAS IgA+ EN ASPIRADO TRAQUEAL
# ================================================================
#
# Estructura de gating:
#
# root
# ├── singlets
# │   └── Bacterias
# │       └── IgA_pos
# │
# └── Beads_scatter
#     └── Beads
#
# Tríadas:
#   SYTO = SYTO/BC sin anti-IgA-APC
#   PRE  = SYTO/BC + anti-IgA-APC antes de separación magnética
#   FR   = fracción post-enriquecimiento magnético
#
# Gating:
#   singlets:
#       FSC-A x FSC-H
#
#   Bacterias:
#       SSC-A x FITC-A/SYTO
#
#   IgA_pos:
#       APC-A >= percentil del control SYTO
#
#   Beads_scatter:
#       FSC-A x SSC-A
#
#   Beads:
#       SSC-A x BV711-A
#
# IMPORTANTE:
#   - No se usa FITC-A para excluir las beads.
#   - FITC-A x BV711-A se utiliza únicamente como diagnóstico.
#   - El recuento absoluto utiliza beads_n.
# ================================================================


# ================================================================
# 0. PAQUETES
# ================================================================

required_cran <- c(
  "data.table",
  "ggplot2"
)

required_bioc <- c(
  "flowCore",
  "flowWorkspace",
  "flowStats"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(flowStats)
  library(data.table)
  library(ggplot2)
})


# ================================================================
# 1. CONFIGURACIÓN
# ================================================================

data_dir <- "C:\\Users\\aaceccato\\Documents\\OneDrive - tauli.cat\\citometria\\Resultados citometria Aspirado"

output_dir <- file.path(
  data_dir,
  "output_IgA_bacterias"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ------------------------------------------------
# Canales exactos de los FCS
# ------------------------------------------------

CHANNEL_FSC_A <- "FSC-A"
CHANNEL_FSC_H <- "FSC-H"
CHANNEL_SSC_A <- "SSC-A"
CHANNEL_SYTO  <- "FITC-A"
CHANNEL_IGA   <- "APC-A"

# Canal usado para identificar las beads.
#
# Debe validarse visualmente en los gráficos QC.
CHANNEL_BEADS <- "BV711-A"

# ------------------------------------------------
# Parámetros IgA
# ------------------------------------------------

IGA_NEGATIVE_QUANTILE <- 0.995

MIN_BACTERIA_EVENTS_FOR_APC_THRESHOLD <- 50
MIN_BACTERIA_EVENTS_FOR_ANALYSIS <- 50

# ------------------------------------------------
# Parámetros de beads
# ------------------------------------------------
#
# Completar con los valores experimentales reales.
#
# Ejemplo:
# BEADS_ADDED_PER_TUBE <- 50000
# SAMPLE_VOLUME_UL <- 500
# ------------------------------------------------

BEADS_ADDED_PER_TUBE <- 15000
SAMPLE_VOLUME_UL <- 15

# ------------------------------------------------
# Metadata
# ------------------------------------------------

metadata_file <- file.path(
  data_dir,
  "metadata.csv"
)

# ------------------------------------------------
# Archivos de gates
# ------------------------------------------------

singlet_gate_file <- file.path(
  output_dir,
  "gate_singlets.rds"
)

bacteria_gate_file <- file.path(
  output_dir,
  "gate_bacterias.rds"
)

bead_scatter_gate_file <- file.path(
  output_dir,
  "gate_beads_scatter_FSC_SSC.rds"
)

bead_gate_file <- file.path(
  output_dir,
  "gate_beads_SSC_BV711.rds"
)

# ------------------------------------------------
# Visualización
# ------------------------------------------------

MAX_EVENTS_PLOT <- 30000

set.seed(20261002)


# ================================================================
# 2. FUNCIONES AUXILIARES
# ================================================================

safe_sample_events <- function(
    fr,
    max_events = MAX_EVENTS_PLOT) {
  
  n <- nrow(exprs(fr))
  
  if (n == 0) {
    return(fr)
  }
  
  if (n <= max_events) {
    return(fr)
  }
  
  idx <- sample.int(
    n = n,
    size = max_events,
    replace = FALSE
  )
  
  fr[idx, ]
}


extract_events <- function(
    gs,
    sample_name,
    population = "root") {
  
  idx <- match(
    sample_name,
    sampleNames(gs)
  )
  
  if (is.na(idx)) {
    stop(
      "No se encontró la muestra en GatingSet: ",
      sample_name
    )
  }
  
  gh_pop_get_data(
    gs[[idx]],
    population,
    returnType = "flowFrame"
  )
}


n_events <- function(
    gs,
    sample_name,
    population = "root") {
  
  fr <- extract_events(
    gs = gs,
    sample_name = sample_name,
    population = population
  )
  
  nrow(exprs(fr))
}


draw_polygon_gate <- function(
    fr,
    x_channel,
    y_channel,
    title,
    xlab = NULL,
    ylab = NULL) {
  
  dat <- exprs(
    safe_sample_events(fr)
  )
  
  if (nrow(dat) == 0) {
    stop(
      "No hay eventos disponibles para dibujar el gate."
    )
  }
  
  plot(
    dat[, x_channel],
    dat[, y_channel],
    pch = ".",
    cex = 0.35,
    xlab = ifelse(
      is.null(xlab),
      x_channel,
      xlab
    ),
    ylab = ifelse(
      is.null(ylab),
      y_channel,
      ylab
    ),
    main = title
  )
  
  message("")
  message(
    "Haz clic alrededor de la población."
  )
  message(
    "Pulsa ESC cuando hayas terminado."
  )
  message(
    "Se requieren al menos 3 vértices."
  )
  message("")
  
  pts <- locator(type = "l")
  
  if (is.null(pts) || length(pts$x) < 3) {
    stop(
      "No se pudo crear el polígono: ",
      "se requieren al menos 3 puntos."
    )
  }
  
  coords <- cbind(
    pts$x,
    pts$y
  )
  
  colnames(coords) <- c(
    x_channel,
    y_channel
  )
  
  lines(
    rbind(
      coords,
      coords[1, ]
    )[, 1],
    rbind(
      coords,
      coords[1, ]
    )[, 2],
    col = "red",
    lwd = 2
  )
  
  coords
}


validate_channels <- function(
    fs,
    channels_required) {
  
  available <- colnames(fs)
  
  missing_channels <- setdiff(
    channels_required,
    available
  )
  
  if (length(missing_channels) > 0) {
    stop(
      paste0(
        "Faltan canales requeridos en los FCS:\n - ",
        paste(
          missing_channels,
          collapse = "\n - "
        ),
        "\n\nCanales disponibles:\n - ",
        paste(
          available,
          collapse = "\n - "
        )
      )
    )
  }
  
  invisible(TRUE)
}


extract_sample_id_from_filename <- function(x) {
  
  file_name <- basename(x)
  
  x_clean <- sub(
    "\\.fcs$",
    "",
    file_name,
    ignore.case = TRUE
  )
  
  pattern_type <- paste0(
    "(?i)([_ -])",
    "(syto|fr\\+|frplus|fr|post|postsort|postsorting|",
    "igapos|iga_positive|pre|pre-sort|pre_sort|",
    "presort|presorting)",
    "([_ -].*)?$"
  )
  
  x_clean <- sub(
    pattern_type,
    "",
    x_clean,
    perl = TRUE
  )
  
  x_clean <- gsub(
    "[_ -]+$",
    "",
    x_clean
  )
  
  trimws(x_clean)
}


infer_file_type <- function(x) {
  
  x_low <- tolower(
    sub(
      "\\.fcs$",
      "",
      basename(x),
      ignore.case = TRUE
    )
  )
  
  if (grepl(
    "(^|[_ -])syto([_ -]|$)",
    x_low,
    perl = TRUE
  )) {
    return("SYTO")
  }
  
  if (grepl(
    "(^|[_ -])(pre|pre-sort|pre_sort|presort|presorting)([_ -]|$)",
    x_low,
    perl = TRUE
  )) {
    return("PRE")
  }
  
  if (grepl(
    "(^|[_ -])(fr\\+|frplus|fr|post|post-sort|post_sort|postsort|postsorting|igapos|iga_positive)([_ -]|$)",
    x_low,
    perl = TRUE
  )) {
    return("FR")
  }
  
  NA_character_
}


remove_population_if_exists <- function(
    gs,
    population_name) {
  
  pop_paths <- gs_get_pop_paths(
    gs,
    path = "auto"
  )
  
  pop_names <- basename(
    pop_paths
  )
  
  if (population_name %in% pop_names) {
    
    message(
      "Eliminando población existente: ",
      population_name
    )
    
    gs_pop_remove(
      gs,
      population_name
    )
    
    recompute(gs)
  }
  
  invisible(NULL)
}


get_exprs_matrix <- function(
    gs,
    sample_name,
    population) {
  
  fr <- extract_events(
    gs = gs,
    sample_name = sample_name,
    population = population
  )
  
  exprs(fr)
}


subsample_matrix <- function(
    mat,
    max_events = MAX_EVENTS_PLOT) {
  
  if (is.null(mat) || nrow(mat) == 0) {
    return(mat)
  }
  
  if (nrow(mat) <= max_events) {
    return(mat)
  }
  
  idx <- sample.int(
    n = nrow(mat),
    size = max_events,
    replace = FALSE
  )
  
  mat[
    idx,
    ,
    drop = FALSE
  ]
}


plot_empty_panel <- function(
    title,
    message_text = "No hay eventos") {
  
  plot.new()
  
  title(
    main = title
  )
  
  text(
    x = 0.5,
    y = 0.5,
    labels = message_text,
    cex = 1.1
  )
}


get_apc_threshold_for_sample <- function(
    sample_name,
    umbral_df) {
  
  idx <- which(
    umbral_df$SYTO == sample_name |
      umbral_df$PRE == sample_name |
      umbral_df$FR == sample_name
  )
  
  if (length(idx) != 1) {
    return(NA_real_)
  }
  
  umbral_df$threshold_APC[idx]
}


get_sample_id_for_fcs <- function(
    sample_name,
    metadata) {
  
  idx <- match(
    sample_name,
    metadata$sample
  )
  
  if (is.na(idx)) {
    return(NA_character_)
  }
  
  as.character(
    metadata$sample_id[idx]
  )
}


get_sample_type_for_fcs <- function(
    sample_name,
    metadata) {
  
  idx <- match(
    sample_name,
    metadata$sample
  )
  
  if (is.na(idx)) {
    return(NA_character_)
  }
  
  as.character(
    metadata$tipo[idx]
  )
}


make_qc_scatter <- function(
    fr,
    x_channel,
    y_channel,
    polygon_coords = NULL,
    title = "",
    output_pdf) {
  
  dat <- exprs(
    safe_sample_events(fr)
  )
  
  if (nrow(dat) == 0) {
    return(invisible(NULL))
  }
  
  df <- data.frame(
    x = dat[, x_channel],
    y = dat[, y_channel]
  )
  
  p <- ggplot(
    df,
    aes(
      x = x,
      y = y
    )
  ) +
    geom_point(
      shape = ".",
      alpha = 0.35
    ) +
    labs(
      title = title,
      x = x_channel,
      y = y_channel
    ) +
    theme_bw(
      base_size = 11
    )
  
  if (!is.null(polygon_coords)) {
    
    poly_df <- data.frame(
      x = c(
        polygon_coords[, 1],
        polygon_coords[1, 1]
      ),
      y = c(
        polygon_coords[, 2],
        polygon_coords[1, 2]
      )
    )
    
    p <- p +
      geom_path(
        data = poly_df,
        aes(
          x = x,
          y = y
        ),
        inherit.aes = FALSE,
        colour = "red",
        linewidth = 0.8
      )
  }
  
  ggsave(
    filename = output_pdf,
    plot = p,
    width = 7,
    height = 6,
    units = "in"
  )
}


make_qc_apc_density <- function(
    fr,
    threshold,
    title,
    output_pdf) {
  
  dat <- exprs(fr)
  
  if (nrow(dat) == 0) {
    return(invisible(NULL))
  }
  
  df <- data.frame(
    APC = dat[, CHANNEL_IGA]
  )
  
  p <- ggplot(
    df,
    aes(x = APC)
  ) +
    geom_density(
      fill = "#6BAED6",
      alpha = 0.5,
      linewidth = 0.4
    ) +
    geom_vline(
      xintercept = threshold,
      colour = "red",
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    labs(
      title = title,
      x = paste0(
        CHANNEL_IGA,
        " transformado logicle"
      ),
      y = "Densidad"
    ) +
    theme_bw(
      base_size = 11
    )
  
  ggsave(
    filename = output_pdf,
    plot = p,
    width = 7,
    height = 5,
    units = "in"
  )
}


# ================================================================
# 3. LECTURA DE FCS
# ================================================================

fcs_files <- list.files(
  path = data_dir,
  pattern = "\\.fcs$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(fcs_files) == 0) {
  stop(
    "No se encontraron archivos .fcs en: ",
    data_dir
  )
}

message(
  "Número de archivos FCS encontrados: ",
  length(fcs_files)
)

fs_raw <- read.flowSet(
  files = fcs_files,
  transformation = FALSE,
  truncate_max_range = FALSE
)

message("")
message(
  "Canales disponibles:"
)
print(
  colnames(fs_raw)
)

channels_required <- c(
  CHANNEL_FSC_A,
  CHANNEL_FSC_H,
  CHANNEL_SSC_A,
  CHANNEL_SYTO,
  CHANNEL_IGA,
  CHANNEL_BEADS
)

validate_channels(
  fs = fs_raw,
  channels_required = channels_required
)

gs_raw <- GatingSet(
  fs_raw
)

message("")
message(
  "Número de muestras leídas: ",
  length(gs_raw)
)

message(
  "Nombres de muestra:"
)
print(
  sampleNames(gs_raw)
)


# ================================================================
# 4. COMPENSACIÓN
# ================================================================
#
# No se aplica compensación automáticamente.
# Confirmar si los FCS están ya compensados antes de continuar.
# ================================================================

first_keywords <- keyword(
  fs_raw[[1]]
)

spill_keywords <- names(first_keywords)[
  grepl(
    "SPILL|COMP",
    names(first_keywords),
    ignore.case = TRUE
  )
]

if (length(spill_keywords) > 0) {
  
  message("")
  message(
    "Posibles keywords de compensación detectados:"
  )
  print(
    spill_keywords
  )
  
} else {
  
  message("")
  message(
    "No se detectaron keywords obvios de compensación."
  )
}

gs <- gs_raw


# ================================================================
# 5. TRANSFORMACIÓN LOGICLE
# ================================================================

fluoro_to_transform <- unique(
  c(
    CHANNEL_SYTO,
    CHANNEL_IGA,
    CHANNEL_BEADS
  )
)

fluoro_to_transform <- intersect(
  fluoro_to_transform,
  colnames(fs_raw)
)

trans_list <- estimateLogicle(
  gs[[1]],
  channels = fluoro_to_transform
)

gs <- transform(
  gs,
  trans_list
)

message("")
message(
  "Canales transformados mediante logicle: ",
  paste(
    fluoro_to_transform,
    collapse = ", "
  )
)


# ================================================================
# 6. METADATA
# ================================================================

sample_names <- sampleNames(gs)

if (file.exists(metadata_file)) {
  
  message("")
  message(
    "Usando metadata.csv: ",
    metadata_file
  )
  
  metadata <- fread(
    metadata_file
  )
  
  required_metadata_columns <- c(
    "sample_id",
    "fcs_file",
    "tipo"
  )
  
  missing_metadata_columns <- setdiff(
    required_metadata_columns,
    colnames(metadata)
  )
  
  if (length(missing_metadata_columns) > 0) {
    stop(
      "Faltan columnas en metadata.csv: ",
      paste(
        missing_metadata_columns,
        collapse = ", "
      )
    )
  }
  
  metadata[
    ,
    tipo := toupper(
      trimws(tipo)
    )
  ]
  
  metadata[
    ,
    fcs_file := basename(fcs_file)
  ]
  
  fcs_lookup <- data.table(
    fcs_file = basename(fcs_files),
    sample = sample_names
  )
  
  metadata <- merge(
    metadata,
    fcs_lookup,
    by = "fcs_file",
    all.x = TRUE,
    sort = FALSE
  )
  
  if (any(is.na(metadata$sample))) {
    
    missing_files <- metadata[
      is.na(sample),
      fcs_file
    ]
    
    stop(
      paste0(
        "Los siguientes archivos de metadata.csv ",
        "no se encontraron entre los FCS:\n - ",
        paste(
          missing_files,
          collapse = "\n - "
        )
      )
    )
  }
  
} else {
  
  message("")
  message(
    "No se encontró metadata.csv."
  )
  message(
    "Se intentará inferir sample_id y tipo."
  )
  
  metadata <- data.table(
    sample = sample_names,
    fcs_file = basename(fcs_files)
  )
  
  metadata[
    ,
    sample_id := vapply(
      fcs_file,
      extract_sample_id_from_filename,
      character(1)
    )
  ]
  
  metadata[
    ,
    tipo := vapply(
      fcs_file,
      infer_file_type,
      character(1)
    )
  ]
}

metadata[
  ,
  tipo := toupper(
    trimws(tipo)
  )
]

valid_types <- c(
  "SYTO",
  "PRE",
  "FR"
)

excluded_files <- metadata[
  is.na(tipo) |
    !tipo %in% valid_types
]

if (nrow(excluded_files) > 0) {
  
  excluded_files[
    ,
    motivo_exclusion :=
      "Tipo no reconocido: no contiene SYTO, PRE o FR"
  ]
  
  fwrite(
    excluded_files,
    file.path(
      output_dir,
      "excluded_files_unrecognized_type.csv"
    )
  )
  
  message("")
  message(
    "Archivos excluidos por tipo no reconocido: ",
    nrow(excluded_files)
  )
}

metadata <- metadata[
  !is.na(tipo) &
    tipo %in% valid_types
]

if (nrow(metadata) == 0) {
  stop(
    "No queda ningún FCS válido después de excluir ",
    "los archivos sin tipo reconocido."
  )
}

if (any(
  is.na(metadata$sample_id) |
  trimws(metadata$sample_id) == ""
)) {
  stop(
    "Hay archivos válidos sin sample_id."
  )
}

duplicates <- metadata[
  ,
  .N,
  by = .(
    sample_id,
    tipo
  )
][N > 1]

if (nrow(duplicates) > 0) {
  
  stop(
    paste0(
      "Hay duplicados por sample_id + tipo.\n",
      "Resolverlos antes de continuar:\n\n",
      paste(
        paste0(
          duplicates$sample_id,
          " / ",
          duplicates$tipo,
          " / n = ",
          duplicates$N
        ),
        collapse = "\n"
      )
    )
  )
}

metadata[
  ,
  tipo := factor(
    tipo,
    levels = c(
      "SYTO",
      "PRE",
      "FR"
    )
  )
]

metadata <- metadata[
  order(
    sample_id,
    tipo
  )
]

fwrite(
  metadata,
  file.path(
    output_dir,
    "metadata_checked.csv"
  )
)


# ================================================================
# 7. TRÍADAS COMPLETAS
# ================================================================

pairs <- dcast(
  metadata,
  sample_id ~ tipo,
  value.var = "sample"
)

for (required_type in c(
  "SYTO",
  "PRE",
  "FR"
)) {
  
  if (!required_type %in% colnames(pairs)) {
    
    pairs[
      ,
      (required_type) := NA_character_
    ]
  }
}

pairs[
  ,
  completa :=
    !is.na(SYTO) &
    !is.na(PRE) &
    !is.na(FR)
]

incomplete_pairs <- pairs[
  completa == FALSE
]

if (nrow(incomplete_pairs) > 0) {
  
  message("")
  message(
    "Tríadas incompletas excluidas:"
  )
  print(
    incomplete_pairs
  )
  
  fwrite(
    incomplete_pairs,
    file.path(
      output_dir,
      "incomplete_triplets.csv"
    )
  )
}

pairs_complete <- pairs[
  completa == TRUE
]

if (nrow(pairs_complete) == 0) {
  stop(
    "No se encontraron tríadas completas SYTO/PRE/FR."
  )
}

samples_keep <- unique(
  c(
    pairs_complete$SYTO,
    pairs_complete$PRE,
    pairs_complete$FR
  )
)

gs <- gs[
  samples_keep
]

metadata <- metadata[
  sample %in% samples_keep
]

setorder(
  metadata,
  sample_id,
  tipo
)

message("")
message(
  "Número de tríadas completas: ",
  nrow(pairs_complete)
)

message(
  "Número de FCS incluidos: ",
  length(gs)
)


# ================================================================
# 8. ELIMINAR GATES PREVIOS
# ================================================================

remove_population_if_exists(
  gs,
  "IgA_pos"
)

remove_population_if_exists(
  gs,
  "Bacterias"
)

remove_population_if_exists(
  gs,
  "singlets"
)

remove_population_if_exists(
  gs,
  "Beads"
)

remove_population_if_exists(
  gs,
  "Beads_scatter"
)


# ================================================================
# 9. GATE DE SINGLETS
# ================================================================

singlet_reference_sample <- pairs_complete$SYTO[1]

fr_root_reference <- extract_events(
  gs = gs,
  sample_name = singlet_reference_sample,
  population = "root"
)

if (file.exists(singlet_gate_file)) {
  
  message("")
  message(
    "Cargando gate de singlets:"
  )
  message(singlet_gate_file)
  
  g_sing <- readRDS(
    singlet_gate_file
  )
  
} else {
  
  message("")
  message(
    "Calculando gate de singlets usando:"
  )
  message(singlet_reference_sample)
  
  g_sing <- flowStats::gate_singlet(
    fr_root_reference,
    area = CHANNEL_FSC_A,
    height = CHANNEL_FSC_H,
    filterId = "singlets"
  )
  
  saveRDS(
    g_sing,
    singlet_gate_file
  )
}

gs_pop_add(
  gs,
  g_sing,
  parent = "root",
  name = "singlets"
)

recompute(gs)

message(
  "Gate singlets añadido."
)


# ================================================================
# 10. GATE DE BACTERIAS
# ================================================================

bacteria_reference_sample <- pairs_complete$SYTO[1]

fr_bacteria_reference <- extract_events(
  gs = gs,
  sample_name = bacteria_reference_sample,
  population = "singlets"
)

if (file.exists(bacteria_gate_file)) {
  
  message("")
  message(
    "Cargando gate bacteriano:"
  )
  message(bacteria_gate_file)
  
  coords_bact <- readRDS(
    bacteria_gate_file
  )
  
} else {
  
  message("")
  message(
    "Dibuja el gate bacteriano en SSC-A x FITC-A."
  )
  message(
    "Muestra de referencia: ",
    bacteria_reference_sample
  )
  
  coords_bact <- draw_polygon_gate(
    fr = fr_bacteria_reference,
    x_channel = CHANNEL_SSC_A,
    y_channel = CHANNEL_SYTO,
    title = paste0(
      "Bacterias: SSC-A x ",
      CHANNEL_SYTO,
      "\n",
      bacteria_reference_sample
    ),
    xlab = CHANNEL_SSC_A,
    ylab = paste0(
      CHANNEL_SYTO,
      " (SYTO; logicle)"
    )
  )
  
  saveRDS(
    coords_bact,
    bacteria_gate_file
  )
}

if (!all(
  c(
    CHANNEL_SSC_A,
    CHANNEL_SYTO
  ) %in%
  colnames(coords_bact)
)) {
  stop(
    "Las coordenadas del gate bacteriano no corresponden a ",
    CHANNEL_SSC_A,
    " y ",
    CHANNEL_SYTO,
    "."
  )
}

g_bact <- polygonGate(
  .gate = coords_bact,
  filterId = "Bacterias"
)

gs_pop_add(
  gs,
  g_bact,
  parent = "singlets",
  name = "Bacterias"
)

recompute(gs)

message(
  "Gate Bacterias añadido."
)


# ================================================================
# 11. GATING DE COUNTING BEADS
# ================================================================
#
# root
# └── Beads_scatter
#     └── Beads
#
# Beads_scatter:
#   FSC-A x SSC-A
#
# Beads:
#   SSC-A x BV711-A
#
# No se utiliza FITC-A para excluir las beads.
# ================================================================

bead_gate_available <- FALSE

coords_beads_scatter <- NULL
coords_beads <- NULL

if (!is.null(CHANNEL_BEADS)) {
  
  bead_reference_sample <- pairs_complete$SYTO[1]
  
  fr_root_beads <- extract_events(
    gs = gs,
    sample_name = bead_reference_sample,
    population = "root"
  )
  
  # ------------------------------------------------
  # 11.1 Beads_scatter: FSC-A x SSC-A
  # ------------------------------------------------
  
  if (file.exists(bead_scatter_gate_file)) {
    
    message("")
    message(
      "Cargando gate Beads_scatter:"
    )
    message(bead_scatter_gate_file)
    
    coords_beads_scatter <- readRDS(
      bead_scatter_gate_file
    )
    
  } else {
    
    message("")
    message(
      "Dibuja una región amplia de beads en FSC-A x SSC-A."
    )
    message(
      "No intentes separar beads y bacterias solo por tamaño."
    )
    
    coords_beads_scatter <- draw_polygon_gate(
      fr = fr_root_beads,
      x_channel = CHANNEL_FSC_A,
      y_channel = CHANNEL_SSC_A,
      title = paste0(
        "Beads_scatter: FSC-A x SSC-A\n",
        bead_reference_sample
      ),
      xlab = CHANNEL_FSC_A,
      ylab = CHANNEL_SSC_A
    )
    
    saveRDS(
      coords_beads_scatter,
      bead_scatter_gate_file
    )
  }
  
  if (!all(
    c(
      CHANNEL_FSC_A,
      CHANNEL_SSC_A
    ) %in%
    colnames(coords_beads_scatter)
  )) {
    stop(
      "El gate Beads_scatter no corresponde a ",
      CHANNEL_FSC_A,
      " y ",
      CHANNEL_SSC_A,
      "."
    )
  }
  
  g_beads_scatter <- polygonGate(
    .gate = coords_beads_scatter,
    filterId = "Beads_scatter"
  )
  
  gs_pop_add(
    gs,
    g_beads_scatter,
    parent = "root",
    name = "Beads_scatter"
  )
  
  recompute(gs)
  
  message(
    "Gate Beads_scatter añadido."
  )
  
  # ------------------------------------------------
  # 11.2 Beads: SSC-A x BV711-A
  # ------------------------------------------------
  
  fr_scatter_beads <- extract_events(
    gs = gs,
    sample_name = bead_reference_sample,
    population = "Beads_scatter"
  )
  
  if (nrow(exprs(fr_scatter_beads)) == 0) {
    stop(
      "Beads_scatter no contiene eventos en la muestra de referencia."
    )
  }
  
  if (file.exists(bead_gate_file)) {
    
    message("")
    message(
      "Cargando gate final de beads:"
    )
    message(bead_gate_file)
    
    coords_beads <- readRDS(
      bead_gate_file
    )
    
  } else {
    
    message("")
    message(
      "Dibuja el gate final de beads en SSC-A x BV711-A."
    )
    message(
      "Selecciona la población bead compacta y fluorescente."
    )
    message(
      "No utilices FITC-A/SYTO para excluir eventos."
    )
    
    coords_beads <- draw_polygon_gate(
      fr = fr_scatter_beads,
      x_channel = CHANNEL_SSC_A,
      y_channel = CHANNEL_BEADS,
      title = paste0(
        "Beads: SSC-A x ",
        CHANNEL_BEADS,
        "\n",
        bead_reference_sample
      ),
      xlab = CHANNEL_SSC_A,
      ylab = paste0(
        CHANNEL_BEADS,
        " (logicle)"
      )
    )
    
    saveRDS(
      coords_beads,
      bead_gate_file
    )
  }
  
  if (!all(
    c(
      CHANNEL_SSC_A,
      CHANNEL_BEADS
    ) %in%
    colnames(coords_beads)
  )) {
    stop(
      "El gate Beads no corresponde a ",
      CHANNEL_SSC_A,
      " y ",
      CHANNEL_BEADS,
      "."
    )
  }
  
  g_beads <- polygonGate(
    .gate = coords_beads,
    filterId = "Beads"
  )
  
  gs_pop_add(
    gs,
    g_beads,
    parent = "Beads_scatter",
    name = "Beads"
  )
  
  recompute(gs)
  
  bead_gate_available <- TRUE
  
  message(
    "Gate Beads añadido bajo Beads_scatter."
  )
}


# ================================================================
# 12. GATE IgA_pos POR TRÍADA
# ================================================================

iga_gates <- vector(
  mode = "list",
  length = length(gs)
)

names(iga_gates) <- sampleNames(gs)

threshold_results <- vector(
  mode = "list",
  length = nrow(pairs_complete)
)

for (k in seq_len(
  nrow(pairs_complete)
)) {
  
  sample_id_k <- pairs_complete$sample_id[k]
  sample_syto <- pairs_complete$SYTO[k]
  sample_pre <- pairs_complete$PRE[k]
  sample_fr <- pairs_complete$FR[k]
  
  fr_ctrl_bact <- extract_events(
    gs = gs,
    sample_name = sample_syto,
    population = "Bacterias"
  )
  
  n_bact_ctrl <- nrow(
    exprs(fr_ctrl_bact)
  )
  
  if (
    n_bact_ctrl <
    MIN_BACTERIA_EVENTS_FOR_APC_THRESHOLD
  ) {
    warning(
      paste0(
        "Muestra ",
        sample_id_k,
        ": solo hay ",
        n_bact_ctrl,
        " eventos bacterianos en SYTO. ",
        "El umbral APC puede ser inestable."
      )
    )
  }
  
  if (n_bact_ctrl == 0) {
    stop(
      "Muestra ",
      sample_id_k,
      ": no hay bacterias en SYTO."
    )
  }
  
  apc_ctrl <- exprs(
    fr_ctrl_bact
  )[, CHANNEL_IGA]
  
  threshold_apc <- unname(
    quantile(
      apc_ctrl,
      probs = IGA_NEGATIVE_QUANTILE,
      na.rm = TRUE,
      names = FALSE
    )
  )
  
  if (!is.finite(threshold_apc)) {
    stop(
      "No se pudo calcular un umbral APC válido para ",
      sample_id_k
    )
  }
  
  gate_iga_k <- do.call(
    rectangleGate,
    c(
      list(
        filterId = "IgA_pos"
      ),
      setNames(
        list(
          c(
            threshold_apc,
            Inf
          )
        ),
        CHANNEL_IGA
      )
    )
  )
  
  iga_gates[[sample_syto]] <- gate_iga_k
  iga_gates[[sample_pre]] <- gate_iga_k
  iga_gates[[sample_fr]] <- gate_iga_k
  
  threshold_results[[k]] <- data.frame(
    sample_id = sample_id_k,
    SYTO = sample_syto,
    PRE = sample_pre,
    FR = sample_fr,
    n_bacterias_control_SYTO = n_bact_ctrl,
    quantile_negative = IGA_NEGATIVE_QUANTILE,
    threshold_APC = threshold_apc,
    stringsAsFactors = FALSE
  )
}

if (any(
  vapply(
    iga_gates,
    is.null,
    logical(1)
  )
)) {
  stop(
    "No se generó un gate IgA para todos los FCS."
  )
}

gs_pop_add(
  gs,
  iga_gates,
  parent = "Bacterias",
  name = "IgA_pos"
)

recompute(gs)

umbral_df <- rbindlist(
  threshold_results
)

fwrite(
  umbral_df,
  file.path(
    output_dir,
    "iga_thresholds_by_triplet.csv"
  )
)

message(
  "Gate IgA_pos añadido."
)


# ================================================================
# 13. RESUMEN DE GATING
# ================================================================

stats_list <- lapply(
  sampleNames(gs),
  function(sample_name) {
    
    n_root <- n_events(
      gs,
      sample_name,
      "root"
    )
    
    n_singlets <- n_events(
      gs,
      sample_name,
      "singlets"
    )
    
    n_bacterias <- n_events(
      gs,
      sample_name,
      "Bacterias"
    )
    
    n_iga_pos <- n_events(
      gs,
      sample_name,
      "IgA_pos"
    )
    
    n_beads_scatter <- NA_real_
    n_beads <- NA_real_
    
    if (bead_gate_available) {
      
      n_beads_scatter <- n_events(
        gs,
        sample_name,
        "Beads_scatter"
      )
      
      n_beads <- n_events(
        gs,
        sample_name,
        "Beads"
      )
    }
    
    data.frame(
      sample = sample_name,
      root_n = n_root,
      singlets_n = n_singlets,
      bacterias_n = n_bacterias,
      iga_pos_n = n_iga_pos,
      beads_scatter_n = n_beads_scatter,
      beads_n = n_beads,
      singlets_pct_root = ifelse(
        n_root > 0,
        100 * n_singlets / n_root,
        NA_real_
      ),
      bacterias_pct_singlets = ifelse(
        n_singlets > 0,
        100 * n_bacterias / n_singlets,
        NA_real_
      ),
      iga_pos_pct_in_bact = ifelse(
        n_bacterias > 0,
        100 * n_iga_pos / n_bacterias,
        NA_real_
      ),
      stringsAsFactors = FALSE
    )
  }
)

stats_df <- rbindlist(
  stats_list
)

stats_df <- merge(
  metadata[
    ,
    .(
      sample,
      sample_id,
      tipo,
      fcs_file
    )
  ],
  stats_df,
  by = "sample",
  all.y = TRUE,
  sort = FALSE
)

setorder(
  stats_df,
  sample_id,
  tipo
)

# ------------------------------------------------
# Cuantificación absoluta
# ------------------------------------------------

stats_df[
  ,
  bacterias_per_tube := NA_real_
]

stats_df[
  ,
  iga_pos_per_tube := NA_real_
]

stats_df[
  ,
  bacterias_per_uL := NA_real_
]

stats_df[
  ,
  iga_pos_per_uL := NA_real_
]

if (
  bead_gate_available &&
  is.finite(BEADS_ADDED_PER_TUBE) &&
  BEADS_ADDED_PER_TUBE > 0
) {
  
  stats_df[
    ,
    bacterias_per_tube :=
      ifelse(
        beads_n > 0,
        bacterias_n /
          beads_n *
          BEADS_ADDED_PER_TUBE,
        NA_real_
      )
  ]
  
  stats_df[
    ,
    iga_pos_per_tube :=
      ifelse(
        beads_n > 0,
        iga_pos_n /
          beads_n *
          BEADS_ADDED_PER_TUBE,
        NA_real_
      )
  ]
  
  if (
    is.finite(SAMPLE_VOLUME_UL) &&
    SAMPLE_VOLUME_UL > 0
  ) {
    
    stats_df[
      ,
      bacterias_per_uL :=
        bacterias_per_tube /
        SAMPLE_VOLUME_UL
    ]
    
    stats_df[
      ,
      iga_pos_per_uL :=
        iga_pos_per_tube /
        SAMPLE_VOLUME_UL
    ]
  }
}

# ------------------------------------------------
# Flags QC
# ------------------------------------------------

stats_df[
  ,
  qc_low_bacteria :=
    bacterias_n <
    MIN_BACTERIA_EVENTS_FOR_ANALYSIS
]

stats_df[
  ,
  qc_zero_beads :=
    ifelse(
      bead_gate_available,
      is.na(beads_n) |
        beads_n == 0,
      NA
    )
]

stats_df[
  ,
  qc_beads_low :=
    ifelse(
      bead_gate_available,
      !is.na(beads_n) &
        beads_n < 100,
      NA
    )
]

stats_df[
  ,
  qc_comment := "OK"
]

stats_df[
  qc_low_bacteria == TRUE,
  qc_comment := paste(
    qc_comment,
    "Pocos eventos bacterianos",
    sep = "; "
  )
]

stats_df[
  qc_zero_beads == TRUE,
  qc_comment := paste(
    qc_comment,
    "No se detectaron beads",
    sep = "; "
  )
]

stats_df[
  qc_beads_low == TRUE,
  qc_comment := paste(
    qc_comment,
    "Pocos eventos bead",
    sep = "; "
  )
]

fwrite(
  stats_df,
  file.path(
    output_dir,
    "gating_summary.csv"
  )
)

message(
  "gating_summary.csv guardado."
)


# ================================================================
# 14. RESUMEN PRE VS FR
# ================================================================

summary_variables <- c(
  "bacterias_n",
  "iga_pos_n",
  "iga_pos_pct_in_bact",
  "beads_n",
  "bacterias_per_uL",
  "iga_pos_per_uL"
)

summary_variables <- intersect(
  summary_variables,
  colnames(stats_df)
)

summary_wide <- dcast(
  stats_df,
  sample_id ~ tipo,
  value.var = summary_variables
)

pre_pct_col <- "iga_pos_pct_in_bact_PRE"
fr_pct_col <- "iga_pos_pct_in_bact_FR"

if (
  pre_pct_col %in% colnames(summary_wide) &&
  fr_pct_col %in% colnames(summary_wide)
) {
  
  summary_wide[
    ,
    iga_pct_difference_FR_minus_PRE :=
      get(fr_pct_col) -
      get(pre_pct_col)
  ]
  
  summary_wide[
    ,
    iga_pct_fold_enrichment_FR_vs_PRE :=
      ifelse(
        get(pre_pct_col) > 0,
        get(fr_pct_col) /
          get(pre_pct_col),
        NA_real_
      )
  ]
}

fwrite(
  summary_wide,
  file.path(
    output_dir,
    "summary_PRE_vs_FR.csv"
  )
)


# ================================================================
# 15. QC INDIVIDUAL
# ================================================================

qc_dir <- file.path(
  output_dir,
  "QC_plots"
)

dir.create(
  qc_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ------------------------------------------------
# 15.1 Bacterias
# ------------------------------------------------

for (sample_name in sampleNames(gs)) {
  
  fr_singlets <- extract_events(
    gs,
    sample_name,
    "singlets"
  )
  
  make_qc_scatter(
    fr = fr_singlets,
    x_channel = CHANNEL_SSC_A,
    y_channel = CHANNEL_SYTO,
    polygon_coords = coords_bact,
    title = paste(
      sample_name,
      "- Bacterias"
    ),
    output_pdf = file.path(
      qc_dir,
      paste0(
        gsub(
          "[^A-Za-z0-9_\\-]",
          "_",
          sample_name
        ),
        "_QC_bacterias.pdf"
      )
    )
  )
}

# ------------------------------------------------
# 15.2 Beads
# ------------------------------------------------

if (bead_gate_available) {
  
  for (sample_name in sampleNames(gs)) {
    
    fr_root <- extract_events(
      gs,
      sample_name,
      "root"
    )
    
    fr_beads_scatter <- extract_events(
      gs,
      sample_name,
      "Beads_scatter"
    )
    
    fr_beads <- extract_events(
      gs,
      sample_name,
      "Beads"
    )
    
    # FSC-A x SSC-A
    make_qc_scatter(
      fr = fr_root,
      x_channel = CHANNEL_FSC_A,
      y_channel = CHANNEL_SSC_A,
      polygon_coords = coords_beads_scatter,
      title = paste(
        sample_name,
        "- Beads_scatter"
      ),
      output_pdf = file.path(
        qc_dir,
        paste0(
          gsub(
            "[^A-Za-z0-9_\\-]",
            "_",
            sample_name
          ),
          "_QC_beads_scatter.pdf"
        )
      )
    )
    
    # SSC-A x BV711-A
    make_qc_scatter(
      fr = fr_beads_scatter,
      x_channel = CHANNEL_SSC_A,
      y_channel = CHANNEL_BEADS,
      polygon_coords = coords_beads,
      title = paste(
        sample_name,
        "- Beads"
      ),
      output_pdf = file.path(
        qc_dir,
        paste0(
          gsub(
            "[^A-Za-z0-9_\\-]",
            "_",
            sample_name
          ),
          "_QC_beads.pdf"
        )
      )
    )
    
    # FITC-A x BV711-A como diagnóstico
    make_qc_scatter(
      fr = fr_beads_scatter,
      x_channel = CHANNEL_SYTO,
      y_channel = CHANNEL_BEADS,
      polygon_coords = NULL,
      title = paste(
        sample_name,
        "- Diagnóstico SYTO x BV711"
      ),
      output_pdf = file.path(
        qc_dir,
        paste0(
          gsub(
            "[^A-Za-z0-9_\\-]",
            "_",
            sample_name
          ),
          "_QC_beads_SYTO_BV711_diagnostic.pdf"
        )
      )
    )
  }
}

# ------------------------------------------------
# 15.3 APC por tríada
# ------------------------------------------------

for (k in seq_len(
  nrow(umbral_df)
)) {
  
  sample_id_k <- umbral_df$sample_id[k]
  threshold_k <- umbral_df$threshold_APC[k]
  
  samples_k <- c(
    SYTO = umbral_df$SYTO[k],
    PRE = umbral_df$PRE[k],
    FR = umbral_df$FR[k]
  )
  
  for (type_k in names(samples_k)) {
    
    sample_k <- samples_k[[type_k]]
    
    fr_bact <- extract_events(
      gs,
      sample_k,
      "Bacterias"
    )
    
    make_qc_apc_density(
      fr = fr_bact,
      threshold = threshold_k,
      title = paste0(
        sample_id_k,
        " - ",
        type_k,
        " - APC dentro de Bacterias"
      ),
      output_pdf = file.path(
        qc_dir,
        paste0(
          gsub(
            "[^A-Za-z0-9_\\-]",
            "_",
            sample_id_k
          ),
          "_",
          type_k,
          "_QC_APC_density.pdf"
        )
      )
    )
  }
}


# ================================================================
# 16. PDF DE GATING COMPLETO
# ================================================================

plot_gating_page <- function(
    gs,
    sample_name,
    metadata,
    umbral_df,
    coords_bact,
    coords_beads = NULL,
    bead_gate_available = FALSE) {
  
  sample_id <- get_sample_id_for_fcs(
    sample_name,
    metadata
  )
  
  sample_type <- get_sample_type_for_fcs(
    sample_name,
    metadata
  )
  
  threshold_apc <- get_apc_threshold_for_sample(
    sample_name,
    umbral_df
  )
  
  root_mat <- get_exprs_matrix(
    gs,
    sample_name,
    "root"
  )
  
  singlet_mat <- get_exprs_matrix(
    gs,
    sample_name,
    "singlets"
  )
  
  bacteria_mat <- get_exprs_matrix(
    gs,
    sample_name,
    "Bacterias"
  )
  
  iga_pos_mat <- get_exprs_matrix(
    gs,
    sample_name,
    "IgA_pos"
  )
  
  beads_mat <- NULL
  
  if (bead_gate_available) {
    
    beads_mat <- get_exprs_matrix(
      gs,
      sample_name,
      "Beads"
    )
  }
  
  root_plot <- subsample_matrix(
    root_mat
  )
  
  singlet_plot <- subsample_matrix(
    singlet_mat
  )
  
  bacteria_plot <- subsample_matrix(
    bacteria_mat
  )
  
  iga_pos_plot <- subsample_matrix(
    iga_pos_mat
  )
  
  beads_plot <- NULL
  
  if (bead_gate_available) {
    beads_plot <- subsample_matrix(
      beads_mat
    )
  }
  
  n_root <- nrow(root_mat)
  n_singlets <- nrow(singlet_mat)
  n_bacteria <- nrow(bacteria_mat)
  n_iga_pos <- nrow(iga_pos_mat)
  n_beads <- NA_integer_
  
  if (bead_gate_available) {
    n_beads <- nrow(beads_mat)
  }
  
  iga_pct <- ifelse(
    n_bacteria > 0,
    100 * n_iga_pos / n_bacteria,
    NA_real_
  )
  
  par(
    mfrow = c(2, 2),
    mar = c(4.1, 4.1, 3.5, 1.2),
    oma = c(0, 0, 3.5, 0)
  )
  
  # Panel 1: singlets
  if (nrow(root_plot) == 0) {
    
    plot_empty_panel(
      "1. root → singlets",
      "No hay eventos root"
    )
    
  } else {
    
    plot(
      root_plot[, CHANNEL_FSC_A],
      root_plot[, CHANNEL_FSC_H],
      pch = ".",
      cex = 0.4,
      col = grDevices::adjustcolor(
        "grey50",
        alpha.f = 0.35
      ),
      xlab = CHANNEL_FSC_A,
      ylab = CHANNEL_FSC_H,
      main = paste0(
        "1. root → singlets\n",
        "root: ",
        format(
          n_root,
          big.mark = ","
        ),
        " | singlets: ",
        format(
          n_singlets,
          big.mark = ","
        )
      )
    )
    
    if (nrow(singlet_plot) > 0) {
      
      points(
        singlet_plot[, CHANNEL_FSC_A],
        singlet_plot[, CHANNEL_FSC_H],
        pch = ".",
        cex = 0.6,
        col = grDevices::adjustcolor(
          "#0072B2",
          alpha.f = 0.7
        )
      )
    }
  }
  
  # Panel 2: bacterias
  if (nrow(singlet_plot) == 0) {
    
    plot_empty_panel(
      "2. singlets → Bacterias",
      "No hay singlets"
    )
    
  } else {
    
    plot(
      singlet_plot[, CHANNEL_SSC_A],
      singlet_plot[, CHANNEL_SYTO],
      pch = ".",
      cex = 0.4,
      col = grDevices::adjustcolor(
        "grey50",
        alpha.f = 0.35
      ),
      xlab = CHANNEL_SSC_A,
      ylab = CHANNEL_SYTO,
      main = paste0(
        "2. singlets → Bacterias\n",
        "Bacterias: ",
        format(
          n_bacteria,
          big.mark = ","
        )
      )
    )
    
    if (nrow(bacteria_plot) > 0) {
      
      points(
        bacteria_plot[, CHANNEL_SSC_A],
        bacteria_plot[, CHANNEL_SYTO],
        pch = ".",
        cex = 0.7,
        col = grDevices::adjustcolor(
          "#009E73",
          alpha.f = 0.75
        )
      )
    }
    
    lines(
      c(
        coords_bact[, CHANNEL_SSC_A],
        coords_bact[1, CHANNEL_SSC_A]
      ),
      c(
        coords_bact[, CHANNEL_SYTO],
        coords_bact[1, CHANNEL_SYTO]
      ),
      col = "red",
      lwd = 2
    )
  }
  
  # Panel 3: beads
  if (!bead_gate_available ||
      nrow(root_plot) == 0) {
    
    plot_empty_panel(
      "3. Beads",
      "Gate de beads no configurado"
    )
    
  } else {
    
    plot(
      root_plot[, CHANNEL_SSC_A],
      root_plot[, CHANNEL_BEADS],
      pch = ".",
      cex = 0.4,
      col = grDevices::adjustcolor(
        "grey50",
        alpha.f = 0.3
      ),
      xlab = CHANNEL_SSC_A,
      ylab = CHANNEL_BEADS,
      main = paste0(
        "3. Beads\n",
        "N = ",
        format(
          n_beads,
          big.mark = ","
        )
      )
    )
    
    if (nrow(beads_plot) > 0) {
      
      points(
        beads_plot[, CHANNEL_SSC_A],
        beads_plot[, CHANNEL_BEADS],
        pch = ".",
        cex = 0.75,
        col = grDevices::adjustcolor(
          "#D55E00",
          alpha.f = 0.85
        )
      )
    }
    
    if (!is.null(coords_beads)) {
      
      lines(
        c(
          coords_beads[, CHANNEL_SSC_A],
          coords_beads[1, CHANNEL_SSC_A]
        ),
        c(
          coords_beads[, CHANNEL_BEADS],
          coords_beads[1, CHANNEL_BEADS]
        ),
        col = "red",
        lwd = 2
      )
    }
  }
  
  # Panel 4: IgA
  if (nrow(bacteria_plot) == 0) {
    
    plot_empty_panel(
      "4. Bacterias → IgA_pos",
      "No hay bacterias"
    )
    
  } else {
    
    plot(
      bacteria_plot[, CHANNEL_SYTO],
      bacteria_plot[, CHANNEL_IGA],
      pch = ".",
      cex = 0.45,
      col = grDevices::adjustcolor(
        "grey50",
        alpha.f = 0.4
      ),
      xlab = CHANNEL_SYTO,
      ylab = CHANNEL_IGA,
      main = paste0(
        "4. Bacterias → IgA_pos\n",
        "IgA+: ",
        format(
          n_iga_pos,
          big.mark = ","
        ),
        " (",
        sprintf(
          "%.2f",
          iga_pct
        ),
        "%)"
      )
    )
    
    if (nrow(iga_pos_plot) > 0) {
      
      points(
        iga_pos_plot[, CHANNEL_SYTO],
        iga_pos_plot[, CHANNEL_IGA],
        pch = ".",
        cex = 0.7,
        col = grDevices::adjustcolor(
          "#D55E00",
          alpha.f = 0.8
        )
      )
    }
    
    if (is.finite(threshold_apc)) {
      
      abline(
        h = threshold_apc,
        col = "#0072B2",
        lwd = 2,
        lty = 2
      )
    }
  }
  
  mtext(
    side = 3,
    outer = TRUE,
    line = 1,
    font = 2,
    cex = 1.1,
    text = paste0(
      "Muestra: ",
      sample_id,
      " | Tipo: ",
      sample_type,
      " | FCS: ",
      sample_name
    )
  )
}


gating_pdf_file <- file.path(
  output_dir,
  "QC_gating_completo_por_muestra.pdf"
)

pdf(
  file = gating_pdf_file,
  width = 12,
  height = 10,
  onefile = TRUE
)

for (sample_name in sampleNames(gs)) {
  
  message(
    "Generando QC completo para: ",
    sample_name
  )
  
  tryCatch(
    
    {
      plot_gating_page(
        gs = gs,
        sample_name = sample_name,
        metadata = metadata,
        umbral_df = umbral_df,
        coords_bact = coords_bact,
        coords_beads = coords_beads,
        bead_gate_available = bead_gate_available
      )
    },
    
    error = function(e) {
      
      par(
        mfrow = c(1, 1),
        mar = c(4, 4, 4, 2)
      )
      
      plot.new()
      
      title(
        main = paste(
          "Error QC:",
          sample_name
        )
      )
      
      text(
        x = 0.5,
        y = 0.5,
        labels = conditionMessage(e)
      )
    }
  )
}

dev.off()

message("")
message(
  "PDF de gating completo creado:"
)
message(gating_pdf_file)


# ================================================================
# 17. PDF ESPECIFICO DE BEADS
# ================================================================

plot_bead_qc_page <- function(
    gs,
    sample_name,
    coords_beads_scatter,
    coords_beads) {
  
  root_mat <- get_exprs_matrix(
    gs,
    sample_name,
    "root"
  )
  
  beads_scatter_mat <- get_exprs_matrix(
    gs,
    sample_name,
    "Beads_scatter"
  )
  
  beads_mat <- get_exprs_matrix(
    gs,
    sample_name,
    "Beads"
  )
  
  root_plot <- subsample_matrix(
    root_mat
  )
  
  beads_scatter_plot <- subsample_matrix(
    beads_scatter_mat
  )
  
  beads_plot <- subsample_matrix(
    beads_mat
  )
  
  par(
    mfrow = c(2, 2),
    mar = c(4.1, 4.1, 3.5, 1.2),
    oma = c(0, 0, 3.5, 0)
  )
  
  # Panel 1: FSC-A x SSC-A
  plot(
    root_plot[, CHANNEL_FSC_A],
    root_plot[, CHANNEL_SSC_A],
    pch = ".",
    cex = 0.4,
    col = grDevices::adjustcolor(
      "grey50",
      alpha.f = 0.3
    ),
    xlab = CHANNEL_FSC_A,
    ylab = CHANNEL_SSC_A,
    main = "1. Beads_scatter"
  )
  
  if (nrow(beads_scatter_plot) > 0) {
    
    points(
      beads_scatter_plot[, CHANNEL_FSC_A],
      beads_scatter_plot[, CHANNEL_SSC_A],
      pch = ".",
      cex = 0.7,
      col = grDevices::adjustcolor(
        "#0072B2",
        alpha.f = 0.8
      )
    )
  }
  
  lines(
    c(
      coords_beads_scatter[, CHANNEL_FSC_A],
      coords_beads_scatter[1, CHANNEL_FSC_A]
    ),
    c(
      coords_beads_scatter[, CHANNEL_SSC_A],
      coords_beads_scatter[1, CHANNEL_SSC_A]
    ),
    col = "red",
    lwd = 2
  )
  
  # Panel 2: SSC-A x BV711-A
  plot(
    beads_scatter_plot[, CHANNEL_SSC_A],
    beads_scatter_plot[, CHANNEL_BEADS],
    pch = ".",
    cex = 0.45,
    col = grDevices::adjustcolor(
      "grey50",
      alpha.f = 0.35
    ),
    xlab = CHANNEL_SSC_A,
    ylab = CHANNEL_BEADS,
    main = "2. Beads"
  )
  
  if (nrow(beads_plot) > 0) {
    
    points(
      beads_plot[, CHANNEL_SSC_A],
      beads_plot[, CHANNEL_BEADS],
      pch = ".",
      cex = 0.75,
      col = grDevices::adjustcolor(
        "#D55E00",
        alpha.f = 0.85
      )
    )
  }
  
  lines(
    c(
      coords_beads[, CHANNEL_SSC_A],
      coords_beads[1, CHANNEL_SSC_A]
    ),
    c(
      coords_beads[, CHANNEL_BEADS],
      coords_beads[1, CHANNEL_BEADS]
    ),
    col = "red",
    lwd = 2
  )
  
  # Panel 3: diagnóstico FITC-A x BV711-A
  plot(
    beads_scatter_plot[, CHANNEL_SYTO],
    beads_scatter_plot[, CHANNEL_BEADS],
    pch = ".",
    cex = 0.45,
    col = grDevices::adjustcolor(
      "grey50",
      alpha.f = 0.35
    ),
    xlab = CHANNEL_SYTO,
    ylab = CHANNEL_BEADS,
    main = "3. Diagnóstico: SYTO x BV711"
  )
  
  if (nrow(beads_plot) > 0) {
    
    points(
      beads_plot[, CHANNEL_SYTO],
      beads_plot[, CHANNEL_BEADS],
      pch = ".",
      cex = 0.75,
      col = grDevices::adjustcolor(
        "#D55E00",
        alpha.f = 0.85
      )
    )
  }
  
  # Panel 4: histograma BV711
  if (nrow(beads_mat) > 0) {
    
    hist(
      beads_mat[, CHANNEL_BEADS],
      breaks = 80,
      col = "#CC79A7",
      border = "white",
      main = "4. BV711 en Beads",
      xlab = CHANNEL_BEADS,
      ylab = "Eventos"
    )
    
  } else {
    
    plot_empty_panel(
      "4. BV711 en Beads",
      "No hay eventos bead"
    )
  }
  
  mtext(
    side = 3,
    outer = TRUE,
    line = 1,
    font = 2,
    cex = 1.1,
    text = paste0(
      "QC beads: ",
      sample_name,
      " | N final = ",
      nrow(beads_mat)
    )
  )
}


beads_qc_pdf <- file.path(
  output_dir,
  "QC_gating_beads_por_muestra.pdf"
)

if (bead_gate_available) {
  
  pdf(
    file = beads_qc_pdf,
    width = 12,
    height = 10,
    onefile = TRUE
  )
  
  for (sample_name in sampleNames(gs)) {
    
    message(
      "Generando QC beads para: ",
      sample_name
    )
    
    tryCatch(
      
      {
        plot_bead_qc_page(
          gs = gs,
          sample_name = sample_name,
          coords_beads_scatter =
            coords_beads_scatter,
          coords_beads =
            coords_beads
        )
      },
      
      error = function(e) {
        
        par(
          mfrow = c(1, 1),
          mar = c(4, 4, 4, 2)
        )
        
        plot.new()
        
        title(
          main = paste(
            "Error QC beads:",
            sample_name
          )
        )
        
        text(
          x = 0.5,
          y = 0.5,
          labels = conditionMessage(e)
        )
      }
    )
  }
  
  dev.off()
  
  message("")
  message(
    "PDF de QC de beads creado:"
  )
  message(beads_qc_pdf)
}


# ================================================================
# 18. GUARDAR OBJETOS FINALES
# ================================================================

saveRDS(
  gs,
  file.path(
    output_dir,
    "GatingSet_final_IgA_bacterias.rds"
  )
)

analysis_objects <- list(
  channels = list(
    FSC_A = CHANNEL_FSC_A,
    FSC_H = CHANNEL_FSC_H,
    SSC_A = CHANNEL_SSC_A,
    SYTO = CHANNEL_SYTO,
    IGA = CHANNEL_IGA,
    BEADS = CHANNEL_BEADS
  ),
  parameters = list(
    IGA_NEGATIVE_QUANTILE =
      IGA_NEGATIVE_QUANTILE,
    MIN_BACTERIA_EVENTS_FOR_APC_THRESHOLD =
      MIN_BACTERIA_EVENTS_FOR_APC_THRESHOLD,
    MIN_BACTERIA_EVENTS_FOR_ANALYSIS =
      MIN_BACTERIA_EVENTS_FOR_ANALYSIS,
    BEADS_ADDED_PER_TUBE =
      BEADS_ADDED_PER_TUBE,
    SAMPLE_VOLUME_UL =
      SAMPLE_VOLUME_UL
  ),
  bacteria_polygon = coords_bact,
  bead_scatter_polygon = coords_beads_scatter,
  bead_polygon = coords_beads,
  thresholds = umbral_df,
  gating_summary = stats_df,
  pre_vs_fr_summary = summary_wide
)

saveRDS(
  analysis_objects,
  file.path(
    output_dir,
    "analysis_objects_IgA_bacterias.rds"
  )
)


# ================================================================
# 19. INFORME FINAL
# ================================================================

message("")
message("================================================")
message("ANALISIS COMPLETADO")
message("================================================")
message("")
message("Archivos principales:")
message(" - metadata_checked.csv")
message(" - excluded_files_unrecognized_type.csv")
message(" - incomplete_triplets.csv")
message(" - iga_thresholds_by_triplet.csv")
message(" - gating_summary.csv")
message(" - summary_PRE_vs_FR.csv")
message(" - GatingSet_final_IgA_bacterias.rds")
message(" - analysis_objects_IgA_bacterias.rds")
message(" - QC_gating_completo_por_muestra.pdf")
message(" - QC_gating_beads_por_muestra.pdf")
message(" - carpeta QC_plots/")
message("")
message(
  "Número de tríadas completas: ",
  nrow(pairs_complete)
)
message(
  "Número de FCS analizados: ",
  nrow(stats_df)
)
message("")
message(
  "Resumen de resultados:"
)

print(
  stats_df[
    ,
    .(
      sample_id,
      tipo,
      bacterias_n,
      iga_pos_n,
      iga_pos_pct_in_bact,
      beads_scatter_n,
      beads_n,
      bacterias_per_tube,
      iga_pos_per_tube,
      bacterias_per_uL,
      iga_pos_per_uL,
      qc_comment
    )
  ]
)