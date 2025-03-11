# Install packages
#install.packages('meteospain')

#Empty workspace
rm(list = ls())

#Import packages
library(ggplot2)
library(maps)
library(dplyr)
library(meteospain)
library(leaflet)

# Define la caja de coordenadas para Barcelona
min_lat <- 41.30  # latitud mínima
max_lat <- 41.49  # latitud máxima
min_lon <-  2.02  # longitud mínima
max_lon <- 2.26   # longitud máxima

# Función para extraer el año del nombre del archivo
# Se asume que el nombre del archivo tiene la forma: pr_CMCC-CM2-SR5_r1i1p1f1_ssp245_20150101-20151231.dat
get_year <- function(filename) {
  bn <- basename(filename)
  bn <- sub("\\.dat$", "", bn)
  # Separa por "_" y extrae el último elemento (la parte de las fechas)
  parts <- strsplit(bn, "_")[[1]]
  date_part <- parts[length(parts)]
  # El año se corresponde con los 4 primeros dígitos de la parte de la fecha
  year_val <- as.numeric(substr(date_part, 1, 4))
  return(year_val)
}

# Función para leer y filtrar un archivo individual
read_and_filter <- function(file) {
  # Lee todas las líneas del archivo
  lines <- readLines(file)
  
  # Extrae las líneas de metadata (se asume que comienzan con "# id;", "# lat;" y "# lon;")
  id_line  <- lines[grep("^#\\s*id;", lines)][1]
  lat_line <- lines[grep("^#\\s*lat;", lines)][1]
  lon_line <- lines[grep("^#\\s*lon;", lines)][1]
  
  # Elimina el "#" y espacios, luego separa por ";"
  id_values  <- unlist(strsplit(sub("^#\\s*", "", id_line), ";"))
  lat_values <- unlist(strsplit(sub("^#\\s*", "", lat_line), ";"))
  lon_values <- unlist(strsplit(sub("^#\\s*", "", lon_line), ";"))
  
  # Quita el primer elemento (la etiqueta)
  id_values  <- id_values[-1]
  lat_values <- as.numeric(lat_values[-1])
  lon_values <- as.numeric(lon_values[-1])
  
  # Determina los índices de los puntos de rejilla que están dentro de Barcelona
  catalonia_idx <- which(lat_values >= min_lat & lat_values <= max_lat &
                           lon_values >= min_lon & lon_values <= max_lon)
  
  # Identifica la primera línea de datos (la primera que no empieza con "#")
  data_start <- which(!grepl("^#", lines))[1]
  
  # Lee la parte de datos (separados por ";" y sin cabecera)
  data_df <- read.table(text = lines[data_start:length(lines)],
                        sep = ";", header = FALSE, stringsAsFactors = FALSE)
  
  # En la tabla, la primera columna es la fecha y el resto corresponden a los puntos de la rejilla.
  # Se conservan la primera columna (fecha) y las columnas de Barcelona.
  cols_to_keep <- c(1, catalonia_idx + 1)
  data_filtered <- data_df[, cols_to_keep, drop = FALSE]
  
  # Opcional: asigna nombres a las columnas ("Fecha" y luego "lat_lon")
  grid_names <- paste0(lat_values[catalonia_idx], "_", lon_values[catalonia_idx])
  colnames(data_filtered) <- c("Fecha", grid_names)
  
  # Convierte la columna de fecha (asumida en formato YYYYMMDD) a tipo Date
  data_filtered$Fecha <- as.Date(as.character(data_filtered$Fecha), format = "%Y%m%d")
  
  return(data_filtered)
}


for (variable in c("Tmax", "Tmin", "Rainfall")){
  for (model in c("CMCC-CM2-SR5", "CNRM-ESM2-1", "EC-Earth3-Veg", "MIROC6", "MPI-ESM1-2-HR", "NorESM2-MM")) {
    
    # Define la carpeta que contiene los archivos .dat
    folder_path <- paste0("/home/jesus/Escritorio/Projecte_BSC_Lluvias/Dades_crues/Rejilla/Historical/", variable, "/", model)
    
    print(folder_path)
    # Lista todos los archivos .dat en la carpeta
    file_list <- list.files(path = folder_path, pattern = "\\.dat$", full.names = TRUE)
    
    # Selecciona únicamente los archivos correspondientes a años entre 2030 y 2060
    selected_files <- file_list[sapply(file_list, function(x) {
      year_val <- get_year(x)
      return(year_val >= 1980 & year_val <= 2010)
    })]
    
    # Lista para almacenar los datos filtrados de cada archivo
    data_list <- list()
    
    # Procesa cada archivo seleccionado
    for (file in selected_files) {
      cat("Procesando archivo:", file, "\n")
      data_temp <- read_and_filter(file)
      data_list[[length(data_list) + 1]] <- data_temp
    }
    
    # Combina todos los datos (se asume que las columnas coinciden entre archivos)
    combined_data <- do.call(rbind, data_list)
    
    # (Opcional) Ordena los datos combinados por fecha
    combined_data <- combined_data[order(combined_data$Fecha), ]
    
    # Define la ruta y nombre del archivo de salida en formato CSV
    output_file <- paste0("/home/jesus/Escritorio/Projecte_BSC_Lluvias/Dades_netes/Historical/", variable, "_Historical_81_10_", model, ".csv")
    
    # Guarda el data frame combinado en un archivo CSV
    write.csv(combined_data, file = output_file, row.names = FALSE)
    
    cat("Datos combinados guardados en:", output_file, "\n")
  }
}





















plot_leaflet_proyecciones <- function(data_filtered) {
  # Extract grid column names (ignoring the "date" column)
  grid_names <- colnames(data_filtered)[-1]
  
  #Split the string to extract lat and lon
  coords_list <- strsplit(grid_names_clean, "_")
  
  # Combine the list into a data frame with columns 'lat' and 'lon'
  coords <- do.call(rbind, coords_list)
  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  colnames(coords) <- c("lat", "lon")
  
  # Convert the coordinate strings to numeric values
  coords$lat <- as.numeric(coords$lat)
  coords$lon <- as.numeric(coords$lon)
  
  # Create the leaflet map centered on Barcelona
  mapa_barcelona <- leaflet() %>%
    addProviderTiles(providers$CartoDB.Positron) %>%
    setView(lng = 2.123006, lat = 41.393496, zoom = 17.45) %>%
    addScaleBar(position = "bottomleft", options = scaleBarOptions(imperial = FALSE))
  
  # Superpose the projection points as circle markers
  mapa_barcelona <- mapa_barcelona %>%
    addCircleMarkers(
      data = coords,
      lng = ~lon,
      lat = ~lat,
      radius = 3,
      color = "blue",
      fillOpacity = 0.8,
      stroke = FALSE,
      popup = ~paste("Lat:", lat, "<br>Lon:", lon)
    )
  
  return(mapa_barcelona)
}

# Ejemplo de uso:
# Suponiendo que 'data_filtered' ya contiene la información procesada,
# ejecuta la función y visualiza el mapa interactivo:
mapa_final <- plot_leaflet_proyecciones(data_filtered)
mapa_final


# 
# # LOAD ANDPROCESS DATAFRAME
# script_dir <- "/home/jesus/Escritorio/Projecte_Barcelona/Dades_crues"
# Meteo_Cat <- read.xlsx(file.path(script_dir, "Meteo_Cat_2019_2024.xlsx")) #Read data
# 
# #Select only stations from Barcelona
# Meteo_BCN <- Meteo_Cat[grepl("Barcelona", Meteo_Cat$station_name), ]
# 
# Meteo_BCN <- Meteo_BCN %>%
#   select(timestamp, station_name, mean_temperature, min_temperature, max_temperature, precipitation, long_station, lat_station) %>% # Select relevant columns
#   rename(
#     Fecha = timestamp,
#     station = station_name,
#     Mean_T = mean_temperature,
#     Min_T = min_temperature,
#     Max_T = max_temperature,
#     Rainfall = precipitation,
#     longitude = long_station,
#     latitude = lat_station
#   )
# 
# Meteo_BCN <- Meteo_BCN %>%
#   group_by(station) %>% # Group by station
#   mutate(
#     Rainfall_21 = rollapply(
#       data = Rainfall,
#       width = 21,
#       FUN = sum,
#       align = "right",
#       fill = NA
#     )
#   ) %>%
#   ungroup() # Remove grouping after computation
# 
# # Filter dates before 2019-03-01
# Meteo_BCN <- Meteo_BCN %>%
#   filter(as.Date(Fecha) >= as.Date("2019-03-01"))
# 
# #Write dataframe
# write.xlsx(Meteo_BCN, "/home/jesus/Escritorio/Projecte_Barcelona/BCN-Model/Simulations_BCN/Dades_netes/Meteo_BCN_2019_2024.xlsx", rowNames = FALSE)
