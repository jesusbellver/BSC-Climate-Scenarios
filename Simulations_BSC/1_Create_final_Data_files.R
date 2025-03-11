# Install packages
#install.packages('meteospain')

#Empty workspace
rm(list = ls())

# Cargar las librerías necesarias
library(dplyr)
library(tidyr)
library(zoo)
library(readr)  # Para parse_double()

# Definir el path base donde se encuentran los archivos
base_path <- "/home/jesus/Escritorio/Projecte_BSC_Lluvias/Dades_netes/SSP585"

# Lista de modelos a procesar
modelos <- c("CMCC-CM2-SR5", "CNRM-ESM2-1", "EC-Earth3-Veg", "MIROC6", "MPI-ESM1-2-HR", "NorESM2-MM")



procesar_archivo <- function(archivo, nombre_variable) {
  df <- read.csv(archivo, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Convertir la columna Fecha a Date (ajustar el formato si es necesario)
  df$Fecha <- as.Date(df$Fecha)
  
  # Transformar de formato ancho a largo: se asume que las columnas (salvo Fecha) tienen nombres "lat_long"
  df_long <- df %>%
    pivot_longer(-Fecha, names_to = "coord", values_to = nombre_variable) %>%
    separate(coord, into = c("latitude", "longitude"), sep = "_") %>%
    mutate(latitude = parse_double(latitude, locale = locale(decimal_mark = ".")),
           longitude = parse_double(longitude, locale = locale(decimal_mark = ".")))
  
  return(df_long)
}

# Procesar cada modelo
for (modelo in modelos) {
  # Construir las rutas de los archivos para cada variable
  archivo_Tmax    <- file.path(base_path, paste0("Tmax_SSP585_31_60_", modelo, ".csv"))
  archivo_Tmin    <- file.path(base_path, paste0("Tmin_SSP585_31_60_", modelo, ".csv"))
  archivo_Rainfall<- file.path(base_path, paste0("Rainfall_SSP585_31_60_", modelo, ".csv"))
  
  # Leer y transformar cada archivo
  df_Tmax <- procesar_archivo(archivo_Tmax, "Tmax")
  df_Tmin <- procesar_archivo(archivo_Tmin, "Tmin")
  df_Rain <- procesar_archivo(archivo_Rainfall, "Rainfall")
  
  # Unir los data frames por Fecha, latitude y longitude
  df_merge <- df_Tmax %>%
    left_join(df_Tmin, by = c("Fecha", "latitude", "longitude")) %>%
    left_join(df_Rain,  by = c("Fecha", "latitude", "longitude"))
  
  # Calcular la media de temperaturas y asignar Tmin y Tmax a Min_T y Max_T
  df_merge <- df_merge %>%
    mutate(Mean_T = (Tmax + Tmin) / 2,
           Min_T = Tmin,
           Max_T = Tmax)
  
  # Calcular la suma acumulada de 21 días para la precipitación para cada ubicación
  # Se agrupa por cada par de coordenadas y se ordena por Fecha
  df_merge <- df_merge %>%
    arrange(latitude, longitude, Fecha) %>%
    group_by(latitude, longitude) %>%
    mutate(Rainfall_21 = rollapply(Rainfall, width = 21, FUN = sum, align = "right", fill = 0)) %>%
    ungroup()
  
  # Seleccionar y reordenar las columnas deseadas
  df_final <- df_merge %>%
    select(Fecha, longitude, latitude, Mean_T, Min_T, Max_T, Rainfall, Rainfall_21)
  
  # Guardar el resultado en un CSV para el modelo actual
  salida <- file.path(base_path, paste0("Meteo_SSP585_31_60_", modelo, ".csv"))
  write.csv(df_final, salida, row.names = FALSE)
  
  cat("Archivo procesado para el modelo:", modelo, "\n")
}
