#!/bin/bash

# Iterar sobre todos los archivos .gb en el directorio actual
for file in *.gb; do
  # Extraer el nombre de la especie, número y la etiqueta del archivo
  filename=$(basename "$file" .gb)
  
  # Separar los componentes del nombre del archivo
  genus=$(echo "$filename" | cut -d'_' -f1)
  specie=$(echo "$filename" | cut -d'_' -f2)
  label=$(echo "$filename" | cut -d'_' -f3,4)

  # Leer el contenido del archivo y reemplazar la línea DEFINITION
  # La etiqueta contiene el formato correcto del archivo (p.ej., 22_Melli.029)
  sed -i "/^DEFINITION/c\DEFINITION  $genus $specie $label" "$file"

  echo "Procesado: $file"
done
