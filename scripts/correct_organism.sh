#!/bin/bash

# Iterar sobre todos los archivos .gb en el directorio actual
for file in *.gb; do
  # Extraer el nombre de la especie y la etiqueta del nombre del archivo
  filename=$(basename "$file" .gb)
  
  # Separar los componentes del nombre del archivo
  genus=$(echo "$filename" | cut -d'_' -f1)
  species=$(echo "$filename" | cut -d'_' -f2)
  label=$(echo "$filename" | cut -d'_' -f3,4)

  # Crear la nueva línea que usaremos para reemplazar
  new_organism="$genus $species $label"
  
  sed -i '/ORGANISM/{N;s/\n//;}' $file # Put the ORGANISM field on a single line.

  # Reemplazar la línea ORGANISM
  sed -i "/^  ORGANISM/c\  ORGANISM  $new_organism" "$file"

  echo "Procesado: $file"
done
