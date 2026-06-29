#!/bin/bash
# Script modificado para compilar el modelo MITgcm
# con opción de usar MPI

# --- Compilación en build ---
echo "Accediendo al directorio 'build'..."
cd ../build/ || { echo "No se pudo acceder al directorio 'build'."; exit 1; }

# Preguntar si se desea limpiar la carpeta build antes de generar el makefile
read -p "¿Desea limpiar la carpeta build antes de generar el makefile? [y/n]: " respuesta_clean
if [[ "$respuesta_clean" == [yY] ]]; then
    echo "Ejecutando make clean..."
    make clean
    echo "Eliminando archivos con rm * ..."
    rm -rf *
fi

GENMAKE2_PATH="../../../tools/genmake2"
BUILD_OPTIONS_PATH="../../../tools/build_options/darwin_arm64_gfortran.opt"

# Preguntar si se desea compilar con MPI
read -p "¿Desea compilar con MPI? [y/n]: " comp_mpi
if [[ "$comp_mpi" == [yY] ]]; then
    echo "Generando makefile con MPI..."
    $GENMAKE2_PATH -mods ../code/ -optfile $BUILD_OPTIONS_PATH -mpi
else
    echo "Generando makefile sin MPI..."
    $GENMAKE2_PATH -mods ../code/ -optfile $BUILD_OPTIONS_PATH
fi

if [ $? -ne 0 ]; then
    echo "Error al generar el makefile."
    exit 1
fi

echo "Ejecutando 'make depend'..."
make depend
if [ $? -ne 0 ]; then
    echo "Error en 'make depend'."
    exit 1
fi

echo "Ejecutando 'make' para compilar el modelo..."

module load intelMPI2022
module load toolsCDFintel2022

make
if [ $? -ne 0 ]; then
    echo "Error durante la compilación."
    exit 1
fi

cd ..

