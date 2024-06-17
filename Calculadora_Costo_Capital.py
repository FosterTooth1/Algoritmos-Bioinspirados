import numpy as np
import csv

# Permutación de los índices de los puntos
permutation = [7, 13, 24, 20, 2, 3, 25, 16, 19, 12, 10, 31, 17, 15, 21, 1, 32, 28, 30, 22, 4, 11, 8, 9, 14, 27, 6, 23, 26, 5, 29, 18]

# Nombre del archivo CSV
archivo_csv = 'distancias.csv'

# Inicializar el costo total
total_cost = 0.0

nombres_ciudades = []

# Abrir el archivo CSV en modo lectura
with open(archivo_csv, newline='', encoding='utf-8') as csvfile:
    # Crear un lector CSV
    csv_reader = csv.reader(csvfile)
    
    # Leer las filas del archivo CSV
    filas_csv = list(csv_reader)
    
    # Obtener el encabezado superior y el encabezado lateral izquierdo
    nombres_ciudades = filas_csv[0][1:]  # Ignoramos la primera celda (es el encabezado lateral)
    
    # Leer las distancias del archivo CSV
    for i in range(len(permutation) - 1):
        start = permutation[i] - 1  # Ajuste de índice
        end = permutation[i + 1] - 1  # Ajuste de índice
        total_cost += float(filas_csv[start + 1][end + 1])  # Sumamos la distancia
        
    # Añadir el costo de regreso al inicio
    start = permutation[-1] - 1  # Último punto en la permutación
    end = permutation[0] - 1  # Primer punto en la permutación
    total_cost += float(filas_csv[start + 1][end + 1])  # Sumamos la distancia

# Ajustar la permutación para imprimir correctamente las ciudades
    permutation_print = [x - 1 for x in permutation]
    
    print("Orden de ciudades recorridas")
    for i in range(len(permutation_print)):
        print(f"{i + 1}. {nombres_ciudades[permutation_print[i]]}")

# Imprimir el costo total
print(f"El costo total de la permutación es: {total_cost}")