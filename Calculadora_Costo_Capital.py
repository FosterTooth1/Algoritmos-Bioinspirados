import numpy as np
import csv

# Permutación de los índices de los puntos
permutation = [31, 13, 8, 17, 9, 24, 3, 2, 25, 6, 7, 18, 27, 30, 22, 4, 26, 5, 19, 29, 28, 20, 12, 32, 16, 11, 14, 15, 21, 10, 23, 1]

# Nombre del archivo CSV
archivo_csv = 'Algoritmos-Bioinspirados\distancias.csv'

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