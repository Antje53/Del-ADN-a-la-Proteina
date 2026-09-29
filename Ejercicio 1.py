# Ejercicio 1. Replicación del ADN
from Bio.Seq import Seq

# ADN original
original = Seq("ATGCCGTTAGCT")
print(f"Secuencia original (5' -> 3'): {original}")

# hebra complementaria
complementaria = original.complement()
print(f"Hebra complementaria (3' -> 5'): {complementaria}")

