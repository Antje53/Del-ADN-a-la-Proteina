from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

# archivo FASTA de prueba en el directorio actual
secuencia_adn = Seq("ATGCCTGAATGC")
registro_fasta = SeqRecord(secuencia_adn, id="Gen_Prueba", description="Hebra codificante")

with open("prueba_ejercicio2.fasta", "w") as f_out:
    SeqIO.write(registro_fasta, f_out, "fasta")

# Leemos el archivo FASTA y transcribimos
for registro in SeqIO.parse("prueba_ejercicio2.fasta", "fasta"):
    print(f"ID del registro: {registro.id}")
    print(f"Secuencia ADN (Codificante 5'->3'): {registro.seq}")
    
    # Transcripción directa
    arnm = registro.seq.transcribe()
    print(f"Transcrito ARNm (5'->3'): {arnm}")
    
    # Cambiando la orientación de la hebra
    hebra_molde = registro.seq.complement()
    print(f"\nHebra molde (3'->5'): {hebra_molde}")
    
    # transcribir la hebra molde directamente
    arnm_invertido = hebra_molde.transcribe()
    print(f"ARNm a partir de la hebra molde invertida: {arnm_invertido}")

# Limpiamos el archivo temporal 
os.remove("prueba_ejercicio2.fasta")