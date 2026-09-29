from Bio import SeqIO
from Bio.Seq import Seq

def dogma_central_pipeline(ruta_fasta):

    
    try:
        # Extraemos la secuencia del archivo FASTA 
        registro = SeqIO.read(ruta_fasta, "fasta")
        adn_codificante = registro.seq
        print(f"-> ID Secuencia: {registro.id}")
        print(f"-> Longitud total: {len(adn_codificante)} pares de bases")
    except FileNotFoundError:
        print(f"\n[!] Error crítico: No se encontró el archivo '{ruta_fasta}'.")
        return
    except Exception as e:
        print(f"\n[!] Error inesperado al procesar el archivo: {e}")
        return
        
    # Replicación
    print("\n[PASO 1] REPLICACIÓN:")
    # La hebra molde es la complementaria a la codificante
    adn_molde = adn_codificante.complement()
    # Mostramos solo los primeros 30 
    print(f"  -> Hebra codificante original (5'->3') [primeras 30 pb]: {adn_codificante[:30]}...")
    print(f"  -> Hebra molde sintetizada (3'->5')    [primeras 30 pb]: {adn_molde[:30]}...")
    
    # Transcripción
    print("\n[PASO 2] TRANSCRIPCIÓN:")
    # Se transcribe asumiendo que adn_codificante es la hebra codificante (T -> U)
    arnm = adn_codificante.transcribe() 
    print("  -> Sintetizando ARNm desde la hebra molde en la ARN polimerasa...")
    print(f"  -> ARNm resultante (5'->3')           [primeras 30 pb]: {arnm[:30]}...")
    
    # Traducción
    print("\n[PASO 3] TRADUCCIÓN:")
    # to_stop=True hace que el ribosoma se detenga al leer el primer codón de paro 
    proteina = arnm.translate(to_stop=True)
    print("  -> Leyendo tripletes (codones) en el ribosoma...")
    print(f"  -> Cadena de aminoácidos resultante  [primeros 10 aa]: {proteina[:10]}...")
    
    print("\nPROCESO COMPLETADO CON ÉXITO.")

if __name__ == "__main__":
    archivo_descargado = r"C:\Users\leono\OneDrive\Desktop\UNI\Cuarto\BIO\ADN a la Proteina\lacZ_Ecoli_NC_000913.3.fasta"
    dogma_central_pipeline(archivo_descargado)