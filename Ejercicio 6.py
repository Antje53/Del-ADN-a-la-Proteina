import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

BASE = Path(__file__).parent
RUTA_POR_DEFECTO = BASE /"lacZ_Ecoli_NC_000913.3.fasta"
N = 30  # nucleótidos 


def info(msg):
    print(f"  -> {msg}")


def corto(seq, n=N):
    return f"{seq[:n]}..." if len(seq) > n else str(seq)


def leer_fasta(ruta):
    print("\n LECTURA DE LA SECUENCIA")
    registro = SeqIO.read(ruta, "fasta")
    seq = registro.seq.upper()
    info(f"Archivo: {ruta.name}")
    info(f"Registro: {registro.description}")
    info(f"Longitud: {len(seq)} nt | Contenido GC: {gc_fraction(seq) * 100:.1f} %")
    invalidos = set(str(seq)) - set("ACGT")
    if invalidos:
        raise ValueError(f"La secuencia contiene caracteres no válidos: {invalidos}")
    return registro, seq


def replicacion(codificante):
    print("\n REPLICACIÓN")
    molde = codificante.complement()  # hebra parental opuesta, 3'->5'
    info(f"  Parental A (codificante) 5' {corto(codificante)} 3'")
    info(f"  Parental B (molde)       3' {corto(molde)} 5'")
    nueva_a = codificante.complement()   # 3'->5', complementaria de A
    nueva_b = molde.complement()         # 5'->3', complementaria de B
    info(f"  Nueva sobre A            3' {corto(nueva_a)} 5'")
    info(f"  Nueva sobre B            5' {corto(nueva_b)} 3'")
    assert nueva_a == molde and nueva_b == codificante
    return molde


def transcripcion(molde, codificante):
    print("\n TRANSCRIPCIÓN")
    info("La ARN polimerasa lee la hebra MOLDE en sentido 3'->5' ...")
    arnm = molde.complement_rna()  # ARN complementario 
    info(f"... y sintetiza el ARNm 5'->3' (U en lugar de T): 5' {corto(arnm)} 3'")
    assert arnm == codificante.transcribe()
    info(f"Longitud del ARNm: {len(arnm)} nt. Comprobación: igual a la hebra "
         "codificante con T->U. [verificado]")
    return arnm


def traduccion(arnm):
    print("\n[PASO 3] TRADUCCIÓN")
    inicio = arnm[:3]
    info(f"Codón de inicio: {inicio} ({'Met, correcto' if inicio == 'AUG' else 'NO es AUG'})")
    if len(arnm) % 3:
        info(f"Aviso: la longitud ({len(arnm)}) no es múltiplo de 3; se recorta.")
        arnm = arnm[: len(arnm) - len(arnm) % 3]
    info(f"El ribosoma lee {len(arnm) // 3} codones de 3 nucleótidos...")
    proteina = arnm.translate(to_stop=True)
    pos_paro = len(proteina) * 3
    paro = arnm[pos_paro:pos_paro + 3]
    info(f"Codón de paro {paro} encontrado en el codón nº {len(proteina) + 1} "
         f"(nucleótidos {pos_paro + 1}-{pos_paro + 3})")
    try:
        arnm.translate(cds=True)
        info("Validación como CDS completa (inicio, sin paros internos, paro final): OK")
    except Exception as e:  # TranslationError
        info(f"Validación como CDS completa: NO superada ({e})")
    info(f"Proteína: {len(proteina)} aminoácidos")
    info(f"  Extremo N: {proteina[:10]}...   Extremo C: ...{proteina[-10:]}")
    return proteina


def guardar(registro, arnm, proteina):
    carpeta = BASE / "resultados"
    carpeta.mkdir(exist_ok=True)
    SeqIO.write(SeqRecord(arnm, id=f"{registro.id}_ARNm", description=""),
                carpeta / "ARNm.fasta", "fasta")
    SeqIO.write(SeqRecord(proteina, id=f"{registro.id}_proteina", description=""),
                carpeta / "proteina.fasta", "fasta")
    print(f"\nResultados guardados en {carpeta.name}/ARNm.fasta y {carpeta.name}/proteina.fasta")


def dogma_central_pipeline(ruta_fasta):
    ruta = Path(ruta_fasta)
    if not ruta.exists():
        sys.exit(f" No se encontró el archivo '{ruta}'.")
    try:
        registro, codificante = leer_fasta(ruta)
    except ValueError as e:
        sys.exit(f" Archivo FASTA no válido: {e}")
    molde = replicacion(codificante)
    arnm = transcripcion(molde, codificante)
    proteina = traduccion(arnm)
    guardar(registro, arnm, proteina)
    print("\n ADN -> ADN (replicación) -> ARNm -> proteína")
    return proteina


if __name__ == "__main__":
    dogma_central_pipeline(sys.argv[1] if len(sys.argv) > 1 else RUTA_POR_DEFECTO)
