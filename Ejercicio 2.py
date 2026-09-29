import sys
from pathlib import Path
from Bio import SeqIO

RUTA_POR_DEFECTO = Path(__file__).parent / "ejercicio2.fasta"


def main():
    ruta = Path(sys.argv[1]) if len(sys.argv) > 1 else RUTA_POR_DEFECTO
    if not ruta.exists():
        sys.exit(f"[!] No se encuentra el archivo FASTA: {ruta}")

    for registro in SeqIO.parse(ruta, "fasta"):
        codificante = registro.seq.upper()            # 5'->3'
        molde_3_5 = codificante.complement()          
        molde_5_3 = codificante.reverse_complement() 

        print(f"Registro: {registro.id}  ({len(codificante)} nt)")
        print(f"  Hebra codificante  5' {codificante} 3'")
        print(f"  Hebra molde        3' {molde_3_5} 5'")

        
        arnm = molde_3_5.complement_rna()
        assert arnm == codificante.transcribe()  # codificante con T->U
        print(f"\n[CORRECTO] ARNm    5' {arnm} 3'  -> proteína: {arnm.translate()}")

        # cambiar la orientación de la hebra
        print("\nExperimento: ¿qué pasa si se usa mal la orientación?")
        error_a = molde_3_5.transcribe()          # molde tratado como si fuese codificante
        error_b = molde_5_3.transcribe()          # molde 5'->3' tratado como codificante
        print(f"  a) Transcribir el molde tal cual (3'->5'):  {error_a}"
              "  -> no es un ARNm: es el complementario del ARNm y está al revés")
        print(f"  b) Molde escrito 5'->3' como si fuera codificante: 5' {error_b} 3'"
              f"  -> proteína: {error_b.translate()}")
        print(f"  c) Forma correcta desde el molde 5'->3' (reverse_complement): "
              f"5' {molde_5_3.reverse_complement_rna()} 3'")



if __name__ == "__main__":
    main()
