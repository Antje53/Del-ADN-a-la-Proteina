import sys
from collections import Counter
import requests

SERVIDOR = "https://rest.ensembl.org"

def simular_splicing():
    exones = {1: "Exón 1", 2: "Exón 2", 3: "Exón 3", 4: "Exón 4", 5: "Exón 5"}
    for combinacion in [(1, 2, 3, 4, 5), (1, 3, 4, 5), (1, 4, 5)]:
        saltados = [e for e in exones if e not in combinacion]
        nombre = "-".join(map(str, combinacion))
        print(f"  ARNm {nombre:<10} exones omitidos: {saltados or 'ninguno'}")
    print()


def consultar_gen(simbolo):
    url = f"{SERVIDOR}/lookup/symbol/homo_sapiens/{simbolo}?expand=1"
    r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=30)
    r.raise_for_status()
    return r.json()


def main():
    simbolo = sys.argv[1] if len(sys.argv) > 1 else "BRCA2"
    comparar_con = sys.argv[2] if len(sys.argv) > 2 else "ENST00000680887"  # Isoforma BRCA2-210

    simular_splicing()
    try:
        gen = consultar_gen(simbolo)
    except requests.RequestException as e:
        sys.exit(f"[!] No se pudo conectar con Ensembl: {e}")

    transcritos = gen.get("Transcript", [])
    print(f"Gen {gen['display_name']} ({gen['id']}), cromosoma {gen['seq_region_name']}, "
          f"{gen.get('assembly_name', '')} ===")
    print(f"Transcritos anotados: {len(transcritos)}")
    print("Por tipo (biotype):")
    for tipo, n in Counter(t["biotype"] for t in transcritos).most_common():
        print(f"  {tipo:<32} {n}")

    codificantes = [t for t in transcritos if t["biotype"] == "protein_coding"]
    # Ordenamos por longitud de la proteína
    codificantes.sort(key=lambda t: -(t.get("Translation") or {}).get("length", 0))
    
    print(f"{'Transcrito':<17}{'Nombre':<12}{'Exones':>7}{'ARN (nt)':>10}{'Proteína (aa)':>15}")
    for t in codificantes:
        aa = (t.get("Translation") or {}).get("length", "-")
        marca = "  <- canónico" if t.get("is_canonical") else ""
        print(f"{t['id']:<17}{t.get('display_name', ''):<12}{len(t.get('Exon', [])):>7}"
              f"{t.get('length', '-'):>10}{aa:>15}{marca}")

    # Comparación exón a exón: canónico vs transcrito elegido
    canonico = next((t for t in transcritos if t.get("is_canonical")), None)
    otro = next((t for t in transcritos if t["id"] == comparar_con), None)
    
    if canonico and otro:
        ex_c = [e["id"] for e in canonico.get("Exon", [])]
        ex_o = [e["id"] for e in otro.get("Exon", [])]
        print(f"\nComparación {canonico['id']} (canónico) vs {otro['id']}:")
        
        for e in ex_c:
            if e not in ex_o:
                print(f"  Exón solo en el canónico: {e} (posición {ex_c.index(e) + 1})")
        for e in ex_o:
            if e not in ex_c:
                print(f"  Exón solo en {otro['id']}: {e} (posición {ex_o.index(e) + 1})")
                
        aa_c = (canonico.get("Translation") or {}).get("length", "-")
        aa_o = (otro.get("Translation") or {}).get("length", "-")
        print(f"  Proteína: {aa_c} aa vs {aa_o} aa")
    else:
        print(f"\n No se pudo realizar la comparación: Transcrito {comparar_con} no encontrado o falta el canónico.")


if __name__ == "__main__":
    main()