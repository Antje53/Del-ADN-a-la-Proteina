from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParam import ProteinAnalysis

KYTE_DOOLITTLE = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
                  "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
                  "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
                  "Y": -1.3, "V": 4.2}

PEPTIDO = "Met-Ile-Ser-Gly-Val-Lys-His"


def tipo(aa):
    v = KYTE_DOOLITTLE[aa]
    return "hidrofóbico" if v > 0 else "hidrofílico"


def main():
    residuos = PEPTIDO.split("-")
    seq = "".join(seq1(r) for r in residuos)
    print(f"Péptido: {PEPTIDO}  ({seq})")
    print(f"  Extremo N:  {residuos[0]}")
    print(f"  Extremo C:  {residuos[-1]}\n")
    print(f"{'Pos':<5}{'aa':<6}{'Kyte-Doolittle':>15}  Tipo")
    for i, (r, a) in enumerate(zip(residuos, seq), 1):
        print(f"{i:<5}{r:<6}{KYTE_DOOLITTLE[a]:>15}  {tipo(a)}")

    gravy = ProteinAnalysis(seq).gravy()
    mutado = seq.replace("V", "K", 1)  # Val (hidrofóbica) -> Lys (hidrofílica)
    gravy_mut = ProteinAnalysis(mutado).gravy()
    print(f"\nHidrofobicidad media (GRAVY) original: {gravy:.2f}")
    print(f"Mutación Val5->Lys ({mutado}): GRAVY = {gravy_mut:.2f}")



if __name__ == "__main__":
    main()
