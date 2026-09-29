from Bio.Seq import Seq

ARNM = Seq("AUGUAUGCUUAA")
MANUAL = "MYA"  # Met-Tyr-Ala


def codones(seq):
    return " ".join(str(seq[i:i + 3]) for i in range(0, len(seq) - len(seq) % 3, 3))


def main():
    print(f"ARNm original: 5' {codones(ARNM)} 3'")
    print(f"  Traducción completa:        {ARNM.translate()} ")
    proteina = ARNM.translate(to_stop=True)
    print(f"  Traducción hasta el paro:   {proteina}")
    assert str(proteina) == MANUAL
    print("  [OK] Coincide con la traducción manual Met-Tyr-Ala.\n")

    # AUG -> GUG en el codón de inicio
    mut_inicio = Seq("G") + ARNM[1:]
    print(f"Mutación AUG->GUG: 5' {codones(mut_inicio)} 3'")
    print(f"  Si se leyera desde la posición 1: {mut_inicio.translate()} "
          "(GUG = Valina: ya no es la Met de inicio)")
    pos = mut_inicio.find("AUG")
    if pos != -1:
        resto = mut_inicio[pos:]
        resto = resto[: len(resto) - len(resto) % 3]
        print(f"  Un ribosoma eucariota buscaría el siguiente AUG (posición {pos + 1}, "
              f"otro marco de lectura): {codones(resto)} -> {resto.translate()} "
              "(péptido distinto y sin codón de paro)")

    # pérdida del codón de paro UAA -> CAA
    mut_paro = ARNM[:9] + Seq("CAA")
    print(f"\nMutación UAA->CAA: 5' {codones(mut_paro)} 3'")
    print(f"  Traducción: {mut_paro.translate()} -> se añade Gln y no hay señal de paro: "
          "el ribosoma seguiría leyendo la región 3' UTR (proteína más larga).")


if __name__ == "__main__":
    main()
