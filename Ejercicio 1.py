import sys
from Bio.Seq import Seq


def replicar(hebra_5_3: Seq):
    parental_a = hebra_5_3               
    parental_b = hebra_5_3.complement()    # 3'->5' 
    nueva_sobre_a = parental_a.complement()  
    nueva_sobre_b = parental_b.complement()  
    return (parental_a, nueva_sobre_a), (parental_b, nueva_sobre_b)


def main():
    secuencia = Seq(sys.argv[1].upper() if len(sys.argv) > 1 else "ATGCCGTTAGCT")
    (pa, na), (pb, nb) = replicar(secuencia)


    print(f"Molécula parental:   5' {pa} 3'")
    print(f"                     3' {pb} 5'\n")
    print("Hija 1:")
    print(f"   parental          5' {pa} 3'")
    print(f"   nueva             3' {na} 5'")
    print("Hija 2:")
    print(f"   nueva             5' {nb} 3'")
    print(f"   parental          3' {pb} 5'\n")

    # Comprobación con el resultado manual
    if len(sys.argv) == 1:
        manual_nueva_1 = "TACGGCAATCGA"   
        manual_nueva_2 = "ATGCCGTTAGCT"   
        assert str(na) == manual_nueva_1 and str(nb) == manual_nueva_2
        print("[OK] Las hebras nuevas coinciden con la resolución manual.")
    
    assert na == pb and nb == pa
    print("[OK] Ambas moléculas hijas son idénticas a la parental ")


if __name__ == "__main__":
    main()
