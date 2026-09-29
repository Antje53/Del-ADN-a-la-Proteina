from Bio.Seq import Seq

arnm_seq = Seq("AUGUAUGCUUAA")
print(f"ARNm original: {arnm_seq}")

# Traducción estándar (muestra el asterisco del codón de parada)
proteina_cruda = arnm_seq.translate()
print(f"Proteína (con señal de paro): {proteina_cruda}")

# Traducción limpia (se detiene en el codón de parada)
proteina_limpia = arnm_seq.translate(to_stop=True)
print(f"Proteína (cadena limpia):     {proteina_limpia}\n")