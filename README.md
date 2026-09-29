# Práctica de Bioinformática: Del ADN a la Proteína

**Institución:** Escuela de Ingeniería Informática — Universidad de Las Palmas de Gran Canaria (ULPGC)  
**Asignatura:** Bioinformática  
**Autoras:** Kimberly Casimiro Torres y Leonoor Antje Barton  

---

## Descripción del Proyecto

Este documento contiene la resolución de los ejercicios prácticos, abarcando el flujo completo del **Dogma Central de la Biología Molecular** (ADN → ARN → Proteína). Cada ejercicio combina la resolución teórica y conceptual con su implementación práctica automatizada mediante scripts en Python utilizando la librería **Biopython** y consultas a bases de datos biológicas e infraestructura REST (Ensembl).

## Estructura del repositorio

```text
.
├── README.md                  
├── Del_ADN_a_la_Proteina.pdf    
├── requirements.txt
├── Ejercicio_1.py               
├── Ejercicio_2.py              
├── Ejercicio_3.py               
├── Ejercicio_4.py              
├── Ejercicio_5.py               
├── Ejercicio_6.py              
├── datos/
    ├── ejercicio2.fasta
    └── lacZ_Ecoli_NC_000913.3.fasta

```

## Requisitos e Instalación

Para ejecutar los scripts de este repositorio se requiere **Python 3.8+** y las librerías indicadas en `requirements.txt` (principalmente `biopython` y `requests`).

```bash
pip install -r requirements.txt

python Ejercicio_1.py
python Ejercicio_2.py 
python Ejercicio_3.py
python Ejercicio_4.py
python Ejercicio_5.py
python Ejercicio_6.py

```

---

## Contenido de los Ejercicios

### Ejercicio 1. Replicación del ADN
* **Objetivo:** Comprender el mecanismo semiconservativo y las enzimas implicadas.
* **Contenido teórico:** 
  * Análisis de la complementariedad de bases a partir de la secuencia padre $5' - \text{ATG CCG TTA GCT} - 3'$.
  * Descripción funcional de las enzimas clave: **Helicasa**, **Primasa**, **ADN polimerasa** y **ADN ligasa**.
  * Reflexión sobre la fidelidad replicativa y las consecuencias de errores no corregidos en el genoma.

|  | Hebra parental (se conserva) | Hebra nueva (sintetizada) |
| --- | --- | --- |
| **Molécula hija 1** | 5'–`ATG CCG TTA GCT`–3' | 3'–`TAC GGC AAT CGA`–5' |
| **Molécula hija 2** | 3'–`TAC GGC AAT CGA`–5' | 5'–`ATG CCG TTA GCT`–3' |

**Reflexión (Error no corregido):** Si la polimerasa incorpora una base errónea y ni la corrección de pruebas ni la reparación (*mismatch repair*) la detectan, en la siguiente ronda de replicación esa base sirve de molde y el error queda fijado como mutación permanente.

**Salida de ejecución:**

```text
Molécula parental:   5' ATGCCGTTAGCT 3'
                     3' TACGGCAATCGA 5'

Hija 1:
   parental          5' ATGCCGTTAGCT 3'
   nueva             3' TACGGCAATCGA 5'
Hija 2:
   nueva             5' ATGCCGTTAGCT 3'
   parental          3' TACGGCAATCGA 5'

[OK] Las hebras nuevas coinciden con la resolución manual.
[OK] Ambas moléculas hijas son idénticas a la parental 

```

### Ejercicio 2. Transcripción del ADN a ARN
* **Objetivo:** Traducir correctamente la información de la cadena molde.
* **Contenido teórico:**
  * Identificación de la cadena molde ($3' \rightarrow 5'$).
  * Síntesis del transcrito de ARNm en dirección $5' \rightarrow 3'$ sustituyendo Timina (T) por Uracilo (U).
  * Diferenciación conceptual entre región promotora y región codificante.

**Salida de ejecución:**

```text
Registro: Ejercicio2  (12 nt)
  Hebra codificante  5' ATGCCTGAATGC 3'
  Hebra molde        3' TACGGACTTACG 5'

[CORRECTO] ARNm    5' AUGCCUGAAUGC 3'  -> proteína: MPEC

Experimento: ¿qué pasa si se usa mal la orientación?
  a) Transcribir el molde tal cual (3'->5'):  UACGGACUUACG  -> no es un ARNm: es el complementario del ARNm y está al revés
  b) Molde escrito 5'->3' como si fuera codificante: 5' GCAUUCAGGCAU 3'  -> proteína: AFRH
  c) Forma correcta desde el molde 5'->3' (reverse_complement): 5' AUGCCUGAAUGC 3'

Conclusión: solo leyendo el molde en sentido 3'->5' (o haciendo el complementario inverso) se obtiene el ARNm que empieza por AUG.

```

### Ejercicio 3. Traducción del ARNm a proteína
* **Objetivo:** Aplicar el código genético y reflexionar sobre mutaciones.
* **Contenido teórico:**
  * Identificación del codón de inicio (`AUG`) y codón de paro (`UAA`) a partir del transcrito `5' - AUG UAU GCU UAA - 3'`.
  * Traducción a cadena polipeptídica: Metionina -- Tirosina -- Alanina ($\text{Met-Tyr-Ala}$).
  * Análisis del impacto de mutaciones en el codón de inicio y pérdida del codón de parada.
**Salida de ejecución:**

```text
ARNm original: 5' AUG UAU GCU UAA 3'
  Traducción completa:        MYA* 
  Traducción hasta el paro:   MYA
  [OK] Coincide con la traducción manual Met-Tyr-Ala.

Mutación AUG->GUG: 5' GUG UAU GCU UAA 3'
  Si se leyera desde la posición 1: VYA* (GUG = Valina: ya no es la Met de inicio)
  Un ribosoma eucariota buscaría el siguiente AUG (posición 5, otro marco de lectura): AUG CUU -> ML (péptido distinto y sin codón de paro)

Mutación UAA->CAA: 5' AUG UAU GCU CAA 3'
  Traducción: MYAQ -> se añade Gln y no hay señal de paro: el ribosoma seguiría leyendo la región 3' UTR (proteína más larga).

```

### Ejercicio 4. Splicing alternativo
* **Objetivo:** Comprender cómo un mismo gen puede generar varias proteínas.
* **Contenido teórico:**
  * Diseño de combinaciones de exones (ej. 1-3-4-5 y 1-4-5) a partir de un gen de 5 exones.
  * Análisis de la diversidad proteica generada sin necesidad de incrementar el número total de genes.

Durante el *splicing* se eliminan los intrones y se unen los exones de manera diferencial. En nuestro análisis con la API REST de Ensembl, hemos analizado el gen **BRCA2** humano. De los transcritos anotados, distintas isoformas proteicas varían según la inclusión de ciertos exones.

**Salida de ejecución:**

```text
  ARNm 1-2-3-4-5  exones omitidos: ninguno
  ARNm 1-3-4-5    exones omitidos: [2]
  ARNm 1-4-5      exones omitidos: [2, 3]

Gen BRCA2 (ENSG00000139618), cromosoma 13, GRCh38 ===
Transcritos anotados: 19
Por tipo (biotype):
  nonsense_mediated_decay          8
  protein_coding                   7
  retained_intron                  4
Transcrito       Nombre       Exones  ARN (nt)  Proteína (aa)
ENST00000544455  BRCA2-206        27     11854            3418
ENST00000380152  BRCA2-201        27     11954            3418  <- canónico
ENST00000680887  BRCA2-210        27     11880            3418
ENST00000700202  BRCA2-214        27     11915            3401
ENST00000713680  BRCA2-219        26     11810            3366
ENST00000530893  BRCA2-204        27     11953            3295
ENST00000713678  BRCA2-217        27     11900            3232

Comparación ENST00000380152 (canónico) vs ENST00000680887:
  Exón solo en el canónico: ENSE00004011581 (posición 1)
  Exón solo en ENST00000680887: ENSE00003854030 (posición 1)
  Proteína: 3418 aa vs 3418 aa

```

### Ejercicio 5. Introducción a las proteínas
* **Objetivo:** Relacionar secuencia, estructura y función.
* **Contenido teórico:**
  * Identificación del extremo amino-terminal ($\text{N}$) y carboxilo-terminal ($\text{C}$) en la secuencia peptídica ($\text{Met-Ile-Ser-Gly-Val-Lys-His}$).
  * Estudio de la estructura primaria y terciaria, y reflexión sobre el efecto de mutaciones que sustituyan un aminoácido hidrofóbico por uno hidrofílico en regiones internas.

**Salida de ejecución:**

```text
Péptido: Met-Ile-Ser-Gly-Val-Lys-His  (MISGVKH)
  Extremo N:  Met
  Extremo C:  His

Pos  aa     Kyte-Doolittle  Tipo
1    Met               1.9  hidrofóbico
2    Ile               4.5  hidrofóbico
3    Ser              -0.8  hidrofílico
4    Gly              -0.4  hidrofílico
5    Val               4.2  hidrofóbico
6    Lys              -3.9  hidrofílico
7    His              -3.2  hidrofílico

Hidrofobicidad media (GRAVY) original: 0.33
Mutación Val5->Lys (MISGKKH): GRAVY = -0.83

```

### Ejercicio 6. Actividad integradora (Del ADN a la proteína)
* **Objetivo:** Recorrer el dogma central completo.
* **Contenido teórico y práctico:**
  * Selección y análisis del gen *lacZ* de *Escherichia coli* (referencia NC\_000913.3).
  * Ejecución paso a paso del pipeline completo: Replicación, Transcripción y Traducción.
  * Reflexión crítica sobre la vulnerabilidad a errores en cada etapa del dogma.

**Salida de ejecución:**

```text
 LECTURA DE LA SECUENCIA
  -> Archivo: lacZ_Ecoli_NC_000913.3.fasta
  -> Registro: NC_000913.3:c366305-363231 Escherichia coli str. K-12 substr. MG1655, complete genome
  -> Longitud: 3075 nt | Contenido GC: 56.3 %

 REPLICACIÓN
  ->   Parental A (codificante) 5' ATGACCATGATTACGGATTCACTGGCCGTC... 3'
  ->   Parental B (molde)       3' TACTGGTACTAATGCCTAAGTGACCGGCAG... 5'
  ->   Nueva sobre A            3' TACTGGTACTAATGCCTAAGTGACCGGCAG... 5'
  ->   Nueva sobre B            5' ATGACCATGATTACGGATTCACTGGCCGTC... 3'

 TRANSCRIPCIÓN
  -> La ARN polimerasa lee la hebra MOLDE en sentido 3'->5' ...
  -> ... y sintetiza el ARNm 5'->3' (U en lugar de T): 5' AUGACCAUGAUUACGGAUUCACUGGCCGUC... 3'
  -> Longitud del ARNm: 3075 nt. Comprobación: igual a la hebra codificante con T->U. [verificado]

[PASO 3] TRADUCCIÓN
  -> Codón de inicio: AUG (Met, correcto)
  -> El ribosoma lee 1025 codones de 3 nucleótidos...
  -> Codón de paro UAA encontrado en el codón nº 1025 (nucleótidos 3073-3075)
  -> Validación como CDS completa (inicio, sin paros internos, paro final): OK
  -> Proteína: 1024 aminoácidos
  ->   Extremo N: MTMITDSLAV...   Extremo C: ...YHYQLVWCQK

Resultados guardados en resultados/ARNm.fasta y resultados/proteina.fasta

 ADN -> ADN (replicación) -> ARNm -> proteína

```

---

## Conclusión

El código genético se conserva gracias a la estricta complementariedad de bases durante estos tres procesos. No obstante, variables como la orientación de la lectura de las hebras, la preservación del marco de lectura y la inclusión o descarte selectivo de exones determinan finalmente qué variante proteica se produce. Un único cambio a nivel nucleotídico tiene el potencial de alterar drásticamente la estructura tridimensional y función de la molécula resultante.

