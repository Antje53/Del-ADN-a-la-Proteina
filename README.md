# Práctica de Bioinformática: Del ADN a la Proteína

**Institución:** Escuela de Ingeniería Informática — Universidad de Las Palmas de Gran Canaria (ULPGC)  
**Asignatura:** Bioinformática  
**Autoras:** Kimberly Casimiro Torres y Leonoor Antje Barton  

---

## Descripción del Proyecto

Este repositorio contiene la resolución de los ejercicios prácticos, abarcando el flujo completo del **Dogma Central de la Biología Molecular** (ADN $\rightarrow$ ARN $\rightarrow$ Proteína). Cada ejercicio combina la resolución teórica y conceptual con su implementación práctica automatizada mediante scripts en Python utilizando la librería **Biopython** y consultas a bases de datos biológicas e infraestructura REST (Ensembl).

---

## Contenido de los Ejercicios

### Ejercicio 1. Replicación del ADN
* **Objetivo:** Comprender el mecanismo semiconservativo y las enzimas implicadas.
* **Contenido teórico:** 
  * Análisis de la complementariedad de bases a partir de la secuencia padre $5' - \text{ATG CCG TTA GCT} - 3'$.
  * Descripción funcional de las enzimas clave: **Helicasa**, **Primasa**, **ADN polimerasa** y **ADN ligasa**.
  * Reflexión sobre la fidelidad replicativa y las consecuencias de errores no corregidos en el genoma.
* **Script asociado:** `Ejercicio 1.py` (generación automática de la hebra complementaria con Biopython).

### Ejercicio 2. Transcripción del ADN a ARN
* **Objetivo:** Traducir correctamente la información de la cadena molde.
* **Contenido teórico:**
  * Identificación de la cadena molde ($3' \rightarrow 5'$).
  * Síntesis del transcrito de ARNm en dirección $5' \rightarrow 3'$ sustituyendo Timina (T) por Uracilo (U).
  * Diferenciación conceptual entre región promotora y región codificante.
* **Script asociado:** `Ejercicio 2.py` (creación de un archivo FASTA temporal, lectura de registros, transcripción directa y análisis de la hebra molde invertida).

### Ejercicio 3. Traducción del ARNm a proteína
* **Objetivo:** Aplicar el código genético y reflexionar sobre mutaciones.
* **Contenido teórico:**
  * Identificación del codón de inicio (`AUG`) y codón de paro (`UAA`) a partir del transcrito `5' - AUG UAU GCU UAA - 3'`.
  * Traducción a cadena polipeptídica: Metionina -- Tirosina -- Alanina ($\text{Met-Tyr-Ala}$).
  * Análisis del impacto de mutaciones en el codón de inicio y pérdida del codón de parada.
* **Script asociado:** `Ejercicio 3.py` (traducción estándar y limpia en Biopython).

### Ejercicio 4. Splicing alternativo
* **Objetivo:** Comprender cómo un mismo gen puede generar varias proteínas.
* **Contenido teórico:**
  * Diseño de combinaciones de exones (ej. 1-3-4-5 y 1-4-5) a partir de un gen de 5 exones.
  * Análisis de la diversidad proteica generada sin necesidad de incrementar el número total de genes.
* **Script asociado:** `Ejercicio 4.py`(consulta automatizada a la API REST de Ensembl para extraer las isoformas y transcritos del gen humano `BRCA2`).

### Ejercicio 5. Introducción a las proteínas
* **Objetivo:** Relacionar secuencia, estructura y función.
* **Contenido teórico:**
  * Identificación del extremo amino-terminal ($\text{N}$) y carboxilo-terminal ($\text{C}$) en la secuencia peptídica ($\text{Met-Ile-Ser-Gly-Val-Lys-His}$).
  * Estudio de la estructura primaria y terciaria, y reflexión sobre el efecto de mutaciones que sustituyan un aminoácido hidrofóbico por uno hidrofílico en regiones internas.

### Ejercicio 6. Actividad integradora (Del ADN a la proteína)
* **Objetivo:** Recorrer el dogma central completo.
* **Contenido teórico y práctico:**
  * Selección y análisis del gen *lacZ* de *Escherichia coli* (referencia NC\_000913.3).
  * Ejecución paso a paso del pipeline completo: Replicación, Transcripción y Traducción.
  * Reflexión crítica sobre la vulnerabilidad a errores en cada etapa del dogma.
* **Script asociado:** `Ejercicio 6.py`(pipeline automatizado con Biopython que informa en tiempo de ejecución de cada proceso biológico simulado).

---

## Requisitos e Instalación

Para ejecutar los scripts de este repositorio se requiere **Python 3.8+** y las siguientes librerías científicas:

```bash
pip install biopython requests
