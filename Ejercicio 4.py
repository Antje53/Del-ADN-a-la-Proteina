import requests

# API REST de Ensembl para buscar el gen BRCA2 en humanos
server = "https://rest.ensembl.org"
ext = "/lookup/symbol/homo_sapiens/BRCA2?expand=1"

respuesta = requests.get(server + ext, headers={"Content-Type": "application/json"})

if respuesta.ok:
    datos_gen = respuesta.json()
    transcritos = datos_gen.get("Transcript", [])
    print(f"Gen encontrado: {datos_gen['display_name']} (ID: {datos_gen['id']})")
    print(f"Número de isoformas (transcritos) encontradas debido al splicing: {len(transcritos)}")
    
    # 3 transcritos como ejemplo
    for i, t in enumerate(transcritos[:3], 1):
        print(f"  {i}. ID: {t['id']}, Longitud: {t['length']} pb")
else:
    print("Error conectando a Ensembl.")
print()