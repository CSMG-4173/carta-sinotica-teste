# -*- coding: utf-8 -*-
"""
download_Data_v01.py

Este script foi criado para baixar dados do METAR, GFS e GOES automaticamente.
Autores: Vitor Tenório e João GM Ribeiro
Data de Criação: 11 jun 2024
Descrição: Este script baixa dados METAR, GFS e GOES para uma data e hora específicas
e armazena os arquivos na pasta 'scripts/data' dentro do repositório GitHub.
"""

# ------------------ IMPORTA BLIBLIOTECA --------------
import requests
import os
import urllib.request
import time as t
from datetime import datetime, timedelta

# ------------------  FUNÇÕES -------------------------
# Função para limpar arquivos antigos
def limpar_arquivos_antigos(diretorio, dias=7):
    agora = datetime.now()
    limite = agora - timedelta(days=dias)

    for filename in os.listdir(diretorio):
        file_path = os.path.join(diretorio, filename)
        if os.path.isfile(file_path):
            file_time = datetime.fromtimestamp(os.path.getmtime(file_path))
            if file_time < limite:
                print(f'Removendo arquivo antigo: {file_path}')
                os.remove(file_path)



# Obter a data e hora atuais
now = datetime.utcnow()

# Extrair dia, mês, ano e hora
dia = now.strftime('%d')
mes = now.strftime('%m')
ano = now.strftime('%Y')
hora = now.strftime('%H')

# Construa o URL do METAR com a hora arredondada
metar_url = f"https://thredds.ucar.edu/thredds/fileServer/noaaport/text/metar/metar_{ano}{mes}{dia}_{hora}00.txt"

# Diretório de download
output_dir = "scripts/data/"

# Nome do arquivo METAR com data e hora
metar_file_name = f"metar_data_{ano}{mes}{dia}_{hora}00.txt"
metar_file_path = os.path.join(output_dir, 'METAR/', metar_file_name)

# Faz o download do arquivo METAR
urllib.request.urlretrieve(metar_url, metar_file_path)

# Start the time counter
start_time = t.time()

print('---------------------------------------')
print('GFS Download (NOMADS) - Script started.')
print('---------------------------------------')

#-----------------------------------------------------------------------------------------------------------
# Download directory
resolution = '25'

# Desired extent
min_lon = '-100.00'
max_lon = '-20.00'
min_lat = '-65.00'
max_lat = '15.00'

# Construir o URL
url = f'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p{resolution}.pl?dir=%2Fgfs.{ano}{mes}{dia}%2F{hora}%2Fatmos&file=gfs.t{hora}z.pgrb2.0p{resolution}.f000&all_var=on&all_lev=on&subregion=&toplat={max_lat}&leftlon={min_lon}&rightlon={max_lon}&bottomlat={min_lat}'
print(url)

# Nome do arquivo a ser baixado
file_name = f'gfs.t{ano}{dia}{hora}z.pgrb2.0p{resolution}.f000'
file_path = os.path.join(output_dir, 'GFS/', file_name)

# Print do nome do arquivo
print("File name: ", file_name)

# Envio da requisição GET para baixar o arquivo
myfile = requests.get(url)

# Download do arquivo
open(file_path, 'wb').write(myfile.content)

# Printando a data e hora do download
print('\n---------------------')
print('Downloading GFS File:')
print('---------------------')
print('Resolution: ' + resolution)
print(f'Date: {dia}/{mes}/{ano}')
print(f'Run: {hora}')

#-----------------------------------------------------------------------------------------------------------

# End the time counter
print('\nTotal Processing Time:', round((t.time() - start_time), 2), 'seconds.')

# Limpar arquivos antigos
limpar_arquivos_antigos(output_dir)

