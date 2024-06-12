# -*- coding: utf-8 -*-
"""
Este script foi criado para baixar dados do METAR, GFS e GOES automaticamente.
Autores: Vitor Tenório e João GM Ribeiro
Data de Criação: 12 jun 2024
Descrição: Este script plota uma carta sinotica com dados METAR, GFS e GOES para uma data e hora específicas
e armazena os arquivos na pasta 'scripts/data/carta' dentro do repositório GitHub.
##Instalando
"""

!pip install matplotlib
!pip install pygrib
!pip install cartopy
!pip install metpy
!pip install netCDF4
!pip install siphon
!apt-get install libproj-dev proj-data proj-bin
!apt-get install libgeos-dev
!pip install cartopy
!apt-get -qq install python-cartopy python3-cartopy
!pip uninstall -y shapely
!pip install geopandas
!apt-get install libeccodes-dev libproj-dev
!pip install importlib-metadata==4.13.0
!pip install cfgrib
!pip install ecCodes
!pip install xarray
!pip install numpy
!pip install wrf-python==1.3.4.1

"""##Importando"""

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="bs4")

import matplotlib.pyplot as plt                                 # Plotting library
import cartopy, cartopy.crs as ccrs                             # Plot maps
import cartopy.io.shapereader as shpreader                      # Import shapefiles
import cartopy.feature as cfeature                              # Common drawing and filtering operations
import os                                                       # Miscellaneous operating system interfaces
import numpy as np                                              # Scientific computing with Python
import requests                                                 # HTTP library for Python
from datetime import timedelta, date, datetime                  # Basic Dates and time types
from metpy.calc import reduce_point_density                     # Provide tools for unit-aware, meteorological calculations
from metpy.io import metar                                      # Parse METAR-formatted data
from metpy.plots import current_weather, sky_cover, StationPlot # Contains functionality for making meteorological plots
import pygrib
import cfgrib
import geopandas as gpd
import glob
import numpy as np
import xarray as xr
import datetime

from datetime import datetime
import urllib.request
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from metpy.units import units
from netCDF4 import num2date
import numpy as np
from scipy.ndimage import gaussian_filter
from siphon.ncss import NCSS
from wrf import (smooth2d)

from datetime import datetime        # Basic Dates and time types
import os                            # Miscellaneous operating system interfaces
import requests                      # HTTP library for Python
import time as t                     # Time access and conversion
from bs4 import BeautifulSoup        # Library for web scraping
from bs4 import UnicodeDammit        # For encoding detection and correction

from google.colab import drive
drive.mount('/content/drive')

"""##Plotando o Mapa

###Função
"""

#-----------------------------------------------------------------------------------------------------------

def plot_maxmin_points(lon, lat, data, extrema, nsize, symbol, color='k',
                       plotValue=True, transform=None):
    """
    This function will find and plot relative maximum and minimum for a 2D grid. The function
    can be used to plot an H for maximum values (e.g., High pressure) and an L for minimum
    values (e.g., low pressue). It is best to used filetered data to obtain  a synoptic scale
    max/min value. The symbol text can be set to a string value and optionally the color of the
    symbol and any plotted value can be set with the parameter color
    lon = plotting longitude values (2D)
    lat = plotting latitude values (2D)
    data = 2D data that you wish to plot the max/min symbol placement
    extrema = Either a value of max for Maximum Values or min for Minimum Values
    nsize = Size of the grid box to filter the max and min values to plot a reasonable number
    symbol = String to be placed at location of max/min value
    color = String matplotlib colorname to plot the symbol (and numerica value, if plotted)
    plot_value = Boolean (True/False) of whether to plot the numeric value of max/min point
    The max/min symbol will be plotted on the current axes within the bounding frame
    (e.g., clip_on=True)
    """
    from scipy.ndimage.filters import maximum_filter, minimum_filter

    if (extrema == 'max'):
        data_ext = maximum_filter(data, nsize, mode='nearest')
    elif (extrema == 'min'):
        data_ext = minimum_filter(data, nsize, mode='nearest')
    else:
        raise ValueError('Value for hilo must be either max or min')

    mxy, mxx = np.where(data_ext == data)

    for i in range(len(mxy)):
        txt1 = ax.annotate(symbol, xy=(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]]), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), color=color, size=12,
                clip_on=True, annotation_clip=True, horizontalalignment='center', verticalalignment='center',
                transform=ccrs.PlateCarree())

        txt2 = ax.annotate('\n' + str(int(data[mxy[i], mxx[i]])), xy=(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]]), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
                color=color, size=6, clip_on=True, annotation_clip=True, fontweight='bold', horizontalalignment='center', verticalalignment='top',
                transform=ccrs.PlateCarree())

#-----------------------------------------------------------------------------------------------------------

"""###Baixando Dados"""

# Obtenha a hora atual
current_time = datetime.utcnow()

# Defina os horários desejados, vale ressaltar que está em UTC
desired_hours = [0, 6, 12, 18]

# Encontre o horário desejado mais próximo
nearest_desired_hour = max(filter(lambda x: x <= current_time.hour, desired_hours))

# Crie a hora arredondada para a URL
rounded_time = current_time.replace(hour=nearest_desired_hour, minute=0, second=0, microsecond=0)

# Formate a hora arredondada para o formato de URL
formatted_time = rounded_time.strftime("%Y%m%d_%H%M")

# Construa o URL do METAR com a hora arredondada
metar_url = f"https://thredds.ucar.edu/thredds/fileServer/noaaport/text/metar/metar_{formatted_time}.txt"

# Agora, você pode usar a variável metar_url para baixar o METAR no horário desejado.
print('A URL do METAR no horário desejado é: ', formatted_time)

# Faz o download do arquivo
urllib.request.urlretrieve(metar_url, "metar_data.txt")

nearest_desired_hour = '{:02d}'.format(nearest_desired_hour)
print(nearest_desired_hour)

#-----------------------------------------------------------------------------------------------------------
#Author: Vitor Tenório
#-----------------------------------------------------------------------------------------------------------

# Start the time counter
start_time = t.time()

print('---------------------------------------')
print('GFS Download (NOMADS) - Script started.')
print('---------------------------------------')

#-----------------------------------------------------------------------------------------------------------
# Download directory
dir = "Samples"; os.makedirs(dir, exist_ok=True)

# Desired date (last 10 days only!): Format - 'YYYYMMDD'
date = datetime.today().strftime('%Y%m%d')

#Vou pegar com resolução de 25KM
resolution = '25'

# Data que foi rodado
hour_run = nearest_desired_hour
nearest_desired_hour = 20240423

# Desired extent
min_lon = '-100.00'
max_lon = '-20.00'
min_lat = '-65.00'
max_lat = '15.00'

# Formate a data arredondada para o formato de URL
formatted_time_link = rounded_time.strftime("%Y%m%d")
def check_gfs_availability(formatted_time_link, nearest_desired_hour, resolution):
    # Construir o URL
    url = f'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p{resolution}.pl?dir=%2Fgfs.{formatted_time_link}%2F{nearest_desired_hour}%2Fatmos&file=gfs.t{nearest_desired_hour}z.pgrb2.0p{resolution}.f000&all_var=on&all_lev=on&subregion=&toplat={max_lat}&leftlon={min_lon}&rightlon={max_lon}&bottomlat={min_lat}'
    print('Por favor confere se o link está correto e funcionando em caso de Erro:')
    print(url)

    # Enviar solicitação GET para obter o conteúdo da página
    response = requests.get(url)

    # Verificar se a resposta foi bem-sucedida (código de status 200)
    if response.status_code == 200:
        # Tenta detectar e corrigir problemas de codificação
        dammit = UnicodeDammit(response.content, ["latin-1", "iso-8859-1", "windows-1251"])

        # Obtém o conteúdo com codificação corrigida
        content = dammit.unicode_markup

        # Cria o objeto BeautifulSoup com o conteúdo corrigido
        soup = BeautifulSoup(content, "html.parser")

        # Verificar se a página contém a mensagem indicando que o arquivo não está presente
        error_message = soup.find(text='Data file is not present:') or soup.find(text='Data file is not present.') or soup.find(text='Data file is not present:')

        if error_message:
            print(f"Os dados do GFS para a hora desejada ({nearest_desired_hour}Z) ainda não estão disponíveis.")
            return False
        else:
            print('---------------------------------------')
            print("Os dados do GFS estão disponíveis para a hora desejada.")
            print('Date: ' + date)
            print('Run: ' + hour_run)
            print('---------------------------------------')
            return True
    else:
        print("Falha na solicitação HTTP:", response.status_code)
        return False

# Verificar disponibilidade dos dados
if not check_gfs_availability(formatted_time_link, nearest_desired_hour, resolution):
    print("Não foi possível baixar os dados do GFS pois inda não estão disponíveis.")
else:
    # Construir o URL de download
    # Construir o URL
    url = f'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p{resolution}.pl?dir=%2Fgfs.{formatted_time_link}%2F{nearest_desired_hour}%2Fatmos&file=gfs.t{nearest_desired_hour}z.pgrb2.0p{resolution}.f000&all_var=on&all_lev=on&subregion=&toplat={max_lat}&leftlon={min_lon}&rightlon={max_lon}&bottomlat={min_lat}'

    # Nome do arquivo a ser baixado
    file_name = f'gfs.t{hour_run}z.pgrb2.0p{resolution}.f000'

    # Print do nome do arquivo
    print("File name: ", file_name)

    # Envio da requisição GET para baixar o arquivo
    myfile = requests.get(url)

    # Download do arquivo
    open(os.path.join(dir, file_name), 'wb').write(myfile.content)

    # Printando a data e hora do download
    print('\n---------------------')
    print('Downloading GFS File:')
    print('---------------------')
    print('Resolution: ' + resolution)
    print('Date: ' + date)
    print('Run: ' + hour_run)

#-----------------------------------------------------------------------------------------------------------

# End the time counter
print('\nTotal Processing Time:', round((t.time() - start_time),2), 'seconds.')

"""#### Salvando p código para baixar o GFS em LOOP"""

'''
#-----------------------------------------------------------------------------------------------------------
# INPE / CPTEC - Training: NWP Data Processing With Python - NWP Download with Python (GFS)
# Author: Diego Souza
#-----------------------------------------------------------------------------------------------------------

# Required modules
from datetime import datetime        # Basic Dates and time types
import os                            # Miscellaneous operating system interfaces
import requests                      # HTTP library for Python
import time as t                     # Time access and conversion
#-----------------------------------------------------------------------------------------------------------

print('---------------------------------------')
print('GFS Download (NOMADS) - Script started.')
print('---------------------------------------')

# Start the time counter
start_time = t.time()

#-----------------------------------------------------------------------------------------------------------

# Download directory
dir = "Samples"; os.makedirs(dir, exist_ok=True)

# Desired date (last 10 days only!): Format - 'YYYYMMDD'
date = datetime.today().strftime('%Y%m%d')

# Desired extent
min_lon = '-100.00'
max_lon = '-20.00'
min_lat = '-65.00'
max_lat = '15.00'

# Desired resolution: '25' or '50' or '1'
resolution = '25'

# Desired run: '0' or '6' or '12' or '18'
hour_run = '00'

# Desired forecast hours
hour_ini = 0  # Init time
hour_end = 24 # End time
hour_int = 6  # Interval

#-----------------------------------------------------------------------------------------------------------

# Link (select "grib filter" and check "Show the URL only for web programming" to verify the URL's):
# https://nomads.ncep.noaa.gov/

def download_gfs(date, iii):

    # Create the URL's based on the resolution
    if (resolution == '25'):
        url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p'+resolution+'.pl?file=gfs.t'+hour_run+'z.pgrb2.0p'+resolution+'.f'+str(hour).zfill(3)+'&all_lev=on&all_var=on&subregion=&leftlon='+min_lon+'&rightlon='+max_lon+'&toplat='+max_lat+'&bottomlat='+min_lat+'&dir=%2Fgfs.'+date+'%2F00%2Fatmos'
        file_name = 'gfs.t'+hour_run+'z.pgrb2.0p'+resolution+'.f'+str(hour).zfill(3)
    elif (resolution == '50'):
        url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p'+resolution+'.pl?file=gfs.t'+hour_run+'z.pgrb2full.0p'+resolution+'.f'+str(hour).zfill(3)+'&all_lev=on&all_var=on&subregion=&leftlon='+min_lon+'&rightlon='+max_lon+'&toplat='+max_lat+'&bottomlat='+min_lat+'&dir=%2Fgfs.'+date+'%2F00%2Fatmos'
        file_name = 'gfs.t'+hour_run+'z.pgrb2.0p'+resolution+'.f'+str(hour).zfill(3)
    elif (resolution == '1'):
        url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_'+resolution+'p00.pl?file=gfs.t'+hour_run+'z.pgrb2.'+resolution+'p00.f'+str(hour).zfill(3)+'&all_lev=on&all_var=on&subregion=&leftlon='+min_lon+'&rightlon='+max_lon+'&toplat='+max_lat+'&bottomlat='+min_lat+'&dir=%2Fgfs.'+date+'%2F00%2Fatmos'
        file_name = 'gfs.t'+hour_run+'z.pgrb2.'+resolution+'p00.f'+str(hour).zfill(3)

    # Print the file name
    print("File name: ", file_name)
    # Sends a GET request to the specified url
    myfile = requests.get(url)

    # Download the file
    open(dir + '//' + file_name, 'wb').write(myfile.content)

#-----------------------------------------------------------------------------------------------------------

# Download loop
for hour in range(hour_ini, hour_end + 1, hour_int):
    print('\n---------------------')
    print('Downloading GFS File:')
    print('---------------------')
    print('Resolution: ' + resolution)
    print('Date: ' + date)
    print('Run: ' + hour_run)
    print('Forecast Hour: f' + str(hour).zfill(3))
    # Call the download function
    download_gfs(date,hour)

#-----------------------------------------------------------------------------------------------------------

# End the time counter
print('\nTotal Processing Time:', round((t.time() - start_time),2), 'seconds.')

'''

"""###Abrindo Dados"""

data = metar.parse_metar_file('metar_data.txt')
data = data.dropna(how='any', subset=['wind_direction', 'wind_speed'])

# Substitua a hora no nome do arquivo
filename = f"/content/Samples/gfs.t{nearest_desired_hour}z.pgrb2.0p25.f000"

# Abra o arquivo GRIB com a hora substituída
grib = pygrib.open(filename)

!rm "/content/Samples/*.idx"
ds_s = xr.open_mfdataset(f"/content/Samples/gfs.t{nearest_desired_hour}z.pgrb2.0p25.f000", engine="cfgrib", concat_dim = 'step', combine='nested',
   backend_kwargs=dict(filter_by_keys={'stepType': 'instant','typeOfLevel':'isobaricInhPa'}))
ds_s

#Convertendo as lat e lon
ds_s.coords['longitude']=(ds_s.coords['longitude'] + 180) % 360 - 180
ds_s = ds_s.sortby(ds_s.longitude)

"""###Plotando o Mapa"""

#selecionando a região de interesse
lon_slice=slice(-100.00, -20.00)
lat_slice=slice(-65.00, 15.00)

# Select the extent [min. lon, min. lat, max. lon, max. lat]
extent = [-100.0, -65.00, -20.00, 15.00]

#-----------------------------------------------------------------------------------------------------------

# Select the variable
prmls = grib.select(name='Pressure reduced to MSL')[0]

# Get information from the file
init  = str(prmls.analDate)      # Init date / time
run   = str(prmls.hour).zfill(2) # Run
ftime = str(prmls.forecastTime)  # Forecast hour
valid = str(prmls.validDate)     # Valid date / time
print('Init: ' + init + ' UTC')
print('Run: ' + run + 'Z')
print('Forecast: +' + ftime)
print('Valid: ' + valid + ' UTC')

# Read the data for a specific region
prmls, lats, lons = prmls.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)

# Convert to hPa
prmls = prmls / 100

# Grab pressure level data
hght_1000 = ds_s['gh'].sel(latitude = lat_slice,longitude = lon_slice, step = ds_s.step[0], isobaricInhPa=1000.)
hght_500  = ds_s['gh'].sel(latitude = lat_slice,longitude = lon_slice, step = ds_s.step[0], isobaricInhPa=500.)

# Calculate and smooth 1000-500 hPa thickness
thickness_1000_500 = gaussian_filter(hght_500 - hght_1000, sigma=3.0)

thickness_1000_500
max_valor_gh = np.max(thickness_1000_500)
min_valor_gh = np.min(thickness_1000_500)
prmls
max_valor_prmls = np.max(prmls)
min_valor_prmls = np.min(prmls)

# Set projection of map display
mapproj = ccrs.LambertConformal(central_longitude=-95, central_latitude=35, standard_parallels=[35])
# Set projection of data
dataproj = ccrs.PlateCarree()

# Altere o DPI da figura resultante. DPI mais alto melhora drasticamente o
#aparência da renderização do texto.
plt.rcParams['savefig.dpi'] = 300

# Crie a figura e um eixo definido para a projeção.
fig = plt.figure(figsize=(20, 10))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# Carregue o shapefile dos estados do Brasil
dir = '/content/drive/MyDrive/Carta_Sinotica/'
estados_shapefile = dir + '/input/BR_UF_2021.shp'
estados_gdf = gpd.read_file(estados_shapefile)

# Plote o mapa dos estados
estados_gdf.boundary.plot(ax=ax, linewidth=0.3, color='grey',alpha=0.5)

# Adicione vários elementos do mapa ao gráfico para torná-lo reconhecível.
ax.add_feature(cfeature.LAND,facecolor='lightgray')
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.COASTLINE, color='grey')
ax.add_feature(cfeature.BORDERS, color='grey')

# Definir limites do gráfico
ax.set_extent((-100, -20, -55, 10))

#
# --------------Aqui está o gráfico real da estação-----------------------------------
#

# Inicie o gráfico da estação especificando os eixos nos quais desenhar, bem como
# lon/lat das estações (com transformação). Também ajustamos o tamanho da fonte para 7 pt.
stationplot = StationPlot(ax, data['longitude'].values, data['latitude'].values,
                          clip_on=True, transform=ccrs.PlateCarree(), fontsize=7)

# Trace a temperatura e o ponto de orvalho na parte superior e inferior esquerda, respectivamente, de
# o ponto central. Cada um usa uma cor diferente.
stationplot.plot_parameter('NW', data['air_temperature'].values, color='red')
stationplot.plot_parameter('SW', data['dew_point_temperature'].values, color='darkgreen')

# Um exemplo mais complexo usa um formatador personalizado para controlar como a pressão ao nível do mar
# valores são plotados. Isso usa os três dígitos finais padrão do valor da pressão em décimos de milibares.
stationplot.plot_parameter('NE', data['air_pressure_at_sea_level'].values,
                           formatter=lambda v: format(10 * v, '.0f')[-3:])

# Trace os símbolos de cobertura de nuvens no local central. Isso usa os códigos feitos acima e
# usa o mapeador `sky_cover` para converter esses valores em códigos de fonte para o fonte do símbolo climático.
stationplot.plot_symbol('C', data['cloud_coverage'].values, sky_cover)

# O mesmo desta vez, mas plote o clima atual à esquerda do centro, usando o
# Mapeador `current_weather` para converter símbolos nos glifos corretos.
stationplot.plot_symbol('W', data['current_wx1_symbol'].values, current_weather)

# Adicione farpas de vento
stationplot.plot_barb(data['eastward_wind'].values, data['northward_wind'].values)

# Também plote o texto real do ID da estação. Em vez de direções cardeais,
# plote mais detalhadamente especificando uma localização de 2 incrementos em x e 0 em y.
# stationplot.plot_text((2, 0), data['station_id'].values)

#--------------------------Plotando Altas e Baixas-------------------------------------------------------------------------
# Select the extent [min. lon, min. lat, max. lon, max. lat]
extent = [-100.0, -65.00, -20.00, 15.00]

#-----------------------------------------------------------------------------------------------------------

# Select the variable
prmls = grib.select(name='Pressure reduced to MSL')[0]

# Get information from the file
init  = str(prmls.analDate)      # Init date / time
run   = str(prmls.hour).zfill(2) # Run
ftime = str(prmls.forecastTime).zfill(2)  # Forecast hour
valid = str(prmls.validDate)     # Valid date / time
print('Init: ' + init + ' UTC')
print('Run: ' + run + 'Z')
print('Valid: ' + valid + ' UTC')

# Read the data for a specific region
prmls, lats, lons = prmls.data(lat1=extent[1],lat2=extent[3],lon1=extent[0]+360,lon2=extent[2]+360)

# Convert to hPa
prmls = prmls / 100

# Define de contour interval
data_min = min_valor_prmls
data_max = max_valor_prmls
interval = 3
levels = np.arange(data_min,data_max,interval)

#---------------------------Plotando as linhas---------------------------------------

# Plot thickness with multiple colors
max_valor_gh = np.max(thickness_1000_500)
min_valor_gh = np.min(thickness_1000_500)

clevs = (np.arange(0, 5400, 40),
         np.array([5400]),
         np.arange(5460, 7000, 40))

# Definindo as cores
azul_claro = (0.4, 0.6, 1)  # Azul claro em formato RGB
azul_escuro = (0, 0, 0.5)  # Azul escuro em formato RGB
vermelho = (1, 0, 0)  # Vermelho em formato RGB

# Usando as cores nas suas funções de plotagem
colors = (azul_claro, azul_escuro, vermelho)

kw_clabels = {'fontsize': 11, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
              'rightside_up': True, 'use_clabeltext': True}

for clevthick, color in zip(clevs, colors):
    cs = ax.contour(lons, lats, thickness_1000_500, levels=clevthick, colors=[color],
                    linewidths=1.0, linestyles='dashed', transform=dataproj)
    plt.clabel(cs, **kw_clabels)

# Smooth the sea level pressure since it tends to be noisy near the
# mountains
prmls = smooth2d(prmls, 3, cenweight=1)

# plotar as linhas
img1 = ax.contour(lons, lats, prmls, colors='black', linewidths=0.5, levels=levels)
ax.clabel(img1, inline=1, inline_spacing=0, fontsize='10',fmt = '%1.0f', colors= 'black')

# Use definition to plot H/L symbols
plot_maxmin_points(lons, lats, prmls, 'max', 50, symbol='A', color='b',  transform=ccrs.PlateCarree())
plot_maxmin_points(lons, lats, prmls, 'min', 25, symbol='B', color='r', transform=ccrs.PlateCarree())

#-----------------------Titulo da figura e Salvamento-------------------------------------------------------------------------------------
# Formate a data e hora para a data
data_para_titulo = rounded_time.strftime('%d/%m/%Y - %H:00h')
plt.title(f'METAR + Centro de Altas e Baixas do dia: {data_para_titulo} UTC', fontsize=12)
data_para_figura = rounded_time.strftime('%d_%m_%Y_%Hh')

# Adicionar linhas de latitude e longitude com rótulos apenas do lado direito (latitude) e na parte inferior (longitude)
gridlines = ax.gridlines(draw_labels=True, linewidth=0, color='black', x_inline=False, y_inline=False)

# Ajustar as localizações das marcações de latitude e longitude
gridlines.ylocator = plt.FixedLocator([-60, -50, -40, -30, -20, -10, 0, 10])
gridlines.xlocator = plt.FixedLocator([-100, -80, -60, -40, -20])

#Créditos
ax.text(0.985, 0.985, 'Desenvolvido por Vitor Tenório\nSupervisado por Michelle Reboita',
        horizontalalignment='right', verticalalignment='top', transform=ax.transAxes,
        fontsize=10, bbox=dict(facecolor='white', alpha=1, edgecolor='none', pad=5),
        zorder=10)

# Ajustar o layout
fig.tight_layout()
plt.savefig( dir + '/output/Carta_Sinotica_'+str(data_para_figura)+'.png', dpi=300, bbox_inches='tight')
plt.show()
