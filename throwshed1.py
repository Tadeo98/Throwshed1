#######################################################################
## THROWSHED ##
#######################################################################

## KNIZNICE
from osgeo import gdal, ogr, osr
import numpy as np
import math
import sys
import matplotlib.pyplot as plt
import os

#######################################################################
## CESTY
dem_path = r"D:\School\STU_SvF_BA\Term10\Diplomovka\Throwshed1\data\dem\dmr.tif" #cesta k dem
point_layer_path = r"D:\School\STU_SvF_BA\Term10\Diplomovka\Throwshed1\data\point\POINT.shp"   #cesta k bodovej vrstve
throwshed_output_folder = r"D:\School\STU_SvF_BA\Term10\Diplomovka\Throwshed1\data\throwshed"  #cesta k priecinku, kde sa ulozi subor
throwshed_file = r"skuska1"   #nazov vystupneho suboru s cistym throwshedom
viewshed_file = r"viewshed" #nazov vystupneho suboru s viewshedom (ak sa ma pouzit)
viewshed_clip_file = r"viewshed_clip"   #nazov vystupneho suboru s orezanym viewshedom (ak sa ma pouzit viewshed)

## NASTAVENIA
use_viewshed = 1 #pouzitie viditelnosti na orezanie throwshedu, nie = 0, ano = 1
band_number = 1 #vybrane pasmo z dem, default = 1
int_compare = 1 #interpolacia DMR vo vypoctovych bodoch, nearest neighbour = 0, linear = 1
keep_point_crs = 0 #vystupna vrstva ma suradnicovy system ako vstupna bodova? ano = 1, nie, nastavim EPSG noveho SS = 0
EPSG = 8353 #EPSG kod (suradnicovy system) vystupnej vrstvy throwshedu

## PREMENNE
h = 1.6 #pociatocna vyska nad povrchom [m]
alfa = 10.0 #uhol hodu/vystrelu [°]
g = -9.8 #gravitacne zrychlenie [m/s^2]
V_0 = 100 #pociatocna rychlost [m/s]
ro = 1.225 #hustota vzduchu [kg/m^3], pri t = 15°C, H = 0m, Fi = 0# (suchy vzduch) sa ro = 1.225 kg/m^3
C_d = 2.5 #koeficient odporu/ťahu objektu vo vzduchu
A = 0.00005 #plocha prierezu šípu [m^2], A = 0.0001 m^2 pri kruhovom priereze s priemerom cca 11 mm, A = 0.00005 m^2 pri kruhovom priereze s priemerom cca 8 mm
m = 0.030 #hmotnost šípu [kg]
dt = 0.001 #casovy interval [s]
dA = 1 #krok v azimute [°]
dr = 1 #krok vzdialenosti, pod ktorou sa bude vzdy interpolovat DMR a porovnavat sa s trajektoriou [m]
h_E = 1.7 # vyska oci strielajuceho pre viewshed, defaultne 1.7
h_T = 1.7 # vyska ciela pre viewshed, defaltne 0


#############################################################################
## VYPOCET


# import DEM
dem_ds = gdal.Open(dem_path)
# vyber pasma
dem_band = dem_ds.GetRasterBand(band_number)
# pridelenie hodnot dem do array
dem_array = dem_band.ReadAsArray()

# import bodovej vrstvy
point_ds = ogr.Open(point_layer_path, 0) #1 = editing, 0 = read only. Datasource
# bodova vrstva
point_layer = point_ds.GetLayer()

# ziskanie poctu bodov vo vrstve
featureCount = point_layer.GetFeatureCount()

# ziskanie konkretneho bodu
point_feature = point_layer.GetFeature(0)
# ziskanie geometrie bodu, X a Y suradnice (odpovedajuce smeru a orientacii osi v QGISe)
point_geom = point_feature.GetGeometryRef()
X_coor_point = point_geom.GetX()
Y_coor_point = point_geom.GetY()


# vytvorenie array so suradnicami X, Y vsetkych buniek (vztiahnute k stredu buniek)
dem_gt = dem_ds.GetGeoTransform()   #pole obsahujuce xorigin, yorigin, sirku a vysku pixela
# dem_array_coor = np.zeros((dem_array.shape[0],dem_array.shape[1],2))
# for i in range(0,dem_array.shape[0]):
#     for j in range(0,dem_array.shape[1]):
#         dem_array_coor[i][j][0] = dem_gt[0]+(j+1/2)*dem_gt[1]
#         dem_array_coor[i][j][1] = dem_gt[3]+(i+1/2)*dem_gt[5]

#VYPOCET SURADNIC KAZDEJ BUNKY ZATIAL NEPOTREBNY

# zistenie vysky bunky, na ktorej sa nachadza bod
dem_cell_column = round(abs((X_coor_point - (dem_gt[0]+dem_gt[1]/2))/dem_gt[1]))
dem_cell_row = round(abs((Y_coor_point - (dem_gt[3]+dem_gt[5]/2))/dem_gt[5]))
dem_cell_height = dem_array[dem_cell_row][dem_cell_column]


# Ziskanie najmensej vysky rastra
if dem_band.GetMinimum() == None:   #niekedy ale nevyhodi hodnotu None
    dem_band.ComputeStatistics(0)   #preto treba spustit tento prikaz na prepocet statistiky
    min_height = dem_band.GetMinimum()  #a uz by malo ziskat ciselnu hodnotu
else:
    min_height = dem_band.GetMinimum()  #ale ak nevyhadzuje None, tak rovno sa prideli hodnota


# Vypocet trajektorie pre mnozinu x, s dragom
# Pociatocny drag
d = -ro*V_0**2*C_d*A/(2*m)    #presny vztah
# d = -C_d/1000*V_0**2    #priblizny vztah, pre sip postacuje

i = 0
r = dr   #radialna vzdialenost (interpolacny skok)
V = V_0 #rychlost
x = [0.0]   #suradnica x v 2D grafe
# y = [h] #suradnica y v 2D grafe, iba vyska nad povrchom
y = [h+dem_cell_height] #suradnica y v 2D grafe
y_r = []    #bude obsahovat vysku kazdych dr metrov

#zlozky pociatocnej rychlosti v smeroch x a y
V_x = V*math.cos(alfa/180*math.pi)
V_y = V*math.sin(alfa/180*math.pi)

#zlozky pociatocneho dragu v smeroch x a y
d_x = d*math.cos(alfa/180*math.pi)
d_y = d*math.sin(alfa/180*math.pi)

#kroky v x a y
dX = V_x*dt+1/2*d_x*dt**2
dY = V_y*dt+1/2*(d_y+g)*dt**2

while True:
    i += 1
    #suradnice
    x = x + [x[-1] + dX]
    y = y + [y[-1] + dY]
    
    #posledna hodnota bude v najmensej vyske celeho rastra, nizsie uz netreba prepocitavat, cyklus sa skonci
    if y[i] < min_height:
        y[i] = min_height
        x[i] = x[i-1] + (y[i-1]-min_height)/math.tan(abs(alfa)/180*math.pi)
        #mala sanca ze krok radialnej vzdialenosti bude presiahnuty, preto pridelenie vysky sipu aj v tomto malo pravdepodobnom bode
        if x[-1] > r:
            y_r.append(((x[-1]-r)*(y[-2]-y[-1]))/(x[-1]-x[-2])+y[-1]) #vyska sipu vo vzdialenosti r
            r += dr
        break

    #v kazdom nasobku nastaveneho skoku radialnej vzialenosti sa vypocita vyska sipu
    if x[-1] > r:
        y_r.append(((x[-1]-r)*(y[-2]-y[-1]))/(x[-1]-x[-2])+y[-1]) #vyska sipu vo vzdialenosti r
        r += dr

    #novy uhol
    alfa = math.atan(dY/dX)/math.pi*180
    
    #nova rychlost
    V = math.sqrt((dX/dt)**2+(dY/dt)**2)
    
    #novy drag
    d = -ro*V**2*C_d*A/(2*m)   #presny vztah
    #d = -C_d/1000*V**2    #priblizny vztah, pre sip postacuje
    
    #zlozky pociatocneho dragu v smeroch x a y
    d_x = d*math.cos(alfa/180*math.pi)
    d_y = d*math.sin(alfa/180*math.pi)
    #zlozky rychlosti v smeroch x a y
    V_x = V_x + d_x*dt
    V_y = V_y + (d_y+g)*dt
    #kroky v x a y
    dX = V_x*dt+1/2*d_x*dt**2
    dY = V_y*dt+1/2*(d_y+g)*dt**2

# print(x[-1])
# print(y[-1])

# print(y_r)

# # vykreslenie
# plt.plot(x, y, 'r')
# plt.xlabel("vzdialenosť [m]")
# plt.ylabel("výška [m]")
# plt.gca().set_aspect('equal', adjustable='box')
# plt.show()


# Otacame pod azimutom a porovnavame hodnoty z y_r s DMR
Azimuth = 0  #azimut, smer sever
X_coor_point_polygon = []
Y_coor_point_polygon = []

while True:
    j = 0
    r = 0
    while True:
        r += dr
        #vypocet suradnic so vzdialenostou dr-nasobku a pod Azimutom
        X_coor_compare_point = X_coor_point + r*math.sin(Azimuth/180*math.pi)
        Y_coor_compare_point = Y_coor_point + r*math.cos(Azimuth/180*math.pi)

        # Interpolacia DMR v porovnavanom bode
        if int_compare == 0:
            #interpolacia DMR v bode porovnania (nearest neighbour)
            dem_int_cell_column = round(abs((X_coor_compare_point - (dem_gt[0]+dem_gt[1]/2))/dem_gt[1]))
            dem_int_cell_row = round(abs((Y_coor_compare_point - (dem_gt[3]+dem_gt[5]/2))/dem_gt[5]))
            dem_int_point_height = dem_array[dem_int_cell_row][dem_int_cell_column]
        elif int_compare == 1:
            #interpolacia DMR v bode porovnania (linear)
            dem_int_cell_column = np.floor(abs((X_coor_compare_point - (dem_gt[0]+dem_gt[1]/2))/dem_gt[1])).astype(np.int32) #najblizsii nizsi stlpec v array od bodu
            dem_int_cell_row = np.floor(abs((Y_coor_compare_point - (dem_gt[3]+dem_gt[5]/2))/dem_gt[5])).astype(np.int32)    #najblizsii nizsi riadok v array od bodu
            X_coor_cell_1 = dem_gt[0] + dem_gt[1]/2 + dem_int_cell_column*dem_gt[1] #X suradnica stredov lavych buniek
            # X_coor_cell_2 = dem_gt[0] + dem_gt[1]/2 + (dem_int_cell_column+1)*dem_gt[1] #X suradnica stredov pravych buniek
            Y_coor_cell_1 = dem_gt[3] + dem_gt[5]/2 + (dem_int_cell_row+1)*dem_gt[5] #Y suradnica stredov dolnych buniek
            # Y_coor_cell_2 = dem_gt[3] + dem_gt[5]/2 + dem_int_cell_row*dem_gt[5] #Y suradnica stredov hornych buniek
            H_1 = dem_array[dem_int_cell_row][dem_int_cell_column]  #H lavej hornej bunky
            H_2 = dem_array[dem_int_cell_row][dem_int_cell_column+1]  #H pravej hornej bunky
            H_3 = dem_array[dem_int_cell_row+1][dem_int_cell_column]  #H lavej dolnej bunky
            H_4 = dem_array[dem_int_cell_row+1][dem_int_cell_column+1]  #H pravej dolnej bunky
            H_int_1 = ((X_coor_compare_point-X_coor_cell_1)*(H_4-H_3))/(abs(dem_gt[1])) + H_3   #Interpolovana vyska na dolnej linii
            H_int_2 = ((X_coor_compare_point-X_coor_cell_1)*(H_2-H_1))/(abs(dem_gt[1])) + H_1   #Interpolovana vyska na hornej linii
            dem_int_point_height = ((Y_coor_compare_point-Y_coor_cell_1)*(H_int_2-H_int_1))/(abs(dem_gt[5])) + H_int_1   #Interpolovana vyska medzi dolnou a hornou liniou
        else:
            print("Hodnota int_compare neznama.")
            exit()

        #porovnanie vysky bunky s vyskou sipu, ak je sip pod DMR, zapise sa suradnica bodu
        if dem_int_point_height >= y_r[j]:
            X_coor_point_polygon.append(X_coor_compare_point)
            Y_coor_point_polygon.append(Y_coor_compare_point)
            break

        #ak by sa stalo ze aj po poslednej vyske sipu je stale sip nad DMR, zapise sa posledna mozna suradnica, aby aj pod tymto azimutom bol predsalen bod (k problemu moze dojst zrejme pri nastaveni velkeho kroku dr)
        if y_r[j] == y_r[-1]:
            X_coor_point_polygon.append(X_coor_compare_point)
            Y_coor_point_polygon.append(Y_coor_compare_point)
            print("Ojedinela situacia. Posledna porovnavana vyska sipu stale vyssie ako DMR. Zapise sa bod polygonu v najvacsej vzdielenosti z vypoctu. Pre vyhnutie sa problemu treba nastavit mensi krok dr")
            break
        j += 1
    Azimuth += dA
    #ukoncenie cyklu s meniacim sa Azimutom
    if Azimuth >= 360:
        break

#vypis suradnic
# print(X_coor_point_polygon)
# print(Y_coor_point_polygon)

#vykreslenie
# plt.plot(X_coor_point_polygon, Y_coor_point_polygon, 'ro')
# plt.show()


#######################################################################
## VYTVORENIE VYSTUPNEJ VRSTVY (VRSTIEV)

# vytvorenie novej geometrie
throwshed_ring = ogr.Geometry(ogr.wkbLinearRing)
# pridanie bodov do geometrie
for i in range(0,len(X_coor_point_polygon)):
    throwshed_ring.AddPoint(X_coor_point_polygon[i], Y_coor_point_polygon[i])
throwshed_ring.AddPoint(X_coor_point_polygon[0], Y_coor_point_polygon[0])   #posledny bod totozny s prvym

# vytvorenie polygonu
throwshed_polygon = ogr.Geometry(ogr.wkbPolygon)
throwshed_polygon.AddGeometry(throwshed_ring)

# ulozenie polygonu do vrstvy
driver = ogr.GetDriverByName("ESRI Shapefile")
throwshed_outds = driver.CreateDataSource(throwshed_output_folder + "\\" + throwshed_file + ".shp")

# definicia referencneho systemu
srs = osr.SpatialReference()
if keep_point_crs == 0:
    srs.ImportFromEPSG(EPSG)    
if keep_point_crs == 1:
    srs = point_layer.GetSpatialRef()
throwshed_outlayer = throwshed_outds.CreateLayer(throwshed_file, srs)

# pridanie polygonu do feature a jeho ulozenie do vystupnej vrstvy
throwshed_feature = ogr.Feature(throwshed_outlayer.GetLayerDefn())
throwshed_feature.SetGeometry(throwshed_polygon)
throwshed_outlayer.CreateFeature(throwshed_feature)

# nakoniec novovytvorena vrstva, datasource aj prvok treba dat rovne None, lebo inak sa nezobrazi spravne v QGISe
throwshed_outds = throwshed_outlayer = throwshed_feature = None

# Vyuzitie viewshedu
if use_viewshed == 1:
    # vytvorenie rastra viditelnosti, ulozi sa ako viewshed.tif do adresara s vystupnym throwshedom
    gdal.ViewshedGenerate(srcBand=dem_band, driverName='GTiff', targetRasterName=throwshed_output_folder + "\\" + viewshed_file + ".tif", creationOptions=None, observerX=X_coor_point, observerY=Y_coor_point, observerHeight=h_E, targetHeight=h_T, visibleVal=1, invisibleVal=0, outOfRangeVal=0, noDataVal=-9999, dfCurvCoeff=0.85714, mode=2, maxDistance=10000)
    # otvorenie viewshed rastra
    viewshed_ds = gdal.Open(throwshed_output_folder + "\\" + viewshed_file + ".tif")
    # orezanie rastra viditelnosti throwshedom
    gdal.Warp(throwshed_output_folder + "\\" + viewshed_clip_file + ".tif", viewshed_ds, cutlineDSName = throwshed_output_folder + "\\" + throwshed_file + ".shp", cropToCutline = True, dstNodata = np.nan)
