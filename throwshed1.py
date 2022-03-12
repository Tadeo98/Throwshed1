#######################################################################
## THROWSHED ##
#######################################################################

## KNIZNICE
from matplotlib.style import use
from osgeo import gdal, ogr, osr
import numpy as np
import math
import os

#######################################################################
## CESTY
dem_path = r"D:\School\STU_SvF_BA\Term10\Diplomovka\Throwshed1\data\dem\dmr.tif" #cesta k dem
point_layer_path = r"D:\School\STU_SvF_BA\Term10\Diplomovka\Throwshed1\data\point\POINT.shp"   #cesta k bodovej vrstve
throwshed_output_folder = r"D:\School\STU_SvF_BA\Term10\Diplomovka\Throwshed1\data\throwshed"  #cesta k priecinku, kde sa ulozi subor
throwshed_file = r"throwshed1"   #nazov vystupneho suboru s cistym throwshedom
viewshed_file = r"viewshed1" #nazov vystupneho suboru s viewshedom (ak sa ma pouzit)
viewshed_clip_file = r"viewshed_clip1"   #nazov vystupneho suboru s orezanym viewshedom (ak sa ma pouzit viewshed)

## NASTAVENIA
use_viewshed = 1 #pouzitie viditelnosti na orezanie throwshedu, nie = 0, ano = 1
band_number = 1 #vybrane pasmo z dem, default = 1
int_compare = 1 #interpolacia DMR vo vypoctovych bodoch, nearest neighbour = 0, linear = 1
keep_point_crs = 0 #vystupna vrstva ma suradnicovy system ako vstupna bodova? ano = 1, nie, nastavim EPSG noveho SS = 0
EPSG = 8353 #EPSG kod (suradnicovy system) vystupnej vrstvy throwshedu

## PREMENNE
h = 1.6 #pociatocna vyska nad povrchom [m]
alfa_min = 45.0 #minimalny uhol hodu/vystrelu [°]
alfa_max = 45.0 #maximalny uhol hodu/vystrelu [°], polozit rovne alfa_min, ak sa striela iba pod jednym uhlom
g = -9.8 #gravitacne zrychlenie [m/s^2]
V_0 = 100 #pociatocna rychlost [m/s]
ro = 1.225 #hustota vzduchu [kg/m^3], pri t = 15°C, H = 0m, Fi = 0# (suchy vzduch) sa ro = 1.225 kg/m^3
C_d = 2.5 #koeficient odporu/ťahu objektu vo vzduchu
A = 0.00005 #plocha prierezu šípu [m^2], A = 0.0001 m^2 pri kruhovom priereze s priemerom cca 11 mm, A = 0.00005 m^2 pri kruhovom priereze s priemerom cca 8 mm
m = 0.030 #hmotnost šípu [kg]
dt = 0.001 #casovy interval [s]
da = 1.0 #krok v uhle vystrelu [°]
dA = 1.0 #krok v azimute [°]
dr = 1.0 #krok vzdialenosti, pod ktorou sa bude vzdy interpolovat DMR a porovnavat sa s trajektoriou [m]
h_E = 1.7 # vyska oci strielajuceho pre viewshed, defaultne 1.7 [m]
h_T = 1.7 # vyska ciela pre viewshed, defaltne 0.0 [m]


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
point_count = point_layer.GetFeatureCount()

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


# VYPOCET TRAJEKTORIE PROJEKTILU (x,y,y_r)
# kontrola, ci su hodnoty uhla vystrelu spravne definovane
if alfa_max < alfa_min:
    print("Minimalny uhol vystrelu ma vacsiu hodnotu ako maximalny. Opravit.")
    exit()
#Vytvorenie listu so vsetkymi hodnotami uhla vystrelu, ktore sa pouziju v cykle
alfa_arange = np.arange(alfa_min, alfa_max, da, dtype = np.float32) #arange() umoznuje vytvorit zoznam hodnot aj s krokom typu float (range to nedokaze)
alfa_list = alfa_arange.tolist() #transformacia typu np.ndarray na list, aby sa dal pouzit append()
alfa_list.append(alfa_max)   #pridelenie aj poslednej hranicnej hodnoty (inak by bola posledna hodnota o cosi mensia ako maximalny uhol vystrelu)

y_r = []    #buduca matica s vyskami projektilu kazdych dr metrov (hodnoty v riadku vedla seba) pod kazdym uhlom vystrelu (niekolko riadkov pod sebou)
j = 0
#cyklenie sa vsetkymi hodnotami alfa
for alfa in alfa_list:
    # Pociatocny drag
    d = -ro*V_0**2*C_d*A/(2*m)    #presny vztah
    # d = -C_d/1000*V_0**2    #priblizny vztah, pre sip postacuje

    r = dr   #radialna vzdialenost (interpolacny skok)
    V = V_0 #rychlost
    x = [0.0]   #suradnica x v 2D grafe
    # y = [h] #suradnica y v 2D grafe, iba vyska nad povrchom
    y = [h+dem_cell_height] #suradnica y v 2D grafe
    y_r1 = []    #bude obsahovat vysku kazdych dr metrov (iba jeden riadok)

    #zlozky pociatocnej rychlosti v smeroch x a y
    V_x = V*math.cos(alfa/180*math.pi)
    V_y = V*math.sin(alfa/180*math.pi)

    #zlozky pociatocneho dragu v smeroch x a y
    d_x = d*math.cos(alfa/180*math.pi)
    d_y = d*math.sin(alfa/180*math.pi)

    #kroky v x a y
    dX = V_x*dt+1/2*d_x*dt**2
    dY = V_y*dt+1/2*(d_y+g)*dt**2

    i = 0
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
                y_r1.append(((x[-1]-r)*(y[-2]-y[-1]))/(x[-1]-x[-2])+y[-1]) #vyska sipu vo vzdialenosti r
                r += dr
            break

        #v kazdom nasobku nastaveneho skoku radialnej vzialenosti sa vypocita vyska sipu
        if x[-1] > r:
            y_r1.append(((x[-1]-r)*(y[-2]-y[-1]))/(x[-1]-x[-2])+y[-1]) #vyska sipu vo vzdialenosti r
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

    y_r = y_r+[y_r1]
    j += 1

    # print(x[-1])
    # print(y[-1])

    # print(y_r)

    # # vykreslenie
    # plt.plot(x, y, 'r')
    # plt.xlabel("vzdialenosť [m]")
    # plt.ylabel("výška [m]")
    # plt.gca().set_aspect('equal', adjustable='box')
    # plt.show()

# definicia vektorov, do ktorych sa budu ukladat suradnice najvzdialenejsich bodov jednotlivych azimutov
X_coor_point_polygon = []
Y_coor_point_polygon = []
# Otacame pod azimutom a porovnavame hodnoty z y_r s DMR
Azimuth = 0  #azimut, smer sever
# cyklus, kde sa meni azimut
while True:
    S = []  #vektor vzdialenosti k najvzdialenejsim bodom pri jednotlivych uhloch vystrelu
    # cyklus kde sa prestriedaju vsetky trajektorie (vsetky uhly vystrelu)
    for i in range(0,len(alfa_list)):
        j = 0
        r = 0
        #cyklus, kde sa meni vzdialenost
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

            #porovnanie vysky bunky s vyskou sipu, ak je sip pod DMR, pripise sa maximalna vzdialenost pre konkretny uhol vystrelu
            if dem_int_point_height >= y_r[i][j]:
                S.append(r)
                break

            #ak by sa stalo ze aj po poslednej vyske sipu je stale sip nad DMR, zapise sa posledna mozna vzdialenost (k problemu moze dojst zrejme pri nastaveni velkeho kroku dr)
            if j == len(y_r[i])-1:
                S.append(r)
                break
                
            j += 1
        #nakoniec sa vyhlada maximalna vzdialenost spomedzi vsetkych pocitanych pre kazdy uhol vystrelu a zapisu sa suradnice najvzdialenejseho bodu pre dany azimut
        if i == range(0,len(alfa_list))[-1]:
            max_r, idx = max((max_r, idx) for (idx, max_r) in enumerate(S))
            X_coor_point_polygon.append(X_coor_point + max_r*math.sin(Azimuth/180*math.pi))
            Y_coor_point_polygon.append(Y_coor_point + max_r*math.cos(Azimuth/180*math.pi))
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
elif keep_point_crs == 1:
    srs = point_layer.GetSpatialRef()
else:
    print("Zle nastavena hodnota keep_point_crs.")
    throwshed_outds = None
    exit()
throwshed_outlayer = throwshed_outds.CreateLayer(throwshed_file, srs)

# pridanie polygonu do feature a jeho ulozenie do vystupnej vrstvy
throwshed_feature = ogr.Feature(throwshed_outlayer.GetLayerDefn())
throwshed_feature.SetGeometry(throwshed_polygon)
throwshed_outlayer.CreateFeature(throwshed_feature)

# Vyuzitie viewshedu
if use_viewshed == 1:
    # treba zavriet polygonovu vrstvu
    throwshed_outds = throwshed_outlayer = throwshed_feature = None
    # vytvorenie rastra viditelnosti, ulozi sa ako viewshed.tif do adresara s vystupnym throwshedom
    gdal.ViewshedGenerate(srcBand=dem_band, driverName='GTiff', targetRasterName=throwshed_output_folder + "\\" + viewshed_file + ".tif", creationOptions=None, observerX=X_coor_point, observerY=Y_coor_point, observerHeight=h_E, targetHeight=h_T, visibleVal=1, invisibleVal=0, outOfRangeVal=0, noDataVal=-9999, dfCurvCoeff=0.85714, mode=2, maxDistance=10000)
    # otvorenie viewshed rastra
    viewshed_ds = gdal.Open(throwshed_output_folder + "\\" + viewshed_file + ".tif")
    # orezanie rastra viditelnosti throwshedom
    gdal.Warp(throwshed_output_folder + "\\" + viewshed_clip_file + ".tif", viewshed_ds, cutlineDSName = throwshed_output_folder + "\\" + throwshed_file + ".shp", cropToCutline = True, dstNodata = np.nan)
    # vymazanie polygonu .shp s throwshedom a rastra .tif s viewshedom
    driver = ogr.GetDriverByName("ESRI Shapefile")
    driver.DeleteDataSource(throwshed_output_folder + "\\" + throwshed_file + ".shp")
    viewshed_ds = None
    os.remove(throwshed_output_folder + "\\" + viewshed_file + ".tif")
elif use_viewshed == 0:
    # najprv treba vytvorit raster a dat mu nastavenia
    throwshed_ds = gdal.GetDriverByName('GTiff').Create(throwshed_output_folder + "\\" + throwshed_file + ".tif", xsize = dem_array.shape[1], ysize = dem_array.shape[0], bands = 1, eType = gdal.GDT_Float32)
    throwshed_ds.SetGeoTransform(dem_gt)
    throwshed_ds.SetProjection(srs.ExportToWkt())   #SS bude nastaveny ako bol aj pri polygonovej vrstve
    throwshed_band = throwshed_ds.GetRasterBand(1)
    throwshed_band.SetNoDataValue(-9999)
    # throwshed polygon sa rasterizuje, [1] - priradenie hodnot do pasma 1, burn_values=[1] - priradenie hodnot buniek = 1
    gdal.RasterizeLayer(throwshed_ds, [1], throwshed_outlayer, burn_values=[1])
    # nakoniec novovytvorena vrstva, datasource aj prvok treba dat rovne None, lebo inak sa nezobrazi spravne v QGISe
    throwshed_outds = throwshed_ds = throwshed_outlayer = throwshed_feature = None
    # vektorova podoba sa vymaze a zostane len rastrova
    driver = ogr.GetDriverByName("ESRI Shapefile")
    driver.DeleteDataSource(throwshed_output_folder + "\\" + throwshed_file + ".shp")
else:
    throwshed_outds = throwshed_outlayer = throwshed_feature = None
    print("Zadana nespravna hodnota pri nastaveni pouzitia viditelnosti.")
    exit()