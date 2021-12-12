#######################################################################
## THROWSHED ##
#######################################################################

## KNIZNICE
from osgeo import gdal, ogr, osr
import numpy as np
import math
import sys
import matplotlib.pyplot as plt

#######################################################################
## PREMENNE A CESTY
dem_path = r"D:\School\STU_SvF_BA\Term10\Diplomovka\Throwshed1\data\dem\dmr.tif" #cesta k dem
point_layer_path = r"D:\School\STU_SvF_BA\Term10\Diplomovka\Throwshed1\data\point\POINT.shp"   #cesta k bodovej vrstve
band_number = 1 #vybrane pasmo z dem, default = 1
throwshed_output_folder = r"D:\School\STU_SvF_BA\Term10\Diplomovka\Throwshed1\data\throwshed"  #cesta k priecinku, kde sa ulozi subor
throwshed_file = r"2_throwshed_nn_1m_1m_rad_25deg"   #nazov vystupneho suboru
EPSG = 8353 #EPSG kod vystupnej vrstvy throwshedu
int_compare = 1 #interpolacia DMR vo vypoctovych bodoch, nearest neighbour = 0, linear = 1
# POTOM ZISTIT PRECO VRSTVA S TAKTO PRIDELENYM REF. SYS. CEZ EPSG (dole) nema priradeny referencny system po importe do QGISu

#konkretny atribut a konkretna jeho hodnota este nevyriesena
#point_id = 1    #vyber bodu z bodovej vrstvy, nie je to id z atributovej tabulky, ale poradie bodu
#point_attribute = "id" 

h = 1.6 #pociatocna vyska nad povrchom [m]
alfa =45.0 #uhol hodu/vystrelu [°]
g = -9.8 #gravitacne zrychlenie [m/s^2]
V_0 = 67.0 #pociatocna rychlost [m/s]
ro = 1.225 #hustota vzduchu [kg/m^3], pri t = 15°C, H = 0m, Fi = 0# (suchy vzduch) sa ro = 1.225 kg/m^3
C_d = 1.0 #koeficient odporu/ťahu objektu vo vzduchu
A = 0.0005 #plocha prierezu šípu [m^2], A = 0.0001 m^2 pri kruhovom priereze s priemerom cca 11 mm
m = 0.030 #hmotnost šípu [kg]
dt = 0.001 #casovy interval [s]
dA = 1 #krok v azimute [°]
dr = 1 #krok vzdialenosti, pod ktorou sa bude vzdy interpolovat DMR a porovnavat sa s trajektoriou [m]


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

# ZISTIT AKO NA TUTO SELEKCIU
# # selekcia bodu s konkretnou hodnotou konkretneho atributu
# point_layer.SetAttributeFilter(point_attribute = point_id)


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


# Vypocet trajektorie pre mnozinu x, s dragom, na rovine
# Pociatocny drag
d = -ro*V_0**2*C_d*A/(2*m)

i = 0
r = dr   #radialna vzdialenost (interpolacny skok)
V = V_0 #rychlost
x = [0.0]   #suradnica x v 2D grafe
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
    
    #v kazdom nasobku nastaveneho skoku radialnej vzialenosti sa vypocita vyska sipu
    if x[-1] > r:
        y_r.append(((x[-1]-r)*(y[-2]-y[-1]))/(x[-1]-x[-2])+y[-1]) #vyska sipu vo vzdialenosti r
        r += dr

    #posledna hodnota bude vo vyske 0 (mala pravdepodobnost, ze niekto bude hladat projektily v zapornej vyske)
    if y[i] < 0:
        y[i] = 0
        x[i] = x[i-1] + y[i-1]/math.sin(abs(alfa)/180*math.pi)
        break

    #novy uhol
    alfa = math.atan(dY/dX)/math.pi*180
    
    #nova rychlost
    V = math.sqrt((dX/dt)**2+(dY/dt)**2)
    
    #novy drag
    d = -ro*V**2*C_d*A/(2*m)
    #d=-C_d*V**2
    
    #zlozky pociatocneho dragu v smeroch x a y
    d_x = d*math.cos(alfa/180*math.pi)
    d_y = d*math.sin(alfa/180*math.pi)
    #zlozky pociatocnej rychlosti v smeroch x a y
    V_x = V_x + d_x*dt
    V_y = V_y + (d_y+g)*dt
    #kroky v x a y
    dX = V_x*dt+1/2*d_x*dt**2
    dY = V_y*dt+1/2*(d_y+g)*dt**2

# print(x[-1])
# print(y[-1])


# print(y_r)


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


# Vytvorenie novej vrstvy s polygonom so suradnicami coor_point_polygon
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
srs = osr.SpatialReference()
srs.ImportFromEPSG(EPSG)    #definicia referencneho systemu
# srs = point_layer.GetSpatialRef()
srs.MorphToESRI()
file = open(throwshed_output_folder + "\\" + throwshed_file + ".prj", 'w')
file.write(srs.ExportToWkt())
file.close()
print("\n", srs.ExportToPrettyWkt(), "\n")
throwshed_outlayer = throwshed_outds.CreateLayer(throwshed_file, srs)

# pridanie polygonu do feature a jeho ulozenie do vystupnej vrstvy
throwshed_feature = ogr.Feature(throwshed_outlayer.GetLayerDefn())
throwshed_feature.SetGeometry(throwshed_polygon)
throwshed_outlayer.CreateFeature(throwshed_feature)

# nakoniec vrstva, datasource aj prvok treba dat rovne None
throwshed_outds = throwshed_outlayer = throwshed_feature = None


