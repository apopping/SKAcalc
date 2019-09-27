from pykml.factory import KML_ElementMaker as KML
from xml.dom.minidom import parseString
import convert_coords as coords
import numpy as np
import csv
import os

def getCoordinates(kml_file):
    #Read KML file as a string
    file = open(kml_file)
    data = file.read()
    file.close()
 
    #Parse that string into a DOM
    dom = parseString(data)
 
    #initialize latitude and longitude lists
    latitudes = []
    longitudes = []
 
    #Iterate through a collection of coordinates elements
    for d in dom.getElementsByTagName('coordinates'):
        #Break them up into latitude and longitude
        coords = d.firstChild.data.split(',')
        longitudes.append(float(coords[0]))
        latitudes.append(float(coords[1]))

    longitudes = np.array(longitudes)    
    latitudes = np.array(latitudes)    
    
    return longitudes,latitudes



#########################


# the file names"
#kml_file='kml_files/arm_p65s30r0.kml'
#config_file = 'config_files/arm_p65s30r0.csv'

kml_name = '/Users/attila/Documents/Work/SKA/Baseline_Configuration/2014_06_09_configs/ska1-sur-vp.wgs84.96x4.kml'
config_file = '/Users/attila/Documents/Work/SKA/Baseline_Configuration/2014_06_09_configs/ska1-sur-vp.wgs84.96x4_config.csv'

kml_name = '/Users/attila/work/ska/sensitivity_calculator/turbogears/calculator/calculator/data/configurations/coords_LOW.kml'
csv_file = '/Users/attila/work/ska/sensitivity_calculator/turbogears/calculator/calculator/data/configurations/coords_LOW.csv'
cfg_file = '/Users/attila/work/ska/sensitivity_calculator/turbogears/calculator/calculator/data/configurations/coords_LOW.cfg'




interval = 1 # 1 for every position


# Read the longitudes and latitudes from kml file

longitudes,latitudes = getCoordinates(kml_name)

# convert the longitudes and latitudes to x and y positions
X,Y = coords.getXYpos(longitudes[0],latitudes[0],longitudes,latitudes)

# save the data
f = open(csv_file,'w')
for i in range(len(X)):
    print i
    if np.mod(i,interval) == 0:
        f.write(str(X[i]) + ',' + str(Y[i])  + '\n')

f.close

# convert the file to a miriad config file:

f  = open(csv_file, "rb")
reader = csv.reader(f)
configuration = list(reader)
f.close()


f = open(cfg_file,'w')
for i in range(len(configuration)):
    # change the order of x and y, as this is the default for miriad
    #f.write(configuration[i][0] + '   ' + configuration[i][1]  + '\n')
    f.write(configuration[i][1] + '   ' + configuration[i][0]  + '\n')
f.flush()
os.fsync(f.fileno())
f.close()
