# concevert between longitude and latitude and meters
import numpy as np



def asRadians(degrees):
    return degrees * np.pi / 180

def getXYpos(lon0,lat0,lon1,lat1):
    """ Calculates X and Y distances in meters.
    """
    deltaLatitude = lat1 - lat0
    deltaLongitude = lon1 - lon0
    latitudeCircumference = 40075160 * np.cos(asRadians(lat0))
    resultX = deltaLongitude * latitudeCircumference / 360
    resultY = deltaLatitude * 40008000 / 360
    return resultX, resultY


def getLonLatpos(lon0,lat0,X,Y):
    """ Calculates Longitude and Latitude in degrees
    """
    latitudeCircumference = 40075160 * np.cos(asRadians(lat0))
    deltaLongitude = X * 360 / latitudeCircumference
    deltaLatitude = Y * 360 /  40008000
    Longitude = deltaLongitude + lon0
    Latitude = deltaLatitude + lat0
    return Longitude, Latitude


def rotateXYpos(X,Y,degrees):
    """ Rotates the xy positions
    """
    z = np.sqrt(X**2 + Y**2)
    degrees_rad = degrees * np.pi / 180
    if X == 0:
        theta = 0
    else:
        theta = np.arctan(Y/X)
    if X < 0:
        theta = theta + np.pi
    theta2 = theta + degrees_rad
    X2 = np.cos(theta2) * z
    Y2 = np.sin(theta2) * z
    return X2,Y2
        
    
