
import numpy as np
import matplotlib.pyplot as plt

#line_table_data = np.load('/data/projects/FIRESIDE/data/line_table_cat_v9.npy', allow_pickle=True)



import numpy as np
import matplotlib.pyplot as plt
import astropy as ap
from astropy.visualization import hist
from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord



from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def line_flux_at_z(lineLum, z):
    d=cosmo.luminosity_distance(z).to('m')
    linelum=lineLum.to("W")
    lineFlux=linelum/(4.*np.pi*d**2)
    return lineFlux

def sensScaled (sensitivity, tobsv, fov, area):
    '''
    fov is the field of view of telescope (area) in units of square degrees.

    sensitivity is the sensitivity of the telescope in units of Watts per square meter. 

    tobsv describes the total time spent observing.

    area is the total area of observation. 
    '''
    fov=fov.to(u.degree**2)
    sensitivity=sensitivity.to(u.W/u.m**2)
    tobsv=tobsv.to(u.hr)
    area=area.to(u.degree**2)

    sensScaled=sensitivity*np.sqrt(1*u.hr/(tobsv*(fov)/(area)))
    return sensScaled

def redshifted_wavelengthFunc(z, waveRest, waveMin, waveMax):

    '''
    z is the redshift of the source observed. Z has no units!

    waveRest is the wavelength of emitted light of the source. The function converts all entered values to the units of 'm'.

    waveMin is the wavelength cutoff min value. The function converts all entered values to the units of 'm'.

    waveMax is the wavelength cutoff max value.The function converts all entered values to the units of 'm'.

    The entered units for this function don't matter, as long as they're covertable to a unit of length, e.g., 'nm' to 'm'.
    '''

    waveRest=waveRest.to(u.micron)
    waveMin=waveMin.to(u.micron)
    waveMax=waveMax.to(u.micron)

    redshifted_wavelength=waveRest*(1+z)
    #detect_galaxies = np.where((redshifted_wavelength>=waveMin) & (redshifted_wavelength<=waveMax),1,0)
    return redshifted_wavelength

#resultRedshift = redshifted_wavelengthFunc(line_table_data[0], 158*u.micron, 25*u.micron, 250*u.micron)
#print(resultRedshift)

def ElementInput(Element):
    
    waveRestData= [
   157.74 * u.micron, 145.53 * u.micron, 121.90 * u.micron, 88.36 * u.micron, 63.18 * u.micron,
   57.32 * u.micron, 51.81 * u.micron, 18.70 * u.micron, 15.60 * u.micron, 12.80 * u.micron,
   25.98 * u.micron, 33.48 * u.micron, 34.82 * u.micron, 10.50 * u.micron, 7.70 * u.micron,
   11.30 * u.micron, 4.05 * u.micron, 12.37 * u.micron, 7.46 * u.micron, 
   14.30 * u.micron, 24.31 * u.micron
]
   
    Assigned_Column_Element = {
                                '[CII] 158':11,
                                '[OI] 145':12, 
                                '[NII] 122':13,
                                '[OIII] 88':14,
                                '[OI] 63':15,
                                '[NIII] 57':16, 
                                '[OIII] 52':17, 
                                '[SIII] 19':18,
                                '[NeIII] 16':19,
                                '[NeII] 13':20,
                                '[OIV] 26':21,
                                '[SIII] 33':22,
                                '[SII] 35':23,
                                '[SIV] 11':24,
                                'PAH 7.7':25,
                                'PAH 11.3':26, 
                                'Bralpha':27,
                                'Hmalpha':28,
                                'Pfalpa':29,
                                'Halpha':30,
                                '[NeV] 14':31, 
                                '[NeV] 24':32 
                                }
   
 

    line_index = Assigned_Column_Element[Element]
    waveRest = waveRestData[line_index-11]
    
    return line_index, waveRest
        
#line_index,waveRest=ElementInput('[CII] 158')
#print(line_index, waveRest)  

def Detected_Galaxies(Element, data_table, waveMin, waveMax, sensitivity, tobsv, fov, area):

    line_index,waveRest=ElementInput(Element)
    NewSensScaled=sensScaled(sensitivity, tobsv, fov, area)
    NewLineFlux=line_flux_at_z(data_table[line_index]*u.solLum, data_table[0])
    redshifted_wavelength=redshifted_wavelengthFunc(data_table[0], waveRest, waveMin, waveMax)


    
     #area = 1.96 #in square degrees
    MaxArea=1.96*u.degree**2
    if area > MaxArea:
        ra_length = MaxArea/np.sqrt(MaxArea)
        factor=int(area/MaxArea)
        dec_angle=MaxArea/np.sqrt(MaxArea)
        
        raise ValueError('Sky area exceedes allowed fov. Area is limitied to 0 - 1.96 sq.degrees')
                         
    elif area <= MaxArea:

        factor=1
        ra_length=area/np.sqrt(area)
        dec_angle=area/np.sqrt(area)

    catIndices = np.where((dec_angle>=data_table[2]*u.degree) & 
                          (ra_length>=data_table[1]*u.degree) & 
                          (NewLineFlux>=NewSensScaled) & 
                          (redshifted_wavelength>=waveMin) & 
                          (redshifted_wavelength<=waveMax))[0]

    
        
    FluxOfGalaxiesDetected = NewLineFlux[catIndices]
    
    Total_Galaxies = factor*len(NewLineFlux[catIndices]) 


        
    return line_index, catIndices, Total_Galaxies, FluxOfGalaxiesDetected

    
    
        
    print('Indices:', catIndices, 'Galaxy detections:', Total_Galaxies, 'Flux of Galaxies Det.:', FluxOfGalaxiesDetected) 


#redshift is really row 1, but that corresponds to index 0  
#['REDSHIFT' , 'RA' ,'DEC','MHALO','MSTAR','QFLAG', 'SFR','MU','ISSB','UMEAN','LIR',
#  1            2     3     4        5       6        7    8     9     10       11     
# 'l[CII] 158','l[OI] 145','l[NII] 122','l[OIII] 88','l[OI] 63','l[NIII] 57','l[OIII] 52','l[SIII] 19',
#  12            13           14           15           16          17         18           19
# 'l[NeIII] 16','l[NeII] 13','l[OIV] 26', 'l[SIII] 33','l[SII] 35','l[SIV] 11','lPAH 7.7','lPAH 11.3','lBralpha','lHmalpha','lPfalpha',
#  20             21          22           23           24            25           26        27         28         29         30
#  'lHalpha','l[NeV] 14','l[NeV] 24',
#  31            32       33
#  '[CII] 158','[OI] 145','[NII] 122','[OIII] 88','[OI] 63','[NIII] 57','[OIII] 52','[SIII] 19',
#   34
#  '[NeIII] 16','[NeII] 13','[OIV] 26','[SIII] 33','[SII] 35','[SIV] 11','PAH 7.7','PAH 11.3','Bralpha','Hmalpha','Pfalpha',
#   'Halpha','[NeV] 14', '[NeV] 24'])
#result = Detected_Galaxies(line_table_data[0], 25*u.micron, 250*u.micron, 10**-19*u.W/u.m**2, 1000*u.hr, 1.6*10**-3*u.degree**2, 1.96*u.degree**2, '[CII] 158')
