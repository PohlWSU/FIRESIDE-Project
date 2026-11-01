
import numpy as np
import matplotlib.pyplot as plt


import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.simplefilter('ignore', category=VerifyWarning)
import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import astropy as ap
from astropy.visualization import hist
from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, QTable, hstack

#AGN_line_table_data = fits.open('/data/projects/FIRESIDE/data/modular_pySIDES_plus_AGN_Q_alpha_0.5(copy).fits')


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
    if not hasattr(waveRest, 'unit') or waveRest.unit is None or waveRest.unit == u.dimensionless_unscaled:
        waveRest = waveRest*u.micron
    else:
        waveRest = waveRest.to(u.micron)
        
    waveMin=waveMin
    waveMax=waveMax

    redshifted_wavelength=waveRest*(1+z)
    #detect_galaxies = np.where((redshifted_wavelength>=waveMin) & (redshifted_wavelength<=waveMax),1,0)
    return redshifted_wavelength

#resultRedshift = redshifted_wavelengthFunc(line_table_data[0], 158*u.micron, 25*u.micron, 250*u.micron)
#print(resultRedshift)
#Load Excel Rest_wave data


def ElementInput(Element, Catalog, rest_wave_table):
    
    if not isinstance(Element, str):
        Element_data = []
        Rest_wave_data = []
        for single_element in Element:
            if (single_element in Catalog[1].columns.names) and (single_element in rest_wave_table.colnames):
                Element_data.append(Catalog[1].data[single_element])
                Rest_wave_data.append(rest_wave_table[single_element])
            elif single_element in Catalog[1].columns.names:
                Element_data.append(Catalog[1].data[single_element])
                print(f"Error:'{single_element}'Rest Wavelength not found")
                Rest_wave_data.append(None)
            elif single_element in rest_wave_table.colnames:
                Rest_wave_data.append(rest_wave_table[single_element])
                print(f"Error:Column'{single_element}'from catalog not found")
                Element_data.append(None)
            else: 
                print(f"Error:Column'{single_element}'from catalog not found")
                print(f"Error:'{single_element}'Rest Wavelength not found")
                Element_data.append(None)
                Rest_wave_data.append(None)
        return Element_data, Rest_wave_data

    else:
        if (Element in Catalog[1].columns.names) and (Element in rest_wave_table.colnames):
            return Catalog[1].data[Element], rest_wave_table[Element]
        elif Element in Catalog[1].columns.names:
            print(f"Error:'{Element}'Rest Wavelength not found")
            return Catalog[1].data[Element], None
        elif Element in rest_wave_table.colnames:
            print(f"Error:Column'{Element}'from catalog not found")
            return None, rest_wave_table[Element]
        else:
            print(f"Error:Column'{Element}'from catalog not found")
            print(f"Error:'{Element}'Rest Wavelength not found")
            return None, None

def Detected_Galaxies(Element, data_table, waveMin, waveMax, sensitivity, tobsv, fov, area, rest_wave_table, upper_f_agn, lower_f_agn, joint):

    z_column_names = ['z', 'redshift', 'Redshift'] #Allowed column headers for z
    correct_z_name = None
    
    for z_name in z_column_names:
        if z_name in data_table[1].columns.names:
            correct_z_name = z_name
            break
    else:
        if correct_z_name is None:
            raise KeyError("Error: Make sure the Redshift column is defined as 'z', 'redshift', or 'Redshift'")


    ra_column_names = ['ra', 'Ra', 'RA', 'right ascension', 'Right ascension', 'Right Ascension'] #Allowed column headers for ra
    correct_ra_name = None
    
    for ra_name in ra_column_names:
        if ra_name in data_table[1].columns.names:
            correct_ra_name = ra_name
            break
    else:
        if correct_ra_name is None:
            raise KeyError("Error: Make sure the Right Ascension column is defined as 'ra', 'Ra', 'RA', 'right ascension', 'Right ascension', or 'Right Ascension'")
        
    dec_column_names = ['dec', 'DEC', 'δ', 'declination', 'Declination' ] #Allowed column headers for dec
    correct_dec_name = None
    
    for dec_name in dec_column_names:
        if dec_name in data_table[1].columns.names:
            correct_dec_name = dec_name
            break
    else:
        if correct_dec_name is None:
            raise KeyError("Error: Make sure the Declination column is defined as 'dec', 'DEC', 'δ', 'declination', or 'Declination'")

    f_agn_column_names = ['f_agn', 'f_AGN', 'F_agn', 'agn_frac', 'AGN_frac', 'AGN_FRAC', 'agn_FRAC', 'F_AGN', 'AGN_fract', 'AGN_FRACT', 'agn_FRACT', 'agn_fract', 'fagn', 'FAGN', 'fAGN', 'Fagn'] #Allowed column headers for f_agn
    correct_f_agn_name = None
    
    for f_agn_name in f_agn_column_names:
        if f_agn_name in data_table[1].columns.names:
            correct_f_agn_name = f_agn_name
            break
    else:
        if correct_f_agn_name is None:
            raise KeyError("Error: Make sure the f_agn column is defined as 'f_agn', 'f_AGN', 'F_agn', 'agn_frac', 'AGN_frac', 'AGN_FRAC', 'agn_FRAC', 'F_AGN', 'AGN_fract', 'AGN_FRACT', 'agn_FRACT', 'agn_fract', 'fagn', 'FAGN', 'fAGN', 'Fagn'")
    #Global calculations for the single & multi line emission loop

    #area of catalog
    ra_min = np.min(data_table[1].data[correct_ra_name])
    ra_max = np.max(data_table[1].data[correct_ra_name])
    dec_min = np.min(data_table[1].data[correct_dec_name])
    dec_max = np.max(data_table[1].data[correct_dec_name])

    #width and height spans
    ra_span = ra_max-ra_min
    dec_span = dec_max-dec_min

    #Restrict MaxArea of catalog
    MaxArea = (ra_span*dec_span)* u.degree**2

    #strips units to compare
    area_value = area.value if hasattr(area, 'value') else area
    MaxArea_value = MaxArea.value if hasattr(MaxArea, 'value') else MaxArea

    if area_value > MaxArea_value:
        #factor if maxArea is exceeded
        ra_length = np.sqrt(MaxArea_value)
        dec_angle = np.sqrt(MaxArea_value)
        factor = area_value/MaxArea_value
                      
    else:
            
        #Within bounds of catalog
        ra_length = np.sqrt(area_value)
        dec_angle = np.sqrt(area_value)
        factor = 1.0
            
    dec_angle = dec_angle.value if hasattr(dec_angle, 'value') else dec_angle
    ra_length = ra_length.value if hasattr(ra_length, 'value') else ra_length
    NewSensScaled=sensScaled(sensitivity, tobsv, fov, area)
    instrument_sens = NewSensScaled.value if hasattr(NewSensScaled, 'value') else NewSensScaled

    w_min = waveMin.value if hasattr(waveMin, 'value') else waveMin
    w_max = waveMax.value if hasattr(waveMax, 'value') else waveMax
    up_f_agn = upper_f_agn.value if hasattr(upper_f_agn, 'value') else upper_f_agn
    low_f_agn = lower_f_agn.value if hasattr(lower_f_agn, 'value') else lower_f_agn

    up_f_agn = float(up_f_agn) if up_f_agn is not None else 1.0
    low_f_agn = float(low_f_agn) if low_f_agn is not None else 0.0

    #caps AGN fraction
    if up_f_agn < 0:
        up_f_agn = 0.0
    elif up_f_agn > 1:
        up_f_agn = 1.0

    if low_f_agn < 0:
        low_f_agn = 0.0
    elif low_f_agn > 1:
        low_f_agn = 1.0
    
    #Actual data filter to grab detected galaxies & their indices
    if not isinstance(Element, str):
        Df_list = []
        dash_list = []
        for multi_elem in Element:
            print(f"Processing element: {multi_elem}")
                        
            Element_data, Rest_wave_data = ElementInput(multi_elem, data_table, rest_wave_table)
            
            if Element_data is None or Rest_wave_data is None:
                print(f"Emission line: {multi_elem} is missing from column and was not processed")
                dash_list.append([])
                Df_list.append(pd.DataFrame(columns=["catIndices"]))
                continue
            
            NewLineFlux = line_flux_at_z(Element_data*u.solLum, data_table[1].data[correct_z_name]) 
            redshifted_wavelength = redshifted_wavelengthFunc(data_table[1].data[correct_z_name], Rest_wave_data, waveMin, waveMax)
        
            line_flux = NewLineFlux.value if hasattr(NewLineFlux, 'value') else NewLineFlux
            z = redshifted_wavelength.value if hasattr(redshifted_wavelength, 'value') else redshifted_wavelength
            
            if correct_f_agn_name is None or correct_f_agn_name not in data_table[1].data.names:
                catIndices = np.where(((dec_angle + dec_min) >= data_table[1].data[correct_dec_name]) &   
                ((ra_length + ra_min) >= data_table[1].data[correct_ra_name]) & 
                (line_flux >= instrument_sens) & 
                (z >= w_min) & 
                (z <= w_max))[0]
            else:     
                catIndices = np.where(((dec_angle + dec_min) >= data_table[1].data[correct_dec_name]) &   
                ((ra_length + ra_min) >= data_table[1].data[correct_ra_name]) & 
                (line_flux >= instrument_sens) & 
                (z >= w_min) & 
                (z <= w_max) & (data_table[1].data[correct_f_agn_name] <= up_f_agn) & 
                (data_table[1].data[correct_f_agn_name] >= low_f_agn))[0]
                
            # Always append the dataframe and list to preserve strict index mapping
            Df = pd.DataFrame({"catIndices": catIndices})
            Df_list.append(Df)
            dash_list.append(catIndices.tolist())

        if joint: 
            if Df_list:
                full_table = Df_list[0]
                for next_table in Df_list[1:]:
                    full_table = pd.merge(full_table, next_table, on='catIndices', how='inner')
                joint_indices = full_table['catIndices'].tolist()
                return joint_indices
        else:
            return dash_list

    #Use .data to pull actual data from catalog since Dash passes FITS as a whole (not just a single data layer)
    else: 
        print(f"Processing element: {Element}")   
        Element_data, Rest_wave_data=ElementInput(Element, data_table, rest_wave_table)
        NewLineFlux=line_flux_at_z(Element_data*u.solLum, data_table[1].data[correct_z_name]) 
        redshifted_wavelength=redshifted_wavelengthFunc(data_table[1].data[correct_z_name], Rest_wave_data, waveMin, waveMax)

        line_flux = NewLineFlux.value if hasattr(NewLineFlux, 'value') else NewLineFlux
        z = redshifted_wavelength.value if hasattr(redshifted_wavelength, 'value') else redshifted_wavelength

    #catIndices pulls the actual detected galaxy & their indices from the data table
        if correct_f_agn_name is None or correct_f_agn_name not in data_table[1].data.names:
            catIndices = np.where(((dec_angle + dec_min) >= data_table[1].data[correct_dec_name]) &   
            ((ra_length + ra_min) >= data_table[1].data[correct_ra_name]) & 
            (line_flux >= instrument_sens) & 
            (z >= w_min) & 
            (z <= w_max))[0]
            
        else: 
        
            catIndices = np.where(((dec_angle + dec_min) >= data_table[1].data[correct_dec_name]) &   
                ((ra_length + ra_min) >= data_table[1].data[correct_ra_name]) & 
                (line_flux >= instrument_sens) & 
                (z >= w_min) & 
                (z <= w_max) & (data_table[1].data[correct_f_agn_name] <= up_f_agn) & 
                (data_table[1].data[correct_f_agn_name] >= low_f_agn))[0]

        Df = pd.DataFrame({"catIndices": catIndices})
        
        return catIndices.tolist()
