#!/usr/bin/env python
# coding: utf-8


# Import a lot of stuff
import numpy as np
import pandas as pd
import sys
import astropy.constants as ac
from tqdm import tqdm
import fileinput
from matej_resolution_functions import * # Import Matej's code

# Specify the file that you want read in and the location
# This code will take it down from 200k to 125k
opacity_file_base = 'Hayley-Opacity-Data-Files/'
opacity_old_file  = 'opacCO.dat'
opacity_new_file  = '125k-opacCO.dat'


# This specifies important stuff for the file\
# The resolution can be pretty picky with Matej's code for reasons that I don't know yet
# For each wavelength there will be a matrix that has N presssures and M temperatures
# These are the hardcoded values that I think are necessary for all the files that I have
new_resolution = 125000
n_wavelengths = 176842 # this is the total number
n_pressure_points = 28
n_temperature_points = 46

# Read in the opacity file
# Make it pandas
df = pd.read_csv(opacity_file_base + opacity_old_file,
                 delim_whitespace= True, dtype=np.float64,
                 names=['P',
                        '500', '600', '700', '800', '900',
                        '1000','1100', '1200', '1300', '1400', '1500', '1600','1700', '1800', '1900',
                        '2000', '2100', '2200', '2300', '2400','2500', '2600', '2700', '2800', '2900',
                        '3000', '3100', '3200', '3300', '3400', '3500', '3600', '3700', '3800', '3900',
                        '4000', '4100', '4200', '4300', '4400', '4500', '4600', '4700', '4800', '4900',
                        '5000'],
                         skiprows=2, nrows=n_wavelengths * (n_pressure_points + 1))

# Get all the wavelengths by parsing every Nth + 1 line
wavelengths = list(df[df.index % (n_pressure_points + 1) == 0].P)

# Take out all the rows that only have the wavelength info
data = df[df.index % (n_pressure_points + 1) != 0]

# Get all the pressures from the dataframe
pressures    = list(data.P[0:n_pressure_points])

# Hardcoded in temperatures that are from reading in the file
temperatures = ['500', '600', '700', '800', '900',
                '1000','1100', '1200', '1300', '1400', '1500', '1600','1700', '1800', '1900',
                '2000', '2100', '2200', '2300', '2400','2500', '2600', '2700', '2800', '2900',
                '3000', '3100', '3200', '3300', '3400', '3500', '3600', '3700', '3800', '3900',
                '4000', '4100', '4200', '4300', '4400', '4500', '4600', '4700', '4800', '4900',
                '5000']

# Drop all the pressure values from the matrixes
data = data.drop(['P'], axis=1)

# Reset the index of the pandas dataframe
data = data.reset_index(drop=True)

# Turn the datframe into a numpy array
# SPEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEED
np_data = data.to_numpy()

# Make a dictionary to hold all the different lists
# There should be N x M lists where each list is # of wavelength points long
data_dictionary = {}


# Fill the dictionary where each key corresponds to the opacities for a pressure-temperature pair
# They corresponding lists are all the wavelengths
# So the dictionary key for a pressure of 10 and a temperature of 500 will be a list, and the order 
# of that list corresponds to all the wavelengths needed
for w in range(n_wavelengths):
    for p in range(len(pressures)):
        for t in range(len(temperatures)):
            key_str = str(pressures[p]) + '_' + str(temperatures[t])
            data_dictionary.setdefault(key_str, []).append(np_data[p + (w * 28)][t])

# Create a new numpy array that will be filled with the new opacities
# The array will only be a series of the NxM matrixes because that is how
# Numpy wants it
# Numpy arrays with varying size rows are not easy to save to text files
data_new = np.zeros((n_wavelengths * (n_pressure_points), n_temperature_points + 1))


# Fill the np array
w = 0
for p in tqdm(range(len(pressures))):
    pressure = pressures[p]

    for t in range(len(temperatures)):
        temperature = temperatures[t]
        
        # Lookup the dictionary list by key
        # The key is just the pressure and temperature
        key_string = str(pressure)+ '_' + str(temperature)
        old_opac  = data_dictionary[key_string]
        old_lambda = wavelengths

        # Take the old wavelengths and the old opacities and calculate them at the lower resolution
        # This takes a really long time
        # I think that this is where all the code gets held up
        new_lambda, new_opac = rebin_spectrum_to_resolution(old_lambda, old_opac, new_resolution, w_unit='cm', type='log')

        # Sometimes Matej's code has the first element be 0
        # I put this in to shame myself
        if new_opac[0] == 0.0:
            new_opac[0] = new_opac[1]
        else:
            pass
        
        # Put the new opacity value in the correct place in the dictionary
        # This might break?
        for w in range(len(new_opac)):
            data_new[p + w*28][t + 1] = new_opac[w]

# Add back in the pressure points for the NxM matrix y axis
for w in range(len(new_opac)):
    for p in range(n_pressure_points):
        data_new[p + w*28][0] = pressures[p]

# Drop the extra rows from the numpy array
# This is necessary because the numpy array was made with a lenght of the
# Old number of wavelengths
data_new = data_new[~np.all(data_new == 0, axis=1)]

# Save the data
np.savetxt(opacity_file_base + opacity_new_file, data_new, fmt='%1.6E')

# Add the initial header from the initial file
i = 0
for linenum,line in enumerate(fileinput.FileInput(opacity_file_base + opacity_new_file, inplace=1)):
    if linenum == 0:
        print ("500.000 600.000 700.000 800.000 900.000 1000.000 1100.000 1200.000 1300.000 1400.000 1500.000 1600.000 1700.000 1800.000 1900.000 2000.000 2100.000 2200.000 2300.000 2400.000 2500.000 2600.000 2700.000 2800.000 2900.000 3000.000 3100.000 3200.000 3300.000 3400.000 3500.000 3600.000 3700.000 3800.000 3900.000 4000.000 4100.000 4200.000 4300.000 4400.000 4500.000 4600.000 4700.000 4800.000 4900.000 5000.000 ")
        print ("1.000000E-01 2.154435E-01 4.641589E-01 1.000000E+00 2.154435E+00 4.641589E+00 1.000000E+01 2.154435E+01 4.641589E+01 1.000000E+02 2.154435E+02 4.641589E+02 1.000000E+03 2.154435E+03 4.641589E+03 1.000000E+04 2.154435E+04 4.641589E+04 1.000000E+05 2.154435E+05 4.641589E+05 1.000000E+06 2.154435E+06 4.641589E+06 1.000000E+07 2.154435E+07 4.641589E+07 1.000000E+08")
        
    if linenum%28 == 0:
        print ("{0:.9E}".format(new_lambda[i]))
        print (line.rstrip())
        i = i + 1
    else:
        print (line.rstrip())