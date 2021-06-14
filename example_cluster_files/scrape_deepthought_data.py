"""
Runs through all the output folders and stitches together the spectra correctly.

author: @arjunsavel
"""

import pdb
from glob import glob
from RT import *
import os
from tqdm import tqdm
from sklearn.model_selection import ParameterGrid
import numpy as np
import pandas as pd
import pdb
from run_RT_deepthought import construct_phases


def scrape_deepthought_data():
	"""
	Scrapes output data. Does *not* clean overlap.
	"""
    def read_file(file):
    	"""
		Reads a spectral output file.

		Inputs
		-------
			:file: (str) path to file to be read.

		Output
		------
			:dat: (pd.DataFrame) Output data. Nonetype if broken.
    	"""
        try:
            dat = pd.read_csv(file, names=['wav', 'depth'], delimiter='\t')
            return dat
        except:
            return None
    

    ##################################### Set up parameter grid ############################
    drag_timescales = ['1e3', '1e4', '1e5', '1e6', '1e7']

    models = ['', '_deep_rcb']
    phases = construct_phases()

    change_elems_list = ['Na', 'Mg', 'Ca_plus', 'Mn', 'K', 'Ti', 'TiO', 'VO', 'H2O', 'Li', 'Cr']


    
    parameter_list = list(ParameterGrid({'drag_timescale': drag_timescales,
               'model': models,
               'change_elems': change_elems_list}))

    ##################################### End parameter grid setup ############################


    os.chdir("RT_3D_Transmission_Code")
    
    dir_dict = {} # dictionary with all the directories as keys and their corr. spectra as values

    for directory in glob('simulation*'):
        if 'transmission.dat' in os.listdir(directory):
            if 'no_doppler' in directory:
                spectrum_num = directory.split('_')[-3]
            else:
                spectrum_num = directory.split('_')[-1]
            dir_dict[directory] = spectrum_num

    unique_spectra = np.unique(list(dir_dict.values()))
    
    # done this way to skip over missing data
    for spectrum in tqdm(unique_spectra):
        parameters = parameter_list[eval(spectrum)]
        spectrum_dirs = [key for key in dir_dict.keys() if dir_dict[key] == spectrum]
        
        # instantiate an RT object
        obj = RT()
        j = 0

        # loop through the directories and link them together
        for i, directory in enumerate(spectrum_dirs):
            data_file = os.path.join(directory, 'transmission.dat')

            # the first file just gets read in, and input file parsed
            if j == 0:
                
                data = read_file(data_file)
                if data is None:
                    continue
                data = data[['wav', 'depth']]
                os.chdir(directory)
                obj.parse_input_file()
                os.chdir('..')
                j += 1

            # subsequent files get appended
            else:
                new_data = read_file(data_file)
                if new_data is None:
                    continue
                data = data.append(new_data[['wav', 'depth']], ignore_index=True, sort=True)

        data[['wav', 'depth']].to_csv('full_transmission.dat')
        

        ########### access grid parameters, assign to RT object, and save to pickle ###########
        
        output = pd.read_csv('full_transmission.dat')
        obj['wav'] = output.wav
        obj['depth'] = output.depth
        
#         model = parameters['model']
#         drag_timescale = parameters['drag_timescale']
#         phase = parameters['phase']

        
#         uniform_los_wind = parameters['uniform_los_wind']
        drag_timescale = parameters['drag_timescale']
        change_elems = parameters['change_elem']
#     phase = parameters['phase']
        model = parameters['model']
        
        obj.condensation = False
        obj.doppler = False
        obj.abundances = change_elems
        obj.scattering = ['H2', 'He', 'H2O',  'CO', 'CO2', 'CH4', 'NH3']
        obj.overlapped_abundances = True
        obj.comments = 'First set of full-spectrum template runs'
    	obj.to_pickle(f'fullspectra/deepthought2/{model}_{drag_timescale}_{change_elems}_doppler_off.pkl')
        
                    
        
if __name__ == "__main__":
    scrape_deepthought_data()
